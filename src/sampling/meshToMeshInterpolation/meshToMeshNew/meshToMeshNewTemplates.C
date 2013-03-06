/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fvMesh.H"
#include "volFields.H"
//#include "ops.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    //- Helper class for list
    template<class Type>
    class ListPlusEqOp
    {
        public:
        void operator()(List<Type>& x, const List<Type> y) const
        {
            if (y.size())
            {
                if (x.size())
                {
                    label sz = x.size();
                    x.setSize(sz + y.size());
                    forAll(y, i)
                    {
                        x[sz++] = y[i];
                    }
                }
                else
                {
                    x = y;
                }
            }
        }
    };

    //- Combine operator for maps/interpolations
    template<class Type, class CombineOp>
    class combineBinaryOp
    {
        const CombineOp& cop_;

        public:

            combineBinaryOp(const CombineOp& cop)
            :
                cop_(cop)
            {}

            void operator()
            (
                Type& x,
                const label faceI,
                const Type& y,
                const scalar weight
            ) const
            {
                cop_(x, weight*y);
            }
    };
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
void Foam::meshToMeshNew::add
(
    UList<Type>& fld,
    const label offset
) const
{
    forAll(fld, i)
    {
        fld[i] += offset;
    }
}


template<class Type, class CombineOp>
void Foam::meshToMeshNew::mapSrcToTgt
(
    const UList<Type>& srcField,
    const CombineOp& cop,
    List<Type>& result
) const
{
    if (result.size() != tgtToSrcCellAddr_.size())
    {
        FatalErrorIn
        (
            "void Foam::meshToMeshNew::mapSrcToTgt"
            "("
                "const UList<Type>&, "
                "const CombineOp&, "
                "List<Type>&"
            ") const"
        )   << "Supplied field size is not equal to target mesh size" << nl
            << "    source mesh   = " << srcToTgtCellAddr_.size() << nl
            << "    target mesh   = " << tgtToSrcCellAddr_.size() << nl
            << "    supplied field = " << result.size()
            << abort(FatalError);
    }

    combineBinaryOp<Type, CombineOp> cbop(cop);

    if (singleMeshProc_ == -1)
    {
        const mapDistribute& map = srcMapPtr_();

        List<Type> work(srcField);
        map.distribute(work);

        forAll(result, cellI)
        {
            const labelList& srcAddress = tgtToSrcCellAddr_[cellI];
            const scalarList& srcWeight = tgtToSrcCellWght_[cellI];

            if (srcAddress.size())
            {
//                result[cellI] = pTraits<Type>::zero;
                result[cellI] *= (1.0 - sum(srcWeight));
                forAll(srcAddress, i)
                {
                    label srcI = srcAddress[i];
                    scalar w = srcWeight[i];
                    cbop(result[cellI], cellI, work[srcI], w);
                }
            }
        }
    }
    else
    {
        forAll(result, cellI)
        {
            const labelList& srcAddress = tgtToSrcCellAddr_[cellI];
            const scalarList& srcWeight = tgtToSrcCellWght_[cellI];

            if (srcAddress.size())
            {
//                result[cellI] = pTraits<Type>::zero;
                result[cellI] *= (1.0 - sum(srcWeight));
                forAll(srcAddress, i)
                {
                    label srcI = srcAddress[i];
                    scalar w = srcWeight[i];
                    cbop(result[cellI], cellI, srcField[srcI], w);
                }
            }
        }
    }
}


template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type> > Foam::meshToMeshNew::mapSrcToTgt
(
    const Field<Type>& srcField,
    const CombineOp& cop
) const
{
    tmp<Field<Type> > tresult
    (
        new Field<Type>
        (
            tgtToSrcCellAddr_.size(),
            pTraits<Type>::zero
        )
    );

    mapSrcToTgt(srcField, cop, tresult());

    return tresult;
}


template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type> > Foam::meshToMeshNew::mapSrcToTgt
(
    const tmp<Field<Type> >& tsrcField,
    const CombineOp& cop
) const
{
    return mapSrcToTgt(tsrcField(), cop);
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::meshToMeshNew::mapSrcToTgt
(
    const Field<Type>& srcField
) const
{
    return mapSrcToTgt(srcField, plusEqOp<Type>());
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::meshToMeshNew::mapSrcToTgt
(
    const tmp<Field<Type> >& tsrcField
) const
{
    return mapSrcToTgt(tsrcField());
}


template<class Type, class CombineOp>
void Foam::meshToMeshNew::mapTgtToSrc
(
    const UList<Type>& tgtField,
    const CombineOp& cop,
    List<Type>& result
) const
{
    if (result.size() != srcToTgtCellAddr_.size())
    {
        FatalErrorIn
        (
            "void Foam::meshToMeshNew::mapTgtToSrc"
            "("
                "const UList<Type>&, "
                "const CombineOp&, "
                "List<Type>&"
            ") const"
        )   << "Supplied field size is not equal to source mesh size" << nl
            << "    source mesh   = " << srcToTgtCellAddr_.size() << nl
            << "    target mesh   = " << tgtToSrcCellAddr_.size() << nl
            << "    supplied field = " << result.size()
            << abort(FatalError);
    }

    combineBinaryOp<Type, CombineOp> cbop(cop);

    if (singleMeshProc_ == -1)
    {
        const mapDistribute& map = tgtMapPtr_();

        List<Type> work(tgtField);
        map.distribute(work);

        forAll(result, cellI)
        {
            const labelList& tgtAddress = srcToTgtCellAddr_[cellI];
            const scalarList& tgtWeight = srcToTgtCellWght_[cellI];

            if (tgtAddress.size())
            {
//                result[cellI] = pTraits<Type>::zero;
                result[cellI] *= (1.0 - sum(tgtWeight));
                forAll(tgtAddress, i)
                {
                    label tgtI = tgtAddress[i];
                    scalar w = tgtWeight[i];
                    cbop(result[cellI], cellI, work[tgtI], w);
                }
            }
        }
    }
    else
    {
        forAll(result, cellI)
        {
            const labelList& tgtAddress = srcToTgtCellAddr_[cellI];
            const scalarList& tgtWeight = srcToTgtCellWght_[cellI];

            if (tgtAddress.size())
            {
//                result[cellI] = pTraits<Type>::zero;
                result[cellI] *= (1.0 - sum(tgtWeight));
                forAll(tgtAddress, i)
                {
                    label tgtI = tgtAddress[i];
                    scalar w = tgtWeight[i];
                    cbop(result[cellI], cellI, tgtField[tgtI], w);
                }
            }
        }
    }
}


template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type> > Foam::meshToMeshNew::mapTgtToSrc
(
    const Field<Type>& tgtField,
    const CombineOp& cop
) const
{
    tmp<Field<Type> > tresult
    (
        new Field<Type>
        (
            srcToTgtCellAddr_.size(),
            pTraits<Type>::zero
        )
    );

    mapTgtToSrc(tgtField, cop, tresult());

    return tresult;
}


template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type> > Foam::meshToMeshNew::mapTgtToSrc
(
    const tmp<Field<Type> >& ttgtField,
    const CombineOp& cop
) const
{
    return mapTgtToSrc(ttgtField(), cop);
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::meshToMeshNew::mapTgtToSrc
(
    const Field<Type>& tgtField
) const
{
    return mapTgtToSrc(tgtField, plusEqOp<Type>());
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::meshToMeshNew::mapTgtToSrc
(
    const tmp<Field<Type> >& ttgtField
) const
{
    return mapTgtToSrc(ttgtField(), plusEqOp<Type>());
}


template<class Type, class CombineOp>
void Foam::meshToMeshNew::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& field,
    const CombineOp& cop,
    GeometricField<Type, fvPatchField, volMesh>& result,
    const bool interpPatches
) const
{
    const fvMesh& mesh = field.mesh();

    if (mesh.name() == srcRegionName_)
    {
        mapSrcToTgt(field, cop, result.internalField());
    }
    else if (mesh.name() == tgtRegionName_)
    {
        mapTgtToSrc(field, cop, result.internalField());
    }
    else
    {
        FatalErrorIn
        (
            "void Foam::meshToMeshNew::interpolate"
            "("
                "const GeometricField<Type, fvPatchField, volMesh>&, "
                "const CombineOp&, "
                "GeometricField<Type, fvPatchField, volMesh>&, "
                "const bool"
            ") const"
        )
            << "Supplied field " << field.name() << " did not originate from "
            << "either the source or target meshes used to create this "
            << "interpolation object"
            << abort(FatalError);
    }

    if (interpPatches)
    {
        switch (method_)
        {
            case imMap:
            {
                result.boundaryField() == field.boundaryField();
                break;
            }
            default:
            {
                notImplemented
                (
                    "void Foam::meshToMeshNew::interpolate"
                    "("
                        "const GeometricField<Type, fvPatchField, volMesh>&, "
                        "const CombineOp&, "
                        "GeometricField<Type, fvPatchField, volMesh>&, "
                        "const bool"
                    ") const - non-conformal patches"
                )

                // do something...
            }
        }
    }
}


template<class Type, class CombineOp>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
Foam::meshToMeshNew::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& field,
    const CombineOp& cop,
    const bool interpPatches
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    const fvMesh& mesh = field.mesh();

    tmp<fieldType> tresult;

    if (mesh.name() == srcRegionName_)
    {
        const fvMesh& tgtMesh =
            mesh.time().lookupObject<fvMesh>(tgtRegionName_);

        tresult =
            new fieldType
            (
                IOobject
                (
                    type() + "::interpolate(" + field.name() + ")",
                    tgtMesh.time().timeName(),
                    tgtMesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                tgtMesh,
                dimensioned<Type>
                (
                    "zero",
                    field.dimensions(),
                    pTraits<Type>::zero
                )
            );

         interpolate(field, cop, tresult(), interpPatches);
    }
    else if (mesh.name() == tgtRegionName_)
    {
        const fvMesh& srcMesh =
            mesh.time().lookupObject<fvMesh>(srcRegionName_);

        tresult =
            new fieldType
            (
                IOobject
                (
                    type() + "::interpolate(" + field.name() + ")",
                    srcMesh.time().timeName(),
                    srcMesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                srcMesh,
                dimensioned<Type>
                (
                    "zero",
                    field.dimensions(),
                    pTraits<Type>::zero
                )
            );

         interpolate(field, cop, tresult(), interpPatches);
    }

    return tresult;
}


template<class Type, class CombineOp>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
Foam::meshToMeshNew::interpolate
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tfield,
    const CombineOp& cop,
    const bool interpPatches
) const
{
    return
        interpolate
        (
            tfield(),
            combineBinaryOp<Type, CombineOp>(cop),
            interpPatches
        );
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
Foam::meshToMeshNew::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& field,
    const bool interpPatches
) const
{
    return interpolate(field, plusEqOp<Type>(), interpPatches);
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
Foam::meshToMeshNew::interpolate
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tfield,
    const bool interpPatches
) const
{
    return interpolate(tfield(), plusEqOp<Type>(), interpPatches);
}


// ************************************************************************* //
