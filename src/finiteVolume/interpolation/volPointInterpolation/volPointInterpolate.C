/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "volPointInterpolation.H"
#include "volFields.H"
#include "pointFields.H"
#include "emptyFvPatch.H"
#include "mapDistribute.H"
#include "coupledPointPatchField.H"
#include "valuePointPatchField.H"
#include "transform.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class CombineOp>
void volPointInterpolation::syncUntransformedData
(
    List<Type>& pointData,
    const CombineOp& cop
) const
{
    // Transfer onto coupled patch
    const globalMeshData& gmd = mesh().globalData();
    const indirectPrimitivePatch& cpp = gmd.coupledPatch();
    const labelList& meshPoints = cpp.meshPoints();

    const mapDistribute& slavesMap = gmd.globalCoPointSlavesMap();
    const labelListList& slaves = gmd.globalCoPointSlaves();

    List<Type> elems(slavesMap.constructSize());
    forAll(meshPoints, i)
    {
        elems[i] = pointData[meshPoints[i]];
    }

    // Pull slave data onto master. No need to update transformed slots.
    slavesMap.distribute(elems, false);

    // Combine master data with slave data
    forAll(slaves, i)
    {
        Type& elem = elems[i];

        const labelList& slavePoints = slaves[i];

        // Combine master with untransformed slave data
        forAll(slavePoints, j)
        {
            cop(elem, elems[slavePoints[j]]);
        }

        // Copy result back to slave slots
        forAll(slavePoints, j)
        {
            elems[slavePoints[j]] = elem;
        }
    }

    // Push slave-slot data back to slaves
    slavesMap.reverseDistribute(elems.size(), elems, false);

    // Extract back onto mesh
    forAll(meshPoints, i)
    {
        pointData[meshPoints[i]] = elems[i];
    }
}


template<class Type>
void volPointInterpolation::pushUntransformedData
(
    List<Type>& pointData
) const
{
    // Transfer onto coupled patch
    const globalMeshData& gmd = mesh().globalData();
    const indirectPrimitivePatch& cpp = gmd.coupledPatch();
    const labelList& meshPoints = cpp.meshPoints();

    const mapDistribute& slavesMap = gmd.globalCoPointSlavesMap();
    const labelListList& slaves = gmd.globalCoPointSlaves();

    List<Type> elems(slavesMap.constructSize());
    forAll(meshPoints, i)
    {
        elems[i] = pointData[meshPoints[i]];
    }

    // Combine master data with slave data
    forAll(slaves, i)
    {
        const labelList& slavePoints = slaves[i];

        // Copy master data to slave slots
        forAll(slavePoints, j)
        {
            elems[slavePoints[j]] = elems[i];
        }
    }

    // Push slave-slot data back to slaves
    slavesMap.reverseDistribute(elems.size(), elems, false);

    // Extract back onto mesh
    forAll(meshPoints, i)
    {
        pointData[meshPoints[i]] = elems[i];
    }
}


template<class Type>
void volPointInterpolation::addSeparated
(
    GeometricField<Type, pointPatchField, pointMesh>& pf
) const
{
    if (debug)
    {
        Pout<< "volPointInterpolation::addSeparated" << endl;
    }

    forAll(pf.boundaryField(), patchI)
    {
        if (pf.boundaryField()[patchI].coupled())
        {
            refCast<coupledPointPatchField<Type> >
                (pf.boundaryField()[patchI]).initSwapAddSeparated
                (
                    Pstream::nonBlocking,
                    pf.internalField()
                );
        }
    }

    // Block for any outstanding requests
    Pstream::waitRequests();

    forAll(pf.boundaryField(), patchI)
    {
        if (pf.boundaryField()[patchI].coupled())
        {
            refCast<coupledPointPatchField<Type> >
                (pf.boundaryField()[patchI]).swapAddSeparated
                (
                    Pstream::nonBlocking,
                    pf.internalField()
                );
        }
    }
}


template<class Type>
void volPointInterpolation::interpolateInternalField
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    GeometricField<Type, pointPatchField, pointMesh>& pf
) const
{
    if (debug)
    {
        Pout<< "volPointInterpolation::interpolateInternalField("
            << "const GeometricField<Type, fvPatchField, volMesh>&, "
            << "GeometricField<Type, pointPatchField, pointMesh>&) : "
            << "interpolating field from cells to points"
            << endl;
    }

    const labelListList& pointCells = vf.mesh().pointCells();

    // Multiply volField by weighting factor matrix to create pointField
    forAll(pointCells, pointi)
    {
        if (!isPatchPoint_[pointi])
        {
            const scalarList& pw = pointWeights_[pointi];
            const labelList& ppc = pointCells[pointi];

            pf[pointi] = pTraits<Type>::zero;

            forAll(ppc, pointCelli)
            {
                pf[pointi] += pw[pointCelli]*vf[ppc[pointCelli]];
            }
        }
    }
}


template<class Type>
tmp<Field<Type> > volPointInterpolation::flatBoundaryField
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = vf.mesh();
    const fvBoundaryMesh& bm = mesh.boundary();

    tmp<Field<Type> > tboundaryVals
    (
        new Field<Type>(mesh.nFaces()-mesh.nInternalFaces())
    );
    Field<Type>& boundaryVals = tboundaryVals();

    forAll(vf.boundaryField(), patchI)
    {
        label bFaceI = bm[patchI].patch().start() - mesh.nInternalFaces();

        if
        (
           !isA<emptyFvPatch>(bm[patchI])
        && !vf.boundaryField()[patchI].coupled()
        )
        {
            SubList<Type>
            (
                boundaryVals,
                vf.boundaryField()[patchI].size(),
                bFaceI
            ).assign(vf.boundaryField()[patchI]);
        }
        else
        {
            const polyPatch& pp = bm[patchI].patch();

            forAll(pp, i)
            {
                boundaryVals[bFaceI++] = pTraits<Type>::zero;
            }
        }
    }

    return tboundaryVals;
}


template<class Type>
void volPointInterpolation::interpolateBoundaryField
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    GeometricField<Type, pointPatchField, pointMesh>& pf,
    const bool overrideFixedValue
) const
{
    const primitivePatch& boundary = boundaryPtr_();

    Field<Type>& pfi = pf.internalField();

    // Get face data in flat list
    tmp<Field<Type> > tboundaryVals(flatBoundaryField(vf));
    const Field<Type>& boundaryVals = tboundaryVals();


    // Do points on 'normal' patches from the surrounding patch faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(boundary.meshPoints(), i)
    {
        label pointI = boundary.meshPoints()[i];

        if (isPatchPoint_[pointI])
        {
            const labelList& pFaces = boundary.pointFaces()[i];
            const scalarList& pWeights = boundaryPointWeights_[i];

            Type& val = pfi[pointI];

            val = pTraits<Type>::zero;
            forAll(pFaces, j)
            {
                if (boundaryIsPatchFace_[pFaces[j]])
                {
                    val += pWeights[j]*boundaryVals[pFaces[j]];
                }
            }
        }
    }

    // Sum collocated contributions
    syncUntransformedData(pfi, plusEqOp<Type>());

    // And add separated contributions
    addSeparated(pf);

    // Push master data to slaves. It is possible (not sure how often) for
    // a coupled point to have its master on a different patch so
    // to make sure just push master data to slaves.
    pushUntransformedData(pfi);



    if (overrideFixedValue)
    {
        forAll(pf.boundaryField(), patchI)
        {
            pointPatchField<Type>& ppf = pf.boundaryField()[patchI];

            if (isA<valuePointPatchField<Type> >(ppf))
            {
                refCast<valuePointPatchField<Type> >(ppf) =
                    ppf.patchInternalField();
            }
        }
    }


    // Override constrained pointPatchField types with the constraint value.
    // This relys on only constrained pointPatchField implementing the evaluate
    // function
    pf.correctBoundaryConditions();

    // Sync any dangling points
    //mesh().globalData().syncPointData(pfi, nopEqOp<Type>());
    pushUntransformedData(pfi);

    // Apply multiple constraints on edge/corner points
    applyCornerConstraints(pf);
}


template<class Type>
void volPointInterpolation::applyCornerConstraints
(
    GeometricField<Type, pointPatchField, pointMesh>& pf
) const
{
    forAll(patchPatchPointConstraintPoints_, pointi)
    {
        pf[patchPatchPointConstraintPoints_[pointi]] = transform
        (
            patchPatchPointConstraintTensors_[pointi],
            pf[patchPatchPointConstraintPoints_[pointi]]
        );
    }
}


template<class Type>
void volPointInterpolation::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    GeometricField<Type, pointPatchField, pointMesh>& pf
) const
{
    if (debug)
    {
        Pout<< "volPointInterpolation::interpolate("
            << "const GeometricField<Type, fvPatchField, volMesh>&, "
            << "GeometricField<Type, pointPatchField, pointMesh>&) : "
            << "interpolating field from cells to points"
            << endl;
    }

    interpolateInternalField(vf, pf);

    // Interpolate to the patches preserving fixed value BCs
    interpolateBoundaryField(vf, pf, false);
}


template<class Type>
tmp<GeometricField<Type, pointPatchField, pointMesh> >
volPointInterpolation::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const wordList& patchFieldTypes
) const
{
    const pointMesh& pm = pointMesh::New(vf.mesh());

    // Construct tmp<pointField>
    tmp<GeometricField<Type, pointPatchField, pointMesh> > tpf
    (
        new GeometricField<Type, pointPatchField, pointMesh>
        (
            IOobject
            (
                "volPointInterpolate(" + vf.name() + ')',
                vf.instance(),
                pm.thisDb()
            ),
            pm,
            vf.dimensions(),
            patchFieldTypes
        )
    );

    interpolateInternalField(vf, tpf());

    // Interpolate to the patches overriding fixed value BCs
    interpolateBoundaryField(vf, tpf(), true);

    return tpf;
}


template<class Type>
tmp<GeometricField<Type, pointPatchField, pointMesh> >
volPointInterpolation::interpolate
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf,
    const wordList& patchFieldTypes
) const
{
    // Construct tmp<pointField>
    tmp<GeometricField<Type, pointPatchField, pointMesh> > tpf =
        interpolate(tvf(), patchFieldTypes);
    tvf.clear();
    return tpf;
}


template<class Type>
tmp<GeometricField<Type, pointPatchField, pointMesh> >
volPointInterpolation::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const pointMesh& pm = pointMesh::New(vf.mesh());

    tmp<GeometricField<Type, pointPatchField, pointMesh> > tpf
    (
        new GeometricField<Type, pointPatchField, pointMesh>
        (
            IOobject
            (
                "volPointInterpolate(" + vf.name() + ')',
                vf.instance(),
                pm.thisDb()
            ),
            pm,
            vf.dimensions()
        )
    );

    interpolateInternalField(vf, tpf());
    interpolateBoundaryField(vf, tpf(), false);

    return tpf;
}


template<class Type>
tmp<GeometricField<Type, pointPatchField, pointMesh> >
volPointInterpolation::interpolate
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf
) const
{
    // Construct tmp<pointField>
    tmp<GeometricField<Type, pointPatchField, pointMesh> > tpf =
        interpolate(tvf());
    tvf.clear();
    return tpf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
