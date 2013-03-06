/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "CoEulerDdtScheme.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<volScalarField> CoEulerDdtScheme<Type>::CorDeltaT() const
{
    const surfaceScalarField cofrDeltaT(CofrDeltaT());

    tmp<volScalarField> tcorDeltaT
    (
        new volScalarField
        (
            IOobject
            (
                "CorDeltaT",
                cofrDeltaT.instance(),
                mesh()
            ),
            mesh(),
            dimensionedScalar("CorDeltaT", cofrDeltaT.dimensions(), 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    volScalarField& corDeltaT = tcorDeltaT();

    const labelUList& owner = mesh().owner();
    const labelUList& neighbour = mesh().neighbour();

    forAll(owner, faceI)
    {
        corDeltaT[owner[faceI]] =
            max(corDeltaT[owner[faceI]], cofrDeltaT[faceI]);

        corDeltaT[neighbour[faceI]] =
            max(corDeltaT[neighbour[faceI]], cofrDeltaT[faceI]);
    }

    volScalarField::GeometricBoundaryField& bcorDeltaT =
        corDeltaT.boundaryField();

    forAll(bcorDeltaT, patchi)
    {
        const fvsPatchScalarField& pcofrDeltaT =
            cofrDeltaT.boundaryField()[patchi];

        const fvPatch& p = pcofrDeltaT.patch();
        const labelUList& faceCells = p.patch().faceCells();

        forAll(pcofrDeltaT, patchFacei)
        {
            corDeltaT[faceCells[patchFacei]] = max
            (
                corDeltaT[faceCells[patchFacei]],
                pcofrDeltaT[patchFacei]
            );
        }
    }

    corDeltaT.correctBoundaryConditions();

    //corDeltaT = max(corDeltaT, max(corDeltaT)/100.0);

    return tcorDeltaT;
}


template<class Type>
tmp<surfaceScalarField> CoEulerDdtScheme<Type>::CofrDeltaT() const
{
    const dimensionedScalar& deltaT = mesh().time().deltaT();

    const surfaceScalarField& phi =
        static_cast<const objectRegistry&>(mesh())
        .lookupObject<surfaceScalarField>(phiName_);

    if (phi.dimensions() == dimensionSet(0, 3, -1, 0, 0))
    {
        surfaceScalarField Co
        (
            mesh().surfaceInterpolation::deltaCoeffs()
           *(mag(phi)/mesh().magSf())
           *deltaT
        );

        return max(Co/maxCo_, scalar(1))/deltaT;
    }
    else if (phi.dimensions() == dimensionSet(1, 0, -1, 0, 0))
    {
        const volScalarField& rho =
            static_cast<const objectRegistry&>(mesh())
           .lookupObject<volScalarField>(rhoName_).oldTime();

        surfaceScalarField Co
        (
            mesh().surfaceInterpolation::deltaCoeffs()
           *(mag(phi)/(fvc::interpolate(rho)*mesh().magSf()))
           *deltaT
        );

        return max(Co/maxCo_, scalar(1))/deltaT;
    }
    else
    {
        FatalErrorIn("CoEulerDdtScheme<Type>::CofrDeltaT() const")
            << "Incorrect dimensions of phi: " << phi.dimensions()
            << abort(FatalError);

        return tmp<surfaceScalarField>(NULL);
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
CoEulerDdtScheme<Type>::fvcDdt
(
    const dimensioned<Type>& dt
)
{
    const volScalarField rDeltaT(CorDeltaT());

    IOobject ddtIOobject
    (
        "ddt("+dt.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    if (mesh().moving())
    {
        tmp<GeometricField<Type, fvPatchField, volMesh> > tdtdt
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                dimensioned<Type>
                (
                    "0",
                    dt.dimensions()/dimTime,
                    pTraits<Type>::zero
                )
            )
        );

        tdtdt().internalField() =
            rDeltaT.internalField()*dt.value()*(1.0 - mesh().V0()/mesh().V());

        return tdtdt;
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                dimensioned<Type>
                (
                    "0",
                    dt.dimensions()/dimTime,
                    pTraits<Type>::zero
                ),
                calculatedFvPatchField<Type>::typeName
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
CoEulerDdtScheme<Type>::fvcDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const volScalarField rDeltaT(CorDeltaT());

    IOobject ddtIOobject
    (
        "ddt("+vf.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*vf.dimensions(),
                rDeltaT.internalField()*
                (
                    vf.internalField()
                  - vf.oldTime().internalField()*mesh().V0()/mesh().V()
                ),
                rDeltaT.boundaryField()*
                (
                    vf.boundaryField() - vf.oldTime().boundaryField()
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*(vf - vf.oldTime())
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
CoEulerDdtScheme<Type>::fvcDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const volScalarField rDeltaT(CorDeltaT());

    IOobject ddtIOobject
    (
        "ddt("+rho.name()+','+vf.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                rDeltaT.internalField()*rho.value()*
                (
                    vf.internalField()
                  - vf.oldTime().internalField()*mesh().V0()/mesh().V()
                ),
                rDeltaT.boundaryField()*rho.value()*
                (
                    vf.boundaryField() - vf.oldTime().boundaryField()
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*rho*(vf - vf.oldTime())
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
CoEulerDdtScheme<Type>::fvcDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const volScalarField rDeltaT(CorDeltaT());

    IOobject ddtIOobject
    (
        "ddt("+rho.name()+','+vf.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                rDeltaT.internalField()*
                (
                    rho.internalField()*vf.internalField()
                  - rho.oldTime().internalField()
                   *vf.oldTime().internalField()*mesh().V0()/mesh().V()
                ),
                rDeltaT.boundaryField()*
                (
                    rho.boundaryField()*vf.boundaryField()
                  - rho.oldTime().boundaryField()
                   *vf.oldTime().boundaryField()
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*(rho*vf - rho.oldTime()*vf.oldTime())
            )
        );
    }
}


template<class Type>
tmp<fvMatrix<Type> >
CoEulerDdtScheme<Type>::fvmDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            vf.dimensions()*dimVol/dimTime
        )
    );

    fvMatrix<Type>& fvm = tfvm();

    scalarField rDeltaT(CorDeltaT()().internalField());

    fvm.diag() = rDeltaT*mesh().V();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT*vf.oldTime().internalField()*mesh().V0();
    }
    else
    {
        fvm.source() = rDeltaT*vf.oldTime().internalField()*mesh().V();
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type> >
CoEulerDdtScheme<Type>::fvmDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm();

    scalarField rDeltaT(CorDeltaT()().internalField());

    fvm.diag() = rDeltaT*rho.value()*mesh().V();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT
            *rho.value()*vf.oldTime().internalField()*mesh().V0();
    }
    else
    {
        fvm.source() = rDeltaT
            *rho.value()*vf.oldTime().internalField()*mesh().V();
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type> >
CoEulerDdtScheme<Type>::fvmDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm();

    scalarField rDeltaT(CorDeltaT()().internalField());

    fvm.diag() = rDeltaT*rho.internalField()*mesh().V();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT
            *rho.oldTime().internalField()
            *vf.oldTime().internalField()*mesh().V0();
    }
    else
    {
        fvm.source() = rDeltaT
            *rho.oldTime().internalField()
            *vf.oldTime().internalField()*mesh().V();
    }

    return tfvm;
}


template<class Type>
tmp<typename CoEulerDdtScheme<Type>::fluxFieldType>
CoEulerDdtScheme<Type>::fvcDdtPhiCorr
(
    const volScalarField& rA,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
    IOobject ddtIOobject
    (
        "ddtPhiCorr(" + rA.name() + ',' + U.name() + ',' + phi.name() + ')',
        mesh().time().timeName(),
        mesh()
    );

    if (mesh().moving())
    {
        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                ddtIOobject,
                mesh(),
                dimensioned<typename flux<Type>::type>
                (
                    "0",
                    rA.dimensions()*phi.dimensions()/dimTime,
                    pTraits<typename flux<Type>::type>::zero
                )
            )
        );
    }
    else
    {
        const volScalarField rDeltaT(CorDeltaT());

        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                ddtIOobject,
                this->fvcDdtPhiCoeff(U.oldTime(), phi.oldTime())*
                (
                    fvc::interpolate(rDeltaT*rA)*phi.oldTime()
                  - (fvc::interpolate(rDeltaT*rA*U.oldTime()) & mesh().Sf())
                )
            )
        );
    }
}


template<class Type>
tmp<typename CoEulerDdtScheme<Type>::fluxFieldType>
CoEulerDdtScheme<Type>::fvcDdtPhiCorr
(
    const volScalarField& rA,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
    IOobject ddtIOobject
    (
        "ddtPhiCorr("
      + rA.name() + ',' + rho.name() + ',' + U.name() + ',' + phi.name() + ')',
        mesh().time().timeName(),
        mesh()
    );

    if (mesh().moving())
    {
        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                ddtIOobject,
                mesh(),
                dimensioned<typename flux<Type>::type>
                (
                    "0",
                    rA.dimensions()*rho.dimensions()*phi.dimensions()/dimTime,
                    pTraits<typename flux<Type>::type>::zero
                )
            )
        );
    }
    else
    {
        const volScalarField rDeltaT(CorDeltaT());

        if
        (
            U.dimensions() == dimVelocity
         && phi.dimensions() == dimVelocity*dimArea
        )
        {
            return tmp<fluxFieldType>
            (
                new fluxFieldType
                (
                    ddtIOobject,
                    this->fvcDdtPhiCoeff(U.oldTime(), phi.oldTime())
                   *(
                        fvc::interpolate(rDeltaT*rA*rho.oldTime())*phi.oldTime()
                      - (fvc::interpolate(rDeltaT*rA*rho.oldTime()*U.oldTime())
                      & mesh().Sf())
                    )
                )
            );
        }
        else if
        (
            U.dimensions() == dimVelocity
         && phi.dimensions() == dimDensity*dimVelocity*dimArea
        )
        {
            return tmp<fluxFieldType>
            (
                new fluxFieldType
                (
                    ddtIOobject,
                    this->fvcDdtPhiCoeff
                    (
                        U.oldTime(),
                        phi.oldTime()/fvc::interpolate(rho.oldTime())
                    )
                   *(
                        fvc::interpolate(rDeltaT*rA*rho.oldTime())
                       *phi.oldTime()/fvc::interpolate(rho.oldTime())
                      - (
                            fvc::interpolate
                            (
                                rDeltaT*rA*rho.oldTime()*U.oldTime()
                            ) & mesh().Sf()
                        )
                    )
                )
            );
        }
        else if
        (
            U.dimensions() == dimDensity*dimVelocity
         && phi.dimensions() == dimDensity*dimVelocity*dimArea
        )
        {
            return tmp<fluxFieldType>
            (
                new fluxFieldType
                (
                    ddtIOobject,
                    this->fvcDdtPhiCoeff
                    (rho.oldTime(), U.oldTime(), phi.oldTime())
                  * (
                        fvc::interpolate(rDeltaT*rA)*phi.oldTime()
                      - (
                            fvc::interpolate(rDeltaT*rA*U.oldTime())&mesh().Sf()
                        )
                    )
                )
            );
        }
        else
        {
            FatalErrorIn
            (
                "CoEulerDdtScheme<Type>::fvcDdtPhiCorr"
            )   << "dimensions of phi are not correct"
                << abort(FatalError);

            return fluxFieldType::null();
        }
    }
}


template<class Type>
tmp<surfaceScalarField> CoEulerDdtScheme<Type>::meshPhi
(
    const GeometricField<Type, fvPatchField, volMesh>&
)
{
    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                "meshPhi",
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensionedScalar("0", dimVolume/dimTime, 0.0)
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
