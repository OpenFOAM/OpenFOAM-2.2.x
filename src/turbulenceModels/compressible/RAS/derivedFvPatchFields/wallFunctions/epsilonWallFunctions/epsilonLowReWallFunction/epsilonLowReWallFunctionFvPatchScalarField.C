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

#include "epsilonLowReWallFunctionFvPatchScalarField.H"
#include "compressible/turbulenceModel/turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

scalar epsilonLowReWallFunctionFvPatchScalarField::yPlusLam
(
    const scalar kappa,
    const scalar E
)
{
    scalar ypl = 11.0;

    for (int i=0; i<10; i++)
    {
        ypl = log(max(E*ypl, 1))/kappa;
    }

    return ypl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

epsilonLowReWallFunctionFvPatchScalarField::
epsilonLowReWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    epsilonWallFunctionFvPatchScalarField(p, iF),
    yPlusLam_(yPlusLam(kappa_, E_))
{}


epsilonLowReWallFunctionFvPatchScalarField::
epsilonLowReWallFunctionFvPatchScalarField
(
    const epsilonLowReWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    epsilonWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    yPlusLam_(ptf.yPlusLam_)
{}


epsilonLowReWallFunctionFvPatchScalarField::
epsilonLowReWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    epsilonWallFunctionFvPatchScalarField(p, iF, dict),
    yPlusLam_(yPlusLam(kappa_, E_))
{}


epsilonLowReWallFunctionFvPatchScalarField::
epsilonLowReWallFunctionFvPatchScalarField
(
    const epsilonLowReWallFunctionFvPatchScalarField& ewfpsf
)
:
    epsilonWallFunctionFvPatchScalarField(ewfpsf),
    yPlusLam_(ewfpsf.yPlusLam_)
{}


epsilonLowReWallFunctionFvPatchScalarField::
epsilonLowReWallFunctionFvPatchScalarField
(
    const epsilonLowReWallFunctionFvPatchScalarField& ewfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    epsilonWallFunctionFvPatchScalarField(ewfpsf, iF),
    yPlusLam_(ewfpsf.yPlusLam_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void epsilonLowReWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchI = patch().index();

    const turbulenceModel& turbulence =
        db().lookupObject<turbulenceModel>("turbulenceModel");
    const scalarField& y = turbulence.y()[patchI];

    volScalarField& G =
        const_cast<volScalarField&>
        (
            db().lookupObject<volScalarField>
            (
                turbulence.GName()
            )
        );

    DimensionedField<scalar, volMesh>& epsilon =
        const_cast<DimensionedField<scalar, volMesh>&>
        (
            dimensionedInternalField()
        );

    const tmp<volScalarField> tk = turbulence.k();
    const volScalarField& k = tk();

    const tmp<volScalarField> tmu = turbulence.mu();
    const scalarField& muw = tmu().boundaryField()[patchI];

    const tmp<volScalarField> tmut = turbulence.mut();
    const volScalarField& mut = tmut();
    const scalarField& mutw = mut.boundaryField()[patchI];

    const scalarField& rhow = turbulence.rho().boundaryField()[patchI];

    const fvPatchVectorField& Uw = turbulence.U().boundaryField()[patchI];
    const scalarField magGradUw(mag(Uw.snGrad()));

    const scalar Cmu25 = pow025(Cmu_);
    const scalar Cmu75 = pow(Cmu_, 0.75);

    // Set epsilon and G
    forAll(mutw, faceI)
    {
        label faceCellI = patch().faceCells()[faceI];

        scalar yPlus = Cmu25*sqrt(k[faceCellI])*y[faceI]/muw[faceI]/rhow[faceI];

        if (yPlus > yPlusLam_)
        {
            epsilon[faceCellI] = Cmu75*pow(k[faceCellI], 1.5)/(kappa_*y[faceI]);
        }
        else
        {
            epsilon[faceCellI] =
                2.0*k[faceCellI]*muw[faceI]/rhow[faceI]/sqr(y[faceI]);
        }

        G[faceCellI] =
            (mutw[faceI] + muw[faceI])
           *magGradUw[faceI]
           *Cmu25*sqrt(k[faceCellI])
           /(kappa_*y[faceI]);
    }

    fixedInternalValueFvPatchField<scalar>::updateCoeffs();

    // TODO: perform averaging for cells sharing more than one boundary face
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    epsilonLowReWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
