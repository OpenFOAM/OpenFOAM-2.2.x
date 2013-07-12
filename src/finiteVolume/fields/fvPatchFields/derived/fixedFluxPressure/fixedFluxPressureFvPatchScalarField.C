/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "fixedFluxPressureFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedFluxPressureFvPatchScalarField::fixedFluxPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    phiHbyAName_("phiHbyA"),
    phiName_("phi"),
    rhoName_("rho"),
    DpName_("Dp"),
    adjoint_(false)
{}


Foam::fixedFluxPressureFvPatchScalarField::fixedFluxPressureFvPatchScalarField
(
    const fixedFluxPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    phiHbyAName_(ptf.phiHbyAName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    DpName_(ptf.DpName_),
    adjoint_(ptf.adjoint_)
{}


Foam::fixedFluxPressureFvPatchScalarField::fixedFluxPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    phiHbyAName_(dict.lookupOrDefault<word>("phiHbyA", "phiHbyA")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    DpName_(dict.lookupOrDefault<word>("Dp", "Dp")),
    adjoint_(dict.lookupOrDefault<Switch>("adjoint", false))
{
    if (dict.found("value") && dict.found("gradient"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
        gradient() = scalarField("gradient", dict, p.size());
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;
    }
}


Foam::fixedFluxPressureFvPatchScalarField::fixedFluxPressureFvPatchScalarField
(
    const fixedFluxPressureFvPatchScalarField& wbppsf
)
:
    fixedGradientFvPatchScalarField(wbppsf),
    phiHbyAName_(wbppsf.phiHbyAName_),
    phiName_(wbppsf.phiName_),
    rhoName_(wbppsf.rhoName_),
    DpName_(wbppsf.DpName_),
    adjoint_(wbppsf.adjoint_)
{}


Foam::fixedFluxPressureFvPatchScalarField::fixedFluxPressureFvPatchScalarField
(
    const fixedFluxPressureFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(wbppsf, iF),
    phiHbyAName_(wbppsf.phiHbyAName_),
    phiName_(wbppsf.phiName_),
    rhoName_(wbppsf.rhoName_),
    DpName_(wbppsf.DpName_),
    adjoint_(wbppsf.adjoint_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedFluxPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const surfaceScalarField& phiHbyA =
        db().lookupObject<surfaceScalarField>(phiHbyAName_);

    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);

    fvsPatchField<scalar> phiHbyAp =
        patch().patchField<surfaceScalarField, scalar>(phiHbyA);

    fvsPatchField<scalar> phip =
        patch().patchField<surfaceScalarField, scalar>(phi);

    /*
    if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchField<scalar>& rhop =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);

        phip /= rhop;
    }

    if (phiHbyA.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchField<scalar>& rhop =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);

        phiHbyAp /= rhop;
    }
    */

    const scalarField *DppPtr = NULL;

    if (db().foundObject<volScalarField>(DpName_))
    {
        DppPtr =
            &patch().lookupPatchField<volScalarField, scalar>(DpName_);
    }
    else if (db().foundObject<surfaceScalarField>(DpName_))
    {
        const surfaceScalarField& Dp =
            db().lookupObject<surfaceScalarField>(DpName_);

        DppPtr =
            &patch().patchField<surfaceScalarField, scalar>(Dp);
    }

    if (adjoint_)
    {
        gradient() = (phip - phiHbyAp)/patch().magSf()/(*DppPtr);
    }
    else
    {
        gradient() = (phiHbyAp - phip)/patch().magSf()/(*DppPtr);
    }

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::fixedFluxPressureFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "phiHbyA", "phiHbyA", phiHbyAName_);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    writeEntryIfDifferent<word>(os, "Dp", "Dp", DpName_);
    if (adjoint_)
    {
        os.writeKeyword("adjoint") << adjoint_ << token::END_STATEMENT << nl;
    }
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedFluxPressureFvPatchScalarField
    );
}

// ************************************************************************* //
