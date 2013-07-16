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

#include "buoyantPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::buoyantPressureFvPatchScalarField::
buoyantPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    rhoName_("rho")
{}


Foam::buoyantPressureFvPatchScalarField::
buoyantPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho"))
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


Foam::buoyantPressureFvPatchScalarField::
buoyantPressureFvPatchScalarField
(
    const buoyantPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    rhoName_(ptf.rhoName_)
{}


Foam::buoyantPressureFvPatchScalarField::
buoyantPressureFvPatchScalarField
(
    const buoyantPressureFvPatchScalarField& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    rhoName_(ptf.rhoName_)
{}


Foam::buoyantPressureFvPatchScalarField::
buoyantPressureFvPatchScalarField
(
    const buoyantPressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    rhoName_(ptf.rhoName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::buoyantPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const uniformDimensionedVectorField& g =
        db().lookupObject<uniformDimensionedVectorField>("g");

    const fvPatchField<scalar>& rho =
        patch().lookupPatchField<volScalarField, scalar>(rhoName_);

    // If the variable name is "p_rgh", "ph_rgh" or "pd"
    // assume it is p? - rho*g.h and set the gradient appropriately.
    // Otherwise assume the variable is the static pressure.
    if
    (
        dimensionedInternalField().name() == "p_rgh"
     || dimensionedInternalField().name() == "ph_rgh"
     || dimensionedInternalField().name() == "pd"
    )
    {
        gradient() = -rho.snGrad()*(g.value() & patch().Cf());
    }
    else
    {
        gradient() = rho*(g.value() & patch().nf());
    }

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::buoyantPressureFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        buoyantPressureFvPatchScalarField
    );
}


// ************************************************************************* //
