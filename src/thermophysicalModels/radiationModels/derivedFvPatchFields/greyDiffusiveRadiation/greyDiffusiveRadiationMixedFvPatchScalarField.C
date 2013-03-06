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

#include "greyDiffusiveRadiationMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

#include "fvDOM.H"
#include "constants.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::greyDiffusiveRadiationMixedFvPatchScalarField::
greyDiffusiveRadiationMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    radiationCoupledBase(p, "undefined", scalarField::null()),
    TName_("T")
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 1.0;
}


Foam::radiation::greyDiffusiveRadiationMixedFvPatchScalarField::
greyDiffusiveRadiationMixedFvPatchScalarField
(
    const greyDiffusiveRadiationMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    radiationCoupledBase
    (
        p,
        ptf.emissivityMethod(),
        ptf.emissivity_
    ),
    TName_(ptf.TName_)
{}


Foam::radiation::greyDiffusiveRadiationMixedFvPatchScalarField::
greyDiffusiveRadiationMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    radiationCoupledBase(p, dict),
    TName_(dict.lookupOrDefault<word>("T", "T"))
{
    if (dict.found("refValue"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // No value given. Restart as fixedValue b.c.

        const scalarField& Tp =
            patch().lookupPatchField<volScalarField, scalar>(TName_);

        //NOTE: Assumes emissivity = 1 as the solidThermo might
        // not be constructed yet
        refValue() =
            4.0*physicoChemical::sigma.value()*pow4(Tp)/pi;
        refGrad() = 0.0;
        valueFraction() = 1.0;

        fvPatchScalarField::operator=(refValue());
    }
}


Foam::radiation::greyDiffusiveRadiationMixedFvPatchScalarField::
greyDiffusiveRadiationMixedFvPatchScalarField
(
    const greyDiffusiveRadiationMixedFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    radiationCoupledBase
    (
        ptf.patch(),
        ptf.emissivityMethod(),
        ptf.emissivity_
    ),
    TName_(ptf.TName_)
{}


Foam::radiation::greyDiffusiveRadiationMixedFvPatchScalarField::
greyDiffusiveRadiationMixedFvPatchScalarField
(
    const greyDiffusiveRadiationMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    radiationCoupledBase
    (
        ptf.patch(),
        ptf.emissivityMethod(),
        ptf.emissivity_
    ),
    TName_(ptf.TName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiation::greyDiffusiveRadiationMixedFvPatchScalarField::
updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    const scalarField& Tp =
        patch().lookupPatchField<volScalarField, scalar>(TName_);

    const radiationModel& radiation =
        db().lookupObject<radiationModel>("radiationProperties");

    const fvDOM& dom(refCast<const fvDOM>(radiation));

    label rayId = -1;
    label lambdaId = -1;
    dom.setRayIdLambdaId(dimensionedInternalField().name(), rayId, lambdaId);

    const label patchI = patch().index();

    if (dom.nLambda() != 1)
    {
        FatalErrorIn
        (
            "Foam::radiation::"
            "greyDiffusiveRadiationMixedFvPatchScalarField::updateCoeffs"
        )   << " a grey boundary condition is used with a non-grey "
            << "absorption model" << nl << exit(FatalError);
    }

    scalarField& Iw = *this;
    const vectorField n(patch().Sf()/patch().magSf());

    radiativeIntensityRay& ray =
        const_cast<radiativeIntensityRay&>(dom.IRay(rayId));

    ray.Qr().boundaryField()[patchI] += Iw*(n & ray.dAve());

    scalarField temissivity = emissivity();

    forAll(Iw, faceI)
    {
        scalar Ir = 0.0;

        for (label rayI=0; rayI < dom.nRay(); rayI++)
        {
            const vector& d = dom.IRay(rayI).d();

            const scalarField& IFace =
                dom.IRay(rayI).ILambda(lambdaId).boundaryField()[patchI];

            if ((-n[faceI] & d) < 0.0)
            {
                // q into the wall
                const vector& dAve = dom.IRay(rayI).dAve();
                Ir += IFace[faceI]*mag(n[faceI] & dAve);
            }
        }

        const vector& d = dom.IRay(rayId).d();

        if ((-n[faceI] & d) > 0.0)
        {
            // direction out of the wall
            refGrad()[faceI] = 0.0;
            valueFraction()[faceI] = 1.0;
            refValue()[faceI] =
                (
                    Ir*(scalar(1.0) - temissivity[faceI])
                  + temissivity[faceI]*physicoChemical::sigma.value()
                  * pow4(Tp[faceI])
                )/pi;

             // Emmited heat flux from this ray direction
            ray.Qem().boundaryField()[patchI][faceI] =
                refValue()[faceI]*(n[faceI] & ray.dAve());
        }
        else
        {
            // direction into the wall
            valueFraction()[faceI] = 0.0;
            refGrad()[faceI] = 0.0;
            refValue()[faceI] = 0.0; //not used

             // Incident heat flux on this ray direction
            ray.Qin().boundaryField()[patchI][faceI] =
                Iw[faceI]*(n[faceI] & ray.dAve());
        }
    }

    // Restore tag
    UPstream::msgType() = oldTag;

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::radiation::greyDiffusiveRadiationMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    radiationCoupledBase::write(os);
    writeEntryIfDifferent<word>(os, "T", "T", TName_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{
    makePatchTypeField
    (
        fvPatchScalarField,
        greyDiffusiveRadiationMixedFvPatchScalarField
    );
}
}


// ************************************************************************* //
