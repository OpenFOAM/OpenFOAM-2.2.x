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

#include "threePhaseMixture.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

//- Calculate and return the laminar viscosity
void Foam::threePhaseMixture::calcNu()
{
    nuModel1_->correct();
    nuModel2_->correct();
    nuModel3_->correct();

    // Average kinematic viscosity calculated from dynamic viscosity
    nu_ = mu()/(alpha1_*rho1_ + alpha2_*rho2_ + alpha3_*rho3_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::threePhaseMixture::threePhaseMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    transportModel(U, phi),

    phase1Name_("phase1"),
    phase2Name_("phase2"),
    phase3Name_("phase3"),

    nuModel1_
    (
        viscosityModel::New
        (
            "nu1",
            subDict(phase1Name_),
            U,
            phi
        )
    ),
    nuModel2_
    (
        viscosityModel::New
        (
            "nu2",
            subDict(phase2Name_),
            U,
            phi
        )
    ),
    nuModel3_
    (
        viscosityModel::New
        (
            "nu3",
            subDict(phase3Name_),
            U,
            phi
        )
    ),

    rho1_(nuModel1_->viscosityProperties().lookup("rho")),
    rho2_(nuModel2_->viscosityProperties().lookup("rho")),
    rho3_(nuModel3_->viscosityProperties().lookup("rho")),

    U_(U),
    phi_(phi),

    alpha1_(U_.db().lookupObject<const volScalarField> ("alpha1")),
    alpha2_(U_.db().lookupObject<const volScalarField> ("alpha2")),
    alpha3_(U_.db().lookupObject<const volScalarField> ("alpha3")),

    nu_
    (
        IOobject
        (
            "nu",
            U_.time().timeName(),
            U_.db()
        ),
        U_.mesh(),
        dimensionedScalar("nu", dimensionSet(0, 2, -1, 0, 0), 0),
        calculatedFvPatchScalarField::typeName
    )
{
    calcNu();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::threePhaseMixture::mu() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            "mu",
            alpha1_*rho1_*nuModel1_->nu()
          + alpha2_*rho2_*nuModel2_->nu()
          + alpha3_*rho3_*nuModel3_->nu()
        )
    );
}


Foam::tmp<Foam::surfaceScalarField> Foam::threePhaseMixture::muf() const
{
    surfaceScalarField alpha1f(fvc::interpolate(alpha1_));
    surfaceScalarField alpha2f(fvc::interpolate(alpha2_));
    surfaceScalarField alpha3f(fvc::interpolate(alpha3_));

    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            "mu",
            alpha1f*rho1_*fvc::interpolate(nuModel1_->nu())
          + alpha2f*rho2_*fvc::interpolate(nuModel2_->nu())
          + alpha3f*rho3_*fvc::interpolate(nuModel3_->nu())
        )
    );
}


Foam::tmp<Foam::surfaceScalarField> Foam::threePhaseMixture::nuf() const
{
    surfaceScalarField alpha1f(fvc::interpolate(alpha1_));
    surfaceScalarField alpha2f(fvc::interpolate(alpha2_));
    surfaceScalarField alpha3f(fvc::interpolate(alpha3_));

    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            "nu",
            (
                alpha1f*rho1_*fvc::interpolate(nuModel1_->nu())
              + alpha2f*rho2_*fvc::interpolate(nuModel2_->nu())
              + alpha3f*rho3_*fvc::interpolate(nuModel3_->nu())
            )/(alpha1f*rho1_ + alpha2f*rho2_ + alpha3f*rho3_)
        )
    );
}


bool Foam::threePhaseMixture::read()
{
    if (transportModel::read())
    {
        if
        (
            nuModel1_().read(*this)
         && nuModel2_().read(*this)
         && nuModel3_().read(*this)
        )
        {
            nuModel1_->viscosityProperties().lookup("rho") >> rho1_;
            nuModel2_->viscosityProperties().lookup("rho") >> rho2_;
            nuModel3_->viscosityProperties().lookup("rho") >> rho3_;

            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
