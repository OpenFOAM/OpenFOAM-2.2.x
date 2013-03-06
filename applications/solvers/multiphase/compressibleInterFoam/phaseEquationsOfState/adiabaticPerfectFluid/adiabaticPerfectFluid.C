/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "adiabaticPerfectFluid.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseEquationsOfState
{
    defineTypeNameAndDebug(adiabaticPerfectFluid, 0);

    addToRunTimeSelectionTable
    (
        phaseEquationOfState,
        adiabaticPerfectFluid,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseEquationsOfState::adiabaticPerfectFluid::adiabaticPerfectFluid
(
    const dictionary& dict
)
:
    phaseEquationOfState(dict),
    p0_("p0", dimPressure, dict.lookup("p0")),
    rho0_("rho0", dimDensity, dict.lookup("rho0")),
    gamma_("gamma", dimless, dict.lookup("gamma")),
    B_("B", dimPressure, dict.lookup("B"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseEquationsOfState::adiabaticPerfectFluid::~adiabaticPerfectFluid()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::phaseEquationsOfState::adiabaticPerfectFluid::rho
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return tmp<Foam::volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "rho",
                p.time().timeName(),
                p.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            rho0_*pow((p + B_)/(p0_ + B_), 1.0/gamma_)
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::phaseEquationsOfState::adiabaticPerfectFluid::psi
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return tmp<Foam::volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "psi",
                p.time().timeName(),
                p.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            (rho0_/(gamma_*(p0_ + B_)))
           *pow((p + B_)/(p0_ + B_), 1.0/gamma_ - 1.0)
        )
    );
}


// ************************************************************************* //
