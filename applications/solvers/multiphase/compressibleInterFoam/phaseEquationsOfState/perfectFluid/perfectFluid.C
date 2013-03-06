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

#include "perfectFluid.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseEquationsOfState
{
    defineTypeNameAndDebug(perfectFluid, 0);

    addToRunTimeSelectionTable
    (
        phaseEquationOfState,
        perfectFluid,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseEquationsOfState::perfectFluid::perfectFluid
(
    const dictionary& dict
)
:
    phaseEquationOfState(dict),
    rho0_("rho0", dimDensity, dict.lookup("rho0")),
    R_("R", dimensionSet(0, 2, -2, -1, 0), dict.lookup("R"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseEquationsOfState::perfectFluid::~perfectFluid()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::phaseEquationsOfState::perfectFluid::rho
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
            rho0_ + psi(p, T)*p
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::phaseEquationsOfState::perfectFluid::psi
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
            1.0/(R_*T)
        )
    );
}


// ************************************************************************* //
