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

#include "linear.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseEquationsOfState
{
    defineTypeNameAndDebug(linear, 0);

    addToRunTimeSelectionTable
    (
        phaseEquationOfState,
        linear,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseEquationsOfState::linear::linear
(
    const dictionary& dict
)
:
    phaseEquationOfState(dict),
    rho0_("rho0", dimDensity, dict.lookup("rho0")),
    psi_("psi", dimDensity/dimPressure, dict.lookup("psi"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseEquationsOfState::linear::~linear()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::phaseEquationsOfState::linear::rho
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
            rho0_ + psi_*p
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::phaseEquationsOfState::linear::psi
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
            p.mesh(),
            psi_
        )
    );
}


// ************************************************************************* //
