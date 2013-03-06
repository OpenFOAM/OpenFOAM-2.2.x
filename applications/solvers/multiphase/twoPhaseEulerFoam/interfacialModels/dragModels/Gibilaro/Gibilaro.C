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

#include "Gibilaro.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Gibilaro, 0);

    addToRunTimeSelectionTable
    (
        dragModel,
        Gibilaro,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Gibilaro::Gibilaro
(
    const dictionary& interfaceDict,
    const volScalarField& alpha1,
    const phaseModel& phase1,
    const phaseModel& phase2
)
:
    dragModel(interfaceDict, alpha1, phase1, phase2)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Gibilaro::~Gibilaro()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::Gibilaro::K
(
    const volScalarField& Ur
) const
{
    volScalarField alpha2(max(scalar(1) - alpha1_, scalar(1.0e-6)));
    volScalarField bp(pow(alpha2, -2.8));
    volScalarField Re(max(alpha2*Ur*phase1_.d()/phase2_.nu(), scalar(1.0e-3)));

    return (17.3/Re + scalar(0.336))*phase2_.rho()*Ur*bp/phase1_.d();
}


// ************************************************************************* //
