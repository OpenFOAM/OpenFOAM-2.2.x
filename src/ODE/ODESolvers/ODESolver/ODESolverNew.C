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

#include "ODESolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::ODESolver> Foam::ODESolver::New
(
    const Foam::word& ODESolverTypeName,
    const Foam::ODE& ode
)
{
    Info<< "Selecting ODE solver " << ODESolverTypeName << endl;

    ODEConstructorTable::iterator cstrIter =
        ODEConstructorTablePtr_->find(ODESolverTypeName);

    if (cstrIter == ODEConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "ODESolver::New(const word& ODESolverTypeName, const ODE& ode)"
        )   << "Unknown ODESolver type "
            << ODESolverTypeName << nl << nl
            << "Valid ODESolvers are : " << endl
            << ODEConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<ODESolver>(cstrIter()(ode));
}


// ************************************************************************* //
