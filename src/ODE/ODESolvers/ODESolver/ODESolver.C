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

#include "ODESolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ODESolver, 0);
    defineRunTimeSelectionTable(ODESolver, ODE);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ODESolver::ODESolver(const ODE& ode)
:
    n_(ode.nEqns()),
    yScale_(n_),
    dydx_(n_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ODESolver::solve
(
    const ODE& ode,
    const scalar xStart,
    const scalar xEnd,
    scalarField& y,
    const scalar eps,
    scalar& hEst
) const
{
    const label MAXSTP = 10000;

    scalar x = xStart;
    scalar h = hEst;
    scalar hNext = 0;
    scalar hPrev = 0;

    for (label nStep=0; nStep<MAXSTP; nStep++)
    {
        ode.derivatives(x, y, dydx_);

        for (label i=0; i<n_; i++)
        {
            yScale_[i] = mag(y[i]) + mag(dydx_[i]*h) + SMALL;
        }

        if ((x + h - xEnd)*(x + h - xStart) > 0.0)
        {
            h = xEnd - x;
            hPrev = hNext;
        }

        hNext = 0;
        scalar hDid;
        solve(ode, x, y, dydx_, eps, yScale_, h, hDid, hNext);

        if ((x - xEnd)*(xEnd - xStart) >= 0.0)
        {
            if (hPrev != 0)
            {
                hEst = hPrev;
            }
            else
            {
                hEst = hNext;
            }

            return;
        }

        h = hNext;
    }

    FatalErrorIn
    (
        "ODESolver::solve"
        "(const ODE& ode, const scalar xStart, const scalar xEnd,"
        "scalarField& yStart, const scalar eps, scalar& hEst) const"
    )   << "Too many integration steps"
        << exit(FatalError);
}

// ************************************************************************* //
