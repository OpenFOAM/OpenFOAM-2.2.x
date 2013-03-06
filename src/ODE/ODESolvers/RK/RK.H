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

Class
    Foam::RK

Description
    Runge-Kutta ODE solver

SourceFiles
    RKCK.C
    RKQS.C

\*---------------------------------------------------------------------------*/

#ifndef RK_H
#define RK_H

#include "ODESolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class RK Declaration
\*---------------------------------------------------------------------------*/

class RK
:
    public ODESolver
{
    // Private data
        static const scalar safety, pGrow, pShrink, errCon;

        static const scalar
            a2, a3, a4, a5, a6,
            b21, b31, b32, b41, b42, b43,
            b51, b52, b53, b54, b61, b62, b63, b64, b65,
            c1, c3, c4, c6,
            dc1, dc3, dc4, dc5, dc6;


        mutable scalarField yTemp_;
        mutable scalarField ak2_;
        mutable scalarField ak3_;
        mutable scalarField ak4_;
        mutable scalarField ak5_;
        mutable scalarField ak6_;

        mutable scalarField yErr_;
        mutable scalarField yTemp2_;

public:

    //- Runtime type information
    TypeName("RK");


    // Constructors

        //- Construct from ODE
        RK(const ODE& ode);


    // Member Functions

        void solve
        (
            const ODE& ode,
            const scalar x,
            const scalarField& y,
            const scalarField& dydx,
            const scalar h,
            scalarField& yout,
            scalarField& yerr
        ) const;


        void solve
        (
            const ODE& ode,
            scalar& x,
            scalarField& y,
            scalarField& dydx,
            const scalar eps,
            const scalarField& yScale,
            const scalar hTry,
            scalar& hDid,
            scalar& hNext
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
