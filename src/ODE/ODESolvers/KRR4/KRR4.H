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
    Foam::KRR4

Description
    Foam::KRR4

SourceFiles
    KRR4CK.C
    KRR4QS.C

\*---------------------------------------------------------------------------*/

#ifndef KRR4_H
#define KRR4_H

#include "ODESolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class KRR4 Declaration
\*---------------------------------------------------------------------------*/

class KRR4
:
    public ODESolver
{
    // Private data

        mutable scalarField yTemp_;
        mutable scalarField dydxTemp_;
        mutable scalarField g1_;
        mutable scalarField g2_;
        mutable scalarField g3_;
        mutable scalarField g4_;
        mutable scalarField yErr_;
        mutable scalarField dfdx_;
        mutable scalarSquareMatrix dfdy_;
        mutable scalarSquareMatrix a_;
        mutable labelList pivotIndices_;

        static const int maxtry = 40;

        static const scalar safety, grow, pgrow, shrink, pshrink, errcon;

        static const scalar
            gamma,
            a21, a31, a32,
            c21, c31, c32, c41, c42, c43,
            b1, b2, b3, b4,
            e1, e2, e3, e4,
            c1X, c2X, c3X, c4X,
            a2X, a3X;


public:

    //- Runtime type information
    TypeName("KRR4");


    // Constructors

        //- Construct from ODE
        KRR4(const ODE& ode);


    // Member Functions

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
