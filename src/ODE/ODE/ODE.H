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
    Foam::ODE

Description
    Abstract base class for the ODE solvers.

\*---------------------------------------------------------------------------*/

#ifndef ODE_H
#define ODE_H

#include "scalarField.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ODE Declaration
\*---------------------------------------------------------------------------*/

class ODE
{

public:

    // Constructors

        //- Construct null
        ODE()
        {}


    //- Destructor
    virtual ~ODE()
    {}


    // Member Functions

        virtual label nEqns() const = 0;

        virtual void derivatives
        (
            const scalar x,
            const scalarField& y,
            scalarField& dydx
        ) const = 0;


        virtual void jacobian
        (
            const scalar x,
            const scalarField& y,
            scalarField& dfdx,
            scalarSquareMatrix& dfdy
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
