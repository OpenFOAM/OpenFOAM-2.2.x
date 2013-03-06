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

Class
    Foam::phaseEquationsOfState::constant

Description
    Constant phase density model.

SourceFiles
    constant.C

\*---------------------------------------------------------------------------*/

#ifndef constant_H
#define constant_H

#include "phaseEquationOfState.H"
#include "dimensionedTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseEquationsOfState
{

/*---------------------------------------------------------------------------*\
                           Class constant Declaration
\*---------------------------------------------------------------------------*/

class constant
:
    public phaseEquationOfState
{
    // Private data

        //- The constant density of the phase
        dimensionedScalar rho_;


public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from components
        constant
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~constant();


    // Member Functions

        tmp<volScalarField> rho
        (
            const volScalarField& p,
            const volScalarField& T
        ) const;

        tmp<volScalarField> psi
        (
            const volScalarField& p,
            const volScalarField& T
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace phaseEquationsOfState
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
