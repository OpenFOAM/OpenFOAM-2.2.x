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
    Foam::phaseEquationsOfState::adiabaticPerfectFluid

Description
    AdiabaticPerfectFluid phase density model.

SourceFiles
    adiabaticPerfectFluid.C

\*---------------------------------------------------------------------------*/

#ifndef adiabaticPerfectFluid_H
#define adiabaticPerfectFluid_H

#include "phaseEquationOfState.H"
#include "dimensionedTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseEquationsOfState
{

/*---------------------------------------------------------------------------*\
                           Class adiabaticPerfectFluid Declaration
\*---------------------------------------------------------------------------*/

class adiabaticPerfectFluid
:
    public phaseEquationOfState
{
    // Private data

        //- Reference pressure
        dimensionedScalar p0_;

        //- Reference density
        dimensionedScalar rho0_;

        //- The isentropic exponent
        dimensionedScalar gamma_;

        //- Pressure offset for a stiffened gas
        dimensionedScalar B_;


public:

    //- Runtime type information
    TypeName("adiabaticPerfectFluid");


    // Constructors

        //- Construct from components
        adiabaticPerfectFluid
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~adiabaticPerfectFluid();


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
