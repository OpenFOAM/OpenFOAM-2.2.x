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
    Foam::phaseEquationOfState

Description
    A2stract base-class for dispersed-phase particle diameter models.

SourceFiles
    phaseEquationOfState.C
    newDiameterModel.C

\*---------------------------------------------------------------------------*/

#ifndef phaseEquationOfState_H
#define phaseEquationOfState_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dictionary.H"
#include "volFieldsFwd.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class phaseEquationOfState Declaration
\*---------------------------------------------------------------------------*/

class phaseEquationOfState
{
protected:

    // Protected data

        const dictionary& dict_;


public:

    //- Runtime type information
    TypeName("phaseEquationOfState");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            phaseEquationOfState,
            dictionary,
            (
                const dictionary& dict
            ),
            (dict)
        );


    // Constructors

        phaseEquationOfState
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~phaseEquationOfState();


    // Selectors

        static autoPtr<phaseEquationOfState> New
        (
            const dictionary& dict
        );


    // Member Functions

        //- Return the phase density
        virtual tmp<volScalarField> rho
        (
            const volScalarField& p,
            const volScalarField& T
        ) const = 0;

        //- Return the phase compressibility
        virtual tmp<volScalarField> psi
        (
            const volScalarField& p,
            const volScalarField& T
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
