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
    Foam::surfaceShearForce

Description
    Film surface shear force

SourceFiles
    surfaceShearForce.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceShearForce_H
#define surfaceShearForce_H

#include "force.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                      Class surfaceShearForce Declaration
\*---------------------------------------------------------------------------*/

class surfaceShearForce
:
    public force
{
private:

    // Private Data

        //- Surface roughness coefficient
        scalar Cf_;



    // Private member functions

        //- Disallow default bitwise copy construct
        surfaceShearForce(const surfaceShearForce&);

        //- Disallow default bitwise assignment
        void operator=(const surfaceShearForce&);


public:

    //- Runtime type information
    TypeName("surfaceShear");


    // Constructors

        //- Construct from surface film model
        surfaceShearForce
        (
            const surfaceFilmModel& owner,
            const dictionary& dict
        );


    //- Destructor
    virtual ~surfaceShearForce();


    // Member Functions

        // Evolution

            //- Correct
            virtual tmp<fvVectorMatrix> correct(volVectorField& U);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
