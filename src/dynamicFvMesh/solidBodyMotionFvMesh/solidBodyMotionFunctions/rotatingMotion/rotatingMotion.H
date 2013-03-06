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
    Foam::solidBodyMotionFunctions::rotatingMotion

Description
    SolidBodyMotionFvMesh 6DoF motion function. Constant
    velocity rotation around CoG.

SourceFiles
    rotatingMotion.C

\*---------------------------------------------------------------------------*/

#ifndef rotatingMotion_H
#define rotatingMotion_H

#include "solidBodyMotionFunction.H"
#include "primitiveFields.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                          Class rotatingMotion Declaration
\*---------------------------------------------------------------------------*/

class rotatingMotion
:
    public solidBodyMotionFunction
{
    // Private data

        //- Centre of gravity
        point CofG_;

        //- Rotational velocity (deg/s)
        vector radialVelocity_;


    // Private Member Functions

        //- Disallow copy construct
        rotatingMotion(const rotatingMotion&);

        //- Disallow default bitwise assignment
        void operator=(const rotatingMotion&);


public:

    //- Runtime type information
    TypeName("rotatingMotion");


    // Constructors

        //- Construct from components
        rotatingMotion
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );


    //- Destructor
    virtual ~rotatingMotion();


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
