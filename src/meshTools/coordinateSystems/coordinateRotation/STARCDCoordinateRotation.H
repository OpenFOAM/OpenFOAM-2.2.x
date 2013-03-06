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
    Foam::STARCDCoordinateRotation

Description
    A coordinateRotation defined by the STAR-CD convention.

    The 3 rotations are defined in the STAR-CD convention
    (around Z, around X' and around Y'').
    The order of the parameter arguments matches this rotation order.

    - the rotation angles are in degrees, unless otherwise explictly specified:

    \verbatim
        coordinateRotation
        {
            type        STARCDRotation;
            degrees     false;
            rotation    (0 0 3.141592654);
        }
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef STARCDCoordinateRotation_H
#define STARCDCoordinateRotation_H

#include "coordinateRotation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class STARCDCoordinateRotation Declaration
\*---------------------------------------------------------------------------*/

class STARCDCoordinateRotation
:
    public coordinateRotation
{
    // Private Member Functions

        //- Calculate transformation tensor
        void calcTransform
        (
            const scalar rotZ,
            const scalar rotX,
            const scalar rotY,
            const bool inDegrees=true
        );


public:

    //- Runtime type information
    TypeName("STARCDRotation");


    // Constructors

        //- Construct null
        STARCDCoordinateRotation();

        //- Construct from rotation vector
        STARCDCoordinateRotation
        (
            const vector& rotZrotXrotY,
            const bool inDegrees=true
        );

        //- Construct from components of rotation vector
        STARCDCoordinateRotation
        (
            const scalar rotZ,
            const scalar rotX,
            const scalar rotY,
            const bool inDegrees=true
        );

        //- Construct from dictionary
        STARCDCoordinateRotation(const dictionary&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
