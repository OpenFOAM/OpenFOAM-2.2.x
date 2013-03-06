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
    Foam::EulerCoordinateRotation

Description
    A coordinateRotation defined in the z-x-y Euler convention.

    The 3 rotations are defined in the Euler convention
    (around Z, around X' and around Z').
    For reference and illustration, see
    http://mathworld.wolfram.com/EulerAngles.html
    Note, however, that it is the reverse transformation
    (local->global) that is defined here.

    - the rotation angles are in degrees, unless otherwise explictly specified:

    \verbatim
        coordinateRotation
        {
            type        EulerRotation
            degrees     false;
            rotation    (0 0 3.141592654);
        }
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef EulerCoordinateRotation_H
#define EulerCoordinateRotation_H

#include "coordinateRotation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class EulerCoordinateRotation Declaration
\*---------------------------------------------------------------------------*/

class EulerCoordinateRotation
:
    public coordinateRotation
{
    // Private Member Functions

        //- Calculate transformation tensor
        void calcTransform
        (
            const scalar phiAngle,
            const scalar thetaAngle,
            const scalar psiAngle,
            const bool inDegrees=true
        );


public:

    //- Runtime type information
    TypeName("EulerRotation");


    // Constructors

        //- Construct null
        EulerCoordinateRotation();

        //- Construct from rotation vector
        EulerCoordinateRotation
        (
            const vector& phiThetaPsi,
            const bool inDegrees=true
        );

        //- Construct from components of rotation vector
        EulerCoordinateRotation
        (
            const scalar phiAngle,
            const scalar thetaAngle,
            const scalar psiAngle,
            const bool inDegrees=true
        );

        //- Construct from dictionary
        EulerCoordinateRotation(const dictionary&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
