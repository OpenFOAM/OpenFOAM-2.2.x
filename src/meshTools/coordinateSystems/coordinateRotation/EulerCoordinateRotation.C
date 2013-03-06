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

\*---------------------------------------------------------------------------*/

#include "EulerCoordinateRotation.H"

#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(EulerCoordinateRotation, 0);
    addToRunTimeSelectionTable
    (
        coordinateRotation,
        EulerCoordinateRotation,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::EulerCoordinateRotation::calcTransform
(
    const scalar phiAngle,
    const scalar thetaAngle,
    const scalar psiAngle,
    const bool inDegrees
)
{
    scalar phi   = phiAngle;
    scalar theta = thetaAngle;
    scalar psi   = psiAngle;

    if (inDegrees)
    {
        phi   *= constant::mathematical::pi/180.0;
        theta *= constant::mathematical::pi/180.0;
        psi   *= constant::mathematical::pi/180.0;
    }

    tensor::operator=
    (
        tensor
        (
            cos(phi)*cos(psi) - sin(phi)*sin(psi)*cos(theta),
            -sin(phi)*cos(psi)*cos(theta) - cos(phi)*sin(psi),
            sin(phi)*sin(theta),

            cos(phi)*sin(psi)*cos(theta) + sin(phi)*cos(psi),
            cos(phi)*cos(psi)*cos(theta) - sin(phi)*sin(psi),
            -cos(phi)*sin(theta),

            sin(psi)*sin(theta),
            cos(psi)*sin(theta),
            cos(theta)
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::EulerCoordinateRotation::EulerCoordinateRotation()
:
    coordinateRotation()
{}


Foam::EulerCoordinateRotation::EulerCoordinateRotation
(
    const vector& phiThetaPsi,
    const bool inDegrees
)
:
    coordinateRotation()
{
    calcTransform
    (
        phiThetaPsi.component(vector::X),
        phiThetaPsi.component(vector::Y),
        phiThetaPsi.component(vector::Z),
        inDegrees
    );
}


Foam::EulerCoordinateRotation::EulerCoordinateRotation
(
    const scalar phiAngle,
    const scalar thetaAngle,
    const scalar psiAngle,
    const bool inDegrees
)
:
    coordinateRotation()
{
    calcTransform(phiAngle, thetaAngle, psiAngle, inDegrees);
}


Foam::EulerCoordinateRotation::EulerCoordinateRotation
(
    const dictionary& dict
)
:
    coordinateRotation()
{
    vector rotation(dict.lookup("rotation"));

    calcTransform
    (
        rotation.component(vector::X),
        rotation.component(vector::Y),
        rotation.component(vector::Z),
        dict.lookupOrDefault("degrees", true)
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
