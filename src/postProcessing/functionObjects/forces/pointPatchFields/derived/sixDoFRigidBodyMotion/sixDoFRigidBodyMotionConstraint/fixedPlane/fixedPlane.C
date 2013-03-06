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

#include "fixedPlane.H"
#include "addToRunTimeSelectionTable.H"
#include "sixDoFRigidBodyMotion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFRigidBodyMotionConstraints
{
    defineTypeNameAndDebug(fixedPlane, 0);

    addToRunTimeSelectionTable
    (
        sixDoFRigidBodyMotionConstraint,
        fixedPlane,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionConstraints::fixedPlane::fixedPlane
(
    const dictionary& sDoFRBMCDict
)
:
    sixDoFRigidBodyMotionConstraint(sDoFRBMCDict),
    fixedPlane_(vector::one)
{
    read(sDoFRBMCDict);
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionConstraints::fixedPlane::~fixedPlane()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::sixDoFRigidBodyMotionConstraints::fixedPlane::constrain
(
    const sixDoFRigidBodyMotion& motion,
    const vector& existingConstraintForce,
    const vector& existingConstraintMoment,
    scalar deltaT,
    vector& constraintPosition,
    vector& constraintForceIncrement,
    vector& constraintMomentIncrement
) const
{
    const point& refPt = fixedPlane_.refPoint();

    const vector& n = fixedPlane_.normal();

    point predictedPosition = motion.predictedPosition
    (
        refPt,
        existingConstraintForce,
        existingConstraintMoment,
        deltaT
    );

    constraintPosition = motion.currentPosition(refPt);

    // Info<< "current position " << constraintPosition << nl
    //     << "next predictedPosition " << predictedPosition
    //     << endl;

    vector error = ((predictedPosition - refPt) & n)*n;

    // Info<< "error " << error << endl;

    constraintForceIncrement =
        -relaxationFactor_*error*motion.mass()/sqr(deltaT);

    constraintMomentIncrement = vector::zero;

    bool converged(mag(error) < tolerance_);

    if (sixDoFRigidBodyMotionConstraint::debug)
    {
        Info<< " error " << error
            << " force " << constraintForceIncrement
            << " moment " << constraintMomentIncrement;

        if (converged)
        {
            Info<< " converged";
        }
        else
        {
            Info<< " not converged";
        }

        Info<< endl;
    }

    return converged;
}


bool Foam::sixDoFRigidBodyMotionConstraints::fixedPlane::read
(
    const dictionary& sDoFRBMCDict
)
{
    sixDoFRigidBodyMotionConstraint::read(sDoFRBMCDict);

    point refPt = sDoFRBMCCoeffs_.lookup("refPoint");

    vector normal = sDoFRBMCCoeffs_.lookup("normal");

    fixedPlane_ = plane(refPt, normal);

    return true;
}


void Foam::sixDoFRigidBodyMotionConstraints::fixedPlane::write
(
    Ostream& os
) const
{
    os.writeKeyword("refPoint")
        << fixedPlane_.refPoint() << token::END_STATEMENT << nl;

    os.writeKeyword("normal")
        << fixedPlane_.normal() << token::END_STATEMENT << nl;
}

// ************************************************************************* //
