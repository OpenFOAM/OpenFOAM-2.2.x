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

#include "fixedAxis.H"
#include "addToRunTimeSelectionTable.H"
#include "sixDoFRigidBodyMotion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFRigidBodyMotionConstraints
{
    defineTypeNameAndDebug(fixedAxis, 0);

    addToRunTimeSelectionTable
    (
        sixDoFRigidBodyMotionConstraint,
        fixedAxis,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionConstraints::fixedAxis::fixedAxis
(
    const dictionary& sDoFRBMCDict
)
:
    sixDoFRigidBodyMotionConstraint(sDoFRBMCDict),
    fixedAxis_()
{
    read(sDoFRBMCDict);
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionConstraints::fixedAxis::~fixedAxis()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::sixDoFRigidBodyMotionConstraints::fixedAxis::constrain
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
    constraintMomentIncrement = vector::zero;

    vector predictedDir = motion.predictedOrientation
    (
        fixedAxis_,
        existingConstraintMoment,
        deltaT
    );

    scalar theta = acos(min(predictedDir & fixedAxis_, 1.0));

    vector rotationAxis = fixedAxis_ ^ predictedDir;

    scalar magRotationAxis = mag(rotationAxis);

    if (magRotationAxis > VSMALL)
    {
        rotationAxis /= magRotationAxis;

        const tensor& Q = motion.orientation();

        // Transform rotationAxis to body local system
        rotationAxis = Q.T() & rotationAxis;

        constraintMomentIncrement =
           -relaxationFactor_
           *(motion.momentOfInertia() & rotationAxis)
           *theta/sqr(deltaT);

        // Transform moment increment to global system
        constraintMomentIncrement = Q & constraintMomentIncrement;

        // Remove any moment that is around the fixedAxis_
        constraintMomentIncrement -=
            (constraintMomentIncrement & fixedAxis_)*fixedAxis_;
    }

    constraintPosition = motion.centreOfMass();

    constraintForceIncrement = vector::zero;

    bool converged(mag(theta) < tolerance_);

    if (sixDoFRigidBodyMotionConstraint::debug)
    {
        Info<< " angle " << theta
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


bool Foam::sixDoFRigidBodyMotionConstraints::fixedAxis::read
(
    const dictionary& sDoFRBMCDict
)
{
    sixDoFRigidBodyMotionConstraint::read(sDoFRBMCDict);

    sDoFRBMCCoeffs_.lookup("axis") >> fixedAxis_;

    scalar magFixedAxis(mag(fixedAxis_));

    if (magFixedAxis > VSMALL)
    {
        fixedAxis_ /= magFixedAxis;
    }
    else
    {
        FatalErrorIn
        (
            "Foam::sixDoFRigidBodyMotionConstraints::fixedAxis::read"
            "("
                "const dictionary& sDoFRBMCDict"
            ")"
        )
            << "axis has zero length"
            << abort(FatalError);
    }

    return true;
}


void Foam::sixDoFRigidBodyMotionConstraints::fixedAxis::write
(
    Ostream& os
) const
{
    os.writeKeyword("axis")
        << fixedAxis_ << token::END_STATEMENT << nl;
}

// ************************************************************************* //
