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

#include "fixedPoint.H"
#include "addToRunTimeSelectionTable.H"
#include "sixDoFRigidBodyMotion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFRigidBodyMotionConstraints
{
    defineTypeNameAndDebug(fixedPoint, 0);

    addToRunTimeSelectionTable
    (
        sixDoFRigidBodyMotionConstraint,
        fixedPoint,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionConstraints::fixedPoint::fixedPoint
(
    const dictionary& sDoFRBMCDict
)
:
    sixDoFRigidBodyMotionConstraint(sDoFRBMCDict),
    fixedPoint_()
{
    read(sDoFRBMCDict);
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionConstraints::fixedPoint::~fixedPoint()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::sixDoFRigidBodyMotionConstraints::fixedPoint::constrain
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
    point predictedPosition = motion.predictedPosition
    (
        fixedPoint_,
        existingConstraintForce,
        existingConstraintMoment,
        deltaT
    );

    constraintPosition = motion.currentPosition(fixedPoint_);

    // Info<< "current position " << constraintPosition << nl
    //     << "next predictedPosition " << predictedPosition
    //     << endl;

    vector error = predictedPosition - fixedPoint_;

    // Info<< "error " << error << endl;

    // Correction force derived from Lagrange multiplier:
    //     G = -lambda*grad(sigma)
    // where
    //     sigma = mag(error) = 0
    // so
    //     grad(sigma) = error/mag(error)
    // Solving for lambda using the SHAKE methodology gives
    //     lambda = mass*mag(error)/sqr(deltaT)
    // This is only strictly applicable (i.e. will converge in one
    // iteration) to constraints at the centre of mass.  Everything
    // else will need to iterate, and may need under-relaxed to be
    // stable.

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


bool Foam::sixDoFRigidBodyMotionConstraints::fixedPoint::read
(
    const dictionary& sDoFRBMCDict
)
{
    sixDoFRigidBodyMotionConstraint::read(sDoFRBMCDict);

    sDoFRBMCCoeffs_.lookup("fixedPoint") >> fixedPoint_;

    return true;
}


void Foam::sixDoFRigidBodyMotionConstraints::fixedPoint::write
(
    Ostream& os
) const
{
    os.writeKeyword("fixedPoint")
        << fixedPoint_ << token::END_STATEMENT << nl;
}

// ************************************************************************* //
