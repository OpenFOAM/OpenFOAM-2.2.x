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

#include "fixedLine.H"
#include "addToRunTimeSelectionTable.H"
#include "sixDoFRigidBodyMotion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFRigidBodyMotionConstraints
{
    defineTypeNameAndDebug(fixedLine, 0);

    addToRunTimeSelectionTable
    (
        sixDoFRigidBodyMotionConstraint,
        fixedLine,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionConstraints::fixedLine::fixedLine
(
    const dictionary& sDoFRBMCDict
)
:
    sixDoFRigidBodyMotionConstraint(sDoFRBMCDict),
    refPt_(),
    dir_()
{
    read(sDoFRBMCDict);
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionConstraints::fixedLine::~fixedLine()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::sixDoFRigidBodyMotionConstraints::fixedLine::constrain
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
        refPt_,
        existingConstraintForce,
        existingConstraintMoment,
        deltaT
    );

    constraintPosition = motion.currentPosition(refPt_);

    // Info<< "current position " << constraintPosition << nl
    //     << "next predictedPosition " << predictedPosition
    //     << endl;

    // Vector from reference point to predicted point
    vector rC = predictedPosition - refPt_;

    vector error = rC - ((rC) & dir_)*dir_;

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


bool Foam::sixDoFRigidBodyMotionConstraints::fixedLine::read
(
    const dictionary& sDoFRBMCDict
)
{
    sixDoFRigidBodyMotionConstraint::read(sDoFRBMCDict);

    sDoFRBMCCoeffs_.lookup("refPoint") >> refPt_;

    sDoFRBMCCoeffs_.lookup("direction") >> dir_;

    scalar magDir(mag(dir_));

    if (magDir > VSMALL)
    {
        dir_ /= magDir;
    }
    else
    {
        FatalErrorIn
        (
            "Foam::sixDoFRigidBodyMotionConstraints::fixedLine::read"
            "("
                "const dictionary& sDoFRBMCDict"
            ")"
        )
            << "line direction has zero length"
            << abort(FatalError);
    }

    return true;
}


void Foam::sixDoFRigidBodyMotionConstraints::fixedLine::write
(
    Ostream& os
) const
{
    os.writeKeyword("refPoint")
        << refPt_ << token::END_STATEMENT << nl;

    os.writeKeyword("direction")
        << dir_ << token::END_STATEMENT << nl;
}

// ************************************************************************* //
