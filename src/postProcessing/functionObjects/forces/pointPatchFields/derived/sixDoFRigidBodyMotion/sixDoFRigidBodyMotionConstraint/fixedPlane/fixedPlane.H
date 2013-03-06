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
    Foam::sixDoFRigidBodyMotionConstraints::fixedPlane

Description
    sixDoFRigidBodyMotionConstraint.  Reference point may only move
    along a plane.

SourceFiles
    fixedPlane.C

\*---------------------------------------------------------------------------*/

#ifndef fixedPlane_H
#define fixedPlane_H

#include "sixDoFRigidBodyMotionConstraint.H"
#include "plane.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace sixDoFRigidBodyMotionConstraints
{

/*---------------------------------------------------------------------------*\
                          Class fixedPlane Declaration
\*---------------------------------------------------------------------------*/

class fixedPlane
:
    public sixDoFRigidBodyMotionConstraint
{
    // Private data

        //- Plane which the transformed reference point is constrained
        //  to move along
        plane fixedPlane_;


public:

    //- Runtime type information
    TypeName("fixedPlane");


    // Constructors

        //- Construct from components
        fixedPlane
        (
            const dictionary& sDoFRBMCDict
        );

        //- Construct and return a clone
        virtual autoPtr<sixDoFRigidBodyMotionConstraint> clone() const
        {
            return autoPtr<sixDoFRigidBodyMotionConstraint>
            (
                new fixedPlane(*this)
            );
        }


    //- Destructor
    virtual ~fixedPlane();


    // Member Functions

        //- Calculate the constraint position, force and moment.
        //  Global reference frame vectors.  Returns boolean stating
        //  whether the constraint been converged to tolerance.
        virtual bool constrain
        (
            const sixDoFRigidBodyMotion& motion,
            const vector& existingConstraintForce,
            const vector& existingConstraintMoment,
            scalar deltaT,
            vector& constraintPosition,
            vector& constraintForceIncrement,
            vector& constraintMomentIncrement
        ) const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& sDoFRBMCCoeff);

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
