/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "sixDoFRigidBodyMotionConstraint.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(sixDoFRigidBodyMotionConstraint, 0);

defineRunTimeSelectionTable(sixDoFRigidBodyMotionConstraint, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionConstraint::sixDoFRigidBodyMotionConstraint
(
    const dictionary& sDoFRBMCDict
)
:
    sDoFRBMCCoeffs_
    (
        sDoFRBMCDict.subDict
        (
            word(sDoFRBMCDict.lookup("sixDoFRigidBodyMotionConstraint"))
          + "Coeffs"
        )
    ),
    tolerance_(readScalar(sDoFRBMCDict.lookup("tolerance"))),
    relaxationFactor_
    (
        sDoFRBMCDict.lookupOrDefault<scalar>("relaxationFactor", 1)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionConstraint::~sixDoFRigidBodyMotionConstraint()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::sixDoFRigidBodyMotionConstraint::read
(
    const dictionary& sDoFRBMCDict
)
{
    tolerance_ = (readScalar(sDoFRBMCDict.lookup("tolerance")));

    relaxationFactor_ = sDoFRBMCDict.lookupOrDefault<scalar>
    (
        "relaxationFactor",
        1
    );

    sDoFRBMCCoeffs_ = sDoFRBMCDict.subDict(type() + "Coeffs");

    return true;
}


void Foam::sixDoFRigidBodyMotionConstraint::write(Ostream& os) const
{
    os.writeKeyword("tolerance")
        << tolerance_ << token::END_STATEMENT << nl;

    os.writeKeyword("relaxationFactor")
        << relaxationFactor_ << token::END_STATEMENT << nl;
}

// ************************************************************************* //
