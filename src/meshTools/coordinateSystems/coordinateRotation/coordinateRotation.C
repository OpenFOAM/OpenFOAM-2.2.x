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

#include "coordinateRotation.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coordinateRotation, 0);
    defineRunTimeSelectionTable(coordinateRotation, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::coordinateRotation::calcTransform
(
    const vector& axis1,
    const vector& axis2,
    const axisOrder& order
)
{
    vector a = axis1 / mag(axis1);
    vector b = axis2;

    // Absorb minor nonorthogonality into axis2
    b = b - (b & a)*a;

    if (mag(b) < SMALL)
    {
        FatalErrorIn("coordinateRotation::calcTransform()")
            << "axis1, axis2 appear co-linear: "
            << axis1 << ", " << axis2 << endl
            << abort(FatalError);
    }

    b = b / mag(b);
    vector c = a ^ b;

    // the global -> local transformation
    tensor Rtr;
    switch (order)
    {
        case e1e2:
            Rtr = tensor(a, b, c);
            break;

        case e2e3:
            Rtr = tensor(c, a, b);
            break;

        case e3e1:
            Rtr = tensor(b, c, a);
            break;

        default:
            FatalErrorIn("coordinateRotation::calcTransform()")
                << "programmer error" << endl
                << abort(FatalError);
            // To satisfy compiler warnings
            Rtr = tensor::zero;
            break;
    }

    // the local -> global transformation
    tensor::operator=( Rtr.T() );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordinateRotation::coordinateRotation()
:
    tensor(sphericalTensor::I)
{}


Foam::coordinateRotation::coordinateRotation
(
    const vector& axis,
    const vector& dir
)
:
    tensor(sphericalTensor::I)
{
    calcTransform(axis, dir, e3e1);
}


Foam::coordinateRotation::coordinateRotation
(
    const dictionary& dict
)
:
    tensor(sphericalTensor::I)
{
    operator=(dict);
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::coordinateRotation> Foam::coordinateRotation::New
(
    const dictionary& dict
)
{
    if (debug)
    {
        Pout<< "coordinateRotation::New(const dictionary&) : "
            << "constructing coordinateRotation"
            << endl;
    }

    // default type is self (alias: "axes")
    word rotType(typeName_());
    dict.readIfPresent("type", rotType);

    // can (must) construct base class directly
    if (rotType == typeName_() || rotType == "axes")
    {
        return autoPtr<coordinateRotation>(new coordinateRotation(dict));
    }


    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(rotType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "coordinateRotation::New(const dictionary&)",
            dict
        )   << "Unknown coordinateRotation type "
            << rotType << nl << nl
            << "Valid coordinateRotation types are :" <<  nl
            << "[default: axes " << typeName_() << "]"
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<coordinateRotation>(cstrIter()(dict));
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::coordinateRotation::clear()
{
    this->tensor::operator=(sphericalTensor::I);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::coordinateRotation::operator=(const dictionary& rhs)
{
    if (debug)
    {
        Pout<< "coordinateRotation::operator=(const dictionary&) : "
            << "assign from " << rhs << endl;
    }

    // allow as embedded sub-dictionary "coordinateRotation"
    const dictionary& dict =
    (
        rhs.found(typeName_())
      ? rhs.subDict(typeName_())
      : rhs
    );

    vector axis1, axis2;
    axisOrder order(e3e1);

    if (dict.readIfPresent("e1", axis1) && dict.readIfPresent("e2", axis2))
    {
        order = e1e2;
    }
    else if (dict.readIfPresent("e2", axis1) && dict.readIfPresent("e3", axis2))
    {
        order = e2e3;
    }
    else if (dict.readIfPresent("e3", axis1) && dict.readIfPresent("e1", axis2))
    {
        order = e3e1;
    }
    else if (dict.found("axis") || dict.found("direction"))
    {
        // let it bomb if only one of axis/direction is defined
        order = e3e1;
        dict.lookup("axis") >> axis1;
        dict.lookup("direction") >> axis2;
    }
    else
    {
        // unspecified axes revert to the global system
        tensor::operator=(sphericalTensor::I);
        return;
    }

    calcTransform(axis1, axis2, order);
}


// ************************************************************************* //
