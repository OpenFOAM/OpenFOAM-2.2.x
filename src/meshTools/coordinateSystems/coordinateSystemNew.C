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

#include "coordinateSystem.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::coordinateSystem> Foam::coordinateSystem::New
(
    const word& name,
    const dictionary& dict
)
{
    if (debug)
    {
        Pout<< "coordinateSystem::New(const word&, const dictionary&) : "
            << "constructing coordinateSystem"
            << endl;
    }

    // construct base class directly, also allow 'cartesian' as an alias
    word coordType(typeName_());
    if
    (
        !dict.readIfPresent("type", coordType)
     || coordType == typeName_()
     || coordType == "cartesian"
    )
    {
        return autoPtr<coordinateSystem>(new coordinateSystem(name, dict));
    }

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(coordType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "coordinateSystem::New(const word&, const dictionary&)",
            dict
        )   << "Unknown coordinateSystem type "
            << coordType << nl << nl
            << "Valid coordinateSystem types are :" << nl
            << "[default: " << typeName_() << "]"
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<coordinateSystem>(cstrIter()(name, dict));
}


Foam::autoPtr<Foam::coordinateSystem> Foam::coordinateSystem::New
(
    const word& coordType,
    const word& name,
    const point& origin,
    const coordinateRotation& cr
)
{
    if (debug)
    {
        Pout<< "coordinateSystem::New(const word&, const word&, "
            << "const point&, const coordinateRotation&) : "
               "constructing coordinateSystem"
            << endl;
    }

    origRotationConstructorTable::iterator cstrIter =
        origRotationConstructorTablePtr_->find(coordType);

    if (cstrIter == origRotationConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "coordinateSystem::New(const word&, const word&, "
            "const point&, const coordinateRotation&) : "
            "constructing coordinateSystem"
        )   << "Unknown coordinateSystem type "
            << coordType << nl << nl
            << "Valid coordinateSystem types are :" << nl
            << origRotationConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<coordinateSystem>(cstrIter()(name, origin, cr));
}


Foam::autoPtr<Foam::coordinateSystem> Foam::coordinateSystem::New
(
    Istream& is
)
{
    const word name(is);
    const dictionary dict(is);

    return New(name, dict);
}


// ************************************************************************* //
