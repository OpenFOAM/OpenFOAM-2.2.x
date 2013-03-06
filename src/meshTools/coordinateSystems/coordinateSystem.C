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

#include "IOstream.H"
#include "coordinateSystem.H"
#include "coordinateSystems.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coordinateSystem, 0);
    defineRunTimeSelectionTable(coordinateSystem, dictionary);
    defineRunTimeSelectionTable(coordinateSystem, origRotation);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordinateSystem::coordinateSystem()
:
    name_(type()),
    note_(),
    origin_(point::zero),
    R_(),
    Rtr_(sphericalTensor::I)
{}


Foam::coordinateSystem::coordinateSystem
(
    const word& name,
    const coordinateSystem& cs
)
:
    name_(name),
    note_(),
    origin_(cs.origin_),
    R_(cs.R_),
    Rtr_(R_.T())
{}


Foam::coordinateSystem::coordinateSystem
(
    const word& name,
    const point& origin,
    const coordinateRotation& cr
)
:
    name_(name),
    note_(),
    origin_(origin),
    R_(cr),
    Rtr_(R_.T())
{}


Foam::coordinateSystem::coordinateSystem
(
    const word& name,
    const point& origin,
    const vector& axis,
    const vector& dirn
)
:
    name_(name),
    note_(),
    origin_(origin),
    R_(axis, dirn),
    Rtr_(R_.T())
{}


Foam::coordinateSystem::coordinateSystem
(
    const word& name,
    const dictionary& dict
)
:
    name_(name),
    note_(),
    origin_(point::zero),
    R_(),
    Rtr_(sphericalTensor::I)
{
    operator=(dict);
}


Foam::coordinateSystem::coordinateSystem(const dictionary& dict)
:
    name_(type()),
    note_(),
    origin_(point::zero),
    R_(),
    Rtr_(sphericalTensor::I)
{
    operator=(dict);
}


Foam::coordinateSystem::coordinateSystem
(
    const dictionary& dict,
    const objectRegistry& obr
)
:
    name_(type()),
    note_(),
    origin_(point::zero),
    R_(),
    Rtr_(sphericalTensor::I)
{
    const entry* entryPtr = dict.lookupEntryPtr(typeName_(), false, false);

    // non-dictionary entry is a lookup into global coordinateSystems
    if (entryPtr && !entryPtr->isDict())
    {
        keyType key(entryPtr->stream());

        const coordinateSystems& lst = coordinateSystems::New(obr);
        const label index = lst.findIndex(key);

        if (debug)
        {
            Info<< "coordinateSystem::coordinateSystem"
                "(const dictionary&, const objectRegistry&):"
                << nl << "using global coordinate system: "
                << key << "=" << index << endl;
        }

        if (index < 0)
        {
            FatalErrorIn
            (
                "coordinateSystem::coordinateSystem"
                "(const dictionary&, const objectRegistry&)"
            )   << "could not find coordinate system: " << key << nl
                << "available coordinate systems: " << lst.toc() << nl << nl
                << exit(FatalError);
        }

        // copy coordinateSystem, but assign the name as the typeName
        // to avoid strange things in writeDict()
        operator=(lst[index]);
        name_ = typeName_();
    }
    else
    {
        operator=(dict);
    }
}


Foam::coordinateSystem::coordinateSystem(Istream& is)
:
    name_(is),
    note_(),
    origin_(point::zero),
    R_(),
    Rtr_(sphericalTensor::I)
{
    dictionary dict(is);
    operator=(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coordinateSystem::~coordinateSystem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dictionary Foam::coordinateSystem::dict(bool ignoreType) const
{
    dictionary dict;

    dict.add("name", name_);

    // only write type for derived types
    if (!ignoreType && type() != typeName_())
    {
        dict.add("type", type());
    }

    // The note entry is optional
    if (note_.size())
    {
        dict.add("note", note_);
    }

    dict.add("origin", origin_);
    dict.add("e1", e1());
    dict.add("e3", e3());

    return dict;
}


Foam::vector Foam::coordinateSystem::localToGlobal
(
    const vector& local,
    bool translate
) const
{
    if (translate)
    {
        return (R_ & local) + origin_;
    }
    else
    {
        return (R_ & local);
    }
}


Foam::tmp<Foam::vectorField> Foam::coordinateSystem::localToGlobal
(
    const vectorField& local,
    bool translate
) const
{
    if (translate)
    {
        return (R_ & local) + origin_;
    }
    else
    {
        return (R_ & local);
    }
}


Foam::vector Foam::coordinateSystem::globalToLocal
(
    const vector& global,
    bool translate
) const
{
    if (translate)
    {
        return (Rtr_ & (global - origin_));
    }
    else
    {
        return (Rtr_ & global);
    }
}


Foam::tmp<Foam::vectorField> Foam::coordinateSystem::globalToLocal
(
    const vectorField& global,
    bool translate
) const
{
    if (translate)
    {
        return (Rtr_ & (global - origin_));
    }
    else
    {
        return (Rtr_ & global);
    }
}


void Foam::coordinateSystem::clear()
{
    note_.clear();
    origin_ = point::zero;
    R_.clear();
    Rtr_ = sphericalTensor::I;
}


void Foam::coordinateSystem::write(Ostream& os) const
{
    os  << type()
        << " origin: " << origin() << " e1: " << e1() << " e3: " << e3();
}


void Foam::coordinateSystem::writeDict(Ostream& os, bool subDict) const
{
    if (subDict)
    {
        os  << indent << name_ << nl
            << indent << token::BEGIN_BLOCK << incrIndent << nl;
    }

    // only write type for derived types
    if (type() != typeName_())
    {
        os.writeKeyword("type") << type() << token::END_STATEMENT << nl;
    }

    // The note entry is optional
    if (note_.size())
    {
        os.writeKeyword("note") << note_ << token::END_STATEMENT << nl;
    }

    os.writeKeyword("origin") << origin_ << token::END_STATEMENT << nl;
    os.writeKeyword("e1") << e1() << token::END_STATEMENT << nl;
    os.writeKeyword("e3") << e3() << token::END_STATEMENT << nl;

    if (subDict)
    {
        os  << decrIndent << indent << token::END_BLOCK << endl;
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::coordinateSystem::operator=(const dictionary& rhs)
{
    if (debug)
    {
        Pout<< "coordinateSystem::operator=(const dictionary&) : "
            << "assign from " << rhs << endl;
    }

    // allow as embedded sub-dictionary "coordinateSystem"
    const dictionary& dict =
    (
        rhs.found(typeName_())
      ? rhs.subDict(typeName_())
      : rhs
    );

    // unspecified origin is (0 0 0)
    origin_ = point::zero;
    dict.readIfPresent("origin", origin_);

    // The note entry is optional
    note_.clear();
    rhs.readIfPresent("note", note_);

    // specify via coordinateRotation sub-dictionary
    if (dict.found("coordinateRotation"))
    {
        R_  = coordinateRotation::New(dict.subDict("coordinateRotation"))();
    }
    else
    {
        // let coordinateRotation constructor extract the axes specification
        R_ = coordinateRotation(dict);
    }

    Rtr_ = R_.T();
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

bool Foam::operator!=(const coordinateSystem& a, const coordinateSystem& b)
{
    return (a.origin() != b.origin() || a.R() != b.R() || a.type() != b.type());
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const coordinateSystem& cs)
{
    cs.write(os);
    os.check("Ostream& operator<<(Ostream&, const coordinateSystem&");
    return os;
}


// ************************************************************************* //
