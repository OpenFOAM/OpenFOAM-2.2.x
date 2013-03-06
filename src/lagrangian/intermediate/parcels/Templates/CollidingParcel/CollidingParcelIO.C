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

#include "CollidingParcel.H"
#include "IOstreams.H"
#include "IOField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::CollidingParcel<ParcelType>::propertyList_ =
    Foam::CollidingParcel<ParcelType>::propertyList();


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::CollidingParcel<ParcelType>::CollidingParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    ParcelType(mesh, is, readFields),
    collisionRecords_()
{
    if (readFields)
    {
        is >> collisionRecords_;
    }

    // Check state of Istream
    is.check
    (
        "CollidingParcel<ParcelType>::Collisions"
        "(const polyMesh&, Istream&, bool)"
    );
}


template<class ParcelType>
template<class CloudType>
void Foam::CollidingParcel<ParcelType>::readFields(CloudType& c)
{
    if (!c.size())
    {
        return;
    }

    ParcelType::readFields(c);

    labelFieldCompactIOField collisionRecordsPairAccessed
    (
        c.fieldIOobject("collisionRecordsPairAccessed", IOobject::MUST_READ)
    );
    c.checkFieldFieldIOobject(c, collisionRecordsPairAccessed);

    labelFieldCompactIOField collisionRecordsPairOrigProcOfOther
    (
        c.fieldIOobject
        (
            "collisionRecordsPairOrigProcOfOther",
            IOobject::MUST_READ
        )
    );
    c.checkFieldFieldIOobject(c, collisionRecordsPairOrigProcOfOther);

    labelFieldCompactIOField collisionRecordsPairOrigIdOfOther
    (
        c.fieldIOobject
        (
            "collisionRecordsPairOrigIdOfOther",
            IOobject::MUST_READ
        )
    );
    c.checkFieldFieldIOobject(c, collisionRecordsPairOrigProcOfOther);

    pairDataFieldCompactIOField collisionRecordsPairData
    (
        c.fieldIOobject("collisionRecordsPairData", IOobject::MUST_READ)
    );
    c.checkFieldFieldIOobject(c, collisionRecordsPairData);

    labelFieldCompactIOField collisionRecordsWallAccessed
    (
        c.fieldIOobject("collisionRecordsWallAccessed", IOobject::MUST_READ)
    );
    c.checkFieldFieldIOobject(c, collisionRecordsWallAccessed);

    vectorFieldCompactIOField collisionRecordsWallPRel
    (
        c.fieldIOobject("collisionRecordsWallPRel", IOobject::MUST_READ)
    );
    c.checkFieldFieldIOobject(c, collisionRecordsWallPRel);

    wallDataFieldCompactIOField collisionRecordsWallData
    (
        c.fieldIOobject("collisionRecordsWallData", IOobject::MUST_READ)
    );
    c.checkFieldFieldIOobject(c, collisionRecordsWallData);

    label i = 0;

    forAllIter(typename CloudType, c, iter)
    {
        CollidingParcel<ParcelType>& p = iter();

        p.collisionRecords_ = collisionRecordList
        (
            collisionRecordsPairAccessed[i],
            collisionRecordsPairOrigProcOfOther[i],
            collisionRecordsPairOrigIdOfOther[i],
            collisionRecordsPairData[i],
            collisionRecordsWallAccessed[i],
            collisionRecordsWallPRel[i],
            collisionRecordsWallData[i]
        );

        i++;
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::CollidingParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);

    label np =  c.size();

    labelFieldCompactIOField collisionRecordsPairAccessed
    (
        c.fieldIOobject("collisionRecordsPairAccessed", IOobject::NO_READ),
        np
    );
    labelFieldCompactIOField collisionRecordsPairOrigProcOfOther
    (
        c.fieldIOobject
        (
            "collisionRecordsPairOrigProcOfOther",
            IOobject::NO_READ
        ),
        np
    );
    labelFieldCompactIOField collisionRecordsPairOrigIdOfOther
    (
        c.fieldIOobject("collisionRecordsPairOrigIdOfOther", IOobject::NO_READ),
        np
    );
    pairDataFieldCompactIOField collisionRecordsPairData
    (
        c.fieldIOobject("collisionRecordsPairData", IOobject::NO_READ),
        np
    );
    labelFieldCompactIOField collisionRecordsWallAccessed
    (
        c.fieldIOobject("collisionRecordsWallAccessed", IOobject::NO_READ),
        np
    );
    vectorFieldCompactIOField collisionRecordsWallPRel
    (
        c.fieldIOobject("collisionRecordsWallPRel", IOobject::NO_READ),
        np
    );
    wallDataFieldCompactIOField collisionRecordsWallData
    (
        c.fieldIOobject("collisionRecordsWallData", IOobject::NO_READ),
        np
    );

    label i = 0;

    forAllConstIter(typename CloudType, c, iter)
    {
        const CollidingParcel<ParcelType>& p = iter();

        collisionRecordsPairAccessed[i] = p.collisionRecords().pairAccessed();
        collisionRecordsPairOrigProcOfOther[i] =
            p.collisionRecords().pairOrigProcOfOther();
        collisionRecordsPairOrigIdOfOther[i] =
            p.collisionRecords().pairOrigIdOfOther();
        collisionRecordsPairData[i] = p.collisionRecords().pairData();
        collisionRecordsWallAccessed[i] = p.collisionRecords().wallAccessed();
        collisionRecordsWallPRel[i] = p.collisionRecords().wallPRel();
        collisionRecordsWallData[i] = p.collisionRecords().wallData();

        i++;
    }

    collisionRecordsPairAccessed.write();
    collisionRecordsPairOrigProcOfOther.write();
    collisionRecordsPairOrigIdOfOther.write();
    collisionRecordsPairData.write();
    collisionRecordsWallAccessed.write();
    collisionRecordsWallPRel.write();
    collisionRecordsWallData.write();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const CollidingParcel<ParcelType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << p.collisionRecords();
    }
    else
    {
        os  << static_cast<const ParcelType&>(p)
            << p.collisionRecords();
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const CollidingParcel<ParcelType>&)"
    );

    return os;
}


// ************************************************************************* //
