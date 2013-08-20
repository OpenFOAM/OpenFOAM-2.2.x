/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "CloudFunctionObject.H"

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class CloudType>
void Foam::CloudFunctionObject<CloudType>::write()
{
    notImplemented("void Foam::CloudFunctionObject<CloudType>::write()");
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CloudFunctionObject<CloudType>::CloudFunctionObject(CloudType& owner)
:
    SubModelBase<CloudType>(owner)
{}


template<class CloudType>
Foam::CloudFunctionObject<CloudType>::CloudFunctionObject
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    SubModelBase<CloudType>(owner, dict, typeName, type, "")
{}


template<class CloudType>
Foam::CloudFunctionObject<CloudType>::CloudFunctionObject
(
    const CloudFunctionObject<CloudType>& ppm
)
:
    SubModelBase<CloudType>(ppm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CloudFunctionObject<CloudType>::~CloudFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::CloudFunctionObject<CloudType>::preEvolve()
{
    // do nothing
}


template<class CloudType>
void Foam::CloudFunctionObject<CloudType>::postEvolve()
{
    if (this->owner().time().outputTime())
    {
        this->write();
    }
}


template<class CloudType>
void Foam::CloudFunctionObject<CloudType>::postMove
(
    const typename CloudType::parcelType&,
    const label,
    const scalar,
    const point&,
    bool&
)
{
    // do nothing
}


template<class CloudType>
void Foam::CloudFunctionObject<CloudType>::postPatch
(
    const typename CloudType::parcelType&,
    const polyPatch&,
    const scalar,
    const tetIndices&,
    bool&
)
{
    // do nothing
}


template<class CloudType>
void Foam::CloudFunctionObject<CloudType>::postFace
(
    const typename CloudType::parcelType&,
    const label,
    bool&
)
{
    // do nothing
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CloudFunctionObjectNew.C"

// ************************************************************************* //
