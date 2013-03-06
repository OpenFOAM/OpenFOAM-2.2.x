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

#include "StochasticCollisionModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::StochasticCollisionModel<CloudType>::StochasticCollisionModel
(
    CloudType& owner
)
:
    SubModelBase<CloudType>(owner)
{}


template<class CloudType>
Foam::StochasticCollisionModel<CloudType>::StochasticCollisionModel
(
    const StochasticCollisionModel<CloudType>& cm
)
:
    SubModelBase<CloudType>(cm)
{}


template<class CloudType>
Foam::StochasticCollisionModel<CloudType>::StochasticCollisionModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    SubModelBase<CloudType>(owner, dict, typeName, type)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::StochasticCollisionModel<CloudType>::~StochasticCollisionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::StochasticCollisionModel<CloudType>::update
(
    const scalar dt,
    cachedRandom& rndGen,
    vector& pos1,
    scalar& m1,
    scalar& d1,
    scalar& N1,
    vector& U1,
    scalar& rho1,
    scalar& T1,
    scalarField& Y1,
    const scalar sigma1,
    const label celli,
    const scalar voli,
    vector& pos2,
    scalar& m2,
    scalar& d2,
    scalar& N2,
    vector& U2,
    scalar& rho2,
    scalar& T2,
    scalarField& Y2,
    const scalar sigma2,
    const label cellj,
    const scalar volj
) const
{
    notImplemented
    (
        "bool Foam::StochasticCollisionModel<CloudType>::update"
        "("
            "const scalar, "
            "cachedRandom&, "
            "vector&, "
            "scalar&, "
            "scalar&, "
            "scalar&, "
            "vector&, "
            "scalar&, "
            "scalar&, "
            "scalarField&, "
            "const scalar, "
            "const label, "
            "const scalar, "
            "vector&, "
            "scalar&, "
            "scalar&, "
            "scalar&, "
            "vector&, "
            "scalar&, "
            "scalar&, "
            "scalarField&, "
            "const scalar, "
            "const label, "
            "const scalar"
        ") const"
    );

    return false;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "StochasticCollisionModelNew.C"

// ************************************************************************* //

