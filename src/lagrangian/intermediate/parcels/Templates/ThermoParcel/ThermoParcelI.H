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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::ThermoParcel<ParcelType>::constantProperties::constantProperties()
:
    ParcelType::constantProperties(),
    T0_(0.0),
    TMin_(0.0),
    TMax_(VGREAT),
    Cp0_(0.0),
    epsilon0_(0.0),
    f0_(0.0)
{}


template<class ParcelType>
inline Foam::ThermoParcel<ParcelType>::constantProperties::constantProperties
(
    const constantProperties& cp
)
:
    ParcelType::constantProperties(cp),
    T0_(cp.T0_),
    TMin_(cp.TMin_),
    TMax_(cp.TMax_),
    Cp0_(cp.Cp0_),
    epsilon0_(cp.epsilon0_),
    f0_(cp.f0_)
{}


template<class ParcelType>
inline Foam::ThermoParcel<ParcelType>::constantProperties::constantProperties
(
    const dictionary& parentDict,
    const bool readFields
)
:
    ParcelType::constantProperties(parentDict, readFields),
    T0_(0.0),
    TMin_(200),
    TMax_(5000),
    Cp0_(0.0),
    epsilon0_(0.0),
    f0_(0.0)
{
    if (readFields)
    {
        if (this->dict().readIfPresent("TMin", TMin_))
        {
            Info<< "    employing parcel TMin of " << TMin_ << endl;
        }
        if (this->dict().readIfPresent("TMax", TMax_))
        {
            Info<< "    employing parcel TMax of " << TMax_ << endl;
        }

        this->dict().lookup("T0") >> T0_;
        this->dict().lookup("Cp0") >> Cp0_;
        this->dict().lookup("epsilon0") >> epsilon0_;
        this->dict().lookup("f0") >> f0_;
    }
}


template<class ParcelType>
inline Foam::ThermoParcel<ParcelType>::constantProperties::constantProperties
(
    const label parcelTypeId,
    const scalar rhoMin,
    const scalar rho0,
    const scalar minParticleMass,
    const scalar youngsModulus,
    const scalar poissonsRatio,
    const scalar T0,
    const scalar TMin,
    const scalar TMax,
    const scalar Cp0,
    const scalar epsilon0,
    const scalar f0
)
:
    ParcelType::constantProperties
    (
        parcelTypeId,
        rhoMin,
        rho0,
        minParticleMass,
        youngsModulus,
        poissonsRatio
    ),
    T0_(T0),
    TMin_(TMin),
    TMax_(TMax),
    Cp0_(Cp0),
    epsilon0_(epsilon0),
    f0_(f0)
{}


template<class ParcelType>
inline Foam::ThermoParcel<ParcelType>::ThermoParcel
(
    const polyMesh& mesh,
    const vector& position,
    const label cellI,
    const label tetFaceI,
    const label tetPtI
)
:
    ParcelType(mesh, position, cellI, tetFaceI, tetPtI),
    T_(0.0),
    Cp_(0.0),
    Tc_(0.0),
    Cpc_(0.0)
{}


template<class ParcelType>
inline Foam::ThermoParcel<ParcelType>::ThermoParcel
(
    const polyMesh& mesh,
    const vector& position,
    const label cellI,
    const label tetFaceI,
    const label tetPtI,
    const label typeId,
    const scalar nParticle0,
    const scalar d0,
    const scalar dTarget0,
    const vector& U0,
    const vector& f0,
    const vector& angularMomentum0,
    const vector& torque0,
    const constantProperties& constProps
)
:
    ParcelType
    (
        mesh,
        position,
        cellI,
        tetFaceI,
        tetPtI,
        typeId,
        nParticle0,
        d0,
        dTarget0,
        U0,
        f0,
        angularMomentum0,
        torque0,
        constProps
    ),
    T_(constProps.T0()),
    Cp_(constProps.Cp0()),
    Tc_(0.0),
    Cpc_(0.0)
{}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

template<class ParcelType>
inline Foam::scalar
Foam::ThermoParcel<ParcelType>::constantProperties::T0() const
{
    return T0_;
}


template<class ParcelType>
inline Foam::scalar
Foam::ThermoParcel<ParcelType>::constantProperties::TMin() const
{
    return TMin_;
}


template<class ParcelType>
inline Foam::scalar
Foam::ThermoParcel<ParcelType>::constantProperties::TMax() const
{
    return TMax_;
}


template<class ParcelType>
inline Foam::scalar&
Foam::ThermoParcel<ParcelType>::constantProperties::TMax()
{
    return TMax_;
}


template<class ParcelType>
inline Foam::scalar
Foam::ThermoParcel<ParcelType>::constantProperties::Cp0() const
{
    return Cp0_;
}


template<class ParcelType>
inline Foam::scalar
Foam::ThermoParcel<ParcelType>::constantProperties::epsilon0() const
{
    return epsilon0_;
}


template<class ParcelType>
inline Foam::scalar
Foam::ThermoParcel<ParcelType>::constantProperties::f0() const
{
    return f0_;
}


// * * * * * * * * * * ThermoParcel Member Functions * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::scalar Foam::ThermoParcel<ParcelType>::T() const
{
    return T_;
}


template<class ParcelType>
inline Foam::scalar Foam::ThermoParcel<ParcelType>::Cp() const
{
    return Cp_;
}


template<class ParcelType>
inline Foam::scalar Foam::ThermoParcel<ParcelType>::hs() const
{
    return Cp_*(T_ - 298.15);
}


template<class ParcelType>
inline Foam::scalar Foam::ThermoParcel<ParcelType>::Tc() const
{
    return Tc_;
}


template<class ParcelType>
inline Foam::scalar Foam::ThermoParcel<ParcelType>::Cpc() const
{
    return Cpc_;
}


template<class ParcelType>
inline Foam::scalar& Foam::ThermoParcel<ParcelType>::T()
{
    return T_;
}


template<class ParcelType>
inline Foam::scalar& Foam::ThermoParcel<ParcelType>::Cp()
{
    return Cp_;
}


// ************************************************************************* //
