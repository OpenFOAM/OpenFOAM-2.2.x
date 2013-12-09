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

#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
inline
Foam::KinematicParcel<ParcelType>::constantProperties::constantProperties()
:
    dict_(dictionary::null),
    parcelTypeId_(-1),
    rhoMin_(0.0),
    rho0_(0.0),
    minParticleMass_(0.0),
    youngsModulus_(0.0),
    poissonsRatio_(0.0)
{}


template<class ParcelType>
inline Foam::KinematicParcel<ParcelType>::constantProperties::constantProperties
(
    const constantProperties& cp
)
:
    dict_(cp.dict_),
    parcelTypeId_(cp.parcelTypeId_),
    rhoMin_(cp.rhoMin_),
    rho0_(cp.rho0_),
    minParticleMass_(cp.minParticleMass_),
    youngsModulus_(cp.youngsModulus_),
    poissonsRatio_(cp.poissonsRatio_)
{}


template<class ParcelType>
inline Foam::KinematicParcel<ParcelType>::constantProperties::constantProperties
(
    const dictionary& parentDict,
    const bool readFields
)
:
    dict_(parentDict.subOrEmptyDict("constantProperties")),
    parcelTypeId_(-1),
    rhoMin_(1e-15),
    rho0_(0.0),
    minParticleMass_(1e-15),
    youngsModulus_(0.0),
    poissonsRatio_(0.0)
{
    if (readFields)
    {
        if (dict_.readIfPresent("parcelTypeId", parcelTypeId_))
        {
            Info<< "    employing parcel parcelTypeId of " << parcelTypeId_
                << endl;
        }
        if (dict_.readIfPresent("rhoMin", rhoMin_))
        {
            Info<< "    employing parcel rhoMin of " << rhoMin_ << endl;
        }
        if (dict_.readIfPresent("minParticleMass", minParticleMass_))
        {
            Info<< "    employing parcel minParticleMass of "
                << minParticleMass_ << endl;
        }

        dict_.lookup("rho0") >> rho0_;
        dict_.lookup("youngsModulus") >> youngsModulus_;
        dict_.lookup("poissonsRatio") >> poissonsRatio_;
    }
}

template<class ParcelType>
inline Foam::KinematicParcel<ParcelType>::constantProperties::constantProperties
(
    const label parcelTypeId,
    const scalar rhoMin,
    const scalar rho0,
    const scalar minParticleMass,
    const scalar youngsModulus,
    const scalar poissonsRatio
)
:
    dict_(dictionary::null),
    parcelTypeId_(parcelTypeId),
    rhoMin_(rhoMin),
    rho0_(rho0),
    minParticleMass_(minParticleMass),
    youngsModulus_(youngsModulus),
    poissonsRatio_(poissonsRatio)
{}


template<class ParcelType>
inline Foam::KinematicParcel<ParcelType>::KinematicParcel
(
    const polyMesh& owner,
    const vector& position,
    const label cellI,
    const label tetFaceI,
    const label tetPtI
)
:
    ParcelType(owner, position, cellI, tetFaceI, tetPtI),
    active_(true),
    typeId_(-1),
    nParticle_(0),
    d_(0.0),
    dTarget_(0.0),
    U_(vector::zero),
    f_(vector::zero),
    angularMomentum_(vector::zero),
    torque_(vector::zero),
    rho_(0.0),
    age_(0.0),
    tTurb_(0.0),
    UTurb_(vector::zero),
    rhoc_(0.0),
    Uc_(vector::zero),
    muc_(0.0)
{}


template<class ParcelType>
inline Foam::KinematicParcel<ParcelType>::KinematicParcel
(
    const polyMesh& owner,
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
    ParcelType(owner, position, cellI, tetFaceI, tetPtI),
    active_(true),
    typeId_(typeId),
    nParticle_(nParticle0),
    d_(d0),
    dTarget_(dTarget0),
    U_(U0),
    f_(f0),
    angularMomentum_(angularMomentum0),
    torque_(torque0),
    rho_(constProps.rho0()),
    age_(0.0),
    tTurb_(0.0),
    UTurb_(vector::zero),
    rhoc_(0.0),
    Uc_(vector::zero),
    muc_(0.0)
{}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

template<class ParcelType>
inline const Foam::dictionary&
Foam::KinematicParcel<ParcelType>::constantProperties::dict() const
{
    return dict_;
}


template<class ParcelType>
inline Foam::label
Foam::KinematicParcel<ParcelType>::constantProperties::parcelTypeId() const
{
    return parcelTypeId_;
}


template<class ParcelType>
inline Foam::scalar
Foam::KinematicParcel<ParcelType>::constantProperties::rhoMin() const
{
    return rhoMin_;
}


template<class ParcelType>
inline Foam::scalar
Foam::KinematicParcel<ParcelType>::constantProperties::rho0() const
{
    return rho0_;
}


template<class ParcelType>
inline Foam::scalar
Foam::KinematicParcel<ParcelType>::constantProperties::minParticleMass() const
{
    return minParticleMass_;
}


template<class ParcelType>
inline Foam::scalar
Foam::KinematicParcel<ParcelType>::constantProperties::youngsModulus() const
{
    return youngsModulus_;
}


template<class ParcelType>
inline Foam::scalar
Foam::KinematicParcel<ParcelType>::constantProperties::poissonsRatio() const
{
    return poissonsRatio_;
}


// * * * * * * * KinematicParcel Member Functions  * * * * * * * //

template<class ParcelType>
inline bool Foam::KinematicParcel<ParcelType>::active() const
{
    return active_;
}


template<class ParcelType>
inline Foam::label Foam::KinematicParcel<ParcelType>::typeId() const
{
    return typeId_;
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::nParticle() const
{
    return nParticle_;
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::d() const
{
    return d_;
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::dTarget() const
{
    return dTarget_;
}


template<class ParcelType>
inline const Foam::vector& Foam::KinematicParcel<ParcelType>::U() const
{
    return U_;
}


template<class ParcelType>
inline const Foam::vector&
Foam::KinematicParcel<ParcelType>::f() const
{
    return f_;
}


template<class ParcelType>
inline const Foam::vector&
Foam::KinematicParcel<ParcelType>::angularMomentum() const
{
    return angularMomentum_;
}


template<class ParcelType>
inline const Foam::vector&
Foam::KinematicParcel<ParcelType>::torque() const
{
    return torque_;
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::rho() const
{
    return rho_;
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::age() const
{
    return age_;
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::tTurb() const
{
    return tTurb_;
}


template<class ParcelType>
inline const Foam::vector& Foam::KinematicParcel<ParcelType>::UTurb() const
{
    return UTurb_;
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::rhoc() const
{
    return rhoc_;
}


template<class ParcelType>
inline const Foam::vector& Foam::KinematicParcel<ParcelType>::Uc() const
{
    return Uc_;
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::muc() const
{
    return muc_;
}


template<class ParcelType>
inline bool& Foam::KinematicParcel<ParcelType>::active()
{
    return active_;
}


template<class ParcelType>
inline Foam::label& Foam::KinematicParcel<ParcelType>::typeId()
{
    return typeId_;
}


template<class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::nParticle()
{
    return nParticle_;
}


template<class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::d()
{
    return d_;
}


template<class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::dTarget()
{
    return dTarget_;
}


template<class ParcelType>
inline Foam::vector& Foam::KinematicParcel<ParcelType>::U()
{
    return U_;
}


template<class ParcelType>
inline Foam::vector& Foam::KinematicParcel<ParcelType>::f()
{
    return f_;
}


template<class ParcelType>
inline Foam::vector& Foam::KinematicParcel<ParcelType>::angularMomentum()
{
    return angularMomentum_;
}


template<class ParcelType>
inline Foam::vector& Foam::KinematicParcel<ParcelType>::torque()
{
    return torque_;
}


template<class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::rho()
{
    return rho_;
}


template<class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::age()
{
    return age_;
}


template<class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::tTurb()
{
    return tTurb_;
}


template<class ParcelType>
inline Foam::vector& Foam::KinematicParcel<ParcelType>::UTurb()
{
    return UTurb_;
}


template<class ParcelType>
inline Foam::label Foam::KinematicParcel<ParcelType>::faceInterpolation() const
{
    // Use volume-based interpolation if dealing with external faces
    if (this->cloud().internalFace(this->face()))
    {
        return this->face();
    }
    else
    {
        return -1;
    }
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::massCell
(
    const label cellI
) const
{
    return rhoc_*this->mesh().cellVolumes()[cellI];
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::mass() const
{
    return rho_*volume();
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::momentOfInertia() const
{
    return 0.1*mass()*sqr(d_);
}


template<class ParcelType>
inline Foam::vector Foam::KinematicParcel<ParcelType>::omega() const
{
    return angularMomentum_/momentOfInertia();
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::volume() const
{
    return volume(d_);
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::volume(const scalar d)
{
    return pi/6.0*pow3(d);
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::areaP() const
{
    return areaP(d_);
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::areaP(const scalar d)
{
    return 0.25*areaS(d);
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::areaS() const
{
    return areaS(d_);
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::areaS(const scalar d)
{
    return pi*d*d;
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::Re
(
    const vector& U,
    const scalar d,
    const scalar rhoc,
    const scalar muc
) const
{
    return rhoc*mag(U - Uc_)*d/(muc + ROOTVSMALL);
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::We
(
    const vector& U,
    const scalar d,
    const scalar rhoc,
    const scalar sigma
) const
{
    return rhoc*magSqr(U - Uc_)*d/(sigma + ROOTVSMALL);
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::Eo
(
    const vector& a,
    const scalar d,
    const scalar sigma
) const
{
    vector dir = U_/(mag(U_) + ROOTVSMALL);
    return mag(a & dir)*(rho_ - rhoc_)*sqr(d)/(sigma + ROOTVSMALL);
}


// ************************************************************************* //
