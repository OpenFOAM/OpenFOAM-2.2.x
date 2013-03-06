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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class thermo>
inline Foam::exponentialSolidTransport<thermo>::exponentialSolidTransport
(
    const thermo& t,
    const scalar kappa0,
    const scalar n0,
    const scalar Tref
)
:
    thermo(t),
    kappa0_(kappa0),
    n0_(n0),
    Tref_(Tref)
{}


template<class thermo>
inline Foam::exponentialSolidTransport<thermo>::exponentialSolidTransport
(
    const word& name,
    const exponentialSolidTransport& ct
)
:
    thermo(name, ct),
    kappa0_(ct.kappa0_),
    n0_(ct.n0_),
    Tref_(ct.Tref_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::exponentialSolidTransport<Thermo> >
Foam::exponentialSolidTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<exponentialSolidTransport<Thermo> >
    (
        new exponentialSolidTransport<Thermo>(dict)
    );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class thermo>
inline Foam::scalar Foam::exponentialSolidTransport<thermo>::kappa
(
    const scalar p, const scalar T
) const
{
    return (kappa0_*pow(T/Tref_, n0_));
}


template<class thermo>
inline Foam::vector Foam::exponentialSolidTransport<thermo>::Kappa
(
    const scalar p, const scalar T
) const
{
    const scalar kappa(kappa0_*pow(T/Tref_, n0_));
    return vector(kappa, kappa, kappa);
}


template<class thermo>
inline Foam::scalar Foam::exponentialSolidTransport<thermo>::
mu(const scalar p, const scalar T) const
{
    notImplemented
    (
        "Foam::scalar Foam::exponentialSolidTransport<thermo>mu::"
        "("
        "    const scalar p, const scalar T"
        ") const"
    );
    return scalar(0);
}


template<class thermo>
inline Foam::scalar Foam::exponentialSolidTransport<thermo>::
alphah(const scalar p, const scalar T) const
{
    return kappa(p, T)/this->Cpv(p, T);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


template<class thermo>
inline Foam::exponentialSolidTransport<thermo>&
Foam::exponentialSolidTransport<thermo>::operator=
(
    const exponentialSolidTransport<thermo>& ct
)
{
    //thermo::operator=(ct);

    kappa0_ = ct.kappa0_;
    n0_ = ct.n0_;
    Tref_ = ct.Tref_;

    return *this;
}


template<class thermo>
inline void Foam::exponentialSolidTransport<thermo>::operator+=
(
    const exponentialSolidTransport<thermo>& ct
)
{
    scalar molr1 = this->nMoles();

    //thermo::operator+=(ct);

    molr1 /= this->nMoles();
    scalar molr2 = ct.nMoles()/this->nMoles();

    kappa0_ = molr1*kappa0_ + molr2*ct.kappa0_;
    n0_ = (molr1*n0_ + molr2*ct.n0_);
    Tref_ = (molr1*Tref_ + molr2*ct.Tref_);
}


template<class thermo>
inline void Foam::exponentialSolidTransport<thermo>::operator-=
(
    const exponentialSolidTransport<thermo>& ct
)
{
    scalar molr1 = this->nMoles();

    //thermo::operator-=(ct);

    molr1 /= this->nMoles();
    scalar molr2 = ct.nMoles()/this->nMoles();

    kappa0_ = (molr1*kappa0_ - molr2*ct.kappa0_);
    n0_ = (molr1*n0_ - molr2*ct.n0_);
    Tref_ = (molr1*Tref_ - molr2*ct.Tref_);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class thermo>
inline Foam::exponentialSolidTransport<thermo> Foam::operator*
(
    const scalar s,
    const exponentialSolidTransport<thermo>& ct
)
{
    return exponentialSolidTransport<thermo>
    (
        s*static_cast<const thermo&>(ct),
        ct.kappa0_,
        ct.n0_,
        ct.Tref_
    );
}

// ************************************************************************* //
