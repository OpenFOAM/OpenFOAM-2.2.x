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

#include "PilchErdman.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PilchErdman<CloudType>::PilchErdman
(
    const dictionary& dict,
    CloudType& owner
)
:
    BreakupModel<CloudType>(dict, owner, typeName),
    B1_(0.375),
    B2_(0.236)
{
    if (!this->defaultCoeffs(true))
    {
        this->coeffDict().lookup("B1") >> B1_;
        this->coeffDict().lookup("B2") >> B2_;
    }
}


template<class CloudType>
Foam::PilchErdman<CloudType>::PilchErdman(const PilchErdman<CloudType>& bum)
:
    BreakupModel<CloudType>(bum),
    B1_(bum.B1_),
    B2_(bum.B2_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PilchErdman<CloudType>::~PilchErdman()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::PilchErdman<CloudType>::update
(
    const scalar dt,
    const vector& g,
    scalar& d,
    scalar& tc,
    scalar& ms,
    scalar& nParticle,
    scalar& KHindex,
    scalar& y,
    scalar& yDot,
    const scalar d0,
    const scalar rho,
    const scalar mu,
    const scalar sigma,
    const vector& U,
    const scalar rhoc,
    const scalar muc,
    const vector& Urel,
    const scalar Urmag,
    const scalar tMom,
    scalar& dChild,
    scalar& massChild
)
{
    scalar semiMass = nParticle*pow3(d);
    scalar We = 0.5*rhoc*sqr(Urmag)*d/sigma;
    scalar Oh = mu/sqrt(rho*d*sigma);

    scalar Wec = 6.0*(1.0 + 1.077*pow(Oh, 1.6));

    if (We > Wec)
    {
        // We > 1335, wave crest stripping
        scalar taubBar = 5.5;

        if (We < 1335)
        {
            if (We > 175.0)
            {
                // sheet stripping
                taubBar = 0.766*pow(2.0*We - 12.0, 0.25);
            }
            else if (We > 22.0)
            {
                // Bag-and-stamen breakup
                taubBar = 14.1*pow(2.0*We - 12.0, -0.25);
            }
            else if (We > 9.0)
            {
                // Bag breakup
                taubBar = 2.45*pow(2.0*We - 12.0, 0.25);
            }
            else if (We > 6.0)
            {
                // Vibrational breakup
                taubBar = 6.0*pow(2.0*We - 12.0, -0.25);
            }
        }

        scalar rho12 = sqrt(rhoc/rho);

        scalar Vd = Urmag*rho12*(B1_*taubBar + B2_*taubBar*taubBar);
        scalar Vd1 = sqr(1.0 - Vd/Urmag);
        Vd1 = max(Vd1, SMALL);
        scalar Ds = 2.0*Wec*sigma/(Vd1*rhoc*sqr(Urmag));
        scalar A = Urmag*rho12/d;

        scalar taub = taubBar/A;

        scalar frac = dt/taub;

        // update the droplet diameter according to the rate eq. (implicitly)
        d = (d + frac*Ds)/(1.0 + frac);

        // correct the number of particles to conserve mass
        nParticle = semiMass/pow3(d);
    }

    return false;
}


// ************************************************************************* //
