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

#include "TrajectoryCollision.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::TrajectoryCollision<CloudType>::TrajectoryCollision
(
    const dictionary& dict,
    CloudType& owner
)
:
    StochasticCollisionModel<CloudType>(dict, owner, typeName),
    cSpace_(readScalar(this->coeffDict().lookup("cSpace"))),
    cTime_(readScalar(this->coeffDict().lookup("cTime"))),
    coalescence_(this->coeffDict().lookup("coalescence"))
{}


template<class CloudType>
Foam::TrajectoryCollision<CloudType>::TrajectoryCollision
(
    const TrajectoryCollision<CloudType>& cm
)
:
    StochasticCollisionModel<CloudType>(cm),
    cSpace_(cm.cSpace_),
    cTime_(cm.cTime_),
    coalescence_(cm.coalescence_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::TrajectoryCollision<CloudType>::~TrajectoryCollision()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::TrajectoryCollision<CloudType>::update
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
    bool coalescence = false;

    vector vRel = U1 - U2;

    vector p = pos2 - pos1;
    scalar dist = mag(p);

    scalar vAlign = vRel & (p/(dist + SMALL));

    if (vAlign > 0)
    {
        scalar sumD = d1 + d2;

        if (vAlign*dt > dist - 0.5*sumD)
        {
            scalar v1Mag = mag(U1);
            scalar v2Mag = mag(U2);
            vector nv1 = U1/v1Mag;
            vector nv2 = U2/v2Mag;

            scalar v1v2 = nv1 & nv2;
            scalar v1p = nv1 & p;
            scalar v2p = nv2 & p;

            scalar det = 1.0 - v1v2*v1v2;

            scalar alpha = 1.0e+20;
            scalar beta = 1.0e+20;

            if (mag(det) > 1.0e-4)
            {
                beta = -(v2p - v1v2*v1p)/det;
                alpha = v1p + v1v2*beta;
            }

            alpha /= v1Mag*dt;
            beta /= v2Mag*dt;

            // is collision possible within this timestep
            if ((alpha>0) && (alpha<1.0) && (beta>0) && (beta<1.0))
            {
                vector p1c = pos1 + alpha*U1*dt;
                vector p2c = pos2 + beta*U2*dt;

                scalar closestDist = mag(p1c-p2c);

                scalar collProb =
                    pow(0.5*sumD/max(0.5*sumD, closestDist), cSpace_)
                  * exp(-cTime_*mag(alpha-beta));

                scalar xx = rndGen.sample01<scalar>();

                // collision occur
                if ((xx < collProb) && (m1 > VSMALL) && (m2 > VSMALL))
                {
                    if (d1 > d2)
                    {
                        coalescence = collideSorted
                        (
                            dt,
                            rndGen,
                            pos1,
                            m1,
                            d1,
                            N1,
                            U1,
                            rho1,
                            T1,
                            Y1,
                            sigma1,
                            celli,
                            voli,
                            pos2,
                            m2,
                            d2,
                            N2,
                            U2,
                            rho2,
                            T2,
                            Y2,
                            sigma2,
                            cellj,
                            volj
                        );
                    }
                    else
                    {
                        coalescence = collideSorted
                        (
                            dt,
                            rndGen,
                            pos2,
                            m2,
                            d2,
                            N2,
                            U2,
                            rho2,
                            T2,
                            Y2,
                            sigma2,
                            cellj,
                            volj,
                            pos1,
                            m1,
                            d1,
                            N1,
                            U1,
                            rho1,
                            T1,
                            Y1,
                            sigma1,
                            celli,
                            voli
                        );
                    }
                }
            }
        }
    }

    return coalescence;
}


template<class CloudType>
bool Foam::TrajectoryCollision<CloudType>::collideSorted
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
    bool coalescence = false;

    vector vRel = U1 - U2;

    scalar mdMin = m2/N2;

    scalar mTot = m1 + m2;

    scalar gamma = d1/max(d2, 1.0e-12);
    scalar f = gamma*gamma*gamma + 2.7*gamma - 2.4*gamma*gamma;

    vector momMax = m1*U1;
    vector momMin = m2*U2;

    // use mass-averaged temperature to calculate We number
    scalar Tm = (T1*m1 + T2*m2)/mTot;

    // and mass averaged fractions ...
    //scalarField Yav((m1*Y1 + m2*Y2)/mTot;

    // interpolate the averaged surface tension
    scalar sigma = sigma1 + (sigma2 - sigma1)*(Tm - T1)/(T2 - T1);

    sigma = max(1.0e-6, sigma);
    scalar Vtot = m1/rho1 + m2/rho2;
    scalar rho = mTot/Vtot;

    scalar dMean = sqrt(d1*d2);
    scalar WeColl = max(1.0e-12, 0.5*rho*magSqr(vRel)*dMean/sigma);

    scalar coalesceProb = min(1.0, 2.4*f/WeColl);

    scalar prob = rndGen.sample01<scalar>();

    // Coalescence
    if ( prob < coalesceProb && coalescence_)
    {
        coalescence = true;
        // How 'many' of the droplets coalesce
        scalar nProb = prob*N2/N1;

        // Conservation of mass, momentum and energy
        scalar m2Org = m2;
        scalar dm = N1*nProb*mdMin;
        m2 -= dm;
        scalar V2 = constant::mathematical::pi*pow3(d2)/6.0;
        N2 = m2/(rho2*V2);

        scalar m1Org = m1;
        m1 += dm;
        T1 = (Tm*mTot - m2*T2)/m1;

        U1 =(momMax + (1.0 - m2/m2Org)*momMin)/m1;

        // update the liquid mass fractions
        Y1 = (m1Org*Y1 + dm*Y2)/m1;
    }
    // Grazing collision (no coalescence)
    else
    {
        scalar gf = sqrt(prob) - sqrt(coalesceProb);
        scalar denom = 1.0 - sqrt(coalesceProb);
        if (denom < 1.0e-5)
        {
            denom = 1.0;
        }
        gf /= denom;

        // if gf negative, this means that coalescence is turned off
        // and these parcels should have coalesced
        gf = max(0.0, gf);

        // gf -> 1 => v1p -> p1().U() ...
        // gf -> 0 => v1p -> momentum/(m1 + m2)

        vector mr = m1*U1 + m2*U2;
        vector v1p = (mr + m2*gf*vRel)/(m1 + m2);
        vector v2p = (mr - m1*gf*vRel)/(m1 + m2);

        if (N1 < N2)
        {
            U1 = v1p;
            U2 = (N1*v2p + (N2 - N1)*U2)/N2;
        }
        else
        {
            U1 = (N2*v1p + (N1 - N2)*U1)/N1;
            U2 = v2p;
        }
    }

    return coalescence;
}


// ************************************************************************* //
