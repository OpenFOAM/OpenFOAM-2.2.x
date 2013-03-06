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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::CollidingParcel<ParcelType>::CollidingParcel
(
    const CollidingParcel<ParcelType>& p
)
:
    ParcelType(p),
    collisionRecords_(p.collisionRecords_)
{}


template<class ParcelType>
Foam::CollidingParcel<ParcelType>::CollidingParcel
(
    const CollidingParcel<ParcelType>& p,
    const polyMesh& mesh
)
:
    ParcelType(p, mesh),
    collisionRecords_(p.collisionRecords_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
bool Foam::CollidingParcel<ParcelType>::move
(
    TrackData& td,
    const scalar trackTime
)
{
    typename TrackData::cloudType::parcelType& p =
        static_cast<typename TrackData::cloudType::parcelType&>(*this);

    td.switchProcessor = false;
    td.keepParticle = true;

    const polyMesh& mesh = td.cloud().pMesh();
    const polyBoundaryMesh& pbMesh = mesh.boundaryMesh();
    const scalarField& V = mesh.cellVolumes();

    switch (td.part())
    {
        case TrackData::tpVelocityHalfStep:
        {
            // First and last leapfrog velocity adjust part, required
            // before and after tracking and force calculation

            p.U() += 0.5*trackTime*p.f()/p.mass();

            p.angularMomentum() += 0.5*trackTime*p.torque();

            break;
        }

        case TrackData::tpLinearTrack:
        {
            scalar tEnd = (1.0 - p.stepFraction())*trackTime;
            const scalar dtMax = tEnd;

            while (td.keepParticle && !td.switchProcessor && tEnd > ROOTVSMALL)
            {
                // Apply correction to position for reduced-D cases
                meshTools::constrainToMeshCentre(mesh, p.position());

                // Set the Lagrangian time-step
                scalar dt = min(dtMax, tEnd);

                // Remember which cell the parcel is in since this
                // will change if a face is hit
                const label cellI = p.cell();

                const scalar magU = mag(p.U());
                if (p.active() && magU > ROOTVSMALL)
                {
                    const scalar d = dt*magU;
                    const scalar maxCo = td.cloud().solution().maxCo();
                    const scalar dCorr = min(d, maxCo*cbrt(V[cellI]));
                    dt *=
                        dCorr/d
                       *p.trackToFace(p.position() + dCorr*p.U()/magU, td);
                }

                tEnd -= dt;
                p.stepFraction() = 1.0 - tEnd/trackTime;

                // Avoid problems with extremely small timesteps
                if (dt > ROOTVSMALL)
                {
                    // Update cell based properties
                    p.setCellValues(td, dt, cellI);

                    if (td.cloud().solution().cellValueSourceCorrection())
                    {
                        p.cellValueSourceCorrection(td, dt, cellI);
                    }

                    p.calc(td, dt, cellI);
                }

                if (p.onBoundary() && td.keepParticle)
                {
                    if (isA<processorPolyPatch>(pbMesh[p.patch(p.face())]))
                    {
                        td.switchProcessor = true;
                    }
                }

                p.age() += dt;
            }

            break;
        }

        case TrackData::tpRotationalTrack:
        {
            notImplemented("TrackData::tpRotationalTrack");

            break;
        }

        default:
        {
            FatalErrorIn
            (
                "CollidingParcel<ParcelType>::move(TrackData&, const scalar)"
            )   << td.part() << " is an invalid part of the tracking method."
                << abort(FatalError);
        }
    }

    return td.keepParticle;
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "CollidingParcelIO.C"

// ************************************************************************* //
