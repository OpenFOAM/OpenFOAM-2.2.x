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

#include "triSurfaceSearch.H"
#include "indexedOctree.H"
#include "boolList.H"
#include "treeDataTriSurface.H"
#include "triSurface.H"
#include "line.H"
#include "cpuTime.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::point Foam::triSurfaceSearch::greatPoint(GREAT, GREAT, GREAT);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from surface. Holds reference!
Foam::triSurfaceSearch::triSurfaceSearch(const triSurface& surface)
:
    surface_(surface),
    treePtr_(NULL)
{
    // Random number generator. Bit dodgy since not exactly random ;-)
    Random rndGen(65431);

    // Slightly extended bb. Slightly off-centred just so on symmetric
    // geometry there are less face/edge aligned items.
    treeBoundBox treeBb
    (
        treeBoundBox(surface_.points(), surface_.meshPoints()).extend
        (
            rndGen,
            1e-4
        )
    );
    treeBb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
    treeBb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

    treePtr_.reset
    (
        new indexedOctree<treeDataTriSurface>
        (
            treeDataTriSurface
            (
                surface_,
                indexedOctree<treeDataTriSurface>::perturbTol()
            ),
            treeBb,
            8,      // maxLevel
            10,     // leafsize
            3.0     // duplicity
        )
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Determine inside/outside for samples
Foam::boolList Foam::triSurfaceSearch::calcInside
(
    const pointField& samples
) const
{
    boolList inside(samples.size());

    forAll(samples, sampleI)
    {
        const point& sample = samples[sampleI];

        if (!tree().bb().contains(sample))
        {
            inside[sampleI] = false;
        }
        else if
        (
            tree().getVolumeType(sample)
         == indexedOctree<treeDataTriSurface>::INSIDE
        )
        {
            inside[sampleI] = true;
        }
        else
        {
            inside[sampleI] = false;
        }
    }
    return inside;
}


Foam::labelList Foam::triSurfaceSearch::calcNearestTri
(
    const pointField& samples,
    const vector& span
) const
{
    labelList nearest(samples.size());

    const scalar nearestDistSqr = 0.25*magSqr(span);

    pointIndexHit hitInfo;

    forAll(samples, sampleI)
    {
        hitInfo = tree().findNearest(samples[sampleI], nearestDistSqr);

        if (hitInfo.hit())
        {
            nearest[sampleI] = hitInfo.index();
        }
        else
        {
            nearest[sampleI] = -1;
        }
    }

    return nearest;
}


// Nearest point on surface
Foam::tmp<Foam::pointField> Foam::triSurfaceSearch::calcNearest
(
    const pointField& samples,
    const vector& span
) const
{
    const scalar nearestDistSqr = 0.25*magSqr(span);

    tmp<pointField> tnearest(new pointField(samples.size()));
    pointField& nearest = tnearest();

    pointIndexHit hitInfo;

    forAll(samples, sampleI)
    {
        hitInfo = tree().findNearest(samples[sampleI], nearestDistSqr);

        if (hitInfo.hit())
        {
            nearest[sampleI] = hitInfo.hitPoint();
        }
        else
        {
            nearest[sampleI] = greatPoint;
        }
    }

    return tnearest;
}


Foam::pointIndexHit Foam::triSurfaceSearch::nearest
(
    const point& pt,
    const vector& span
)
const
{
    const scalar nearestDistSqr = 0.25*magSqr(span);

    return tree().findNearest(pt, nearestDistSqr);
}


void Foam::triSurfaceSearch::findLineAll
(
    const point& start,
    const point& end,
    List<pointIndexHit>& hits
)
const
{
    // See if any intersection between pt and end
    pointIndexHit inter = tree().findLine(start, end);

    if (inter.hit())
    {
        hits.setSize(1);
        hits[0] = inter;

        const vector dirVec(end-start);
        const scalar magSqrDirVec(magSqr(dirVec));
        const vector smallVec
        (
            indexedOctree<treeDataTriSurface>::perturbTol()*dirVec
          + vector(ROOTVSMALL,ROOTVSMALL,ROOTVSMALL)
        );


        // Initial perturbation amount
        vector perturbVec(smallVec);

        while (true)
        {
            // Start tracking from last hit.
            point pt = hits.last().hitPoint() + perturbVec;

            if (((pt-start)&dirVec) > magSqrDirVec)
            {
                return;
            }

            // See if any intersection between pt and end
            pointIndexHit inter = tree().findLine(pt, end);

            if (!inter.hit())
            {
                return;
            }

            // Check if already found this intersection
            bool duplicateHit = false;
            forAllReverse(hits, i)
            {
                if (hits[i].index() == inter.index())
                {
                    duplicateHit = true;
                    break;
                }
            }


            if (duplicateHit)
            {
                // Hit same triangle again. Increase perturbVec and try again.
                perturbVec *= 2;
            }
            else
            {
                // Proper hit
                label sz = hits.size();
                hits.setSize(sz+1);
                hits[sz] = inter;
                // Restore perturbVec
                perturbVec = smallVec;
            }
        }
    }
    else
    {
        hits.clear();
    }
}


// ************************************************************************* //
