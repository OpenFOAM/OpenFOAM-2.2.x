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

#include "treeDataEdge.H"
#include "indexedOctree.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(treeDataEdge, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::treeBoundBox Foam::treeDataEdge::calcBb(const label edgeI) const
{
    const edge& e = edges_[edgeI];
    const point& p0 = points_[e[0]];
    const point& p1 = points_[e[1]];

    return treeBoundBox(min(p0, p1), max(p0, p1));
}


void Foam::treeDataEdge::update()
{
    if (cacheBb_)
    {
        bbs_.setSize(edgeLabels_.size());

        forAll(edgeLabels_, i)
        {
            bbs_[i] = calcBb(edgeLabels_[i]);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::treeDataEdge::treeDataEdge
(
    const bool cacheBb,
    const edgeList& edges,
    const pointField& points,
    const labelUList& edgeLabels
)
:
    edges_(edges),
    points_(points),
    edgeLabels_(edgeLabels),
    cacheBb_(cacheBb)
{
    update();
}


Foam::treeDataEdge::treeDataEdge
(
    const bool cacheBb,
    const edgeList& edges,
    const pointField& points,
    const Xfer<labelList>& edgeLabels
)
:
    edges_(edges),
    points_(points),
    edgeLabels_(edgeLabels),
    cacheBb_(cacheBb)
{
    update();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::pointField Foam::treeDataEdge::shapePoints() const
{
    pointField eMids(edgeLabels_.size());

    forAll(edgeLabels_, i)
    {
        const edge& e = edges_[edgeLabels_[i]];

        eMids[i] = e.centre(points_);
    }
    return eMids;
}


//- Get type (inside,outside,mixed,unknown) of point w.r.t. surface.
//  Only makes sense for closed surfaces.
Foam::label Foam::treeDataEdge::getVolumeType
(
    const indexedOctree<treeDataEdge>& oc,
    const point& sample
) const
{
    return indexedOctree<treeDataEdge>::UNKNOWN;
}


// Check if any point on shape is inside cubeBb.
bool Foam::treeDataEdge::overlaps
(
    const label index,
    const treeBoundBox& cubeBb
) const
{
    const edge& e = edges_[edgeLabels_[index]];

    const point& start = points_[e.start()];
    const point& end = points_[e.end()];

    point intersect;

    return cubeBb.intersects(start, end, intersect);
}


// Check if any point on shape is inside sphere.
bool Foam::treeDataEdge::overlaps
(
    const label index,
    const point& centre,
    const scalar radiusSqr
) const
{
    const edge& e = edges_[edgeLabels_[index]];

    const pointHit nearHit = e.line(points_).nearestDist(centre);

    const scalar distSqr = sqr(nearHit.distance());

    if (distSqr <= radiusSqr)
    {
        return true;
    }

    return false;
}


// Calculate nearest point to sample. Updates (if any) nearestDistSqr, minIndex,
// nearestPoint.
void Foam::treeDataEdge::findNearest
(
    const labelUList& indices,
    const point& sample,

    scalar& nearestDistSqr,
    label& minIndex,
    point& nearestPoint
) const
{
    forAll(indices, i)
    {
        const label index = indices[i];

        const edge& e = edges_[edgeLabels_[index]];

        pointHit nearHit = e.line(points_).nearestDist(sample);

        scalar distSqr = sqr(nearHit.distance());

        if (distSqr < nearestDistSqr)
        {
            nearestDistSqr = distSqr;
            minIndex = index;
            nearestPoint = nearHit.rawPoint();
        }
    }
}


//- Calculates nearest (to line) point in shape.
//  Returns point and distance (squared)
void Foam::treeDataEdge::findNearest
(
    const labelUList& indices,
    const linePointRef& ln,

    treeBoundBox& tightest,
    label& minIndex,
    point& linePoint,
    point& nearestPoint
) const
{
    // Best so far
    scalar nearestDistSqr = magSqr(linePoint - nearestPoint);

    forAll(indices, i)
    {
        const label index = indices[i];

        const edge& e = edges_[edgeLabels_[index]];

        // Note: could do bb test ? Worthwhile?

        // Nearest point on line
        point ePoint, lnPt;
        scalar dist = e.line(points_).nearestDist(ln, ePoint, lnPt);
        scalar distSqr = sqr(dist);

        if (distSqr < nearestDistSqr)
        {
            nearestDistSqr = distSqr;
            minIndex = index;
            linePoint = lnPt;
            nearestPoint = ePoint;

            {
                point& minPt = tightest.min();
                minPt = min(ln.start(), ln.end());
                minPt.x() -= dist;
                minPt.y() -= dist;
                minPt.z() -= dist;
            }
            {
                point& maxPt = tightest.max();
                maxPt = max(ln.start(), ln.end());
                maxPt.x() += dist;
                maxPt.y() += dist;
                maxPt.z() += dist;
            }
        }
    }
}


// ************************************************************************* //
