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

#include "treeDataTriSurface.H"
#include "triSurfaceTools.H"
#include "triangleFuncs.H"
#include "indexedOctree.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(treeDataTriSurface, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// // Fast distance to triangle calculation. From
// // "Distance Between Point and Triangle in 3D"
// // David Eberly, Magic Software Inc. Aug. 2003.
// // Works on function Q giving distance to point and tries to minimize this.
// Foam::scalar Foam::treeDataTriSurface::nearestCoords
// (
//     const point& base,
//     const point& E0,
//     const point& E1,
//     const scalar a,
//     const scalar b,
//     const scalar c,
//     const point& P,
//     scalar& s,
//     scalar& t
// )
// {
//     // distance vector
//     const vector D(base - P);

//     // Precalculate distance factors.
//     const scalar d = E0 & D;
//     const scalar e = E1 & D;

//     // Do classification
//     const scalar det = a*c - b*b;

//     s = b*e - c*d;
//     t = b*d - a*e;

//     if (s+t < det)
//     {
//         if (s < 0)
//         {
//             if (t < 0)
//             {
//                 //region 4
//                 if (e > 0)
//                 {
//                     //min on edge t = 0
//                     t = 0;
//                     s = (d >= 0 ? 0 : (-d >= a ? 1 : -d/a));
//                 }
//                 else
//                 {
//                     //min on edge s=0
//                     s = 0;
//                     t = (e >= 0 ? 0 : (-e >= c ? 1 : -e/c));
//                 }
//             }
//             else
//             {
//                 //region 3. Min on edge s = 0
//                 s = 0;
//                 t = (e >= 0 ? 0 : (-e >= c ? 1 : -e/c));
//             }
//         }
//         else if (t < 0)
//         {
//             //region 5
//             t = 0;
//             s = (d >= 0 ? 0 : (-d >= a ? 1 : -d/a));
//         }
//         else
//         {
//             //region 0
//             const scalar invDet = 1/det;
//             s *= invDet;
//             t *= invDet;
//         }
//     }
//     else
//     {
//         if (s < 0)
//         {
//             //region 2
//             const scalar tmp0 = b + d;
//             const scalar tmp1 = c + e;
//             if (tmp1 > tmp0)
//             {
//                 //min on edge s+t=1
//                 const scalar numer = tmp1 - tmp0;
//                 const scalar denom = a-2*b+c;
//                 s = (numer >= denom ? 1 : numer/denom);
//                 t = 1 - s;
//             }
//             else
//             {
//                 //min on edge s=0
//                 s = 0;
//                 t = (tmp1 <= 0 ? 1 : (e >= 0 ? 0 : - e/c));
//             }
//         }
//         else if (t < 0)
//         {
//             //region 6
//             const scalar tmp0 = b + d;
//             const scalar tmp1 = c + e;
//             if (tmp1 > tmp0)
//             {
//                 //min on edge s+t=1
//                 const scalar numer = tmp1 - tmp0;
//                 const scalar denom = a-2*b+c;
//                 s = (numer >= denom ? 1 : numer/denom);
//                 t = 1 - s;
//             }
//             else
//             {
//                 //min on edge t=0
//                 t = 0;
//                 s = (tmp1 <= 0 ? 1 : (d >= 0 ? 0 : - d/a));
//             }
//         }
//         else
//         {
//             //region 1
//             const scalar numer = c+e-(b+d);
//             if (numer <= 0)
//             {
//                 s = 0;
//             }
//             else
//             {
//                 const scalar denom = a-2*b+c;
//                 s = (numer >= denom ? 1 : numer/denom);
//             }
//         }
//         t = 1 - s;
//     }


//     // Calculate distance.
//     // Note: abs should not be needed but truncation error causes problems
//     // with points very close to one of the triangle vertices.
//     // (seen up to -9e-15). Alternatively add some small value.

//     const scalar f = D & D;
//     return Foam::mag(a*s*s + 2*b*s*t + c*t*t + 2*d*s + 2*e*t + f);
// }


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::treeDataTriSurface::treeDataTriSurface
(
    const triSurface& surface,
    const scalar planarTol
)
:
    surface_(surface),
    planarTol_(planarTol)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::pointField Foam::treeDataTriSurface::shapePoints() const
{
    const pointField& points = surface_.points();

    pointField centres(surface_.size());

    forAll(surface_, triI)
    {
        centres[triI] = surface_[triI].centre(points);
    }
    return centres;
}


//- Get type of sample (inside/outside/mixed) w.r.t. surface.
Foam::label Foam::treeDataTriSurface::getVolumeType
(
    const indexedOctree<treeDataTriSurface>& tree,
    const point& sample
) const
{
    // Find nearest point
    const treeBoundBox& treeBb = tree.bb();

    pointIndexHit pHit = tree.findNearest
    (
        sample,
        max
        (
            Foam::sqr(GREAT),
            Foam::magSqr(treeBb.span())
        )
    );

    if (!pHit.hit())
    {
        FatalErrorIn("treeDataTriSurface::getVolumeType(..)")
            << "treeBb:" << treeBb
            << " sample:" << sample
            << " pHit:" << pHit
            << abort(FatalError);
    }

    triSurfaceTools::sideType t = triSurfaceTools::surfaceSide
    (
        surface_,
        sample,
        pHit.index()
    );

    if (t == triSurfaceTools::UNKNOWN)
    {
        return indexedOctree<treeDataTriSurface>::UNKNOWN;
    }
    else if (t == triSurfaceTools::INSIDE)
    {
        return indexedOctree<treeDataTriSurface>::INSIDE;
    }
    else if (t == triSurfaceTools::OUTSIDE)
    {
        return indexedOctree<treeDataTriSurface>::OUTSIDE;
    }
    else
    {
        FatalErrorIn("treeDataTriSurface::getVolumeType(..)")
            << "problem" << abort(FatalError);
        return indexedOctree<treeDataTriSurface>::UNKNOWN;
    }
}


// Check if any point on triangle is inside cubeBb.
bool Foam::treeDataTriSurface::overlaps
(
    const label index,
    const treeBoundBox& cubeBb
) const
{
    const pointField& points = surface_.points();
    const labelledTri& f = surface_[index];

    treeBoundBox triBb(points, surface_[index]);

    //- For testing: robust one
    //return cubeBb.overlaps(triBb);

    //- Exact test of triangle intersecting bb

    // Quick rejection. If whole bounding box of tri is outside cubeBb then
    // there will be no intersection.
    if (!cubeBb.overlaps(triBb))
    {
        return false;
    }

    // Check if one or more triangle point inside
    if (cubeBb.containsAny(points, f))
    {
        return true;
    }

    // Triangle points
    const point& p0 = points[f[0]];
    const point& p1 = points[f[1]];
    const point& p2 = points[f[2]];


    // Now we have the difficult case: all points are outside but connecting
    // edges might go through cube. Use fast intersection of bounding box.

    //return triangleFuncs::intersectBbExact(p0, p1, p2, cubeBb);
    return triangleFuncs::intersectBb(p0, p1, p2, cubeBb);
}


// Calculate nearest point to sample. Updates (if any) nearestDistSqr, minIndex,
// nearestPoint.
void Foam::treeDataTriSurface::findNearest
(
    const labelUList& indices,
    const point& sample,

    scalar& nearestDistSqr,
    label& minIndex,
    point& nearestPoint
) const
{
    const pointField& points = surface_.points();

    forAll(indices, i)
    {
        label index = indices[i];
        const triSurface::FaceType& f = surface_[index];

        ////- Possible optimization: do quick rejection of triangle if bounding
        ////  sphere does not intersect triangle bounding box. From simplistic
        ////  test was not found to speed up things.
        //
        //// Triangle bounding box.
        //point triBbMin = min(p0, min(p1, p2));
        //point triBbMax = max(p0, max(p1, p2));
        //
        //if
        //(
        //    indexedOctree<treeDataTriSurface>::intersects
        //    (
        //        triBbMin,
        //        triBbMax,
        //        nearestDistSqr,
        //        sample
        //    )
        //)
        {
            // // Get spanning vectors of triangle
            // vector base(p1);
            // vector E0(p0 - p1);
            // vector E1(p2 - p1);

            // scalar a(E0& E0);
            // scalar b(E0& E1);
            // scalar c(E1& E1);

            // // Get nearest point in s,t coordinates (s is along E0, t
            // // is along E1)
            // scalar s;
            // scalar t;

            // scalar distSqr = nearestCoords
            // (
            //     base,
            //     E0,
            //     E1,
            //     a,
            //     b,
            //     c,
            //     sample,

            //     s,
            //     t
            // );

            pointHit pHit = f.nearestPoint(sample, points);

            scalar distSqr = sqr(pHit.distance());

            if (distSqr < nearestDistSqr)
            {
                nearestDistSqr = distSqr;
                minIndex = index;
                nearestPoint = pHit.rawPoint();
            }
        }
    }
}


// Calculate nearest point to line. Updates (if any) nearestDistSqr, minIndex,
// nearestPoint.
void Foam::treeDataTriSurface::findNearest
(
    const labelUList& indices,
    const linePointRef& ln,

    treeBoundBox& tightest,
    label& minIndex,
    point& linePoint,
    point& nearestPoint
) const
{
    notImplemented
    (
        "treeDataTriSurface::findNearest(const labelUList&"
        ", const linePointRef&, treeBoundBox&, label&, point&, point&) const"
    );
}


bool Foam::treeDataTriSurface::intersects
(
    const label index,
    const point& start,
    const point& end,
    point& intersectionPoint
) const
{
    const pointField& points = surface_.points();

    const triSurface::FaceType& f = surface_[index];

    // Do quick rejection test
    treeBoundBox triBb(points, f);

    const direction startBits(triBb.posBits(start));
    const direction endBits(triBb.posBits(end));

    if ((startBits & endBits) != 0)
    {
        // start and end in same block outside of triBb.
        return false;
    }

    const vector dir(end - start);

    // Use relative tolerance (from octree) to determine intersection.
    pointHit inter = f.intersection
    (
        start,
        dir,
        points,
        intersection::HALF_RAY,
        planarTol_
    );

    if (inter.hit() && inter.distance() <= 1)
    {
        // Note: no extra test on whether intersection is in front of us
        // since using half_ray.
        intersectionPoint = inter.hitPoint();

        return true;
    }
    else
    {
        return false;
    }


    //- Using exact intersections
    //pointHit info = f.tri(points).intersectionExact(start, end);
    //
    //if (info.hit())
    //{
    //    intersectionPoint = info.hitPoint();
    //}
    //return info.hit();
}


// ************************************************************************* //
