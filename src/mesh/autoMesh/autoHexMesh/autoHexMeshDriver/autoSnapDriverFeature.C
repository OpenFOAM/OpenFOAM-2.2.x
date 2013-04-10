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

#include "autoSnapDriver.H"
#include "polyTopoChange.H"
#include "OFstream.H"
#include "syncTools.H"
#include "fvMesh.H"
#include "OFstream.H"
#include "motionSmoother.H"
#include "refinementSurfaces.H"
#include "refinementFeatures.H"
#include "unitConversion.H"
#include "plane.H"
#include "featureEdgeMesh.H"
#include "treeDataPoint.H"
#include "indexedOctree.H"
#include "snapParameters.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    class listTransform
    {
    public:

        void operator()
        (
            const vectorTensorTransform& vt,
            const bool forward,
            List<List<point> >& fld
        ) const
        {
            const tensor T
            (
                forward
              ? vt.R()
              : vt.R().T()
            );

            forAll(fld, i)
            {
                List<point>& elems = fld[i];
                forAll(elems, elemI)
                {
                    elems[elemI] = transform(T, elems[elemI]);
                }
            }
        }
    };

    template<class T>
    class listPlusEqOp
    {
    public:

        void operator()
        (
            List<T>& x,
            const List<T>& y
        ) const
        {
            label sz = x.size();
            x.setSize(sz+y.size());
            forAll(y, i)
            {
                x[sz++] = y[i];
            }
        }
    };
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::autoSnapDriver::isFeaturePoint
(
    const scalar featureCos,
    const indirectPrimitivePatch& pp,
    const PackedBoolList& isFeatureEdge,
    const label pointI
) const
{
    const pointField& points = pp.localPoints();
    const edgeList& edges = pp.edges();
    const labelList& pEdges = pp.pointEdges()[pointI];

    label nFeatEdges = 0;

    forAll(pEdges, i)
    {
        if (isFeatureEdge[pEdges[i]])
        {
            nFeatEdges++;

            for (label j = i+1; j < pEdges.size(); j++)
            {
                if (isFeatureEdge[pEdges[j]])
                {
                    const edge& eI = edges[pEdges[i]];
                    const edge& eJ = edges[pEdges[j]];

                    const point& p = points[pointI];
                    const point& pI = points[eI.otherVertex(pointI)];
                    const point& pJ = points[eJ.otherVertex(pointI)];

                    vector vI = p-pI;
                    scalar vIMag = mag(vI);

                    vector vJ = pJ-p;
                    scalar vJMag = mag(vJ);

                    if
                    (
                        vIMag > SMALL
                     && vJMag > SMALL
                     && ((vI/vIMag & vJ/vJMag) < featureCos)
                    )
                    {
                        return true;
                    }
                }
            }
        }
    }

    if (nFeatEdges == 1)
    {
        // End of feature-edge string
        return true;
    }

    return false;
}


void Foam::autoSnapDriver::smoothAndConstrain
(
    const indirectPrimitivePatch& pp,
    const List<pointConstraint>& constraints,
    vectorField& disp
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    for (label avgIter = 0; avgIter < 20; avgIter++)
    {
        // Calculate average displacement of neighbours
        // - unconstrained (i.e. surface) points use average of all
        //   neighbouring points
        // - from testing it has been observed that it is not beneficial
        //   to have edge constrained points use average of all edge or point
        //   constrained neighbours since they're already attracted to
        //   the nearest point on the feature.
        //   Having them attract to point-constrained neighbours does not
        //   make sense either since there is usually just one of them so
        //   it severely distorts it.
        // - same for feature points. They are already attracted to the
        //   nearest feature point.

        vectorField dispSum(pp.nPoints(), vector::zero);
        labelList dispCount(pp.nPoints(), 0);

        const labelListList& pointEdges = pp.pointEdges();
        const edgeList& edges = pp.edges();

        forAll(pointEdges, pointI)
        {
            const labelList& pEdges = pointEdges[pointI];

            label nConstraints = constraints[pointI].first();

            if (nConstraints <= 1)
            {
                forAll(pEdges, i)
                {
                    label nbrPointI = edges[pEdges[i]].otherVertex(pointI);
                    if (constraints[nbrPointI].first() >= nConstraints)
                    {
                        dispSum[pointI] += disp[nbrPointI];
                        dispCount[pointI]++;
                    }
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            dispSum,
            plusEqOp<point>(),
            vector::zero,
            mapDistribute::transform()
        );
        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            dispCount,
            plusEqOp<label>(),
            0,
            mapDistribute::transform()
        );

        // Constraints
        forAll(constraints, pointI)
        {
            if (dispCount[pointI] > 0)
            {
                // Mix my displacement with neighbours' displacement
                disp[pointI] =
                    0.5
                   *(disp[pointI] + dispSum[pointI]/dispCount[pointI]);
            }
        }
    }
}
//XXXXXX
//TODO: make proper parallel so coupled edges don't have double influence
void Foam::autoSnapDriver::smoothAndConstrain2
(
    const bool applyConstraints,
    const indirectPrimitivePatch& pp,
    const List<pointConstraint>& constraints,
    vectorField& disp
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    for (label avgIter = 0; avgIter < 20; avgIter++)
    {
        vectorField dispSum(pp.nPoints(), vector::zero);
        labelList dispCount(pp.nPoints(), 0);

        const labelListList& pointEdges = pp.pointEdges();
        const edgeList& edges = pp.edges();

        forAll(pointEdges, pointI)
        {
            const labelList& pEdges = pointEdges[pointI];

            forAll(pEdges, i)
            {
                label nbrPointI = edges[pEdges[i]].otherVertex(pointI);
                dispSum[pointI] += disp[nbrPointI];
                dispCount[pointI]++;
            }
        }

        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            dispSum,
            plusEqOp<point>(),
            vector::zero,
            mapDistribute::transform()
        );
        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            dispCount,
            plusEqOp<label>(),
            0,
            mapDistribute::transform()
        );

        // Constraints
        forAll(constraints, pointI)
        {
            if (dispCount[pointI] > 0)// && constraints[pointI].first() <= 1)
            {
                // Mix my displacement with neighbours' displacement
                disp[pointI] =
                    0.5
                   *(disp[pointI] + dispSum[pointI]/dispCount[pointI]);

                if (applyConstraints)
                {
                    disp[pointI] = transform
                    (
                        constraints[pointI].constraintTransformation(),
                        disp[pointI]
                    );
                }
            }
        }
    }
}
//XXXXXX


void Foam::autoSnapDriver::calcNearestFace
(
    const label iter,
    const indirectPrimitivePatch& pp,
    vectorField& faceDisp,
    vectorField& faceSurfaceNormal,
    labelList& faceSurfaceGlobalRegion,
    vectorField& faceRotation
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const refinementSurfaces& surfaces = meshRefiner_.surfaces();

    // Displacement and orientation per pp face.
    faceDisp.setSize(pp.size());
    faceDisp = vector::zero;
    faceSurfaceNormal.setSize(pp.size());
    faceSurfaceNormal = vector::zero;
    faceSurfaceGlobalRegion.setSize(pp.size());
    faceSurfaceGlobalRegion = -1;

    // Divide surfaces into zoned and unzoned
    labelList zonedSurfaces = surfaces.getNamedSurfaces();
    labelList unzonedSurfaces = surfaces.getUnnamedSurfaces();

    // Per pp face the current surface snapped to
    labelList snapSurf(pp.size(), -1);


    // Do zoned surfaces
    // ~~~~~~~~~~~~~~~~~
    // Zoned faces only attract to corresponding surface

    // Extract faces per zone
    const wordList& faceZoneNames = surfaces.faceZoneNames();

    forAll(zonedSurfaces, i)
    {
        label zoneSurfI = zonedSurfaces[i];

        // Get indices of faces on pp that are also in zone
        label zoneI = mesh.faceZones().findZoneID(faceZoneNames[zoneSurfI]);
        if (zoneI == -1)
        {
            FatalErrorIn
            (
                "autoSnapDriver::calcNearestFace(..)"
            )   << "Problem. Cannot find zone " << faceZoneNames[zoneSurfI]
                << exit(FatalError);
        }
        const faceZone& fZone = mesh.faceZones()[zoneI];
        PackedBoolList isZonedFace(mesh.nFaces());
        forAll(fZone, i)
        {
            isZonedFace[fZone[i]] = 1;
        }

        DynamicList<label> ppFaces(fZone.size());
        DynamicList<label> meshFaces(fZone.size());
        forAll(pp.addressing(), i)
        {
            if (isZonedFace[pp.addressing()[i]])
            {
                snapSurf[i] = zoneSurfI;
                ppFaces.append(i);
                meshFaces.append(pp.addressing()[i]);
            }
        }

        //Pout<< "For faceZone " << fZone.name()
        //    << " found " << ppFaces.size() << " out of " << pp.size()
        //    << endl;

        pointField fc
        (
            indirectPrimitivePatch
            (
                IndirectList<face>(mesh.faces(), meshFaces),
                mesh.points()
            ).faceCentres()
        );

        List<pointIndexHit> hitInfo;
        labelList hitSurface;
        labelList hitRegion;
        vectorField hitNormal;
        surfaces.findNearestRegion
        (
            labelList(1, zoneSurfI),
            fc,
            sqr(scalarField(fc.size(), GREAT)),// sqr of attract dist
            hitSurface,
            hitInfo,
            hitRegion,
            hitNormal
        );

        forAll(hitInfo, hitI)
        {
            if (hitInfo[hitI].hit())
            {
                label faceI = ppFaces[hitI];
                faceDisp[faceI] = hitInfo[hitI].hitPoint() - fc[hitI];
                faceSurfaceNormal[faceI] = hitNormal[hitI];
                faceSurfaceGlobalRegion[faceI] = surfaces.globalRegion
                (
                    hitSurface[hitI],
                    hitRegion[hitI]
                );
            }
        }
    }


    // Do unzoned surfaces
    // ~~~~~~~~~~~~~~~~~~~
    // Unzoned faces attract to any unzoned surface

    DynamicList<label> ppFaces(pp.size());
    DynamicList<label> meshFaces(pp.size());
    forAll(pp.addressing(), i)
    {
        if (snapSurf[i] == -1)
        {
            ppFaces.append(i);
            meshFaces.append(pp.addressing()[i]);
        }
    }
    //Pout<< "Found " << ppFaces.size() << " unzoned faces out of "
    //   << pp.size() << endl;

    pointField fc
    (
        indirectPrimitivePatch
        (
            IndirectList<face>(mesh.faces(), meshFaces),
            mesh.points()
        ).faceCentres()
    );

    List<pointIndexHit> hitInfo;
    labelList hitSurface;
    labelList hitRegion;
    vectorField hitNormal;
    surfaces.findNearestRegion
    (
        unzonedSurfaces,
        fc,
        sqr(scalarField(fc.size(), GREAT)),// sqr of attract dist
        hitSurface,
        hitInfo,
        hitRegion,
        hitNormal
    );

    forAll(hitInfo, hitI)
    {
        if (hitInfo[hitI].hit())
        {
            label faceI = ppFaces[hitI];
            faceDisp[faceI] = hitInfo[hitI].hitPoint() - fc[hitI];
            faceSurfaceNormal[faceI] = hitNormal[hitI];
            faceSurfaceGlobalRegion[faceI] = surfaces.globalRegion
            (
                hitSurface[hitI],
                hitRegion[hitI]
            );
        }
    }


    // Determine rotation
    // ~~~~~~~~~~~~~~~~~~

    // Determine rotation axis
    faceRotation.setSize(pp.size());
    faceRotation = vector::zero;

    forAll(faceRotation, faceI)
    {
        // Note: extend to >180 degrees checking
        faceRotation[faceI] =
            pp.faceNormals()[faceI]
          ^ faceSurfaceNormal[faceI];
    }

    if (debug&meshRefinement::OBJINTERSECTIONS)
    {
        dumpMove
        (
            mesh.time().path()
          / "faceDisp_" + name(iter) + ".obj",
            pp.faceCentres(),
            pp.faceCentres() + faceDisp
        );
        dumpMove
        (
            mesh.time().path()
          / "faceRotation_" + name(iter) + ".obj",
            pp.faceCentres(),
            pp.faceCentres() + faceRotation
        );
    }
}


// Collect (possibly remote) per point data of all surrounding faces
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// - faceSurfaceNormal
// - faceDisp
// - faceCentres&faceNormal
void Foam::autoSnapDriver::calcNearestFacePointProperties
(
    const label iter,
    const indirectPrimitivePatch& pp,

    const vectorField& faceDisp,
    const vectorField& faceSurfaceNormal,
    const labelList& faceSurfaceGlobalRegion,

    List<List<point> >& pointFaceSurfNormals,
    List<List<point> >& pointFaceDisp,
    List<List<point> >& pointFaceCentres,
    List<labelList>&    pointFacePatchID
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    // For now just get all surrounding face data. Expensive - should just
    // store and sync data on coupled points only
    // (see e.g PatchToolsNormals.C)

    pointFaceSurfNormals.setSize(pp.nPoints());
    pointFaceDisp.setSize(pp.nPoints());
    pointFaceCentres.setSize(pp.nPoints());
    pointFacePatchID.setSize(pp.nPoints());

    // Fill local data
    forAll(pp.pointFaces(), pointI)
    {
        const labelList& pFaces = pp.pointFaces()[pointI];
        List<point>& pNormals = pointFaceSurfNormals[pointI];
        pNormals.setSize(pFaces.size());
        List<point>& pDisp = pointFaceDisp[pointI];
        pDisp.setSize(pFaces.size());
        List<point>& pFc = pointFaceCentres[pointI];
        pFc.setSize(pFaces.size());
        labelList& pFid = pointFacePatchID[pointI];
        pFid.setSize(pFaces.size());

        forAll(pFaces, i)
        {
            label faceI = pFaces[i];
            pNormals[i] = faceSurfaceNormal[faceI];
            pDisp[i] = faceDisp[faceI];
            pFc[i] = pp.faceCentres()[faceI];
            //label meshFaceI = pp.addressing()[faceI];
            //pFid[i] = mesh.boundaryMesh().whichPatch(meshFaceI);
            pFid[i] = globalToMasterPatch_[faceSurfaceGlobalRegion[faceI]];
        }
    }


    // Collect additionally 'normal' boundary faces for boundaryPoints of pp
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // points on the boundary of pp should pick up non-pp normals
    // as well for the feature-reconstruction to behave correctly.
    // (the movement is already constrained outside correctly so it
    //  is only that the unconstrained attraction vector is calculated
    //  correctly)
    {
        const polyBoundaryMesh& pbm = mesh.boundaryMesh();
        labelList patchID(pbm.patchID());

        // Unmark all non-coupled boundary faces
        forAll(pbm, patchI)
        {
            const polyPatch& pp = pbm[patchI];

            if (pp.coupled() || isA<emptyPolyPatch>(pp))
            {
                forAll(pp, i)
                {
                    label meshFaceI = pp.start()+i;
                    patchID[meshFaceI-mesh.nInternalFaces()] = -1;
                }
            }
        }

        // Remove any meshed faces
        forAll(pp.addressing(), i)
        {
            label meshFaceI = pp.addressing()[i];
            patchID[meshFaceI-mesh.nInternalFaces()] = -1;
        }

        // See if pp point uses any non-meshed boundary faces

        const labelList& boundaryPoints = pp.boundaryPoints();
        forAll(boundaryPoints, i)
        {
            label pointI = boundaryPoints[i];
            label meshPointI = pp.meshPoints()[pointI];
            const point& pt = mesh.points()[meshPointI];
            const labelList& pFaces = mesh.pointFaces()[meshPointI];

            List<point>& pNormals = pointFaceSurfNormals[pointI];
            List<point>& pDisp = pointFaceDisp[pointI];
            List<point>& pFc = pointFaceCentres[pointI];
            labelList& pFid = pointFacePatchID[pointI];

            forAll(pFaces, i)
            {
                label meshFaceI = pFaces[i];
                if (!mesh.isInternalFace(meshFaceI))
                {
                    label patchI = patchID[meshFaceI-mesh.nInternalFaces()];

                    if (patchI != -1)
                    {
                        vector fn = mesh.faceAreas()[meshFaceI];
                        pNormals.append(fn/mag(fn));
                        pDisp.append(mesh.faceCentres()[meshFaceI]-pt);
                        pFc.append(mesh.faceCentres()[meshFaceI]);
                        pFid.append(patchI);
                    }
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        pointFaceSurfNormals,
        listPlusEqOp<point>(),
        List<point>(),
        listTransform()
    );
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        pointFaceDisp,
        listPlusEqOp<point>(),
        List<point>(),
        listTransform()
    );
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        pointFaceCentres,
        listPlusEqOp<point>(),
        List<point>(),
        listTransform()
    );
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        pointFacePatchID,
        listPlusEqOp<label>(),
        List<label>()
    );
}


// Gets passed in offset to nearest point on feature edge. Calculates
// if the point has a different number of faces on either side of the feature
// and if so attracts the point to that non-dominant plane.
void Foam::autoSnapDriver::correctAttraction
(
    const DynamicList<point>& surfacePoints,
    const DynamicList<label>& surfaceCounts,
    const point& edgePt,
    const vector& edgeNormal,       // normalised normal
    const point& pt,

    vector& edgeOffset              // offset from pt to point on edge
) const
{
    // Tangential component along edge
    scalar tang = ((pt-edgePt)&edgeNormal);

    labelList order;
    Foam::sortedOrder(surfaceCounts, order);

    if (order[0] < order[1])
    {
        // There is a non-dominant plane. Use the point on the plane to
        // attract to.
        vector attractD = surfacePoints[order[0]]-edgePt;
        // Tangential component along edge
        scalar tang2 = (attractD&edgeNormal);
        // Normal component
        attractD -= tang2*edgeNormal;
        // Calculate fraction of normal distances
        scalar magAttractD = mag(attractD);
        scalar fraction = magAttractD/(magAttractD+mag(edgeOffset));

        point linePt =
            edgePt
          + ((1.0-fraction)*tang2 + fraction*tang)*edgeNormal;
        edgeOffset = linePt-pt;
    }
}


Foam::pointIndexHit Foam::autoSnapDriver::findMultiPatchPoint
(
    const point& pt,
    const labelList& patchIDs,
    const List<point>& faceCentres
) const
{
    // Determine if multiple patchIDs
    if (patchIDs.size())
    {
        label patch0 = patchIDs[0];

        for (label i = 1; i < patchIDs.size(); i++)
        {
            if (patchIDs[i] != patch0)
            {
                return pointIndexHit(true, pt, labelMax);
            }
        }
    }
    return pointIndexHit(false, vector::zero, labelMax);
}


void Foam::autoSnapDriver::binFeatureFace
(
    const label iter,
    const scalar featureCos,

    const indirectPrimitivePatch& pp,
    const scalar snapDist,

    const point& fc,
    const vector& faceSurfaceNormal,
    const vector& faceDisp,

    DynamicList<point>& surfacePoints,
    DynamicList<vector>& surfaceNormals,
    DynamicList<label>& surfaceCounts
) const
{
    // What to do with very far attraction? For now just ignore the face
    if (magSqr(faceDisp) < sqr(snapDist) && mag(faceSurfaceNormal) > VSMALL)
    {
        const point pt = fc + faceDisp;

        bool same = false;
        forAll(surfaceNormals, j)
        {
            scalar cosAngle = (faceSurfaceNormal&surfaceNormals[j]);

            if
            (
                (cosAngle >= featureCos)
             || (cosAngle < (-1+0.001)) // triangle baffles
            )
            {
                same = true;
                surfaceCounts[j]++;
                break;
            }
        }

        if (!same)
        {
            // Now check if the planes go through the same edge or point

            if (surfacePoints.size() <= 1)
            {
                surfacePoints.append(pt);
                surfaceNormals.append(faceSurfaceNormal);
                surfaceCounts.append(1);
            }
            else if (surfacePoints.size() == 2)
            {
                plane pl0(surfacePoints[0], surfaceNormals[0]);
                plane pl1(surfacePoints[1], surfaceNormals[1]);
                plane::ray r(pl0.planeIntersect(pl1));
                vector featureNormal = r.dir() / mag(r.dir());

                if (mag(faceSurfaceNormal&featureNormal) >= 0.001)
                {
                    // Definitely makes a feature point
                    surfacePoints.append(pt);
                    surfaceNormals.append(faceSurfaceNormal);
                    surfaceCounts.append(1);
                }
            }
            else if (surfacePoints.size() == 3)
            {
                // Have already feature point. See if this new plane is the
                // same point or not.
                plane pl0(surfacePoints[0], surfaceNormals[0]);
                plane pl1(surfacePoints[1], surfaceNormals[1]);
                plane pl2(surfacePoints[2], surfaceNormals[2]);
                point p012(pl0.planePlaneIntersect(pl1, pl2));

                plane::ray r(pl0.planeIntersect(pl1));
                vector featureNormal = r.dir() / mag(r.dir());
                if (mag(faceSurfaceNormal&featureNormal) >= 0.001)
                {
                    plane pl3(pt, faceSurfaceNormal);
                    point p013(pl0.planePlaneIntersect(pl1, pl3));

                    if (mag(p012-p013) > 0.0001)    //TBD
                    {
                        // Different feature point
                        surfacePoints.append(pt);
                        surfaceNormals.append(faceSurfaceNormal);
                        surfaceCounts.append(1);
                    }
                }
            }
        }
    }
}


// Check the faces surrounding a point. Bin them according to normal.
void Foam::autoSnapDriver::binFeatureFaces
(
    const label iter,
    const scalar featureCos,

    const indirectPrimitivePatch& pp,
    const scalarField& snapDist,

    const label pointI,

    const List<List<point> >& pointFaceSurfNormals,
    const List<List<point> >& pointFaceDisp,
    const List<List<point> >& pointFaceCentres,

    DynamicList<point>& surfacePoints,
    DynamicList<vector>& surfaceNormals,
    DynamicList<label>& surfaceCounts
) const
{
    const List<point>& pfSurfNormals = pointFaceSurfNormals[pointI];
    const List<point>& pfDisp = pointFaceDisp[pointI];
    const List<point>& pfCentres = pointFaceCentres[pointI];

    // Collect all different directions
    forAll(pfSurfNormals, i)
    {
        binFeatureFace
        (
            iter,
            featureCos,

            pp,
            snapDist[pointI],

            pfCentres[i],
            pfSurfNormals[i],
            pfDisp[i],

            surfacePoints,
            surfaceNormals,
            surfaceCounts
        );
    }
}


void Foam::autoSnapDriver::featureAttractionUsingReconstruction
(
    const label iter,
    const scalar featureCos,

    const indirectPrimitivePatch& pp,
    const scalarField& snapDist,
    const label pointI,

    const List<List<point> >& pointFaceSurfNormals,
    const List<List<point> >& pointFaceDisp,
    const List<List<point> >& pointFaceCentres,
    const labelListList& pointFacePatchID,

    vector& patchAttraction,
    pointConstraint& patchConstraint
) const
{
    patchAttraction = vector::zero;
    patchConstraint = pointConstraint();

    // Collect all different directions
    DynamicList<point> surfacePoints(4);
    DynamicList<vector> surfaceNormals(4);
    DynamicList<label> surfaceCounts(4);

    binFeatureFaces
    (
        iter,
        featureCos,

        pp,
        snapDist,
        pointI,

        pointFaceSurfNormals,
        pointFaceDisp,
        pointFaceCentres,

        surfacePoints,
        surfaceNormals,
        surfaceCounts
    );

    const point& pt = pp.localPoints()[pointI];

    // Check the number of directions
    if (surfaceNormals.size() == 1)
    {
        // Normal distance to plane
        vector d =
            ((surfacePoints[0]-pt) & surfaceNormals[0])
           *surfaceNormals[0];

        // Trim to snap distance
        if (magSqr(d) > sqr(snapDist[pointI]))
        {
            d *= Foam::sqrt(sqr(snapDist[pointI])/magSqr(d));
        }

        patchAttraction = d;

        // Store constraints
        patchConstraint.applyConstraint(surfaceNormals[0]);
    }
    else if (surfaceNormals.size() == 2)
    {
        plane pl0(surfacePoints[0], surfaceNormals[0]);
        plane pl1(surfacePoints[1], surfaceNormals[1]);
        plane::ray r(pl0.planeIntersect(pl1));
        vector n = r.dir() / mag(r.dir());

        // Get nearest point on infinite ray
        vector d = r.refPoint()-pt;
        d -= (d&n)*n;

        //- This does not help much but distorts a perfectly aligned mesh
        //  so disabled for now.
        //// Correct for attraction to non-dominant face
        //correctAttraction
        //(
        //    surfacePoints,
        //    surfaceCounts,
        //    r.refPoint(),
        //    n,                  // normalised normal
        //    pt,
        //
        //    d                   // perpendicular offset vector
        //);

        // Trim to snap distance
        if (magSqr(d) > sqr(snapDist[pointI]))
        {
            d *= Foam::sqrt(sqr(snapDist[pointI])/magSqr(d));
        }

        patchAttraction = d;

        // Store constraints
        patchConstraint.applyConstraint(surfaceNormals[0]);
        patchConstraint.applyConstraint(surfaceNormals[1]);
    }
    else if (surfaceNormals.size() == 3)
    {
        // Calculate point from the faces.
        plane pl0(surfacePoints[0], surfaceNormals[0]);
        plane pl1(surfacePoints[1], surfaceNormals[1]);
        plane pl2(surfacePoints[2], surfaceNormals[2]);
        point cornerPt(pl0.planePlaneIntersect(pl1, pl2));
        vector d = cornerPt - pt;

        // Trim to snap distance
        if (magSqr(d) > sqr(snapDist[pointI]))
        {
            d *= Foam::sqrt(sqr(snapDist[pointI])/magSqr(d));
        }

        patchAttraction = d;

        // Store constraints
        patchConstraint.applyConstraint(surfaceNormals[0]);
        patchConstraint.applyConstraint(surfaceNormals[1]);
        patchConstraint.applyConstraint(surfaceNormals[2]);

        //Pout<< "# Feature point " << pt << nl;
        //meshTools::writeOBJ(Pout, pt);
        //meshTools::writeOBJ(Pout, surfacePoints[0]);
        //meshTools::writeOBJ(Pout, surfacePoints[1]);
        //meshTools::writeOBJ(Pout, surfacePoints[2]);
        //Pout<< "l 1 2" << nl
        //    << "l 1 3" << nl
        //    << "l 1 4" << nl;
    }
}


// Special version that calculates attraction in one go
void Foam::autoSnapDriver::featureAttractionUsingReconstruction
(
    const label iter,
    const scalar featureCos,

    const indirectPrimitivePatch& pp,
    const scalarField& snapDist,

    const List<List<point> >& pointFaceSurfNormals,
    const List<List<point> >& pointFaceDisp,
    const List<List<point> >& pointFaceCentres,
    const labelListList& pointFacePatchID,

    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    autoPtr<OFstream> feStr;
    label feVertI = 0;
    autoPtr<OFstream> fpStr;
    label fpVertI = 0;

    if (debug&meshRefinement::OBJINTERSECTIONS)
    {
        feStr.reset
        (
            new OFstream
            (
                meshRefiner_.mesh().time().path()
              / "implicitFeatureEdge_" + name(iter) + ".obj"
            )
        );
        Pout<< "Dumping implicit feature-edge direction to "
            << feStr().name() << endl;

        fpStr.reset
        (
            new OFstream
            (
                meshRefiner_.mesh().time().path()
              / "implicitFeaturePoint_" + name(iter) + ".obj"
            )
        );
        Pout<< "Dumping implicit feature-point direction to "
            << fpStr().name() << endl;
    }


    forAll(pp.localPoints(), pointI)
    {
        vector attraction = vector::zero;
        pointConstraint constraint;

        featureAttractionUsingReconstruction
        (
            iter,
            featureCos,

            pp,
            snapDist,
            pointI,

            pointFaceSurfNormals,
            pointFaceDisp,
            pointFaceCentres,
            pointFacePatchID,

            attraction,
            constraint
        );

        if
        (
            (constraint.first() > patchConstraints[pointI].first())
         || (magSqr(attraction) < magSqr(patchAttraction[pointI]))
        )
        {
            patchAttraction[pointI] = attraction;
            patchConstraints[pointI] = constraint;

            const point& pt = pp.localPoints()[pointI];

            if (patchConstraints[pointI].first() == 2 && feStr.valid())
            {
                meshTools::writeOBJ(feStr(), pt);
                feVertI++;
                meshTools::writeOBJ(feStr(), pt+patchAttraction[pointI]);
                feVertI++;
                feStr() << "l " << feVertI-1 << ' ' << feVertI << nl;
            }
            else if (patchConstraints[pointI].first() == 3 && fpStr.valid())
            {
                meshTools::writeOBJ(fpStr(), pt);
                fpVertI++;
                meshTools::writeOBJ(fpStr(), pt+patchAttraction[pointI]);
                fpVertI++;
                fpStr() << "l " << fpVertI-1 << ' ' << fpVertI << nl;
            }
        }
    }
}


void Foam::autoSnapDriver::stringFeatureEdges
(
    const label iter,
    const scalar featureCos,

    const indirectPrimitivePatch& pp,
    const scalarField& snapDist,

    const vectorField& rawPatchAttraction,
    const List<pointConstraint>& rawPatchConstraints,

    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    // Snap edges to feature edges
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Walk existing edges and snap remaining ones (that are marked as
    // feature edges in rawPatchConstraints)

    // What this does is fill in any faces where not all points
    // on the face are being attracted:
    /*
           +
          / \
         /   \
      ---+    +---
         \   /
          \ /
           +
    */
    // so the top and bottom will never get attracted since the nearest
    // back from the feature edge will always be one of the left or right
    // points since the face is diamond like. So here we walk the feature edges
    // and add any non-attracted points.


    while (true)
    {
        label nChanged = 0;

        const labelListList& pointEdges = pp.pointEdges();
        forAll(pointEdges, pointI)
        {
            if (patchConstraints[pointI].first() == 2)
            {
                const point& pt = pp.localPoints()[pointI];
                const labelList& pEdges = pointEdges[pointI];
                const vector& featVec = patchConstraints[pointI].second();

                // Detect whether there are edges in both directions.
                // (direction along the feature edge that is)
                bool hasPos = false;
                bool hasNeg = false;

                forAll(pEdges, pEdgeI)
                {
                    const edge& e = pp.edges()[pEdges[pEdgeI]];
                    label nbrPointI = e.otherVertex(pointI);

                    if (patchConstraints[nbrPointI].first() > 1)
                    {
                        const point& nbrPt = pp.localPoints()[nbrPointI];
                        const point featPt =
                            nbrPt + patchAttraction[nbrPointI];
                        const scalar cosAngle = (featVec & (featPt-pt));

                        if (cosAngle > 0)
                        {
                            hasPos = true;
                        }
                        else
                        {
                            hasNeg = true;
                        }
                    }
                }

                if (!hasPos || !hasNeg)
                {
                    //Pout<< "**Detected feature string end at  "
                    //    << pp.localPoints()[pointI] << endl;

                    // No string. Assign best choice on either side
                    label bestPosPointI = -1;
                    scalar minPosDistSqr = GREAT;
                    label bestNegPointI = -1;
                    scalar minNegDistSqr = GREAT;

                    forAll(pEdges, pEdgeI)
                    {
                        const edge& e = pp.edges()[pEdges[pEdgeI]];
                        label nbrPointI = e.otherVertex(pointI);

                        if
                        (
                            patchConstraints[nbrPointI].first() <= 1
                         && rawPatchConstraints[nbrPointI].first() > 1
                        )
                        {
                            const vector& nbrFeatVec =
                                rawPatchConstraints[pointI].second();

                            if (mag(featVec&nbrFeatVec) > featureCos)
                            {
                                // nbrPointI attracted to sameish feature
                                // Note: also check on position.

                                scalar d2 = magSqr
                                (
                                    rawPatchAttraction[nbrPointI]
                                );

                                const point featPt =
                                    pp.localPoints()[nbrPointI]
                                  + rawPatchAttraction[nbrPointI];
                                const scalar cosAngle =
                                    (featVec & (featPt-pt));

                                if (cosAngle > 0)
                                {
                                    if (!hasPos && d2 < minPosDistSqr)
                                    {
                                        minPosDistSqr = d2;
                                        bestPosPointI = nbrPointI;
                                    }
                                }
                                else
                                {
                                    if (!hasNeg && d2 < minNegDistSqr)
                                    {
                                        minNegDistSqr = d2;
                                        bestNegPointI = nbrPointI;
                                    }
                                }
                            }
                        }
                    }

                    if (bestPosPointI != -1)
                    {
                        // Use reconstructed-feature attraction. Use only
                        // part of it since not sure...
                        //const point& bestPt =
                        //    pp.localPoints()[bestPosPointI];
                        //Pout<< "**Overriding point " << bestPt
                        //    << " on reconstructed feature edge at "
                        //    << rawPatchAttraction[bestPosPointI]+bestPt
                        //    << " to attracted-to-feature-edge." << endl;
                        patchAttraction[bestPosPointI] =
                            0.5*rawPatchAttraction[bestPosPointI];
                        patchConstraints[bestPosPointI] =
                            rawPatchConstraints[bestPosPointI];

                        nChanged++;
                    }
                    if (bestNegPointI != -1)
                    {
                        // Use reconstructed-feature attraction. Use only
                        // part of it since not sure...
                        //const point& bestPt =
                        //    pp.localPoints()[bestNegPointI];
                        //Pout<< "**Overriding point " << bestPt
                        //    << " on reconstructed feature edge at "
                        //    << rawPatchAttraction[bestNegPointI]+bestPt
                        //    << " to attracted-to-feature-edge." << endl;
                        patchAttraction[bestNegPointI] =
                            0.5*rawPatchAttraction[bestNegPointI];
                        patchConstraints[bestNegPointI] =
                            rawPatchConstraints[bestNegPointI];

                        nChanged++;
                    }
                }
            }
        }


        reduce(nChanged, sumOp<label>());
        Info<< "Stringing feature edges : changed " << nChanged << " points"
            << endl;
        if (nChanged == 0)
        {
            break;
        }
    }
}


Foam::pointIndexHit Foam::autoSnapDriver::findNearFeatureEdge
(
    const indirectPrimitivePatch& pp,
    const scalarField& snapDist,
    const label pointI,
    const point& estimatedPt,

    label& featI,
    List<List<DynamicList<point> > >& edgeAttractors,
    List<List<DynamicList<pointConstraint> > >& edgeConstraints,
    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    const refinementFeatures& features = meshRefiner_.features();

    labelList nearEdgeFeat;
    List<pointIndexHit> nearEdgeInfo;
    features.findNearestEdge
    (
        pointField(1, estimatedPt),
        scalarField(1, sqr(snapDist[pointI])),
        nearEdgeFeat,
        nearEdgeInfo
    );

    const pointIndexHit& nearInfo = nearEdgeInfo[0];
    featI = nearEdgeFeat[0];

    if (nearInfo.hit())
    {
        // So we have a point on the feature edge. Use this
        // instead of our estimate from planes.
        edgeAttractors[featI][nearInfo.index()].append
        (
            nearInfo.hitPoint()
        );
        pointConstraint c;
        const edge e = features[featI].edges()[nearInfo.index()];
        vector eVec = e.vec(features[featI].points());
        eVec /= mag(eVec)+VSMALL;
        c.first() = 2;
        c.second() = eVec;
        edgeConstraints[featI][nearInfo.index()].append(c);

        // Store for later use
        patchAttraction[pointI] =
            nearInfo.hitPoint()-pp.localPoints()[pointI];
        patchConstraints[pointI] = c;
    }
    return nearInfo;
}
Foam::labelPair Foam::autoSnapDriver::findNearFeaturePoint
(
    const indirectPrimitivePatch& pp,
    const scalarField& snapDist,
    const label pointI,
    const point& estimatedPt,

    // Feature-point to pp point
    List<labelList>& pointAttractor,
    List<List<pointConstraint> >& pointConstraints,
    // Feature-edge to pp point
    List<List<DynamicList<point> > >& edgeAttractors,
    List<List<DynamicList<pointConstraint> > >& edgeConstraints,
    // pp point to nearest feature
    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    const refinementFeatures& features = meshRefiner_.features();

    labelList nearFeat;
    labelList nearIndex;
    features.findNearestPoint
    (
        pointField(1, estimatedPt),
        scalarField(1, sqr(snapDist[pointI])),
        nearFeat,
        nearIndex
    );

    label featI = nearFeat[0];
    label featPointI = -1;

    if (featI != -1)
    {
        const point& pt = pp.localPoints()[pointI];

        const treeDataPoint& shapes =
            features.pointTrees()[featI].shapes();
        featPointI = shapes.pointLabels()[nearIndex[0]];
        const point& featPt = shapes.points()[featPointI];
        scalar distSqr = magSqr(featPt-pt);

        // Check if already attracted
        label oldPointI = pointAttractor[featI][featPointI];

        if (oldPointI != -1)
        {
            // Check distance
            if (distSqr >= magSqr(featPt-pp.localPoints()[oldPointI]))
            {
                // oldPointI nearest. Keep.
                featI = -1;
                featPointI = -1;
            }
            else
            {
                // Current pointI nearer.
                pointAttractor[featI][featPointI] = pointI;
                pointConstraints[featI][featPointI].first() = 3;
                pointConstraints[featI][featPointI].second() = vector::zero;

                // Store for later use
                patchAttraction[pointI] = featPt-pt;
                patchConstraints[pointI] =
                    pointConstraints[featI][featPointI];

                // Reset oldPointI to nearest on feature edge
                patchAttraction[oldPointI] = vector::zero;
                patchConstraints[oldPointI] = pointConstraint();

                label edgeFeatI;
                findNearFeatureEdge
                (
                    pp,
                    snapDist,
                    oldPointI,
                    pp.localPoints()[oldPointI],

                    edgeFeatI,
                    edgeAttractors,
                    edgeConstraints,
                    patchAttraction,
                    patchConstraints
                );
            }
        }
        else
        {
            // Current pointI nearer.
            pointAttractor[featI][featPointI] = pointI;
            pointConstraints[featI][featPointI].first() = 3;
            pointConstraints[featI][featPointI].second() = vector::zero;

            // Store for later use
            patchAttraction[pointI] = featPt-pt;
            patchConstraints[pointI] = pointConstraints[featI][featPointI];
        }
    }

    return labelPair(featI, featPointI);
}


// Determines for every pp point - that is on multiple faces that form
// a feature - the nearest feature edge/point.
void Foam::autoSnapDriver::determineFeatures
(
    const label iter,
    const scalar featureCos,
    const bool multiRegionFeatureSnap,

    const indirectPrimitivePatch& pp,
    const scalarField& snapDist,

    const List<List<point> >& pointFaceSurfNormals,
    const List<List<point> >& pointFaceDisp,
    const List<List<point> >& pointFaceCentres,
    const labelListList& pointFacePatchID,

    // Feature-point to pp point
    List<labelList>& pointAttractor,
    List<List<pointConstraint> >& pointConstraints,
    // Feature-edge to pp point
    List<List<DynamicList<point> > >& edgeAttractors,
    List<List<DynamicList<pointConstraint> > >& edgeConstraints,
    // pp point to nearest feature
    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    autoPtr<OFstream> featureEdgeStr;
    label featureEdgeVertI = 0;
    autoPtr<OFstream> missedEdgeStr;
    label missedVertI = 0;
    autoPtr<OFstream> featurePointStr;
    label featurePointVertI = 0;

    if (debug&meshRefinement::OBJINTERSECTIONS)
    {
        featureEdgeStr.reset
        (
            new OFstream
            (
                meshRefiner_.mesh().time().path()
              / "featureEdge_" + name(iter) + ".obj"
            )
        );
        Pout<< "Dumping feature-edge sampling to "
            << featureEdgeStr().name() << endl;

        missedEdgeStr.reset
        (
            new OFstream
            (
                meshRefiner_.mesh().time().path()
              / "missedFeatureEdge_" + name(iter) + ".obj"
            )
        );
        Pout<< "Dumping feature-edges that are too far away to "
            << missedEdgeStr().name() << endl;

        featurePointStr.reset
        (
            new OFstream
            (
                meshRefiner_.mesh().time().path()
              / "featurePoint_" + name(iter) + ".obj"
            )
        );
        Pout<< "Dumping feature-point sampling to "
            << featurePointStr().name() << endl;
    }

    const refinementFeatures& features = meshRefiner_.features();

    forAll(pp.localPoints(), pointI)
    {
        const point& pt = pp.localPoints()[pointI];

        vector attraction = vector::zero;
        pointConstraint constraint;

        featureAttractionUsingReconstruction
        (
            iter,
            featureCos,

            pp,
            snapDist,
            pointI,

            pointFaceSurfNormals,
            pointFaceDisp,
            pointFaceCentres,
            pointFacePatchID,

            attraction,
            constraint
        );

        if
        (
            (constraint.first() > patchConstraints[pointI].first())
         || (magSqr(attraction) < magSqr(patchAttraction[pointI]))
        )
        {
            patchAttraction[pointI] = attraction;
            patchConstraints[pointI] = constraint;

            // Check the number of directions
            if (patchConstraints[pointI].first() == 1)
            {
                // Flat surface. Check for different patchIDs
                if (multiRegionFeatureSnap)
                {
                    pointIndexHit multiPatchPt
                    (
                        findMultiPatchPoint
                        (
                            pt,
                            pointFacePatchID[pointI],
                            pointFaceCentres[pointI]
                        )
                    );
                    if (multiPatchPt.hit())
                    {
                        // Behave like when having two surface normals so
                        // attract to nearest feature edge (with a guess for
                        // the multipatch point as starting point)
                        label featI = -1;
                        pointIndexHit nearInfo = findNearFeatureEdge
                        (
                            pp,
                            snapDist,
                            pointI,
                            multiPatchPt.hitPoint(),        //estimatedPt

                            featI,
                            edgeAttractors,
                            edgeConstraints,

                            patchAttraction,
                            patchConstraints
                        );

                        if (nearInfo.hit())
                        {
                            // Dump
                            if (featureEdgeStr.valid())
                            {
                                meshTools::writeOBJ(featureEdgeStr(), pt);
                                featureEdgeVertI++;
                                meshTools::writeOBJ
                                (
                                    featureEdgeStr(),
                                    nearInfo.hitPoint()
                                );
                                featureEdgeVertI++;
                                featureEdgeStr()
                                    << "l " << featureEdgeVertI-1 << ' '
                                    << featureEdgeVertI << nl;
                            }
                        }
                        else
                        {
                            if (missedEdgeStr.valid())
                            {
                                meshTools::writeOBJ(missedEdgeStr(), pt);
                                missedVertI++;
                                meshTools::writeOBJ
                                (
                                    missedEdgeStr(),
                                    nearInfo.missPoint()
                                );
                                missedVertI++;
                                missedEdgeStr()
                                    << "l " << missedVertI-1 << ' '
                                    << missedVertI << nl;
                            }
                        }
                    }
                }
            }
            else if (patchConstraints[pointI].first() == 2)
            {
                // Mark point on the nearest feature edge. Note that we
                // only search within the surrounding since the plane
                // reconstruction might find a feature where there isn't one.
                const point estimatedPt(pt + patchAttraction[pointI]);

                // Determine nearest point on feature edge. Store constraint
                // (calculated from feature edge, alternative would be to
                //  use constraint calculated from both surfaceNormals)
                label featI = -1;
                pointIndexHit nearInfo = findNearFeatureEdge
                (
                    pp,
                    snapDist,
                    pointI,
                    estimatedPt,

                    featI,
                    edgeAttractors,
                    edgeConstraints,

                    patchAttraction,
                    patchConstraints
                );

                if (nearInfo.hit())
                {
                    // Dump
                    if (featureEdgeStr.valid())
                    {
                        meshTools::writeOBJ(featureEdgeStr(), pt);
                        featureEdgeVertI++;
                        meshTools::writeOBJ
                        (
                            featureEdgeStr(),
                            nearInfo.hitPoint()
                        );
                        featureEdgeVertI++;
                        featureEdgeStr()
                            << "l " << featureEdgeVertI-1 << ' '
                            << featureEdgeVertI << nl;
                    }
                }
                else
                {
                    if (missedEdgeStr.valid())
                    {
                        meshTools::writeOBJ(missedEdgeStr(), pt);
                        missedVertI++;
                        meshTools::writeOBJ
                        (
                            missedEdgeStr(),
                            nearInfo.missPoint()
                        );
                        missedVertI++;
                        missedEdgeStr()
                            << "l " << missedVertI-1 << ' '
                            << missedVertI << nl;
                    }
                }
            }
            else if (patchConstraints[pointI].first() == 3)
            {
                // Mark point on the nearest feature point.
                const point estimatedPt(pt + patchAttraction[pointI]);

                labelPair nearInfo = findNearFeaturePoint
                (
                    pp,
                    snapDist,
                    pointI,
                    estimatedPt,

                    // Feature-point to pp point
                    pointAttractor,
                    pointConstraints,
                    // Feature-edge to pp point
                    edgeAttractors,
                    edgeConstraints,
                    // pp point to nearest feature
                    patchAttraction,
                    patchConstraints
                );

                if (nearInfo.first() != -1)
                {
                    // Dump
                    if (featurePointStr.valid())
                    {
                        const treeDataPoint& shapes =
                            features.pointTrees()[nearInfo.first()].shapes();
                        const point& featPt =
                            shapes.points()[nearInfo.second()];

                        meshTools::writeOBJ(featurePointStr(), pt);
                        featurePointVertI++;
                        meshTools::writeOBJ(featurePointStr(), featPt);
                        featurePointVertI++;
                        featurePointStr()
                            << "l " << featurePointVertI-1 << ' '
                            << featurePointVertI << nl;
                    }
                }
            }
        }
    }
}


void Foam::autoSnapDriver::featureAttractionUsingFeatureEdges
(
    const label iter,
    const scalar featureCos,
    const bool multiRegionFeatureSnap,

    const indirectPrimitivePatch& pp,
    const scalarField& snapDist,

    const List<List<point> >& pointFaceSurfNormals,
    const List<List<point> >& pointFaceDisp,
    const List<List<point> >& pointFaceCentres,
    const labelListList& pointFacePatchID,

    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    const refinementFeatures& features = meshRefiner_.features();

    // Collect ordered attractions on feature edges
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Per feature, per feature-edge a list of attraction points and their
    // originating vertex.
    List<List<DynamicList<point> > > edgeAttractors(features.size());
    List<List<DynamicList<pointConstraint> > > edgeConstraints
    (
        features.size()
    );
    forAll(features, featI)
    {
        label nFeatEdges = features[featI].edges().size();
        edgeAttractors[featI].setSize(nFeatEdges);
        edgeConstraints[featI].setSize(nFeatEdges);
    }

    // Per feature, per feature-point the pp point that is attracted to it.
    // This list is only used to subset the feature-points that are actually
    // used.
    List<labelList> pointAttractor(features.size());
    List<List<pointConstraint> > pointConstraints(features.size());
    forAll(features, featI)
    {
        label nFeatPoints = features[featI].points().size();
        pointAttractor[featI].setSize(nFeatPoints, -1);
        pointConstraints[featI].setSize(nFeatPoints);
    }

    // Reverse: from pp point to nearest feature
    vectorField rawPatchAttraction(pp.nPoints(), vector::zero);
    List<pointConstraint> rawPatchConstraints(pp.nPoints());

    determineFeatures
    (
        iter,
        featureCos,
        multiRegionFeatureSnap,

        pp,
        snapDist,

        pointFaceSurfNormals,
        pointFaceDisp,
        pointFaceCentres,
        pointFacePatchID,

        // Feature-point to pp point
        pointAttractor,
        pointConstraints,
        // Feature-edge to pp point
        edgeAttractors,
        edgeConstraints,
        // pp point to nearest feature
        rawPatchAttraction,
        rawPatchConstraints
    );



    // Baffle handling
    // ~~~~~~~~~~~~~~~
    // Override pointAttractor, edgeAttractor, rawPatchAttration etc. to
    // implement 'baffle' handling.
    // Baffle: the mesh pp point originates from a loose standing
    // baffle.
    // Sampling the surface with the surrounding face-centres only picks up
    // a single triangle normal so above determineFeatures will not have
    // detected anything. So explicitly pick up feature edges on the pp
    // (after duplicating points & smoothing so will already have been
    // expanded) and match these to the features.
    {
        const fvMesh& mesh = meshRefiner_.mesh();

        // Calculate edge-faces
        List<List<point> > edgeFaceNormals(pp.nEdges());

        // Fill local data
        forAll(pp.edgeFaces(), edgeI)
        {
            const labelList& eFaces = pp.edgeFaces()[edgeI];
            List<point>& eFc = edgeFaceNormals[edgeI];
            eFc.setSize(eFaces.size());
            forAll(eFaces, i)
            {
                label faceI = eFaces[i];
                eFc[i] = pp.faceNormals()[faceI];
            }
        }

        {
            // Precalculate mesh edges for pp.edges.
            const labelList meshEdges
            (
                pp.meshEdges(mesh.edges(), mesh.pointEdges())
            );
            syncTools::syncEdgeList
            (
                mesh,
                meshEdges,
                edgeFaceNormals,
                listPlusEqOp<point>(),
                List<point>(),
                listTransform()
            );
        }

        // Detect baffle edges. Assume initial mesh will have 0,90 or 180
        // (baffle) degree angles so smoothing should make 0,90
        // to be less than 90.
        const scalar baffleFeatureCos = Foam::cos(degToRad(91));


        autoPtr<OFstream> baffleEdgeStr;
        label baffleEdgeVertI = 0;
        if (debug&meshRefinement::OBJINTERSECTIONS)
        {
            baffleEdgeStr.reset
            (
                new OFstream
                (
                    meshRefiner_.mesh().time().path()
                  / "baffleEdge_" + name(iter) + ".obj"
                )
            );
            Info<< nl << "Dumping baffle-edges to "
                << baffleEdgeStr().name() << endl;
        }


        // Is edge on baffle
        PackedBoolList isBaffleEdge(pp.nEdges());
        // Is point on
        //  0 : baffle-edge (0)
        //  1 : baffle-feature-point (1)
        // -1 : rest
        labelList pointStatus(pp.nPoints(), -1);

        forAll(edgeFaceNormals, edgeI)
        {
            const List<point>& efn = edgeFaceNormals[edgeI];

            if (efn.size() == 2 && (efn[0]&efn[1]) < baffleFeatureCos)
            {
                isBaffleEdge[edgeI] = true;
                const edge& e = pp.edges()[edgeI];
                pointStatus[e[0]] = 0;
                pointStatus[e[1]] = 0;

                if (baffleEdgeStr.valid())
                {
                    const point& p0 = pp.localPoints()[e[0]];
                    const point& p1 = pp.localPoints()[e[1]];
                    meshTools::writeOBJ(baffleEdgeStr(), p0);
                    baffleEdgeVertI++;
                    meshTools::writeOBJ(baffleEdgeStr(), p1);
                    baffleEdgeVertI++;
                    baffleEdgeStr() << "l " << baffleEdgeVertI-1
                        << ' ' << baffleEdgeVertI << nl;
                }
            }
        }

        forAll(pp.pointEdges(), pointI)
        {
            if
            (
                isFeaturePoint
                (
                    featureCos,
                    pp,
                    isBaffleEdge,
                    pointI
                )
            )
            {
                //Pout<< "Detected feature point:" << pp.localPoints()[pointI]
                //    << endl;
                //-TEMPORARILY DISABLED:
                //pointStatus[pointI] = 1;
            }
        }

        forAll(pointStatus, pointI)
        {
            const point& pt = pp.localPoints()[pointI];

            if (pointStatus[pointI] == 0)   // baffle edge
            {
                label featI;
                const pointIndexHit nearInfo = findNearFeatureEdge
                (
                    pp,
                    snapDist,
                    pointI,
                    pt,

                    featI,
                    edgeAttractors,
                    edgeConstraints,
                    rawPatchAttraction,
                    rawPatchConstraints
                );

                if (!nearInfo.hit())
                {
                    //Pout<< "*** Failed to find close edge to point " << pt
                    //    << endl;
                }
            }
            else if (pointStatus[pointI] == 1)   // baffle point
            {
                labelList nearFeat;
                labelList nearIndex;
                features.findNearestPoint
                (
                    pointField(1, pt),
                    scalarField(1, sqr(snapDist[pointI])),
                    nearFeat,
                    nearIndex
                );

                label featI = nearFeat[0];

                if (featI != -1)
                {
                    const treeDataPoint& shapes =
                        features.pointTrees()[featI].shapes();
                    label featPointI = shapes.pointLabels()[nearIndex[0]];
                    const point& featPt = shapes.points()[featPointI];
                    scalar distSqr = magSqr(featPt-pt);

                    // Check if already attracted
                    label oldPointI = pointAttractor[featI][featPointI];

                    if
                    (
                        oldPointI == -1
                     || (
                            distSqr
                          < magSqr(featPt-pp.localPoints()[oldPointI])
                        )
                    )
                    {
                        pointAttractor[featI][featPointI] = pointI;
                        pointConstraints[featI][featPointI].first() = 3;
                        pointConstraints[featI][featPointI].second() =
                            vector::zero;

                        // Store for later use
                        rawPatchAttraction[pointI] = featPt-pt;
                        rawPatchConstraints[pointI] =
                            pointConstraints[featI][featPointI];

                        if (oldPointI != -1)
                        {
                            // The current point is closer so wins. Reset
                            // the old point to attract to nearest edge
                            // instead.
                            label edgeFeatI;
                            findNearFeatureEdge
                            (
                                pp,
                                snapDist,
                                oldPointI,
                                pp.localPoints()[oldPointI],

                                edgeFeatI,
                                edgeAttractors,
                                edgeConstraints,
                                rawPatchAttraction,
                                rawPatchConstraints
                            );
                        }
                    }
                    else
                    {
                        // Make it fall through to check below
                        featI = -1;
                    }
                }

                // Not found a feature point or another point is already
                // closer to that feature
                if (featI == -1)
                {
                    //Pout<< "*** Falling back to finding nearest feature edge"
                    //    << " for baffle-feature-point " << pt
                    //    << endl;

                    label featI;
                    findNearFeatureEdge
                    (
                        pp,
                        snapDist,
                        pointI,
                        pt,                     // starting point

                        featI,
                        edgeAttractors,
                        edgeConstraints,
                        rawPatchAttraction,
                        rawPatchConstraints
                    );
                }
            }
        }
    }


    //
    // Reverse lookup
    // ~~~~~~~~~~~~~~
    //


    // Find nearest mesh point to feature edge
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Reverse lookup : go through all edgeAttractors and find the
    // nearest point on pp

    // Get search domain and extend it a bit
    treeBoundBox bb(pp.localPoints());
    {
        // Random number generator. Bit dodgy since not exactly random ;-)
        Random rndGen(65431);

        // Slightly extended bb. Slightly off-centred just so on symmetric
        // geometry there are less face/edge aligned items.
        bb = bb.extend(rndGen, 1e-4);
        bb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
        bb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
    }

    indexedOctree<treeDataPoint> ppTree
    (
        treeDataPoint(pp.localPoints()),
        bb,                             // overall search domain
        8,                              // maxLevel
        10,                             // leafsize
        3.0                             // duplicity
    );

    // Per mesh point the point on nearest feature edge.
    patchAttraction.setSize(pp.nPoints());
    patchAttraction = vector::zero;
    patchConstraints.setSize(pp.nPoints());
    patchConstraints = pointConstraint();

    forAll(edgeAttractors, featI)
    {
        const List<DynamicList<point> >& edgeAttr = edgeAttractors[featI];
        const List<DynamicList<pointConstraint> >& edgeConstr =
            edgeConstraints[featI];

        forAll(edgeAttr, featEdgeI)
        {
            const DynamicList<point>& attr = edgeAttr[featEdgeI];
            forAll(attr, i)
            {
                // Find nearest pp point
                const point& featPt = attr[i];
                pointIndexHit nearInfo = ppTree.findNearest
                (
                    featPt,
                    sqr(GREAT)
                );

                if (nearInfo.hit())
                {
                    label pointI = nearInfo.index();
                    const point attraction = featPt-pp.localPoints()[pointI];

                    // Check if this point is already being attracted. If so
                    // override it only if nearer.
                    if
                    (
                        patchConstraints[pointI].first() <= 1
                     || magSqr(attraction) < magSqr(patchAttraction[pointI])
                    )
                    {
                        patchAttraction[pointI] = attraction;
                        patchConstraints[pointI] = edgeConstr[featEdgeI][i];
                    }
                }
                else
                {
                    WarningIn
                    (
                        "autoSnapDriver::featureAttractionUsingFeatureEdges"
                        "(..)"
                    )   << "Did not find pp point near " << featPt
                        << endl;
                }
            }
        }
    }


    // Different procs might have different patchAttraction,patchConstraints
    // however these only contain geometric information, no topology
    // so as long as we synchronise after overriding with feature points
    // there is no problem, just possibly a small error.


    // Find nearest mesh point to feature point
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // (overrides attraction to feature edge)
    forAll(pointAttractor, featI)
    {
        const labelList& pointAttr = pointAttractor[featI];
        const List<pointConstraint>& pointConstr = pointConstraints[featI];

        forAll(pointAttr, featPointI)
        {
            if (pointAttr[featPointI] != -1)
            {
                const point& featPt = features[featI].points()
                [
                    featPointI
                ];

                // Find nearest pp point
                pointIndexHit nearInfo = ppTree.findNearest
                (
                    featPt,
                    sqr(GREAT)
                );

                if (nearInfo.hit())
                {
                    label pointI = nearInfo.index();
                    const point& pt = pp.localPoints()[pointI];
                    const point attraction = featPt-pt;

                    // - already attracted to feature edge : point always wins
                    // - already attracted to feature point: nearest wins

                    if (patchConstraints[pointI].first() <= 1)
                    {
                        patchAttraction[pointI] = attraction;
                        patchConstraints[pointI] = pointConstr[featPointI];
                    }
                    else if (patchConstraints[pointI].first() == 2)
                    {
                        patchAttraction[pointI] = attraction;
                        patchConstraints[pointI] = pointConstr[featPointI];
                    }
                    else if (patchConstraints[pointI].first() == 3)
                    {
                        // Only if nearer
                        if
                        (
                            magSqr(attraction)
                          < magSqr(patchAttraction[pointI])
                        )
                        {
                            patchAttraction[pointI] = attraction;
                            patchConstraints[pointI] =
                                pointConstr[featPointI];
                        }
                    }
                }
            }
        }
    }



    //MEJ: any faces that have multi-patch points only keep the multi-patch
    //     points. The other points on the face will be dragged along
    //     (hopefully)
    if (multiRegionFeatureSnap)
    {
        autoPtr<OFstream> multiPatchStr;
        if (debug&meshRefinement::OBJINTERSECTIONS)
        {
            multiPatchStr.reset
            (
                new OFstream
                (
                    meshRefiner_.mesh().time().path()
                  / "multiPatch_" + name(iter) + ".obj"
                )
            );
            Pout<< "Dumping removed constraints due to same-face"
                << " multi-patch points to "
                << multiPatchStr().name() << endl;
        }


        // 1. Mark points on multiple patches
        PackedBoolList isMultiPatchPoint(pp.size());

        forAll(pointFacePatchID, pointI)
        {
            pointIndexHit multiPatchPt = findMultiPatchPoint
            (
                pp.localPoints()[pointI],
                pointFacePatchID[pointI],
                pointFaceCentres[pointI]
            );
            isMultiPatchPoint[pointI] = multiPatchPt.hit();
        }

        // 2. Make sure multi-patch points are also attracted
        forAll(isMultiPatchPoint, pointI)
        {
            if (isMultiPatchPoint[pointI])
            {
                if
                (
                    patchConstraints[pointI].first() <= 1
                 && rawPatchConstraints[pointI].first() > 1
                )
                {
                    patchAttraction[pointI] = rawPatchAttraction[pointI];
                    patchConstraints[pointI] = rawPatchConstraints[pointI];

                    if (multiPatchStr.valid())
                    {
                        Pout<< "Adding constraint on multiPatchPoint:"
                            << pp.localPoints()[pointI]
                            << " constraint:" << patchConstraints[pointI]
                            << " attraction:" << patchAttraction[pointI]
                            << endl;
                    }
                }
            }
        }

        // Up to here it is all parallel ok.


        // 3. Knock out any attraction on faces with multi-patch points
        label nChanged = 0;
        forAll(pp.localFaces(), faceI)
        {
            const face& f = pp.localFaces()[faceI];

            label nMultiPatchPoints = 0;
            forAll(f, fp)
            {
                label pointI = f[fp];
                if
                (
                    isMultiPatchPoint[pointI]
                 && patchConstraints[pointI].first() > 1
                )
                {
                    nMultiPatchPoints++;
                }
            }

            if (nMultiPatchPoints > 0)
            {
                forAll(f, fp)
                {
                    label pointI = f[fp];
                    if
                    (
                       !isMultiPatchPoint[pointI]
                     && patchConstraints[pointI].first() > 1
                    )
                    {
                        //Pout<< "Knocking out constraint"
                        //    << " on non-multiPatchPoint:"
                        //    << pp.localPoints()[pointI] << endl;
                        patchAttraction[pointI] = vector::zero;
                        patchConstraints[pointI] = pointConstraint();
                        nChanged++;

                        if (multiPatchStr.valid())
                        {
                            meshTools::writeOBJ
                            (
                                multiPatchStr(),
                                pp.localPoints()[pointI]
                            );
                        }
                    }
                }
            }
        }

        reduce(nChanged, sumOp<label>());
        Info<< "Removing constraints near multi-patch points : changed "
            << nChanged << " points" << endl;
    }



    // Dump
    if (debug&meshRefinement::OBJINTERSECTIONS)
    {
        OFstream featureEdgeStr
        (
            meshRefiner_.mesh().time().path()
          / "edgeAttractors_" + name(iter) + ".obj"
        );
        label featureEdgeVertI = 0;
        Pout<< "Dumping feature-edge attraction to "
            << featureEdgeStr.name() << endl;

        OFstream featurePointStr
        (
            meshRefiner_.mesh().time().path()
          / "pointAttractors_" + name(iter) + ".obj"
        );
        label featurePointVertI = 0;
        Pout<< "Dumping feature-point attraction to "
            << featurePointStr.name() << endl;

        forAll(patchConstraints, pointI)
        {
            const point& pt = pp.localPoints()[pointI];

            if (patchConstraints[pointI].first() == 2)
            {
                meshTools::writeOBJ(featureEdgeStr, pt);
                featureEdgeVertI++;
                meshTools::writeOBJ
                (
                    featureEdgeStr,
                    pt+patchAttraction[pointI]
                );
                featureEdgeVertI++;
                featureEdgeStr << "l " << featureEdgeVertI-1
                    << ' ' << featureEdgeVertI << nl;
            }
            else if (patchConstraints[pointI].first() == 3)
            {
                meshTools::writeOBJ(featurePointStr, pt);
                featurePointVertI++;
                meshTools::writeOBJ
                (
                    featurePointStr,
                    pt+patchAttraction[pointI]
                );
                featurePointVertI++;
                featurePointStr << "l " << featurePointVertI-1
                    << ' ' << featurePointVertI << nl;
            }
        }
    }



    // Snap edges to feature edges
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Walk existing edges and snap remaining ones (that are marked as
    // feature edges in rawPatchConstraints)

    stringFeatureEdges
    (
        iter,
        featureCos,

        pp,
        snapDist,

        rawPatchAttraction,
        rawPatchConstraints,

        patchAttraction,
        patchConstraints
    );


    // Avoid diagonal attraction
    // ~~~~~~~~~~~~~~~~~~~~~~~~~
    // Attract one of the non-diagonal points.

    forAll(pp.localFaces(), faceI)
    {
        const face& f = pp.localFaces()[faceI];
        // For now just detect any attraction. Improve this to look at
        // actual attraction position and only if would form a problem add
        // the non-diagonal point
        if (f.size() == 4)
        {
            label nAttract = 0;
            label firstAttract = -1;
            forAll(f, fp)
            {
                label pointI = f[fp];
                if (patchConstraints[pointI].first() == 2)
                {
                    nAttract++;
                    if (firstAttract == -1)
                    {
                        firstAttract = fp;
                    }
                }
            }
            if (nAttract == 2)
            {
                label nextAttract = f.fcIndex(f.fcIndex(firstAttract));
                label pointI = f[nextAttract];

                if (patchConstraints[pointI].first() == 2)
                {
                    // Found two diagonal points that being attracted.
                    // For now just attract my one to the average of those.
                    const label i0 = f[firstAttract];
                    const point pt0 =
                        pp.localPoints()[i0]+patchAttraction[i0];
                    const label i1 = f[nextAttract];
                    const point pt1 =
                        pp.localPoints()[i1]+patchAttraction[i1];
                    const point mid = 0.5*(pt0+pt1);


                    const scalar cosAngle = mag
                    (
                        patchConstraints[i0].second()
                      & patchConstraints[i1].second()
                    );

                    //Pout<< "Found diagonal attraction at indices:"
                    //    << firstAttract
                    //    << " and " << nextAttract
                    //    << " with cosAngle:" << cosAngle
                    //    << " mid:" << mid << endl;

                    if (cosAngle > featureCos)
                    {
                        // Add the nearest of the other two points as
                        // attractor
                        label minFp = -1;
                        scalar minDistSqr = GREAT;
                        forAll(f, fp)
                        {
                            label pointI = f[fp];
                            if (patchConstraints[pointI].first() <= 1)
                            {
                                const point& pt = pp.localPoints()[pointI];
                                scalar distSqr = magSqr(mid-pt);
                                if (distSqr < minDistSqr)
                                {
                                    distSqr = minDistSqr;
                                    minFp = fp;
                                }
                            }
                        }
                        if (minFp != -1)
                        {
                            label minPointI = f[minFp];
                            patchAttraction[minPointI] =
                                mid-pp.localPoints()[minPointI];
                            patchConstraints[minPointI] =
                                patchConstraints[f[firstAttract]];
                        }
                    }
                    else
                    {
                        //Pout<< "Diagonal attractors at" << nl
                        //    << "    pt0:" << pt0
                        //    << "    constraint:"
                        //    << patchConstraints[i0].second() << nl
                        //    << "    pt1:" << pt1
                        //    << "    constraint:"
                        //    << patchConstraints[i1].second() << nl
                        //    << "    make too large an angle:"
                        //    <<  mag
                        //        (
                        //            patchConstraints[i0].second()
                        //          & patchConstraints[i1].second()
                        //        )
                        //    << endl;
                    }
                }
            }
        }
    }


    if (debug&meshRefinement::OBJINTERSECTIONS)
    {
        dumpMove
        (
            meshRefiner_.mesh().time().path()
          / "patchAttraction_" + name(iter) + ".obj",
            pp.localPoints(),
            pp.localPoints() + patchAttraction
        );
    }
}


// Correct for squeezing of face
void Foam::autoSnapDriver::preventFaceSqueeze
(
    const label iter,
    const scalar featureCos,

    const indirectPrimitivePatch& pp,
    const scalarField& snapDist,

    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    pointField points;
    face singleF;
    forAll(pp.localFaces(), faceI)
    {
        const face& f = pp.localFaces()[faceI];

        if (f.size() != points.size())
        {
            points.setSize(f.size());
            singleF.setSize(f.size());
            for (label i = 0; i < f.size(); i++)
            {
                singleF[i] = i;
            }
        }
        label nConstraints = 0;
        forAll(f, fp)
        {
            label pointI = f[fp];
            const point& pt = pp.localPoints()[pointI];

            if (patchConstraints[pointI].first() > 1)
            {
                points[fp] = pt + patchAttraction[pointI];
                nConstraints++;
            }
            else
            {
                points[fp] = pt;
            }
        }

        if (nConstraints == f.size())
        {
            scalar oldArea = f.mag(pp.localPoints());
            scalar newArea = singleF.mag(points);
            if (newArea < 0.1*oldArea)
            {
                // For now remove the point with largest distance
                label maxFp = -1;
                scalar maxS = -1;
                forAll(f, fp)
                {
                    scalar s = magSqr(patchAttraction[f[fp]]);
                    if (s > maxS)
                    {
                        maxS = s;
                        maxFp = fp;
                    }
                }
                if (maxFp != -1)
                {
                    label pointI = f[maxFp];
                    // Lower attraction on pointI
                    patchAttraction[pointI] *= 0.5;
                }
            }
        }
    }
}


Foam::vectorField Foam::autoSnapDriver::calcNearestSurfaceFeature
(
    const snapParameters& snapParams,
    const label iter,
    const scalar featureCos,
    const scalar featureAttract,
    const scalarField& snapDist,
    const vectorField& nearestDisp,
    motionSmoother& meshMover
) const
{
    const Switch implicitFeatureAttraction = snapParams.implicitFeatureSnap();
    const Switch explicitFeatureAttraction = snapParams.explicitFeatureSnap();
    const Switch multiRegionFeatureSnap = snapParams.multiRegionFeatureSnap();

    Info<< "Overriding displacement on features :" << nl
        << "   implicit features    : " << implicitFeatureAttraction << nl
        << "   explicit features    : " << explicitFeatureAttraction << nl
        << "   multi-patch features : " << multiRegionFeatureSnap << nl
        << endl;


    const indirectPrimitivePatch& pp = meshMover.patch();
    const pointField& localPoints = pp.localPoints();
    const fvMesh& mesh = meshRefiner_.mesh();


    // Displacement and orientation per pp face
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // vector from point on surface back to face centre
    vectorField faceDisp(pp.size(), vector::zero);
    // normal of surface at point on surface
    vectorField faceSurfaceNormal(pp.size(), vector::zero);
    labelList faceSurfaceGlobalRegion(pp.size(), -1);
    vectorField faceRotation(pp.size(), vector::zero);

    calcNearestFace
    (
        iter,
        pp,
        faceDisp,
        faceSurfaceNormal,
        faceSurfaceGlobalRegion,
        faceRotation
    );


    //// Displacement and orientation per pp point
    //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //vectorField pointDisp(pp.nPoints(), vector::zero);
    //vectorField pointSurfaceNormal(pp.nPoints(), vector::zero);
    //vectorField pointRotation(pp.nPoints(), vector::zero);
    //calcNearest
    //(
    //    iter,
    //    pp,
    //    pointDisp,
    //    pointSurfaceNormal,
    //    pointRotation
    //);


    // Collect (possibly remote) per point data of all surrounding faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // - faceSurfaceNormal
    // - faceDisp
    // - faceCentres&faceNormal
    List<List<point> > pointFaceSurfNormals(pp.nPoints());
    List<List<point> > pointFaceDisp(pp.nPoints());
    List<List<point> > pointFaceCentres(pp.nPoints());
    List<labelList>    pointFacePatchID(pp.nPoints());

    calcNearestFacePointProperties
    (
        iter,
        pp,

        faceDisp,
        faceSurfaceNormal,
        faceSurfaceGlobalRegion,

        pointFaceSurfNormals,
        pointFaceDisp,
        pointFaceCentres,
        pointFacePatchID
    );


    // Start off with nearest point on surface
    vectorField patchDisp = nearestDisp;


    // Main calculation
    // ~~~~~~~~~~~~~~~~
    // This is the main intelligence which calculates per point the vector to
    // attract it to the nearest surface. There are lots of possibilities
    // here.

    // Nearest feature
    vectorField patchAttraction(localPoints.size(), vector::zero);
    // Constraints at feature
    List<pointConstraint> patchConstraints(localPoints.size());


    if (implicitFeatureAttraction)
    {
        // Sample faces around each point and see if nearest surface normal
        // differs. Reconstruct a feature edge/point if possible and snap to
        // it.
        featureAttractionUsingReconstruction
        (
            iter,
            featureCos,

            pp,
            snapDist,

            pointFaceSurfNormals,
            pointFaceDisp,
            pointFaceCentres,
            pointFacePatchID,

            patchAttraction,
            patchConstraints
        );
    }

    if (explicitFeatureAttraction)
    {
        // Sample faces around each point and see if nearest surface normal
        // differs. For those find the nearest real feature edge/point and
        // store the correspondence. Then loop over feature edge/point
        // and attract those nearest mesh point. (the first phase just is
        // a subsetting of candidate points, the second makes sure that only
        // one mesh point gets attracted per feature)
        featureAttractionUsingFeatureEdges
        (
            iter,
            featureCos,
            multiRegionFeatureSnap,

            pp,
            snapDist,

            pointFaceSurfNormals,
            pointFaceDisp,
            pointFaceCentres,
            pointFacePatchID,

            patchAttraction,
            patchConstraints
        );
    }

    preventFaceSqueeze
    (
        iter,
        featureCos,

        pp,
        snapDist,

        patchAttraction,
        patchConstraints
    );


    Info<< "Attraction:" << endl
        << "     linear   : max:" << gMax(patchDisp)
        << " avg:" << gAverage(patchDisp)
        << endl
        << "     feature  : max:" << gMax(patchAttraction)
        << " avg:" << gAverage(patchAttraction)
        << endl;


    // So now we have:
    // - patchDisp          : point movement to go to nearest point on surface
    //                       (either direct or through interpolation of
    //                        face nearest)
    // - patchAttraction    : direct attraction to features
    // - patchConstraints   : type of features

    // Use any combination of patchDisp and direct feature attraction.


    // Mix with direct feature attraction
    forAll(patchConstraints, pointI)
    {
        if (patchConstraints[pointI].first() > 1)
        {
            patchDisp[pointI] =
                (1.0-featureAttract)*patchDisp[pointI]
              + featureAttract*patchAttraction[pointI];
        }
    }



    // Count
    {
        label nPlanar = 0;
        label nEdge = 0;
        label nPoint = 0;

        forAll(patchConstraints, pointI)
        {
            if (patchConstraints[pointI].first() == 1)
            {
                nPlanar++;
            }
            else if (patchConstraints[pointI].first() == 2)
            {
                nEdge++;
            }
            else if (patchConstraints[pointI].first() == 3)
            {
                nPoint++;
            }
        }

        label nTotPoints = returnReduce(pp.nPoints(), sumOp<label>());
        reduce(nPlanar, sumOp<label>());
        reduce(nEdge, sumOp<label>());
        reduce(nPoint, sumOp<label>());
        Info<< "Feature analysis : total points:"
            << nTotPoints
            << " attraction to :" << nl
            << "    feature point   : " << nPoint << nl
            << "    feature edge    : " << nEdge << nl
            << "    nearest surface : " << nPlanar << nl
            << "    rest            : " << nTotPoints-nPoint-nEdge-nPlanar
            << nl
            << endl;
    }


    // Now we have the displacement per patch point to move onto the surface
    // Split into tangential and normal direction.
    // - start off with all non-constrained points following the constrained
    //   ones since point normals not relevant.
    // - finish with only tangential component smoothed.
    // Note: tangential is most
    // likely to come purely from face-centre snapping, not face rotation.
    // Note: could use the constraints here (constraintTransformation())
    //       but this is not necessarily accurate and we're smoothing to
    //       get out of problems.

    if (featureAttract < 1-0.001)
    {
        // 1. Smoothed all displacement
        vectorField smoothedPatchDisp = patchDisp;
        smoothAndConstrain
        (
            pp,
            patchConstraints,
            smoothedPatchDisp
        );


        // 2. Smoothed tangential component
        vectorField tangPatchDisp = patchDisp;
        tangPatchDisp -= (pp.pointNormals() & patchDisp) * pp.pointNormals();
        smoothAndConstrain
        (
            pp,
            patchConstraints,
            tangPatchDisp
        );

        // Re-add normal component
        tangPatchDisp += (pp.pointNormals() & patchDisp) * pp.pointNormals();

        if (debug&meshRefinement::OBJINTERSECTIONS)
        {
            dumpMove
            (
                mesh.time().path()
              / "tangPatchDispConstrained_" + name(iter) + ".obj",
                pp.localPoints(),
                pp.localPoints() + tangPatchDisp
            );
        }

        patchDisp =
             (1.0-featureAttract)*smoothedPatchDisp
           + featureAttract*tangPatchDisp;
    }


    const scalar relax = featureAttract;
    patchDisp *= relax;


    // Points on zones in one domain but only present as point on other
    // will not do condition 2 on all. Sync explicitly.
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        patchDisp,
        minMagSqrEqOp<point>(),         // combine op
        vector(GREAT, GREAT, GREAT)     // null value (note: cant use VGREAT)
    );

    return patchDisp;
}


// ************************************************************************* //
