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

#include "meshRefinement.H"
#include "trackedParticle.H"
#include "syncTools.H"
#include "Time.H"
#include "refinementSurfaces.H"
#include "refinementFeatures.H"
#include "shellSurfaces.H"
#include "faceSet.H"
#include "decompositionMethod.H"
#include "fvMeshDistribute.H"
#include "polyTopoChange.H"
#include "mapDistributePolyMesh.H"
#include "featureEdgeMesh.H"
#include "Cloud.H"
//#include "globalIndex.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Get faces (on the new mesh) that have in some way been affected by the
// mesh change. Picks up all faces but those that are between two
// unrefined faces. (Note that of an unchanged face the edge still might be
// split but that does not change any face centre or cell centre.
Foam::labelList Foam::meshRefinement::getChangedFaces
(
    const mapPolyMesh& map,
    const labelList& oldCellsToRefine
)
{
    const polyMesh& mesh = map.mesh();

    labelList changedFaces;

    // For reporting: number of masterFaces changed
    label nMasterChanged = 0;
    PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh));

    {
        // Mark any face on a cell which has been added or changed
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // Note that refining a face changes the face centre (for a warped face)
        // which changes the cell centre. This again changes the cellcentre-
        // cellcentre edge across all faces of the cell.
        // Note: this does not happen for unwarped faces but unfortunately
        // we don't have this information.

        const labelList& faceOwner = mesh.faceOwner();
        const labelList& faceNeighbour = mesh.faceNeighbour();
        const cellList& cells = mesh.cells();
        const label nInternalFaces = mesh.nInternalFaces();

        // Mark refined cells on old mesh
        PackedBoolList oldRefineCell(map.nOldCells());

        forAll(oldCellsToRefine, i)
        {
            oldRefineCell.set(oldCellsToRefine[i], 1u);
        }

        // Mark refined faces
        PackedBoolList refinedInternalFace(nInternalFaces);

        // 1. Internal faces

        for (label faceI = 0; faceI < nInternalFaces; faceI++)
        {
            label oldOwn = map.cellMap()[faceOwner[faceI]];
            label oldNei = map.cellMap()[faceNeighbour[faceI]];

            if
            (
                oldOwn >= 0
             && oldRefineCell.get(oldOwn) == 0u
             && oldNei >= 0
             && oldRefineCell.get(oldNei) == 0u
            )
            {
                // Unaffected face since both neighbours were not refined.
            }
            else
            {
                refinedInternalFace.set(faceI, 1u);
            }
        }


        // 2. Boundary faces

        boolList refinedBoundaryFace(mesh.nFaces()-nInternalFaces, false);

        forAll(mesh.boundaryMesh(), patchI)
        {
            const polyPatch& pp = mesh.boundaryMesh()[patchI];

            label faceI = pp.start();

            forAll(pp, i)
            {
                label oldOwn = map.cellMap()[faceOwner[faceI]];

                if (oldOwn >= 0 && oldRefineCell.get(oldOwn) == 0u)
                {
                    // owner did exist and wasn't refined.
                }
                else
                {
                    refinedBoundaryFace[faceI-nInternalFaces] = true;
                }
                faceI++;
            }
        }

        // Synchronise refined face status
        syncTools::syncBoundaryFaceList
        (
            mesh,
            refinedBoundaryFace,
            orEqOp<bool>()
        );


        // 3. Mark all faces affected by refinement. Refined faces are in
        //    - refinedInternalFace
        //    - refinedBoundaryFace
        boolList changedFace(mesh.nFaces(), false);

        forAll(refinedInternalFace, faceI)
        {
            if (refinedInternalFace.get(faceI) == 1u)
            {
                const cell& ownFaces = cells[faceOwner[faceI]];
                forAll(ownFaces, ownI)
                {
                    changedFace[ownFaces[ownI]] = true;
                }
                const cell& neiFaces = cells[faceNeighbour[faceI]];
                forAll(neiFaces, neiI)
                {
                    changedFace[neiFaces[neiI]] = true;
                }
            }
        }

        forAll(refinedBoundaryFace, i)
        {
            if (refinedBoundaryFace[i])
            {
                const cell& ownFaces = cells[faceOwner[i+nInternalFaces]];
                forAll(ownFaces, ownI)
                {
                    changedFace[ownFaces[ownI]] = true;
                }
            }
        }

        syncTools::syncFaceList
        (
            mesh,
            changedFace,
            orEqOp<bool>()
        );


        // Now we have in changedFace marked all affected faces. Pack.
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        changedFaces = findIndices(changedFace, true);

        // Count changed master faces.
        nMasterChanged = 0;

        forAll(changedFace, faceI)
        {
            if (changedFace[faceI] && isMasterFace[faceI])
            {
                nMasterChanged++;
            }
        }

    }

    if (debug&meshRefinement::MESH)
    {
        Pout<< "getChangedFaces : Detected "
            << " local:" << changedFaces.size()
            << " global:" << returnReduce(nMasterChanged, sumOp<label>())
            << " changed faces out of " << mesh.globalData().nTotalFaces()
            << endl;

        faceSet changedFacesSet(mesh, "changedFaces", changedFaces);
        Pout<< "getChangedFaces : Writing " << changedFaces.size()
            << " changed faces to faceSet " << changedFacesSet.name()
            << endl;
        changedFacesSet.write();
    }

    return changedFaces;
}


// Mark cell for refinement (if not already marked). Return false if
// refinelimit hit. Keeps running count (in nRefine) of cells marked for
// refinement
bool Foam::meshRefinement::markForRefine
(
    const label markValue,
    const label nAllowRefine,

    label& cellValue,
    label& nRefine
)
{
    if (cellValue == -1)
    {
        cellValue = markValue;
        nRefine++;
    }

    return nRefine <= nAllowRefine;
}


// Calculates list of cells to refine based on intersection with feature edge.
Foam::label Foam::meshRefinement::markFeatureRefinement
(
    const point& keepPoint,
    const label nAllowRefine,

    labelList& refineCell,
    label& nRefine
) const
{
    // We want to refine all cells containing a feature edge.
    // - don't want to search over whole mesh
    // - don't want to build octree for whole mesh
    // - so use tracking to follow the feature edges
    //
    // 1. find non-manifold points on feature edges (i.e. start of feature edge
    //    or 'knot')
    // 2. seed particle starting at keepPoint going to this non-manifold point
    // 3. track particles to their non-manifold point
    //
    // 4. track particles across their connected feature edges, marking all
    //    visited cells with their level (through trackingData)
    // 5. do 4 until all edges have been visited.


    // Find all start cells of features. Is done by tracking from keepPoint.
    Cloud<trackedParticle> cloud(mesh_, IDLList<trackedParticle>());


    // Features are identical on all processors. Number them so we know
    // what to seed. Do this on only the processor that
    // holds the keepPoint.

    label cellI = -1;
    label tetFaceI = -1;
    label tetPtI = -1;

    mesh_.findCellFacePt(keepPoint, cellI, tetFaceI, tetPtI);

    if (cellI != -1)
    {
        // I am the processor that holds the keepPoint

        forAll(features_, featI)
        {
            const featureEdgeMesh& featureMesh = features_[featI];
            const label featureLevel = features_.levels()[featI][0];
            const labelListList& pointEdges = featureMesh.pointEdges();

            // Find regions on edgeMesh
            labelList edgeRegion;
            label nRegions = featureMesh.regions(edgeRegion);


            PackedBoolList regionVisited(nRegions);


            // 1. Seed all 'knots' in edgeMesh


            forAll(pointEdges, pointI)
            {
                if (pointEdges[pointI].size() != 2)
                {
                    if (debug&meshRefinement::FEATURESEEDS)
                    {
                        Pout<< "Adding particle from point:" << pointI
                            << " coord:" << featureMesh.points()[pointI]
                            << " since number of emanating edges:"
                            << pointEdges[pointI].size()
                            << endl;
                    }

                    // Non-manifold point. Create particle.
                    cloud.addParticle
                    (
                        new trackedParticle
                        (
                            mesh_,
                            keepPoint,
                            cellI,
                            tetFaceI,
                            tetPtI,
                            featureMesh.points()[pointI],   // endpos
                            featureLevel,                   // level
                            featI,                          // featureMesh
                            pointI                          // end point
                        )
                    );

                    // Mark
                    if (pointEdges[pointI].size() > 0)
                    {
                        label e0 = pointEdges[pointI][0];
                        label regionI = edgeRegion[e0];
                        regionVisited[regionI] = 1u;
                    }
                }
            }


            // 2. Any regions that have not been visited at all? These can
            //    only be circular regions!
            forAll(featureMesh.edges(), edgeI)
            {
                if (regionVisited.set(edgeRegion[edgeI], 1u))
                {
                    const edge& e = featureMesh.edges()[edgeI];
                    label pointI = e.start();
                    if (debug&meshRefinement::FEATURESEEDS)
                    {
                        Pout<< "Adding particle from point:" << pointI
                            << " coord:" << featureMesh.points()[pointI]
                            << " on circular region:" << edgeRegion[edgeI]
                            << endl;
                    }

                    // Non-manifold point. Create particle.
                    cloud.addParticle
                    (
                        new trackedParticle
                        (
                            mesh_,
                            keepPoint,
                            cellI,
                            tetFaceI,
                            tetPtI,
                            featureMesh.points()[pointI],   // endpos
                            featureLevel,                   // level
                            featI,                          // featureMesh
                            pointI                          // end point
                        )
                    );
                }
            }
        }
    }


    // Largest refinement level of any feature passed through
    labelList maxFeatureLevel(mesh_.nCells(), -1);

    // Database to pass into trackedParticle::move
    trackedParticle::trackingData td(cloud, maxFeatureLevel);

    // Track all particles to their end position (= starting feature point)
    // Note that the particle might have started on a different processor
    // so this will transfer across nicely until we can start tracking proper.
    cloud.move(td, GREAT);

    // Reset level
    maxFeatureLevel = -1;

    // Whether edge has been visited.
    List<PackedBoolList> featureEdgeVisited(features_.size());

    forAll(features_, featI)
    {
        featureEdgeVisited[featI].setSize(features_[featI].edges().size());
        featureEdgeVisited[featI] = 0u;
    }

    while (true)
    {
        label nParticles = 0;

        // Make particle follow edge.
        forAllIter(Cloud<trackedParticle>, cloud, iter)
        {
            trackedParticle& tp = iter();

            label featI = tp.i();
            label pointI = tp.j();

            const featureEdgeMesh& featureMesh = features_[featI];
            const labelList& pEdges = featureMesh.pointEdges()[pointI];

            // Particle now at pointI. Check connected edges to see which one
            // we have to visit now.

            bool keepParticle = false;

            forAll(pEdges, i)
            {
                label edgeI = pEdges[i];

                if (featureEdgeVisited[featI].set(edgeI, 1u))
                {
                    // Unvisited edge. Make the particle go to the other point
                    // on the edge.

                    const edge& e = featureMesh.edges()[edgeI];
                    label otherPointI = e.otherVertex(pointI);

                    tp.end() = featureMesh.points()[otherPointI];
                    tp.j() = otherPointI;
                    keepParticle = true;
                    break;
                }
            }

            if (!keepParticle)
            {
                // Particle at 'knot' where another particle already has been
                // seeded. Delete particle.
                cloud.deleteParticle(tp);
            }
            else
            {
                // Keep particle
                nParticles++;
            }
        }

        reduce(nParticles, sumOp<label>());
        if (nParticles == 0)
        {
            break;
        }

        // Track all particles to their end position.
        cloud.move(td, GREAT);
    }


    // See which cells to refine. maxFeatureLevel will hold highest level
    // of any feature edge that passed through.

    const labelList& cellLevel = meshCutter_.cellLevel();

    label oldNRefine = nRefine;

    forAll(maxFeatureLevel, cellI)
    {
        if (maxFeatureLevel[cellI] > cellLevel[cellI])
        {
            // Mark
            if
            (
               !markForRefine
                (
                    0,                      // surface (n/a)
                    nAllowRefine,
                    refineCell[cellI],
                    nRefine
                )
            )
            {
                // Reached limit
                break;
            }
        }
    }

    if
    (
        returnReduce(nRefine, sumOp<label>())
      > returnReduce(nAllowRefine, sumOp<label>())
    )
    {
        Info<< "Reached refinement limit." << endl;
    }

    return returnReduce(nRefine-oldNRefine,  sumOp<label>());
}


// Mark cells for distance-to-feature based refinement.
Foam::label Foam::meshRefinement::markInternalDistanceToFeatureRefinement
(
    const label nAllowRefine,

    labelList& refineCell,
    label& nRefine
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const pointField& cellCentres = mesh_.cellCentres();

    // Detect if there are any distance shells
    if (features_.maxDistance() <= 0.0)
    {
        return 0;
    }

    label oldNRefine = nRefine;

    // Collect cells to test
    pointField testCc(cellLevel.size()-nRefine);
    labelList testLevels(cellLevel.size()-nRefine);
    label testI = 0;

    forAll(cellLevel, cellI)
    {
        if (refineCell[cellI] == -1)
        {
            testCc[testI] = cellCentres[cellI];
            testLevels[testI] = cellLevel[cellI];
            testI++;
        }
    }

    // Do test to see whether cells is inside/outside shell with higher level
    labelList maxLevel;
    features_.findHigherLevel(testCc, testLevels, maxLevel);

    // Mark for refinement. Note that we didn't store the original cellID so
    // now just reloop in same order.
    testI = 0;
    forAll(cellLevel, cellI)
    {
        if (refineCell[cellI] == -1)
        {
            if (maxLevel[testI] > testLevels[testI])
            {
                bool reachedLimit = !markForRefine
                (
                    maxLevel[testI],    // mark with any positive value
                    nAllowRefine,
                    refineCell[cellI],
                    nRefine
                );

                if (reachedLimit)
                {
                    if (debug)
                    {
                        Pout<< "Stopped refining internal cells"
                            << " since reaching my cell limit of "
                            << mesh_.nCells()+7*nRefine << endl;
                    }
                    break;
                }
            }
            testI++;
        }
    }

    if
    (
        returnReduce(nRefine, sumOp<label>())
      > returnReduce(nAllowRefine, sumOp<label>())
    )
    {
        Info<< "Reached refinement limit." << endl;
    }

    return returnReduce(nRefine-oldNRefine, sumOp<label>());
}


// Mark cells for non-surface intersection based refinement.
Foam::label Foam::meshRefinement::markInternalRefinement
(
    const label nAllowRefine,

    labelList& refineCell,
    label& nRefine
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const pointField& cellCentres = mesh_.cellCentres();

    label oldNRefine = nRefine;

    // Collect cells to test
    pointField testCc(cellLevel.size()-nRefine);
    labelList testLevels(cellLevel.size()-nRefine);
    label testI = 0;

    forAll(cellLevel, cellI)
    {
        if (refineCell[cellI] == -1)
        {
            testCc[testI] = cellCentres[cellI];
            testLevels[testI] = cellLevel[cellI];
            testI++;
        }
    }

    // Do test to see whether cells is inside/outside shell with higher level
    labelList maxLevel;
    shells_.findHigherLevel(testCc, testLevels, maxLevel);

    // Mark for refinement. Note that we didn't store the original cellID so
    // now just reloop in same order.
    testI = 0;
    forAll(cellLevel, cellI)
    {
        if (refineCell[cellI] == -1)
        {
            if (maxLevel[testI] > testLevels[testI])
            {
                bool reachedLimit = !markForRefine
                (
                    maxLevel[testI],    // mark with any positive value
                    nAllowRefine,
                    refineCell[cellI],
                    nRefine
                );

                if (reachedLimit)
                {
                    if (debug)
                    {
                        Pout<< "Stopped refining internal cells"
                            << " since reaching my cell limit of "
                            << mesh_.nCells()+7*nRefine << endl;
                    }
                    break;
                }
            }
            testI++;
        }
    }

    if
    (
        returnReduce(nRefine, sumOp<label>())
      > returnReduce(nAllowRefine, sumOp<label>())
    )
    {
        Info<< "Reached refinement limit." << endl;
    }

    return returnReduce(nRefine-oldNRefine, sumOp<label>());
}


// Collect faces that are intersected and whose neighbours aren't yet marked
// for refinement.
Foam::labelList Foam::meshRefinement::getRefineCandidateFaces
(
    const labelList& refineCell
) const
{
    labelList testFaces(mesh_.nFaces());

    label nTest = 0;

    forAll(surfaceIndex_, faceI)
    {
        if (surfaceIndex_[faceI] != -1)
        {
            label own = mesh_.faceOwner()[faceI];

            if (mesh_.isInternalFace(faceI))
            {
                label nei = mesh_.faceNeighbour()[faceI];

                if (refineCell[own] == -1 || refineCell[nei] == -1)
                {
                    testFaces[nTest++] = faceI;
                }
            }
            else
            {
                if (refineCell[own] == -1)
                {
                    testFaces[nTest++] = faceI;
                }
            }
        }
    }
    testFaces.setSize(nTest);

    return testFaces;
}


// Mark cells for surface intersection based refinement.
Foam::label Foam::meshRefinement::markSurfaceRefinement
(
    const label nAllowRefine,
    const labelList& neiLevel,
    const pointField& neiCc,

    labelList& refineCell,
    label& nRefine
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const pointField& cellCentres = mesh_.cellCentres();

    label oldNRefine = nRefine;

    // Use cached surfaceIndex_ to detect if any intersection. If so
    // re-intersect to determine level wanted.

    // Collect candidate faces
    // ~~~~~~~~~~~~~~~~~~~~~~~

    labelList testFaces(getRefineCandidateFaces(refineCell));

    // Collect segments
    // ~~~~~~~~~~~~~~~~

    pointField start(testFaces.size());
    pointField end(testFaces.size());
    labelList minLevel(testFaces.size());

    forAll(testFaces, i)
    {
        label faceI = testFaces[i];

        label own = mesh_.faceOwner()[faceI];

        if (mesh_.isInternalFace(faceI))
        {
            label nei = mesh_.faceNeighbour()[faceI];

            start[i] = cellCentres[own];
            end[i] = cellCentres[nei];
            minLevel[i] = min(cellLevel[own], cellLevel[nei]);
        }
        else
        {
            label bFaceI = faceI - mesh_.nInternalFaces();

            start[i] = cellCentres[own];
            end[i] = neiCc[bFaceI];
            minLevel[i] = min(cellLevel[own], neiLevel[bFaceI]);
        }
    }

    // Extend segments a bit
    {
        const vectorField smallVec(Foam::sqrt(SMALL)*(end-start));
        start -= smallVec;
        end += smallVec;
    }


    // Do test for higher intersections
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList surfaceHit;
    labelList surfaceMinLevel;
    surfaces_.findHigherIntersection
    (
        start,
        end,
        minLevel,

        surfaceHit,
        surfaceMinLevel
    );


    // Mark cells for refinement
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(testFaces, i)
    {
        label faceI = testFaces[i];

        label surfI = surfaceHit[i];

        if (surfI != -1)
        {
            // Found intersection with surface with higher wanted
            // refinement. Check if the level field on that surface
            // specifies an even higher level. Note:this is weird. Should
            // do the check with the surfaceMinLevel whilst intersecting the
            // surfaces?

            label own = mesh_.faceOwner()[faceI];

            if (surfaceMinLevel[i] > cellLevel[own])
            {
                // Owner needs refining
                if
                (
                   !markForRefine
                    (
                        surfI,
                        nAllowRefine,
                        refineCell[own],
                        nRefine
                    )
                )
                {
                    break;
                }
            }

            if (mesh_.isInternalFace(faceI))
            {
                label nei = mesh_.faceNeighbour()[faceI];
                if (surfaceMinLevel[i] > cellLevel[nei])
                {
                    // Neighbour needs refining
                    if
                    (
                       !markForRefine
                        (
                            surfI,
                            nAllowRefine,
                            refineCell[nei],
                            nRefine
                        )
                    )
                    {
                        break;
                    }
                }
            }
        }
    }

    if
    (
        returnReduce(nRefine, sumOp<label>())
      > returnReduce(nAllowRefine, sumOp<label>())
    )
    {
        Info<< "Reached refinement limit." << endl;
    }

    return returnReduce(nRefine-oldNRefine, sumOp<label>());
}


// Checks if multiple intersections of a cell (by a surface with a higher
// max than the cell level) and if so if the normals at these intersections
// make a large angle.
// Returns false if the nRefine limit has been reached, true otherwise.
bool Foam::meshRefinement::checkCurvature
(
    const scalar curvature,
    const label nAllowRefine,

    const label surfaceLevel,   // current intersection max level
    const vector& surfaceNormal,// current intersection normal

    const label cellI,

    label& cellMaxLevel,        // cached max surface level for this cell
    vector& cellMaxNormal,      // cached surface normal for this cell

    labelList& refineCell,
    label& nRefine
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();

    // Test if surface applicable
    if (surfaceLevel > cellLevel[cellI])
    {
        if (cellMaxLevel == -1)
        {
            // First visit of cell. Store
            cellMaxLevel = surfaceLevel;
            cellMaxNormal = surfaceNormal;
        }
        else
        {
            // Second or more visit. Check curvature.
            if ((cellMaxNormal & surfaceNormal) < curvature)
            {
                return markForRefine
                (
                    surfaceLevel,   // mark with any non-neg number.
                    nAllowRefine,
                    refineCell[cellI],
                    nRefine
                );
            }

            // Set normal to that of highest surface. Not really necessary
            // over here but we reuse cellMax info when doing coupled faces.
            if (surfaceLevel > cellMaxLevel)
            {
                cellMaxLevel = surfaceLevel;
                cellMaxNormal = surfaceNormal;
            }
        }
    }

    // Did not reach refinement limit.
    return true;
}


// Mark cells for surface curvature based refinement.
Foam::label Foam::meshRefinement::markSurfaceCurvatureRefinement
(
    const scalar curvature,
    const label nAllowRefine,
    const labelList& neiLevel,
    const pointField& neiCc,

    labelList& refineCell,
    label& nRefine
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const pointField& cellCentres = mesh_.cellCentres();

    label oldNRefine = nRefine;

    // 1. local test: any cell on more than one surface gets refined
    // (if its current level is < max of the surface max level)

    // 2. 'global' test: any cell on only one surface with a neighbour
    // on a different surface gets refined (if its current level etc.)


    // Collect candidate faces (i.e. intersecting any surface and
    // owner/neighbour not yet refined.
    labelList testFaces(getRefineCandidateFaces(refineCell));

    // Collect segments
    pointField start(testFaces.size());
    pointField end(testFaces.size());
    labelList minLevel(testFaces.size());

    forAll(testFaces, i)
    {
        label faceI = testFaces[i];

        label own = mesh_.faceOwner()[faceI];

        if (mesh_.isInternalFace(faceI))
        {
            label nei = mesh_.faceNeighbour()[faceI];

            start[i] = cellCentres[own];
            end[i] = cellCentres[nei];
            minLevel[i] = min(cellLevel[own], cellLevel[nei]);
        }
        else
        {
            label bFaceI = faceI - mesh_.nInternalFaces();

            start[i] = cellCentres[own];
            end[i] = neiCc[bFaceI];
            minLevel[i] = min(cellLevel[own], neiLevel[bFaceI]);
        }
    }

    // Extend segments a bit
    {
        const vectorField smallVec(Foam::sqrt(SMALL)*(end-start));
        start -= smallVec;
        end += smallVec;
    }


    // Test for all intersections (with surfaces of higher max level than
    // minLevel) and cache per cell the max surface level and the local normal
    // on that surface.
    labelList cellMaxLevel(mesh_.nCells(), -1);
    vectorField cellMaxNormal(mesh_.nCells(), vector::zero);

    {
        // Per segment the normals of the surfaces
        List<vectorList> surfaceNormal;
        // Per segment the list of levels of the surfaces
        labelListList surfaceLevel;

        surfaces_.findAllHigherIntersections
        (
            start,
            end,
            minLevel,           // max level of surface has to be bigger
                                // than min level of neighbouring cells
            surfaceNormal,
            surfaceLevel
        );
        // Clear out unnecessary data
        start.clear();
        end.clear();
        minLevel.clear();

        // Extract per cell information on the surface with the highest max
        forAll(testFaces, i)
        {
            label faceI = testFaces[i];
            label own = mesh_.faceOwner()[faceI];

            const vectorList& fNormals = surfaceNormal[i];
            const labelList& fLevels = surfaceLevel[i];

            forAll(fLevels, hitI)
            {
                checkCurvature
                (
                    curvature,
                    nAllowRefine,

                    fLevels[hitI],
                    fNormals[hitI],

                    own,
                    cellMaxLevel[own],
                    cellMaxNormal[own],

                    refineCell,
                    nRefine
                );
            }

            if (mesh_.isInternalFace(faceI))
            {
                label nei = mesh_.faceNeighbour()[faceI];

                forAll(fLevels, hitI)
                {
                    checkCurvature
                    (
                        curvature,
                        nAllowRefine,

                        fLevels[hitI],
                        fNormals[hitI],

                        nei,
                        cellMaxLevel[nei],
                        cellMaxNormal[nei],

                        refineCell,
                        nRefine
                    );
                }
            }
        }
    }

    // 2. Find out a measure of surface curvature and region edges.
    // Send over surface region and surface normal to neighbour cell.

    labelList neiBndMaxLevel(mesh_.nFaces()-mesh_.nInternalFaces());
    vectorField neiBndMaxNormal(mesh_.nFaces()-mesh_.nInternalFaces());

    for (label faceI = mesh_.nInternalFaces(); faceI < mesh_.nFaces(); faceI++)
    {
        label bFaceI = faceI-mesh_.nInternalFaces();
        label own = mesh_.faceOwner()[faceI];

        neiBndMaxLevel[bFaceI] = cellMaxLevel[own];
        neiBndMaxNormal[bFaceI] = cellMaxNormal[own];
    }
    syncTools::swapBoundaryFaceList(mesh_, neiBndMaxLevel);
    syncTools::swapBoundaryFaceList(mesh_, neiBndMaxNormal);

    // Loop over all faces. Could only be checkFaces.. except if they're coupled

    // Internal faces
    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        label own = mesh_.faceOwner()[faceI];
        label nei = mesh_.faceNeighbour()[faceI];

        if (cellMaxLevel[own] != -1 && cellMaxLevel[nei] != -1)
        {
            // Have valid data on both sides. Check curvature.
            if ((cellMaxNormal[own] & cellMaxNormal[nei]) < curvature)
            {
                // See which side to refine
                if (cellLevel[own] < cellMaxLevel[own])
                {
                    if
                    (
                        !markForRefine
                        (
                            cellMaxLevel[own],
                            nAllowRefine,
                            refineCell[own],
                            nRefine
                        )
                    )
                    {
                        if (debug)
                        {
                            Pout<< "Stopped refining since reaching my cell"
                                << " limit of " << mesh_.nCells()+7*nRefine
                                << endl;
                        }
                        break;
                    }
                }

                if (cellLevel[nei] < cellMaxLevel[nei])
                {
                    if
                    (
                        !markForRefine
                        (
                            cellMaxLevel[nei],
                            nAllowRefine,
                            refineCell[nei],
                            nRefine
                        )
                    )
                    {
                        if (debug)
                        {
                            Pout<< "Stopped refining since reaching my cell"
                                << " limit of " << mesh_.nCells()+7*nRefine
                                << endl;
                        }
                        break;
                    }
                }
            }
        }
    }
    // Boundary faces
    for (label faceI = mesh_.nInternalFaces(); faceI < mesh_.nFaces(); faceI++)
    {
        label own = mesh_.faceOwner()[faceI];
        label bFaceI = faceI - mesh_.nInternalFaces();

        if (cellLevel[own] < cellMaxLevel[own] && neiBndMaxLevel[bFaceI] != -1)
        {
            // Have valid data on both sides. Check curvature.
            if ((cellMaxNormal[own] & neiBndMaxNormal[bFaceI]) < curvature)
            {
                if
                (
                    !markForRefine
                    (
                        cellMaxLevel[own],
                        nAllowRefine,
                        refineCell[own],
                        nRefine
                    )
                )
                {
                    if (debug)
                    {
                        Pout<< "Stopped refining since reaching my cell"
                            << " limit of " << mesh_.nCells()+7*nRefine
                            << endl;
                    }
                    break;
                }
            }
        }
    }

    if
    (
        returnReduce(nRefine, sumOp<label>())
      > returnReduce(nAllowRefine, sumOp<label>())
    )
    {
        Info<< "Reached refinement limit." << endl;
    }

    return returnReduce(nRefine-oldNRefine, sumOp<label>());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Calculate list of cells to refine. Gets for any edge (start - end)
// whether it intersects the surface. Does more accurate test and checks
// the wanted level on the surface intersected.
// Does approximate precalculation of how many cells can be refined before
// hitting overall limit maxGlobalCells.
Foam::labelList Foam::meshRefinement::refineCandidates
(
    const point& keepPoint,
    const scalar curvature,

    const bool featureRefinement,
    const bool featureDistanceRefinement,
    const bool internalRefinement,
    const bool surfaceRefinement,
    const bool curvatureRefinement,
    const label maxGlobalCells,
    const label maxLocalCells
) const
{
    label totNCells = mesh_.globalData().nTotalCells();

    labelList cellsToRefine;

    if (totNCells >= maxGlobalCells)
    {
        Info<< "No cells marked for refinement since reached limit "
            << maxGlobalCells << '.' << endl;
    }
    else
    {
        // Every cell I refine adds 7 cells. Estimate the number of cells
        // I am allowed to refine.
        // Assume perfect distribution so can only refine as much the fraction
        // of the mesh I hold. This prediction step prevents us having to do
        // lots of reduces to keep count of the total number of cells selected
        // for refinement.

        //scalar fraction = scalar(mesh_.nCells())/totNCells;
        //label myMaxCells = label(maxGlobalCells*fraction);
        //label nAllowRefine = (myMaxCells - mesh_.nCells())/7;
        ////label nAllowRefine = (maxLocalCells - mesh_.nCells())/7;
        //
        //Pout<< "refineCandidates:" << nl
        //    << "    total cells:" << totNCells << nl
        //    << "    local cells:" << mesh_.nCells() << nl
        //    << "    local fraction:" << fraction << nl
        //    << "    local allowable cells:" << myMaxCells << nl
        //    << "    local allowable refinement:" << nAllowRefine << nl
        //    << endl;

        //- Disable refinement shortcut. nAllowRefine is per processor limit.
        label nAllowRefine = labelMax / Pstream::nProcs();

        // Marked for refinement (>= 0) or not (-1). Actual value is the
        // index of the surface it intersects.
        labelList refineCell(mesh_.nCells(), -1);
        label nRefine = 0;


        // Swap neighbouring cell centres and cell level
        labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());
        pointField neiCc(mesh_.nFaces()-mesh_.nInternalFaces());
        calcNeighbourData(neiLevel, neiCc);



        // Cells pierced by feature lines
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (featureRefinement)
        {
            label nFeatures = markFeatureRefinement
            (
                keepPoint,
                nAllowRefine,

                refineCell,
                nRefine
            );

            Info<< "Marked for refinement due to explicit features             "
                << ": " << nFeatures << " cells."  << endl;
        }

        // Inside distance-to-feature shells
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (featureDistanceRefinement)
        {
            label nShell = markInternalDistanceToFeatureRefinement
            (
                nAllowRefine,

                refineCell,
                nRefine
            );
            Info<< "Marked for refinement due to distance to explicit features "
                ": " << nShell << " cells."  << endl;
        }

        // Inside refinement shells
        // ~~~~~~~~~~~~~~~~~~~~~~~~

        if (internalRefinement)
        {
            label nShell = markInternalRefinement
            (
                nAllowRefine,

                refineCell,
                nRefine
            );
            Info<< "Marked for refinement due to refinement shells             "
                << ": " << nShell << " cells."  << endl;
        }

        // Refinement based on intersection of surface
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (surfaceRefinement)
        {
            label nSurf = markSurfaceRefinement
            (
                nAllowRefine,
                neiLevel,
                neiCc,

                refineCell,
                nRefine
            );
            Info<< "Marked for refinement due to surface intersection          "
                << ": " << nSurf << " cells."  << endl;
        }

        // Refinement based on curvature of surface
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if
        (
            curvatureRefinement
         && (curvature >= -1 && curvature <= 1)
         && (surfaces_.minLevel() != surfaces_.maxLevel())
        )
        {
            label nCurv = markSurfaceCurvatureRefinement
            (
                curvature,
                nAllowRefine,
                neiLevel,
                neiCc,

                refineCell,
                nRefine
            );
            Info<< "Marked for refinement due to curvature/regions             "
                << ": " << nCurv << " cells."  << endl;
        }

        // Pack cells-to-refine
        // ~~~~~~~~~~~~~~~~~~~~

        cellsToRefine.setSize(nRefine);
        nRefine = 0;

        forAll(refineCell, cellI)
        {
            if (refineCell[cellI] != -1)
            {
                cellsToRefine[nRefine++] = cellI;
            }
        }
    }

    return cellsToRefine;
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::refine
(
    const labelList& cellsToRefine
)
{
    // Mesh changing engine.
    polyTopoChange meshMod(mesh_);

    // Play refinement commands into mesh changer.
    meshCutter_.setRefinement(cellsToRefine, meshMod);

    // Create mesh (no inflation), return map from old to new mesh.
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false);

    // Update fields
    mesh_.updateMesh(map);

    // Optionally inflate mesh
    if (map().hasMotionPoints())
    {
        mesh_.movePoints(map().preMotionPoints());
    }
    else
    {
        // Delete mesh volumes.
        mesh_.clearOut();
    }

    // Reset the instance for if in overwrite mode
    mesh_.setInstance(timeName());

    // Update intersection info
    updateMesh(map, getChangedFaces(map, cellsToRefine));

    return map;
}


// Do refinement of consistent set of cells followed by truncation and
// load balancing.
Foam::autoPtr<Foam::mapDistributePolyMesh>
Foam::meshRefinement::refineAndBalance
(
    const string& msg,
    decompositionMethod& decomposer,
    fvMeshDistribute& distributor,
    const labelList& cellsToRefine,
    const scalar maxLoadUnbalance
)
{
    // Do all refinement
    refine(cellsToRefine);

    if (debug&meshRefinement::MESH)
    {
        Pout<< "Writing refined but unbalanced " << msg
            << " mesh to time " << timeName() << endl;
        write
        (
            debug,
            mesh_.time().path()
           /timeName()
        );
        Pout<< "Dumped debug data in = "
            << mesh_.time().cpuTimeIncrement() << " s" << endl;

        // test all is still synced across proc patches
        checkData();
    }

    Info<< "Refined mesh in = "
        << mesh_.time().cpuTimeIncrement() << " s" << endl;
    printMeshInfo(debug, "After refinement " + msg);


    // Load balancing
    // ~~~~~~~~~~~~~~

    autoPtr<mapDistributePolyMesh> distMap;

    if (Pstream::nProcs() > 1)
    {
        scalar nIdealCells =
            mesh_.globalData().nTotalCells()
          / Pstream::nProcs();

        scalar unbalance = returnReduce
        (
            mag(1.0-mesh_.nCells()/nIdealCells),
            maxOp<scalar>()
        );

        if (unbalance <= maxLoadUnbalance)
        {
            Info<< "Skipping balancing since max unbalance " << unbalance
                << " is less than allowable " << maxLoadUnbalance
                << endl;
        }
        else
        {
            scalarField cellWeights(mesh_.nCells(), 1);

            distMap = balance
            (
                false,  //keepZoneFaces
                false,  //keepBaffles
                cellWeights,
                decomposer,
                distributor
            );

            Info<< "Balanced mesh in = "
                << mesh_.time().cpuTimeIncrement() << " s" << endl;

            printMeshInfo(debug, "After balancing " + msg);


            if (debug&meshRefinement::MESH)
            {
                Pout<< "Writing balanced " << msg
                    << " mesh to time " << timeName() << endl;
                write
                (
                    debug,
                    mesh_.time().path()/timeName()
                );
                Pout<< "Dumped debug data in = "
                    << mesh_.time().cpuTimeIncrement() << " s" << endl;

                // test all is still synced across proc patches
                checkData();
            }
        }
    }

    return distMap;
}


// Do load balancing followed by refinement of consistent set of cells.
Foam::autoPtr<Foam::mapDistributePolyMesh>
Foam::meshRefinement::balanceAndRefine
(
    const string& msg,
    decompositionMethod& decomposer,
    fvMeshDistribute& distributor,
    const labelList& initCellsToRefine,
    const scalar maxLoadUnbalance
)
{
    labelList cellsToRefine(initCellsToRefine);

    //{
    //    globalIndex globalCells(mesh_.nCells());
    //
    //    Info<< "** Distribution before balancing/refining:" << endl;
    //    for (label procI = 0; procI < Pstream::nProcs(); procI++)
    //    {
    //        Info<< "    " << procI << '\t'
    //            << globalCells.localSize(procI) << endl;
    //    }
    //    Info<< endl;
    //}
    //{
    //    globalIndex globalCells(cellsToRefine.size());
    //
    //    Info<< "** Cells to be refined:" << endl;
    //    for (label procI = 0; procI < Pstream::nProcs(); procI++)
    //    {
    //        Info<< "    " << procI << '\t'
    //            << globalCells.localSize(procI) << endl;
    //    }
    //    Info<< endl;
    //}


    // Load balancing
    // ~~~~~~~~~~~~~~

    autoPtr<mapDistributePolyMesh> distMap;

    if (Pstream::nProcs() > 1)
    {
        // First check if we need to balance at all. Precalculate number of
        // cells after refinement and see what maximum difference is.
        scalar nNewCells = scalar(mesh_.nCells() + 7*cellsToRefine.size());
        scalar nIdealNewCells =
            returnReduce(nNewCells, sumOp<scalar>())
          / Pstream::nProcs();
        scalar unbalance = returnReduce
        (
            mag(1.0-nNewCells/nIdealNewCells),
            maxOp<scalar>()
        );

        if (unbalance <= maxLoadUnbalance)
        {
            Info<< "Skipping balancing since max unbalance " << unbalance
                << " is less than allowable " << maxLoadUnbalance
                << endl;
        }
        else
        {
            scalarField cellWeights(mesh_.nCells(), 1);
            forAll(cellsToRefine, i)
            {
                cellWeights[cellsToRefine[i]] += 7;
            }

            distMap = balance
            (
                false,  //keepZoneFaces
                false,  //keepBaffles
                cellWeights,
                decomposer,
                distributor
            );

            // Update cells to refine
            distMap().distributeCellIndices(cellsToRefine);

            Info<< "Balanced mesh in = "
                << mesh_.time().cpuTimeIncrement() << " s" << endl;
        }

        //{
        //    globalIndex globalCells(mesh_.nCells());
        //
        //    Info<< "** Distribution after balancing:" << endl;
        //    for (label procI = 0; procI < Pstream::nProcs(); procI++)
        //    {
        //        Info<< "    " << procI << '\t'
        //            << globalCells.localSize(procI) << endl;
        //    }
        //    Info<< endl;
        //}

        printMeshInfo(debug, "After balancing " + msg);

        if (debug&meshRefinement::MESH)
        {
            Pout<< "Writing balanced " << msg
                << " mesh to time " << timeName() << endl;
            write
            (
                debug,
                mesh_.time().path()/timeName()
            );
            Pout<< "Dumped debug data in = "
                << mesh_.time().cpuTimeIncrement() << " s" << endl;

            // test all is still synced across proc patches
            checkData();
        }
    }


    // Refinement
    // ~~~~~~~~~~

    refine(cellsToRefine);

    if (debug&meshRefinement::MESH)
    {
        Pout<< "Writing refined " << msg
            << " mesh to time " << timeName() << endl;
        write
        (
            debug,
            mesh_.time().path()
           /timeName()
        );
        Pout<< "Dumped debug data in = "
            << mesh_.time().cpuTimeIncrement() << " s" << endl;

        // test all is still synced across proc patches
        checkData();
    }

    Info<< "Refined mesh in = "
        << mesh_.time().cpuTimeIncrement() << " s" << endl;

    //{
    //    globalIndex globalCells(mesh_.nCells());
    //
    //    Info<< "** After refinement distribution:" << endl;
    //    for (label procI = 0; procI < Pstream::nProcs(); procI++)
    //    {
    //        Info<< "    " << procI << '\t'
    //            << globalCells.localSize(procI) << endl;
    //    }
    //    Info<< endl;
    //}

    printMeshInfo(debug, "After refinement " + msg);

    return distMap;
}


// ************************************************************************* //
