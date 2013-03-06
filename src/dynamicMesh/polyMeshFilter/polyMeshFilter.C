/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

#include "polyMeshFilter.H"
#include "polyMesh.H"
#include "fvMesh.H"
#include "unitConversion.H"
#include "edgeCollapser.H"
#include "syncTools.H"
#include "polyTopoChange.H"
#include "globalIndex.H"
#include "PackedBoolList.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(polyMeshFilter, 0);
}


Foam::autoPtr<Foam::fvMesh> Foam::polyMeshFilter::copyMesh(const fvMesh& mesh)
{
    polyTopoChange originalMeshToNewMesh(mesh);

    autoPtr<fvMesh> meshCopy;
    autoPtr<mapPolyMesh> mapPtr = originalMeshToNewMesh.makeMesh
    (
        meshCopy,
        IOobject
        (
            mesh.name(),
            mesh.polyMesh::instance(),
            mesh.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        true // parallel sync
    );

    const mapPolyMesh& map = mapPtr();

    // Update fields
    meshCopy().updateMesh(map);
    if (map.hasMotionPoints())
    {
        meshCopy().movePoints(map.preMotionPoints());
    }

    return meshCopy;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::polyMeshFilter::updatePointErrorCount
(
    const PackedBoolList& isErrorPoint,
    const labelList& oldToNewMesh,
    labelList& pointErrorCount
) const
{
    forAll(mesh_.points(), pI)
    {
        if (isErrorPoint[oldToNewMesh[pI]])
        {
            pointErrorCount[pI]++;
        }
    }
}


void Foam::polyMeshFilter::checkMeshEdgesAndRelaxEdges
(
    const polyMesh& newMesh,
    const labelList& oldToNewMesh,
    const PackedBoolList& isErrorPoint,
    const labelList& pointErrorCount
)
{
    const edgeList& edges = mesh_.edges();

    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];
        label newStart = oldToNewMesh[e[0]];
        label newEnd = oldToNewMesh[e[1]];

        if
        (
            pointErrorCount[e[0]] >= maxPointErrorCount_
         || pointErrorCount[e[1]] >= maxPointErrorCount_
        )
        {
            minEdgeLen_[edgeI] = 0;
        }

        if
        (
            (newStart >= 0 && isErrorPoint[newStart])
         || (newEnd >= 0 && isErrorPoint[newEnd])
        )
        {
            minEdgeLen_[edgeI] *= edgeReductionFactor_;
        }
    }

    syncTools::syncEdgeList(mesh_, minEdgeLen_, minEqOp<scalar>(), scalar(0.0));

    for (label smoothIter = 0; smoothIter < maxSmoothIters_; ++smoothIter)
    {
        // Smooth minEdgeLen
        forAll(mesh_.edges(), edgeI)
        {
            const edge& e = mesh_.edges()[edgeI];

            scalar sumMinEdgeLen = 0;
            label nEdges = 0;

            forAll(e, pointI)
            {
                const labelList& pEdges = mesh_.pointEdges()[e[pointI]];

                forAll(pEdges, pEdgeI)
                {
                    const label pEdge = pEdges[pEdgeI];
                    sumMinEdgeLen += minEdgeLen_[pEdge];
                    nEdges++;
                }
            }

            minEdgeLen_[edgeI] = min
            (
                minEdgeLen_[edgeI],
                sumMinEdgeLen/nEdges
            );
        }

        syncTools::syncEdgeList
        (
            mesh_,
            minEdgeLen_,
            minEqOp<scalar>(),
            scalar(0.0)
        );
    }
}


void Foam::polyMeshFilter::checkMeshFacesAndRelaxEdges
(
    const polyMesh& newMesh,
    const labelList& oldToNewMesh,
    const PackedBoolList& isErrorPoint,
    const labelList& pointErrorCount
)
{
    const faceList& faces = mesh_.faces();

    forAll(faces, faceI)
    {
        const face& f = faces[faceI];

        forAll(f, fpI)
        {
            const label ptIndex = oldToNewMesh[f[fpI]];

            if (pointErrorCount[f[fpI]] >= maxPointErrorCount_)
            {
                faceFilterFactor_[faceI] = 0;
            }

            if (isErrorPoint[ptIndex])
            {
                faceFilterFactor_[faceI] *= faceReductionFactor_;

                break;
            }
        }
    }

    syncTools::syncFaceList(mesh_, faceFilterFactor_, minEqOp<scalar>());

    for (label smoothIter = 0; smoothIter < maxSmoothIters_; ++smoothIter)
    {
        // Smooth faceFilterFactor
        forAll(faces, faceI)
        {
            const labelList& fEdges = mesh_.faceEdges()[faceI];

            scalar sumFaceFilterFactors = 0;
            label nFaces = 0;

            // This is important: Only smooth around faces that share an
            // edge with a bad face
            bool skipFace = true;

            forAll(fEdges, fEdgeI)
            {
                const labelList& eFaces = mesh_.edgeFaces()[fEdges[fEdgeI]];

                forAll(eFaces, eFaceI)
                {
                    const label eFace = eFaces[eFaceI];

                    const face& f = faces[eFace];

                    forAll(f, fpI)
                    {
                        const label ptIndex = oldToNewMesh[f[fpI]];

                        if (isErrorPoint[ptIndex])
                        {
                            skipFace = false;
                            break;
                        }
                    }

                    if (eFace != faceI)
                    {
                        sumFaceFilterFactors += faceFilterFactor_[eFace];
                        nFaces++;
                    }
                }
            }

            if (skipFace)
            {
                continue;
            }

            faceFilterFactor_[faceI] = min
            (
                faceFilterFactor_[faceI],
                sumFaceFilterFactors/nFaces
            );
        }

        // Face filter factor needs to be synchronised!
        syncTools::syncFaceList(mesh_, faceFilterFactor_, minEqOp<scalar>());
    }
}


Foam::labelList Foam::polyMeshFilter::findBoundaryPoints
(
    const polyMesh& mesh//,
//    labelIOList& boundaryIOPts
) const
{
    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

    labelList boundaryPoint(mesh.nPoints(), -1);

    // Get all processor boundary points and the processor patch label
    // that they are on.
    forAll(bMesh, patchI)
    {
        const polyPatch& patch = bMesh[patchI];

        if (!isA<coupledPolyPatch>(patch))
        {
            forAll(patch, fI)
            {
                const face& f = patch[fI];

                forAll(f, fp)
                {
                    boundaryPoint[f[fp]] = 0; //boundaryIOPts[f[fp]];
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        boundaryPoint,
        maxEqOp<label>(),
        labelMin
    );

    return boundaryPoint;
}


void Foam::polyMeshFilter::printScalarFieldStats
(
    const string desc,
    const scalarField& fld
) const
{
    Info<< incrIndent << indent << desc
        << ": min = " << returnReduce(min(fld), minOp<scalar>())
        << " av = "
        << returnReduce(sum(fld), sumOp<scalar>())
          /returnReduce(fld.size(), sumOp<label>())
        << " max = " << returnReduce(max(fld), maxOp<scalar>())
        << decrIndent << endl;
}


void Foam::polyMeshFilter::mapOldMeshEdgeFieldToNewMesh
(
    const polyMesh& newMesh,
    const labelList& pointMap,
    scalarField& newMeshMinEdgeLen
) const
{
    scalarField tmp(newMesh.nEdges());

    const edgeList& newEdges = newMesh.edges();

    forAll(newEdges, newEdgeI)
    {
        const edge& newEdge = newEdges[newEdgeI];
        const label pStart = newEdge.start();
        const label pEnd = newEdge.end();

        tmp[newEdgeI] = min
        (
            newMeshMinEdgeLen[pointMap[pStart]],
            newMeshMinEdgeLen[pointMap[pEnd]]
        );
    }

    newMeshMinEdgeLen.transfer(tmp);

    syncTools::syncEdgeList
    (
        newMesh,
        newMeshMinEdgeLen,
        maxEqOp<scalar>(),
        scalar(0.0)
    );
}


void Foam::polyMeshFilter::mapOldMeshFaceFieldToNewMesh
(
    const polyMesh& newMesh,
    const labelList& faceMap,
    scalarField& newMeshFaceFilterFactor
) const
{
    scalarField tmp(newMesh.nFaces());

    forAll(faceMap, newFaceI)
    {
        const label oldFaceI = faceMap[newFaceI];

        tmp[newFaceI] = newMeshFaceFilterFactor[oldFaceI];
    }

    newMeshFaceFilterFactor.transfer(tmp);

    syncTools::syncFaceList
    (
        newMesh,
        newMeshFaceFilterFactor,
        maxEqOp<scalar>()
    );
}


void Foam::polyMeshFilter::updateOldToNewPointMap
(
    const labelList& currToNew,
    labelList& origToCurrentPointMap
) const
{
    forAll(origToCurrentPointMap, origPointI)
    {
        label oldPointI = origToCurrentPointMap[origPointI];

        if (oldPointI < currToNew.size())
        {
            label newPointI = currToNew[oldPointI];

            if (newPointI != -1)
            {
                origToCurrentPointMap[origPointI] = newPointI;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyMeshFilter::polyMeshFilter(const fvMesh& mesh)
:
    mesh_(mesh),
    newMeshPtr_(),
    dict_
    (
        IOobject
        (
            "collapseDict",
            mesh.time().system(),
            mesh.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    controlMeshQuality_
    (
        dict_.lookupOrDefault<Switch>("controlMeshQuality", false)
    ),
    collapseEdgesCoeffDict_(dict_.subDict("collapseEdgesCoeffs")),
    collapseFacesCoeffDict_(dict_.subOrEmptyDict("collapseFacesCoeffs")),
    meshQualityCoeffDict_(dict_.subOrEmptyDict("controlMeshQualityCoeffs")),
    minLen_(readScalar(collapseEdgesCoeffDict_.lookup("minimumEdgeLength"))),
    maxCos_
    (
        ::cos
        (
            degToRad
            (
                readScalar(collapseEdgesCoeffDict_.lookup("maximumMergeAngle"))
            )
        )
    ),
    edgeReductionFactor_
    (
        meshQualityCoeffDict_.lookupOrDefault<scalar>("edgeReductionFactor", -1)
    ),
    maxIterations_
    (
        meshQualityCoeffDict_.lookupOrAddDefault<label>("maximumIterations", 1)
    ),
    maxSmoothIters_
    (
        meshQualityCoeffDict_.lookupOrAddDefault<label>
        (
            "maximumSmoothingIterations",
            0
        )
    ),
    initialFaceLengthFactor_
    (
        collapseFacesCoeffDict_.lookupOrAddDefault<scalar>
        (
            "initialFaceLengthFactor",
            -1
        )
    ),
    faceReductionFactor_
    (
        meshQualityCoeffDict_.lookupOrAddDefault<scalar>
        (
            "faceReductionFactor",
            -1
        )
    ),
    maxPointErrorCount_
    (
        meshQualityCoeffDict_.lookupOrAddDefault<label>("maxPointErrorCount", 0)
    ),
    minEdgeLen_(),
    faceFilterFactor_()
{
    Info<< "Merging:" << nl
        << "    edges with length less than " << minLen_ << " meters" << nl
        << "    edges split by a point with edges in line to within "
        << radToDeg(::acos(maxCos_)) << " degrees" << nl
        << "    Minimum edge length reduction factor = "
        << edgeReductionFactor_ << nl
        << endl;

    Info<< "Collapse faces with reduction factor = " << faceReductionFactor_
        << endl;

    Info<< "Selectively disabling wanted collapses until resulting quality"
        << " satisfies constraints in system/meshQualityDict" << nl
        << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::polyMeshFilter::~polyMeshFilter()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::label Foam::polyMeshFilter::filter(const label nOriginalBadFaces)
{
    label nBadFaces = labelMax;
    label nOuterIterations = 0;

    minEdgeLen_.resize(mesh_.nEdges(), minLen_);
    faceFilterFactor_.resize(mesh_.nFaces(), initialFaceLengthFactor_);

    // Maintain the number of times a point has been part of a bad face
    labelList pointErrorCount(mesh_.nPoints(), 0);

    // Main loop
    // ~~~~~~~~~
    // It tries and do some collapses, checks the resulting mesh and
    // 'freezes' some edges (by marking in minEdgeLen) and tries again.
    // This will iterate ultimately to the situation where every edge is
    // frozen and nothing gets collapsed.
    while
    (
        nOuterIterations < maxIterations_
     && nBadFaces > nOriginalBadFaces
    )
    {
        Info<< nl << "Outer Iteration = " << nOuterIterations++ << nl
            << endl;

        printScalarFieldStats("Edge Filter Factor", minEdgeLen_);
        printScalarFieldStats("Face Filter Factor", faceFilterFactor_);

        // Reset the new mesh to the old mesh
        newMeshPtr_ = copyMesh(mesh_);
        fvMesh& newMesh = newMeshPtr_();

        scalarField newMeshFaceFilterFactor = faceFilterFactor_;

        labelList origToCurrentPointMap(identity(newMesh.nPoints()));
        {

            label nInnerIterations = 0;
            label nPrevLocalCollapse = labelMax;

            Info<< incrIndent;

            while (true)
            {
                Info<< nl << indent << "Inner iteration = "
                    << nInnerIterations++ << nl << incrIndent << endl;

                // Per edge collapse status
                PackedBoolList collapseEdge(newMesh.nEdges());

                Map<point> collapsePointToLocation(newMesh.nPoints());

                // Mark points on boundary
                const labelList boundaryPoint = findBoundaryPoints(newMesh);

                edgeCollapser collapser(newMesh, collapseFacesCoeffDict_);

                {
                    // Collapse faces
                    labelPair nCollapsedPtEdge = collapser.markSmallSliverFaces
                    (
                        newMeshFaceFilterFactor,
                        boundaryPoint,
                        collapseEdge,
                        collapsePointToLocation
                    );

                    label nCollapsed = 0;
                    forAll(nCollapsedPtEdge, collapseTypeI)
                    {
                        nCollapsed += nCollapsedPtEdge[collapseTypeI];
                    }

                    reduce(nCollapsed, sumOp<label>());

                    Info<< indent
                        << "Collapsing " << nCollapsed << " faces"
                        << " (to point = "
                        << returnReduce
                           (
                               nCollapsedPtEdge.first(),
                               sumOp<label>()
                           )
                        << ", to edge = "
                        << returnReduce
                           (
                               nCollapsedPtEdge.second(),
                               sumOp<label>()
                           )
                        << ")" << endl;

                    if (nCollapsed == 0)
                    {
                        Info<< decrIndent;
                        Info<< decrIndent;
                        break;
                    }
                }

                // Merge edge collapses into consistent collapse-network.
                // Make sure no cells get collapsed.
                List<pointEdgeCollapse> allPointInfo;
                const globalIndex globalPoints(newMesh.nPoints());

                collapser.consistentCollapse
                (
                    globalPoints,
                    boundaryPoint,
                    collapsePointToLocation,
                    collapseEdge,
                    allPointInfo
                );

                label nLocalCollapse = collapseEdge.count();

                reduce(nLocalCollapse, sumOp<label>());
                Info<< nl << indent << "Collapsing " << nLocalCollapse
                    << " edges after synchronisation and PointEdgeWave" << endl;

                if (nLocalCollapse >= nPrevLocalCollapse)
                {
                    Info<< decrIndent;
                    Info<< decrIndent;
                    break;
                }
                else
                {
                    nPrevLocalCollapse = nLocalCollapse;
                }

                {
                    // Apply collapses to current mesh
                    polyTopoChange newMeshMod(newMesh);

                    // Insert mesh refinement into polyTopoChange.
                    collapser.setRefinement(allPointInfo, newMeshMod);

                    Info<< indent << "Apply changes to the current mesh"
                        << decrIndent << endl;

                    // Apply changes to current mesh
                    autoPtr<mapPolyMesh> newMapPtr = newMeshMod.changeMesh
                    (
                        newMesh,
                        false
                    );
                    const mapPolyMesh& newMap = newMapPtr();

                    // Update fields
                    newMesh.updateMesh(newMap);
                    if (newMap.hasMotionPoints())
                    {
                        newMesh.movePoints(newMap.preMotionPoints());
                    }

                    // Relabel the boundary points
        //            labelList newBoundaryPoints(newMesh.nPoints(), -1);
        //            forAll(newBoundaryPoints, pI)
        //            {
        //                const label newToOldptI= map.pointMap()[pI];
        //                newBoundaryPoints[pI] = boundaryIOPts[newToOldptI];
        //            }
        //            boundaryIOPts = newBoundaryPoints;


                    mapOldMeshFaceFieldToNewMesh
                    (
                        newMesh,
                        newMap.faceMap(),
                        newMeshFaceFilterFactor
                    );

                    updateOldToNewPointMap
                    (
                        newMap.reversePointMap(),
                        origToCurrentPointMap
                    );
                }
            }
        }


        scalarField newMeshMinEdgeLen = minEdgeLen_;

        label nInnerIterations = 0;
        label nPrevLocalCollapse = labelMax;

        while (true)
        {
            Info<< nl << indent << "Inner iteration = "
                << nInnerIterations++ << nl << incrIndent << endl;

            // Per edge collapse status
            PackedBoolList collapseEdge(newMesh.nEdges());

            Map<point> collapsePointToLocation(newMesh.nPoints());

            // Mark points on boundary
            const labelList boundaryPoint = findBoundaryPoints
            (
                newMesh//,
//                    boundaryIOPts
            );

            edgeCollapser collapser(newMesh, collapseFacesCoeffDict_);

            // Work out which edges to collapse
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // This is by looking at minEdgeLen (to avoid frozen edges)
            // and marking in collapseEdge.
            label nSmallCollapsed = collapser.markSmallEdges
            (
                newMeshMinEdgeLen,
                boundaryPoint,
                collapseEdge,
                collapsePointToLocation
            );

            reduce(nSmallCollapsed, sumOp<label>());
            Info<< indent << "Collapsing " << nSmallCollapsed
                << " small edges" << endl;

            // Merge inline edges
            label nMerged = collapser.markMergeEdges
            (
                maxCos_,
                boundaryPoint,
                collapseEdge,
                collapsePointToLocation
            );

            reduce(nMerged, sumOp<label>());
            Info<< indent << "Collapsing " << nMerged << " in line edges"
                << endl;

            if (nMerged + nSmallCollapsed == 0)
            {
                Info<< decrIndent;
                break;
            }

            // Merge edge collapses into consistent collapse-network.
            // Make sure no cells get collapsed.
            List<pointEdgeCollapse> allPointInfo;
            const globalIndex globalPoints(newMesh.nPoints());

            collapser.consistentCollapse
            (
                globalPoints,
                boundaryPoint,
                collapsePointToLocation,
                collapseEdge,
                allPointInfo
            );

            label nLocalCollapse = collapseEdge.count();

            reduce(nLocalCollapse, sumOp<label>());
            Info<< nl << indent << "Collapsing " << nLocalCollapse
                << " edges after synchronisation and PointEdgeWave" << endl;

            if (nLocalCollapse >= nPrevLocalCollapse)
            {
                Info<< decrIndent;
                break;
            }
            else
            {
                nPrevLocalCollapse = nLocalCollapse;
            }

            // Apply collapses to current mesh
            polyTopoChange newMeshMod(newMesh);

            // Insert mesh refinement into polyTopoChange.
            collapser.setRefinement(allPointInfo, newMeshMod);

            Info<< indent << "Apply changes to the current mesh"
                << decrIndent << endl;

            // Apply changes to current mesh
            autoPtr<mapPolyMesh> newMapPtr = newMeshMod.changeMesh
            (
                newMesh,
                false
            );
            const mapPolyMesh& newMap = newMapPtr();

            // Update fields
            newMesh.updateMesh(newMap);
            if (newMap.hasMotionPoints())
            {
                newMesh.movePoints(newMap.preMotionPoints());
            }

            // Relabel the boundary points
//            labelList newBoundaryPoints(newMesh.nPoints(), -1);
//            forAll(newBoundaryPoints, pI)
//            {
//                const label newToOldptI= map.pointMap()[pI];
//                newBoundaryPoints[pI] = boundaryIOPts[newToOldptI];
//            }
//            boundaryIOPts = newBoundaryPoints;

            // Synchronise the factors
            mapOldMeshEdgeFieldToNewMesh
            (
                newMesh,
                newMap.pointMap(),
                newMeshMinEdgeLen
            );

            updateOldToNewPointMap
            (
                newMap.reversePointMap(),
                origToCurrentPointMap
            );
        }


        // Mesh check
        // ~~~~~~~~~~~~~~~~~~
        // Do not allow collapses in regions of error.
        // Updates minEdgeLen, nRelaxedEdges

        if (controlMeshQuality_)
        {
            PackedBoolList isErrorPoint(newMesh.nPoints());
            nBadFaces = edgeCollapser::checkMeshQuality
            (
                newMesh,
                meshQualityCoeffDict_,
                isErrorPoint
            );

            Info<< nl << "    Number of bad faces     : " << nBadFaces << nl
                << "    Number of marked points : "
                << returnReduce(isErrorPoint.count(), sumOp<unsigned int>())
                << endl;

            updatePointErrorCount
            (
                isErrorPoint,
                origToCurrentPointMap,
                pointErrorCount
            );

            checkMeshEdgesAndRelaxEdges
            (
                newMesh,
                origToCurrentPointMap,
                isErrorPoint,
                pointErrorCount
            );

            checkMeshFacesAndRelaxEdges
            (
                newMesh,
                origToCurrentPointMap,
                isErrorPoint,
                pointErrorCount
            );
        }
        else
        {
            return -1;
        }
    }

    return nBadFaces;
}


Foam::label Foam::polyMeshFilter::filterEdges
(
    const label nOriginalBadFaces
)
{
    label nBadFaces = labelMax/2;
    label nPreviousBadFaces = labelMax;
    label nOuterIterations = 0;

    minEdgeLen_.resize(mesh_.nEdges(), minLen_);
    faceFilterFactor_.resize(0);

    labelList pointErrorCount(mesh_.nPoints(), 0);

    // Main loop
    // ~~~~~~~~~
    // It tries and do some collapses, checks the resulting mesh and
    // 'freezes' some edges (by marking in minEdgeLen) and tries again.
    // This will iterate ultimately to the situation where every edge is
    // frozen and nothing gets collapsed.
    while
    (
        nOuterIterations < maxIterations_
     && nBadFaces > nOriginalBadFaces
     && nBadFaces < nPreviousBadFaces
    )
    {
        Info<< nl << "Outer Iteration = " << nOuterIterations++ << nl
            << endl;

        printScalarFieldStats("Edge Filter Factor", minEdgeLen_);

        nPreviousBadFaces = nBadFaces;

        // Reset the new mesh to the old mesh
        newMeshPtr_ = copyMesh(mesh_);
        fvMesh& newMesh = newMeshPtr_();

        scalarField newMeshMinEdgeLen = minEdgeLen_;

        labelList origToCurrentPointMap(identity(newMesh.nPoints()));

        label nInnerIterations = 0;
        label nPrevLocalCollapse = labelMax;

        Info<< incrIndent;

        while (true)
        {
            Info<< nl << indent << "Inner iteration = "
                << nInnerIterations++ << nl << incrIndent << endl;

            // Per edge collapse status
            PackedBoolList collapseEdge(newMesh.nEdges());

            Map<point> collapsePointToLocation(newMesh.nPoints());

            // Mark points on boundary
            const labelList boundaryPoint = findBoundaryPoints(newMesh);

            edgeCollapser collapser(newMesh, collapseFacesCoeffDict_);

            // Work out which edges to collapse
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // This is by looking at minEdgeLen (to avoid frozen edges)
            // and marking in collapseEdge.
            label nSmallCollapsed = collapser.markSmallEdges
            (
                newMeshMinEdgeLen,
                boundaryPoint,
                collapseEdge,
                collapsePointToLocation
            );

            reduce(nSmallCollapsed, sumOp<label>());
            Info<< indent << "Collapsing " << nSmallCollapsed
                << " small edges" << endl;

            // Merge inline edges
            label nMerged = collapser.markMergeEdges
            (
                maxCos_,
                boundaryPoint,
                collapseEdge,
                collapsePointToLocation
            );

            reduce(nMerged, sumOp<label>());
            Info<< indent << "Collapsing " << nMerged << " in line edges"
                << endl;

            if (nMerged + nSmallCollapsed == 0)
            {
                Info<< decrIndent;
                Info<< decrIndent;
                break;
            }

            // Merge edge collapses into consistent collapse-network.
            // Make sure no cells get collapsed.
            List<pointEdgeCollapse> allPointInfo;
            const globalIndex globalPoints(newMesh.nPoints());

            collapser.consistentCollapse
            (
                globalPoints,
                boundaryPoint,
                collapsePointToLocation,
                collapseEdge,
                allPointInfo
            );

            label nLocalCollapse = collapseEdge.count();

            reduce(nLocalCollapse, sumOp<label>());
            Info<< nl << indent << "Collapsing " << nLocalCollapse
                << " edges after synchronisation and PointEdgeWave" << endl;

            if (nLocalCollapse >= nPrevLocalCollapse)
            {
                Info<< decrIndent;
                Info<< decrIndent;
                break;
            }
            else
            {
                nPrevLocalCollapse = nLocalCollapse;
            }

            // Apply collapses to current mesh
            polyTopoChange newMeshMod(newMesh);

            // Insert mesh refinement into polyTopoChange.
            collapser.setRefinement(allPointInfo, newMeshMod);

            Info<< indent << "Apply changes to the current mesh"
                << decrIndent << endl;

            // Apply changes to current mesh
            autoPtr<mapPolyMesh> newMapPtr = newMeshMod.changeMesh
            (
                newMesh,
                false
            );
            const mapPolyMesh& newMap = newMapPtr();

            // Update fields
            newMesh.updateMesh(newMap);
            if (newMap.hasMotionPoints())
            {
                newMesh.movePoints(newMap.preMotionPoints());
            }

            // Relabel the boundary points
//            labelList newBoundaryPoints(newMesh.nPoints(), -1);
//            forAll(newBoundaryPoints, pI)
//            {
//                const label newToOldptI= map.pointMap()[pI];
//                newBoundaryPoints[pI] = boundaryIOPts[newToOldptI];
//            }
//            boundaryIOPts = newBoundaryPoints;

            // Synchronise the factors
            mapOldMeshEdgeFieldToNewMesh
            (
                newMesh,
                newMap.pointMap(),
                newMeshMinEdgeLen
            );

            updateOldToNewPointMap
            (
                newMap.reversePointMap(),
                origToCurrentPointMap
            );
        }

        // Mesh check
        // ~~~~~~~~~~~~~~~~~~
        // Do not allow collapses in regions of error.
        // Updates minEdgeLen, nRelaxedEdges

        if (controlMeshQuality_)
        {
            PackedBoolList isErrorPoint(newMesh.nPoints());
            nBadFaces = edgeCollapser::checkMeshQuality
            (
                newMesh,
                meshQualityCoeffDict_,
                isErrorPoint
            );

            Info<< nl << "    Number of bad faces     : " << nBadFaces << nl
                << "    Number of marked points : "
                << returnReduce(isErrorPoint.count(), sumOp<unsigned int>())
                << endl;

            updatePointErrorCount
            (
                isErrorPoint,
                origToCurrentPointMap,
                pointErrorCount
            );

            checkMeshEdgesAndRelaxEdges
            (
                newMesh,
                origToCurrentPointMap,
                isErrorPoint,
                pointErrorCount
            );
        }
        else
        {
            return -1;
        }
    }

    return nBadFaces;
}


Foam::label Foam::polyMeshFilter::filterIndirectPatchFaces()
{
    newMeshPtr_ = copyMesh(mesh_);
    fvMesh& newMesh = newMeshPtr_();

    label nIterations = 0;
    label nBadFaces = 0;

    while (true)
    {
        Info<< nl << indent << "Iteration = "
            << nIterations++ << nl << incrIndent << endl;

        // Per edge collapse status
        PackedBoolList collapseEdge(newMesh.nEdges());

        Map<point> collapsePointToLocation(newMesh.nPoints());

        labelList boundaryPoint(newMesh.nPoints());

        edgeCollapser collapser(newMesh, collapseFacesCoeffDict_);

        collapser.markIndirectPatchFaces
        (
            collapseEdge,
            collapsePointToLocation
        );

        // Merge edge collapses into consistent collapse-network.
        // Make sure no cells get collapsed.
        List<pointEdgeCollapse> allPointInfo;
        const globalIndex globalPoints(newMesh.nPoints());

        collapser.consistentCollapse
        (
            globalPoints,
            boundaryPoint,
            collapsePointToLocation,
            collapseEdge,
            allPointInfo
        );

        label nLocalCollapse = collapseEdge.count();

        reduce(nLocalCollapse, sumOp<label>());
        Info<< nl << indent << "Collapsing " << nLocalCollapse
            << " edges after synchronisation and PointEdgeWave" << endl;

        if (nLocalCollapse == 0)
        {
            break;
        }

        // Apply collapses to current mesh
        polyTopoChange newMeshMod(newMesh);

        // Insert mesh refinement into polyTopoChange.
        collapser.setRefinement(allPointInfo, newMeshMod);

        Info<< indent << "Apply changes to the current mesh"
            << decrIndent << endl;

        // Apply changes to current mesh
        autoPtr<mapPolyMesh> newMapPtr = newMeshMod.changeMesh
        (
            newMesh,
            false
        );
        const mapPolyMesh& newMap = newMapPtr();

        // Update fields
        newMesh.updateMesh(newMap);
        if (newMap.hasMotionPoints())
        {
            newMesh.movePoints(newMap.preMotionPoints());
        }

        // Mesh check
        // ~~~~~~~~~~~~~~~~~~
        // Do not allow collapses in regions of error.
        // Updates minEdgeLen, nRelaxedEdges

        if (controlMeshQuality_)
        {
            PackedBoolList isErrorPoint(newMesh.nPoints());
            nBadFaces = edgeCollapser::checkMeshQuality
            (
                newMesh,
                meshQualityCoeffDict_,
                isErrorPoint
            );

            Info<< nl << "    Number of bad faces     : " << nBadFaces << nl
                << "    Number of marked points : "
                << returnReduce(isErrorPoint.count(), sumOp<unsigned int>())
                << endl;
        }
    }

    return nBadFaces;
}


const Foam::autoPtr<Foam::fvMesh>& Foam::polyMeshFilter::filteredMesh() const
{
    return newMeshPtr_;
}


// ************************************************************************* //
