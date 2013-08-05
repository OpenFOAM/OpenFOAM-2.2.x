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

#include "autoRefineDriver.H"
#include "meshRefinement.H"
#include "fvMesh.H"
#include "Time.H"
#include "cellSet.H"
#include "syncTools.H"
#include "refinementParameters.H"
#include "refinementSurfaces.H"
#include "refinementFeatures.H"
#include "shellSurfaces.H"
#include "mapDistributePolyMesh.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(autoRefineDriver, 0);

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::autoRefineDriver::autoRefineDriver
(
    meshRefinement& meshRefiner,
    decompositionMethod& decomposer,
    fvMeshDistribute& distributor,
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch
)
:
    meshRefiner_(meshRefiner),
    decomposer_(decomposer),
    distributor_(distributor),
    globalToMasterPatch_(globalToMasterPatch),
    globalToSlavePatch_(globalToSlavePatch)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::autoRefineDriver::featureEdgeRefine
(
    const refinementParameters& refineParams,
    const label maxIter,
    const label minRefine
)
{
    const fvMesh& mesh = meshRefiner_.mesh();

    label iter = 0;

    if (meshRefiner_.features().size() && maxIter > 0)
    {
        for (; iter < maxIter; iter++)
        {
            Info<< nl
                << "Feature refinement iteration " << iter << nl
                << "------------------------------" << nl
                << endl;

            labelList candidateCells
            (
                meshRefiner_.refineCandidates
                (
                    refineParams.keepPoints()[0],    // For now only use one.
                    refineParams.curvature(),

                    true,               // featureRefinement
                    false,              // featureDistanceRefinement
                    false,              // internalRefinement
                    false,              // surfaceRefinement
                    false,              // curvatureRefinement
                    refineParams.maxGlobalCells(),
                    refineParams.maxLocalCells()
                )
            );
            labelList cellsToRefine
            (
                meshRefiner_.meshCutter().consistentRefinement
                (
                    candidateCells,
                    true
                )
            );
            Info<< "Determined cells to refine in = "
                << mesh.time().cpuTimeIncrement() << " s" << endl;



            label nCellsToRefine = cellsToRefine.size();
            reduce(nCellsToRefine, sumOp<label>());

            Info<< "Selected for feature refinement : " << nCellsToRefine
                << " cells (out of " << mesh.globalData().nTotalCells()
                << ')' << endl;

            if (nCellsToRefine <= minRefine)
            {
                Info<< "Stopping refining since too few cells selected."
                    << nl << endl;
                break;
            }


            if (debug > 0)
            {
                const_cast<Time&>(mesh.time())++;
            }


            if
            (
                returnReduce
                (
                    (mesh.nCells() >= refineParams.maxLocalCells()),
                    orOp<bool>()
                )
            )
            {
                meshRefiner_.balanceAndRefine
                (
                    "feature refinement iteration " + name(iter),
                    decomposer_,
                    distributor_,
                    cellsToRefine,
                    refineParams.maxLoadUnbalance()
                );
            }
            else
            {
                meshRefiner_.refineAndBalance
                (
                    "feature refinement iteration " + name(iter),
                    decomposer_,
                    distributor_,
                    cellsToRefine,
                    refineParams.maxLoadUnbalance()
                );
            }
        }
    }
    return iter;
}


Foam::label Foam::autoRefineDriver::surfaceOnlyRefine
(
    const refinementParameters& refineParams,
    const label maxIter
)
{
    const fvMesh& mesh = meshRefiner_.mesh();

    // Determine the maximum refinement level over all surfaces. This
    // determines the minumum number of surface refinement iterations.
    label overallMaxLevel = max(meshRefiner_.surfaces().maxLevel());

    label iter;
    for (iter = 0; iter < maxIter; iter++)
    {
        Info<< nl
            << "Surface refinement iteration " << iter << nl
            << "------------------------------" << nl
            << endl;


        // Determine cells to refine
        // ~~~~~~~~~~~~~~~~~~~~~~~~~
        // Only look at surface intersections (minLevel and surface curvature),
        // do not do internal refinement (refinementShells)

        labelList candidateCells
        (
            meshRefiner_.refineCandidates
            (
                refineParams.keepPoints()[0],
                refineParams.curvature(),

                false,              // featureRefinement
                false,              // featureDistanceRefinement
                false,              // internalRefinement
                true,               // surfaceRefinement
                true,               // curvatureRefinement
                refineParams.maxGlobalCells(),
                refineParams.maxLocalCells()
            )
        );
        labelList cellsToRefine
        (
            meshRefiner_.meshCutter().consistentRefinement
            (
                candidateCells,
                true
            )
        );
        Info<< "Determined cells to refine in = "
            << mesh.time().cpuTimeIncrement() << " s" << endl;


        label nCellsToRefine = cellsToRefine.size();
        reduce(nCellsToRefine, sumOp<label>());

        Info<< "Selected for refinement : " << nCellsToRefine
            << " cells (out of " << mesh.globalData().nTotalCells()
            << ')' << endl;

        // Stop when no cells to refine or have done minimum necessary
        // iterations and not enough cells to refine.
        if
        (
            nCellsToRefine == 0
         || (
                iter >= overallMaxLevel
             && nCellsToRefine <= refineParams.minRefineCells()
            )
        )
        {
            Info<< "Stopping refining since too few cells selected."
                << nl << endl;
            break;
        }


        if (debug)
        {
            const_cast<Time&>(mesh.time())++;
        }


        if
        (
            returnReduce
            (
                (mesh.nCells() >= refineParams.maxLocalCells()),
                orOp<bool>()
            )
        )
        {
            meshRefiner_.balanceAndRefine
            (
                "surface refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams.maxLoadUnbalance()
            );
        }
        else
        {
            meshRefiner_.refineAndBalance
            (
                "surface refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams.maxLoadUnbalance()
            );
        }
    }
    return iter;
}


void Foam::autoRefineDriver::removeInsideCells
(
    const refinementParameters& refineParams,
    const label nBufferLayers
)
{
    Info<< nl
        << "Removing mesh beyond surface intersections" << nl
        << "------------------------------------------" << nl
        << endl;

    const fvMesh& mesh = meshRefiner_.mesh();

    if (debug)
    {
       const_cast<Time&>(mesh.time())++;
    }

    meshRefiner_.splitMesh
    (
        nBufferLayers,                  // nBufferLayers
        globalToMasterPatch_,
        globalToSlavePatch_,
        refineParams.keepPoints()[0]
    );

    if (debug&meshRefinement::MESH)
    {
        Pout<< "Writing subsetted mesh to time "
            << meshRefiner_.timeName() << '.' << endl;
        meshRefiner_.write(debug, mesh.time().path()/meshRefiner_.timeName());
        Pout<< "Dumped mesh in = "
            << mesh.time().cpuTimeIncrement() << " s\n" << nl << endl;
    }
}


Foam::label Foam::autoRefineDriver::shellRefine
(
    const refinementParameters& refineParams,
    const label maxIter
)
{
    const fvMesh& mesh = meshRefiner_.mesh();

    // Mark current boundary faces with 0. Have meshRefiner maintain them.
    meshRefiner_.userFaceData().setSize(1);

    // mark list to remove any refined faces
    meshRefiner_.userFaceData()[0].first() = meshRefinement::REMOVE;
    meshRefiner_.userFaceData()[0].second() = createWithValues<labelList>
    (
        mesh.nFaces(),
        -1,
        meshRefiner_.intersectedFaces(),
        0
    );

    // Determine the maximum refinement level over all volume refinement
    // regions. This determines the minumum number of shell refinement
    // iterations.
    label overallMaxShellLevel = meshRefiner_.shells().maxLevel();

    label iter;
    for (iter = 0; iter < maxIter; iter++)
    {
        Info<< nl
            << "Shell refinement iteration " << iter << nl
            << "----------------------------" << nl
            << endl;

        labelList candidateCells
        (
            meshRefiner_.refineCandidates
            (
                refineParams.keepPoints()[0],
                refineParams.curvature(),

                false,              // featureRefinement
                true,               // featureDistanceRefinement
                true,               // internalRefinement
                false,              // surfaceRefinement
                false,              // curvatureRefinement
                refineParams.maxGlobalCells(),
                refineParams.maxLocalCells()
            )
        );

        if (debug&meshRefinement::MESH)
        {
            Pout<< "Dumping " << candidateCells.size()
                << " cells to cellSet candidateCellsFromShells." << endl;

            cellSet c(mesh, "candidateCellsFromShells", candidateCells);
            c.instance() = meshRefiner_.timeName();
            c.write();
        }

        // Problem choosing starting faces for bufferlayers (bFaces)
        //  - we can't use the current intersected boundary faces
        //    (intersectedFaces) since this grows indefinitely
        //  - if we use 0 faces we don't satisfy bufferLayers from the
        //    surface.
        //  - possibly we want to have bFaces only the initial set of faces
        //    and maintain the list while doing the refinement.
        labelList bFaces
        (
            findIndices(meshRefiner_.userFaceData()[0].second(), 0)
        );

        //Info<< "Collected boundary faces : "
        //    << returnReduce(bFaces.size(), sumOp<label>()) << endl;

        labelList cellsToRefine;

        if (refineParams.nBufferLayers() <= 2)
        {
            cellsToRefine = meshRefiner_.meshCutter().consistentSlowRefinement
            (
                refineParams.nBufferLayers(),
                candidateCells,                     // cells to refine
                bFaces,                             // faces for nBufferLayers
                1,                                  // point difference
                meshRefiner_.intersectedPoints()    // points to check
            );
        }
        else
        {
            cellsToRefine = meshRefiner_.meshCutter().consistentSlowRefinement2
            (
                refineParams.nBufferLayers(),
                candidateCells,                 // cells to refine
                bFaces                          // faces for nBufferLayers
            );
        }

        Info<< "Determined cells to refine in = "
            << mesh.time().cpuTimeIncrement() << " s" << endl;


        label nCellsToRefine = cellsToRefine.size();
        reduce(nCellsToRefine, sumOp<label>());

        Info<< "Selected for internal refinement : " << nCellsToRefine
            << " cells (out of " << mesh.globalData().nTotalCells()
            << ')' << endl;

        // Stop when no cells to refine or have done minimum necessary
        // iterations and not enough cells to refine.
        if
        (
            nCellsToRefine == 0
         || (
                iter >= overallMaxShellLevel
             && nCellsToRefine <= refineParams.minRefineCells()
            )
        )
        {
            Info<< "Stopping refining since too few cells selected."
                << nl << endl;
            break;
        }


        if (debug)
        {
            const_cast<Time&>(mesh.time())++;
        }

        if
        (
            returnReduce
            (
                (mesh.nCells() >= refineParams.maxLocalCells()),
                orOp<bool>()
            )
        )
        {
            meshRefiner_.balanceAndRefine
            (
                "shell refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams.maxLoadUnbalance()
            );
        }
        else
        {
            meshRefiner_.refineAndBalance
            (
                "shell refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams.maxLoadUnbalance()
            );
        }
    }
    meshRefiner_.userFaceData().clear();

    return iter;
}


void Foam::autoRefineDriver::baffleAndSplitMesh
(
    const refinementParameters& refineParams,
    const bool handleSnapProblems,
    const dictionary& motionDict
)
{
    Info<< nl
        << "Splitting mesh at surface intersections" << nl
        << "---------------------------------------" << nl
        << endl;

    const fvMesh& mesh = meshRefiner_.mesh();

    // Introduce baffles at surface intersections. Note:
    // meshRefiment::surfaceIndex() will
    // be like boundary face from now on so not coupled anymore.
    meshRefiner_.baffleAndSplitMesh
    (
        handleSnapProblems,             // detect&remove potential snap problem
        false,                          // perpendicular edge connected cells
        scalarField(0),                 // per region perpendicular angle
        !handleSnapProblems,            // merge free standing baffles?
        motionDict,
        const_cast<Time&>(mesh.time()),
        globalToMasterPatch_,
        globalToSlavePatch_,
        refineParams.keepPoints()[0]
    );
}


void Foam::autoRefineDriver::zonify
(
    const refinementParameters& refineParams
)
{
    // Mesh is at its finest. Do zoning
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // This puts all faces with intersection across a zoneable surface
    // into that surface's faceZone. All cells inside faceZone get given the
    // same cellZone.

    if (meshRefiner_.surfaces().getNamedSurfaces().size())
    {
        Info<< nl
            << "Introducing zones for interfaces" << nl
            << "--------------------------------" << nl
            << endl;

        const fvMesh& mesh = meshRefiner_.mesh();

        if (debug)
        {
            const_cast<Time&>(mesh.time())++;
        }

        meshRefiner_.zonify
        (
            refineParams.keepPoints()[0],
            refineParams.allowFreeStandingZoneFaces()
        );

        if (debug&meshRefinement::MESH)
        {
            Pout<< "Writing zoned mesh to time "
                << meshRefiner_.timeName() << '.' << endl;
            meshRefiner_.write
            (
                debug,
                mesh.time().path()/meshRefiner_.timeName()
            );
        }

        // Check that all faces are synced
        meshRefinement::checkCoupledFaceZones(mesh);
    }
}


void Foam::autoRefineDriver::splitAndMergeBaffles
(
    const refinementParameters& refineParams,
    const bool handleSnapProblems,
    const dictionary& motionDict
)
{
    Info<< nl
        << "Handling cells with snap problems" << nl
        << "---------------------------------" << nl
        << endl;

    const fvMesh& mesh = meshRefiner_.mesh();

    // Introduce baffles and split mesh
    if (debug)
    {
        const_cast<Time&>(mesh.time())++;
    }

    const scalarField& perpAngle = meshRefiner_.surfaces().perpendicularAngle();

    meshRefiner_.baffleAndSplitMesh
    (
        handleSnapProblems,
        handleSnapProblems,                 // remove perp edge connected cells
        perpAngle,                          // perp angle
        false,                              // merge free standing baffles?
        //true,                               // merge free standing baffles?
        motionDict,
        const_cast<Time&>(mesh.time()),
        globalToMasterPatch_,
        globalToSlavePatch_,
        refineParams.keepPoints()[0]
    );

    if (debug)
    {
        const_cast<Time&>(mesh.time())++;
    }

    // Duplicate points on baffles that are on more than one cell
    // region. This will help snapping pull them to separate surfaces.
    meshRefiner_.dupNonManifoldPoints();


    // Merge all baffles that are still remaining after duplicating points.
    List<labelPair> couples
    (
        meshRefiner_.getDuplicateFaces   // get all baffles
        (
            identity(mesh.nFaces()-mesh.nInternalFaces())
          + mesh.nInternalFaces()
        )
    );

    label nCouples = returnReduce(couples.size(), sumOp<label>());

    Info<< "Detected unsplittable baffles : "
        << nCouples << endl;

    if (nCouples > 0)
    {
        // Actually merge baffles. Note: not exactly parallellized. Should
        // convert baffle faces into processor faces if they resulted
        // from them.
        meshRefiner_.mergeBaffles(couples);

        if (debug)
        {
            // Debug:test all is still synced across proc patches
            meshRefiner_.checkData();
        }
        Info<< "Merged free-standing baffles in = "
            << mesh.time().cpuTimeIncrement() << " s." << endl;
    }

    if (debug&meshRefinement::MESH)
    {
        Pout<< "Writing handleProblemCells mesh to time "
            << meshRefiner_.timeName() << '.' << endl;
        meshRefiner_.write(debug, mesh.time().path()/meshRefiner_.timeName());
    }
}


void Foam::autoRefineDriver::mergePatchFaces
(
    const refinementParameters& refineParams,
    const dictionary& motionDict
)
{
    Info<< nl
        << "Merge refined boundary faces" << nl
        << "----------------------------" << nl
        << endl;

    const fvMesh& mesh = meshRefiner_.mesh();

    meshRefiner_.mergePatchFacesUndo
    (
        Foam::cos(degToRad(45.0)),
        Foam::cos(degToRad(45.0)),
        meshRefiner_.meshedPatches(),
        motionDict,
        labelList(mesh.nFaces(), -1)
    );

    if (debug)
    {
        meshRefiner_.checkData();
    }

    meshRefiner_.mergeEdgesUndo(Foam::cos(degToRad(45.0)), motionDict);

    if (debug)
    {
        meshRefiner_.checkData();
    }
}


void Foam::autoRefineDriver::doRefine
(
    const dictionary& refineDict,
    const refinementParameters& refineParams,
    const bool prepareForSnapping,
    const dictionary& motionDict
)
{
    Info<< nl
        << "Refinement phase" << nl
        << "----------------" << nl
        << endl;

    const fvMesh& mesh = meshRefiner_.mesh();

    // Check that all the keep points are inside the mesh.
    refineParams.findCells(mesh);

    // Refine around feature edges
    featureEdgeRefine
    (
        refineParams,
        100,    // maxIter
        0       // min cells to refine
    );

    // Refine based on surface
    surfaceOnlyRefine
    (
        refineParams,
        100     // maxIter
    );

    // Remove cells (a certain distance) beyond surface intersections
    removeInsideCells
    (
        refineParams,
        1       // nBufferLayers
    );

    // Internal mesh refinement
    shellRefine
    (
        refineParams,
        100    // maxIter
    );

    // Introduce baffles at surface intersections. Remove sections unreachable
    // from keepPoint.
    baffleAndSplitMesh(refineParams, prepareForSnapping, motionDict);

    // Mesh is at its finest. Do optional zoning.
    zonify(refineParams);

    // Pull baffles apart
    splitAndMergeBaffles(refineParams, prepareForSnapping, motionDict);

    // Do something about cells with refined faces on the boundary
    if (prepareForSnapping)
    {
        mergePatchFaces(refineParams, motionDict);
    }


    if (Pstream::parRun())
    {
        Info<< nl
            << "Doing final balancing" << nl
            << "---------------------" << nl
            << endl;

        //if (debug)
        //{
        //    const_cast<Time&>(mesh.time())++;
        //}

        // Do final balancing. Keep zoned faces on one processor since the
        // snap phase will convert them to baffles and this only works for
        // internal faces.
        meshRefiner_.balance
        (
            true,
            false,
            scalarField(mesh.nCells(), 1), // dummy weights
            decomposer_,
            distributor_
        );
    }
}


// ************************************************************************* //
