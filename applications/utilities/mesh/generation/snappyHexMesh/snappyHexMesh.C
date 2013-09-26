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

Application
    snappyHexMesh

Description
    Automatic split hex mesher. Refines and snaps to surface.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "autoRefineDriver.H"
#include "autoSnapDriver.H"
#include "autoLayerDriver.H"
#include "searchableSurfaces.H"
#include "refinementSurfaces.H"
#include "refinementFeatures.H"
#include "shellSurfaces.H"
#include "decompositionMethod.H"
#include "noDecomp.H"
#include "fvMeshDistribute.H"
#include "wallPolyPatch.H"
#include "refinementParameters.H"
#include "snapParameters.H"
#include "layerParameters.H"
#include "vtkSetWriter.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Check writing tolerance before doing any serious work
scalar getMergeDistance(const polyMesh& mesh, const scalar mergeTol)
{
    const boundBox& meshBb = mesh.bounds();
    scalar mergeDist = mergeTol * meshBb.mag();

    Info<< nl
        << "Overall mesh bounding box  : " << meshBb << nl
        << "Relative tolerance         : " << mergeTol << nl
        << "Absolute matching distance : " << mergeDist << nl
        << endl;

    // check writing tolerance
    if (mesh.time().writeFormat() == IOstream::ASCII)
    {
        const scalar writeTol = std::pow
        (
            scalar(10.0),
            -scalar(IOstream::defaultPrecision())
        );

        if (mergeTol < writeTol)
        {
            FatalErrorIn("getMergeDistance(const polyMesh&, const dictionary&)")
                << "Your current settings specify ASCII writing with "
                << IOstream::defaultPrecision() << " digits precision." << nl
                << "Your merging tolerance (" << mergeTol
                << ") is finer than this." << nl
                << "Change to binary writeFormat, "
                << "or increase the writePrecision" << endl
                << "or adjust the merge tolerance (mergeTol)."
                << exit(FatalError);
        }
    }

    return mergeDist;
}


// Write mesh and additional information
void writeMesh
(
    const string& msg,
    const meshRefinement& meshRefiner,
    const bool writeLevel,
    const label debug
)
{
    const fvMesh& mesh = meshRefiner.mesh();

    meshRefiner.printMeshInfo(debug, msg);
    Info<< "Writing mesh to time " << meshRefiner.timeName() << endl;

    label flag = meshRefinement::MESH;
    if (writeLevel)
    {
        flag |= meshRefinement::SCALARLEVELS;
    }
    if (debug & meshRefinement::OBJINTERSECTIONS)
    {
        flag |= meshRefinement::OBJINTERSECTIONS;
    }
    meshRefiner.write(flag, mesh.time().path()/meshRefiner.timeName());
    Info<< "Wrote mesh in = "
        << mesh.time().cpuTimeIncrement() << " s." << endl;
}


int main(int argc, char *argv[])
{
#   include "addOverwriteOption.H"
    Foam::argList::addBoolOption
    (
        "checkGeometry",
        "check all surface geometry for quality"
    );
    Foam::argList::addBoolOption
    (
        "writeLevel",
        "write pointLevel and cellLevel postprocessing files"
    );

#   include "setRootCase.H"
#   include "createTime.H"
    runTime.functionObjects().off();
#   include "createMesh.H"

    Info<< "Read mesh in = "
        << runTime.cpuTimeIncrement() << " s" << endl;

    const bool overwrite = args.optionFound("overwrite");
    const bool checkGeometry = args.optionFound("checkGeometry");
    bool writeLevel = args.optionFound("writeLevel");

    // Check patches and faceZones are synchronised
    mesh.boundaryMesh().checkParallelSync(true);
    meshRefinement::checkCoupledFaceZones(mesh);


    // Read meshing dictionary
    IOdictionary meshDict
    (
       IOobject
       (
            "snappyHexMeshDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
       )
    );

    // all surface geometry
    const dictionary& geometryDict = meshDict.subDict("geometry");

    // refinement parameters
    const dictionary& refineDict = meshDict.subDict("castellatedMeshControls");

    // mesh motion and mesh quality parameters
    const dictionary& motionDict = meshDict.subDict("meshQualityControls");

    // snap-to-surface parameters
    const dictionary& snapDict = meshDict.subDict("snapControls");

    // layer addition parameters
    const dictionary& layerDict = meshDict.subDict("addLayersControls");

    // absolute merge distance
    const scalar mergeDist = getMergeDistance
    (
        mesh,
        readScalar(meshDict.lookup("mergeTolerance"))
    );



    // Read decomposePar dictionary
    dictionary decomposeDict;
    {
        if (Pstream::parRun())
        {
            decomposeDict = IOdictionary
            (
                IOobject
                (
                    "decomposeParDict",
                    runTime.system(),
                    mesh,
                    IOobject::MUST_READ_IF_MODIFIED,
                    IOobject::NO_WRITE
                )
            );
        }
        else
        {
            decomposeDict.add("method", "none");
            decomposeDict.add("numberOfSubdomains", 1);
        }
    }


    // Debug
    // ~~~~~

    const label debug = meshDict.lookupOrDefault<label>("debug", 0);
    if (debug > 0)
    {
        meshRefinement::debug   = debug;
        autoRefineDriver::debug = debug;
        autoSnapDriver::debug   = debug;
        autoLayerDriver::debug  = debug;
    }

    writeLevel = meshDict.lookupOrDefault<bool>("writeLevel", writeLevel);

    // Read geometry
    // ~~~~~~~~~~~~~

    searchableSurfaces allGeometry
    (
        IOobject
        (
            "abc",                      // dummy name
            mesh.time().constant(),     // instance
            //mesh.time().findInstance("triSurface", word::null),// instance
            "triSurface",               // local
            mesh.time(),                // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        geometryDict
    );


    // Read refinement surfaces
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    Info<< "Reading refinement surfaces." << endl;
    refinementSurfaces surfaces
    (
        allGeometry,
        refineDict.subDict("refinementSurfaces")
    );
    Info<< "Read refinement surfaces in = "
        << mesh.time().cpuTimeIncrement() << " s" << nl << endl;


    // Checking only?

    if (checkGeometry)
    {
        // Extract patchInfo
        List<wordList> patchTypes(allGeometry.size());

        const PtrList<dictionary>& patchInfo = surfaces.patchInfo();
        const labelList& surfaceGeometry = surfaces.surfaces();
        forAll(surfaceGeometry, surfI)
        {
            label geomI = surfaceGeometry[surfI];
            const wordList& regNames = allGeometry.regionNames()[geomI];

            patchTypes[geomI].setSize(regNames.size());
            forAll(regNames, regionI)
            {
                label globalRegionI = surfaces.globalRegion(surfI, regionI);

                if (patchInfo.set(globalRegionI))
                {
                    patchTypes[geomI][regionI] =
                        word(patchInfo[globalRegionI].lookup("type"));
                }
                else
                {
                    patchTypes[geomI][regionI] = wallPolyPatch::typeName;
                }
            }
        }

        // Write some stats
        allGeometry.writeStats(patchTypes, Info);
        // Check topology
        allGeometry.checkTopology(true);
        // Check geometry
        allGeometry.checkGeometry
        (
            100.0,      // max size ratio
            1e-9,       // intersection tolerance
            autoPtr<writer<scalar> >(new vtkSetWriter<scalar>()),
            0.01,       // min triangle quality
            true
        );

        return 0;
    }



    // Read refinement shells
    // ~~~~~~~~~~~~~~~~~~~~~~

    Info<< "Reading refinement shells." << endl;
    shellSurfaces shells
    (
        allGeometry,
        refineDict.subDict("refinementRegions")
    );
    Info<< "Read refinement shells in = "
        << mesh.time().cpuTimeIncrement() << " s" << nl << endl;


    Info<< "Setting refinement level of surface to be consistent"
        << " with shells." << endl;
    surfaces.setMinLevelFields(shells);
    Info<< "Checked shell refinement in = "
        << mesh.time().cpuTimeIncrement() << " s" << nl << endl;


    // Read feature meshes
    // ~~~~~~~~~~~~~~~~~~~

    Info<< "Reading features." << endl;
    refinementFeatures features
    (
        mesh,
        refineDict.lookup("features")
    );
    Info<< "Read features in = "
        << mesh.time().cpuTimeIncrement() << " s" << nl << endl;



    // Refinement engine
    // ~~~~~~~~~~~~~~~~~

    Info<< nl
        << "Determining initial surface intersections" << nl
        << "-----------------------------------------" << nl
        << endl;

    // Main refinement engine
    meshRefinement meshRefiner
    (
        mesh,
        mergeDist,          // tolerance used in sorting coordinates
        overwrite,          // overwrite mesh files?
        surfaces,           // for surface intersection refinement
        features,           // for feature edges/point based refinement
        shells              // for volume (inside/outside) refinement
    );
    Info<< "Calculated surface intersections in = "
        << mesh.time().cpuTimeIncrement() << " s" << nl << endl;

    // Some stats
    meshRefiner.printMeshInfo(debug, "Initial mesh");

    meshRefiner.write
    (
        debug & meshRefinement::OBJINTERSECTIONS,
        mesh.time().path()/meshRefiner.timeName()
    );


    // Add all the surface regions as patches
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //- Global surface region to patch (non faceZone surface) or patches
    //  (faceZone surfaces)
    labelList globalToMasterPatch;
    labelList globalToSlavePatch;
    {
        Info<< nl
            << "Adding patches for surface regions" << nl
            << "----------------------------------" << nl
            << endl;

        // From global region number to mesh patch.
        globalToMasterPatch.setSize(surfaces.nRegions(), -1);
        globalToSlavePatch.setSize(surfaces.nRegions(), -1);

        Info<< "Patch\tType\tRegion" << nl
            << "-----\t----\t------"
            << endl;

        const labelList& surfaceGeometry = surfaces.surfaces();
        const PtrList<dictionary>& surfacePatchInfo = surfaces.patchInfo();

        forAll(surfaceGeometry, surfI)
        {
            label geomI = surfaceGeometry[surfI];

            const wordList& regNames = allGeometry.regionNames()[geomI];

            Info<< surfaces.names()[surfI] << ':' << nl << nl;

            if (surfaces.faceZoneNames()[surfI].empty())
            {
                // 'Normal' surface
                forAll(regNames, i)
                {
                    label globalRegionI = surfaces.globalRegion(surfI, i);

                    label patchI;

                    if (surfacePatchInfo.set(globalRegionI))
                    {
                        patchI = meshRefiner.addMeshedPatch
                        (
                            regNames[i],
                            surfacePatchInfo[globalRegionI]
                        );
                    }
                    else
                    {
                        dictionary patchInfo;
                        patchInfo.set("type", wallPolyPatch::typeName);

                        patchI = meshRefiner.addMeshedPatch
                        (
                            regNames[i],
                            patchInfo
                        );
                    }

                    Info<< patchI << '\t' << mesh.boundaryMesh()[patchI].type()
                        << '\t' << regNames[i] << nl;

                    globalToMasterPatch[globalRegionI] = patchI;
                    globalToSlavePatch[globalRegionI] = patchI;
                }
            }
            else
            {
                // Zoned surface
                forAll(regNames, i)
                {
                    label globalRegionI = surfaces.globalRegion(surfI, i);

                    // Add master side patch
                    {
                        label patchI;

                        if (surfacePatchInfo.set(globalRegionI))
                        {
                            patchI = meshRefiner.addMeshedPatch
                            (
                                regNames[i],
                                surfacePatchInfo[globalRegionI]
                            );
                        }
                        else
                        {
                            dictionary patchInfo;
                            patchInfo.set("type", wallPolyPatch::typeName);

                            patchI = meshRefiner.addMeshedPatch
                            (
                                regNames[i],
                                patchInfo
                            );
                        }

                        Info<< patchI << '\t'
                            << mesh.boundaryMesh()[patchI].type()
                            << '\t' << regNames[i] << nl;

                        globalToMasterPatch[globalRegionI] = patchI;
                    }
                    // Add slave side patch
                    {
                        const word slaveName = regNames[i] + "_slave";
                        label patchI;

                        if (surfacePatchInfo.set(globalRegionI))
                        {
                            patchI = meshRefiner.addMeshedPatch
                            (
                                slaveName,
                                surfacePatchInfo[globalRegionI]
                            );
                        }
                        else
                        {
                            dictionary patchInfo;
                            patchInfo.set("type", wallPolyPatch::typeName);

                            patchI = meshRefiner.addMeshedPatch
                            (
                                slaveName,
                                patchInfo
                            );
                        }

                        Info<< patchI << '\t'
                            << mesh.boundaryMesh()[patchI].type()
                            << '\t' << slaveName << nl;

                        globalToSlavePatch[globalRegionI] = patchI;
                    }
                }
            }

            Info<< nl;
        }
        Info<< "Added patches in = "
            << mesh.time().cpuTimeIncrement() << " s" << nl << endl;
    }


    // Parallel
    // ~~~~~~~~

    // Decomposition
    autoPtr<decompositionMethod> decomposerPtr
    (
        decompositionMethod::New
        (
            decomposeDict
        )
    );
    decompositionMethod& decomposer = decomposerPtr();

    if (Pstream::parRun() && !decomposer.parallelAware())
    {
        FatalErrorIn(args.executable())
            << "You have selected decomposition method "
            << decomposer.typeName
            << " which is not parallel aware." << endl
            << "Please select one that is (hierarchical, ptscotch)"
            << exit(FatalError);
    }

    // Mesh distribution engine (uses tolerance to reconstruct meshes)
    fvMeshDistribute distributor(mesh, mergeDist);





    // Now do the real work -refinement -snapping -layers
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const Switch wantRefine(meshDict.lookup("castellatedMesh"));
    const Switch wantSnap(meshDict.lookup("snap"));
    const Switch wantLayers(meshDict.lookup("addLayers"));

    if (wantRefine)
    {
        cpuTime timer;

        autoRefineDriver refineDriver
        (
            meshRefiner,
            decomposer,
            distributor,
            globalToMasterPatch,
            globalToSlavePatch
        );

        // Refinement parameters
        refinementParameters refineParams(refineDict);

        if (!overwrite && !debug)
        {
            const_cast<Time&>(mesh.time())++;
        }

        refineDriver.doRefine(refineDict, refineParams, wantSnap, motionDict);

        writeMesh
        (
            "Refined mesh",
            meshRefiner,
            writeLevel,
            debug
        );

        Info<< "Mesh refined in = "
            << timer.cpuTimeIncrement() << " s." << endl;
    }

    if (wantSnap)
    {
        cpuTime timer;

        autoSnapDriver snapDriver
        (
            meshRefiner,
            globalToMasterPatch,
            globalToSlavePatch
        );

        // Snap parameters
        snapParameters snapParams(snapDict);
        // Temporary hack to get access to resolveFeatureAngle
        scalar curvature;
        {
            refinementParameters refineParams(refineDict);
            curvature = refineParams.curvature();
        }

        if (!overwrite && !debug)
        {
            const_cast<Time&>(mesh.time())++;
        }

        snapDriver.doSnap(snapDict, motionDict, curvature, snapParams);

        writeMesh
        (
            "Snapped mesh",
            meshRefiner,
            writeLevel,
            debug
        );

        Info<< "Mesh snapped in = "
            << timer.cpuTimeIncrement() << " s." << endl;
    }

    if (wantLayers)
    {
        cpuTime timer;

        autoLayerDriver layerDriver
        (
            meshRefiner,
            globalToMasterPatch,
            globalToSlavePatch
        );

        // Layer addition parameters
        layerParameters layerParams(layerDict, mesh.boundaryMesh());

        //!!! Temporary hack to get access to maxLocalCells
        bool preBalance;
        {
            refinementParameters refineParams(refineDict);

            preBalance = returnReduce
            (
                (mesh.nCells() >= refineParams.maxLocalCells()),
                orOp<bool>()
            );
        }


        if (!overwrite &&  !debug)
        {
            const_cast<Time&>(mesh.time())++;
        }

        layerDriver.doLayers
        (
            layerDict,
            motionDict,
            layerParams,
            preBalance,
            decomposer,
            distributor
        );

        writeMesh
        (
            "Layer mesh",
            meshRefiner,
            writeLevel,
            debug
        );

        Info<< "Layers added in = "
            << timer.cpuTimeIncrement() << " s." << endl;
    }


    Info<< "Finished meshing in = "
        << runTime.elapsedCpuTime() << " s." << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
