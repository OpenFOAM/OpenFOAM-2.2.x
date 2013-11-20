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

#include "meshToMeshNew.H"
#include "OFstream.H"
#include "Time.H"
#include "globalIndex.H"
#include "mergePoints.H"
#include "treeBoundBox.H"
#include "tetOverlapVolume.H"
#include "indexedOctree.H"
#include "treeDataCell.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(meshToMeshNew, 0);

    template<>
    const char* Foam::NamedEnum
    <
        Foam::meshToMeshNew::interpolationMethod,
        2
    >::names[] =
    {
        "map",
        "cellVolumeWeight"
    };

    const NamedEnum<meshToMeshNew::interpolationMethod, 2>
        meshToMeshNew::interpolationMethodNames_;
}

Foam::scalar Foam::meshToMeshNew::tolerance_ = 1e-6;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::meshToMeshNew::writeConnectivity
(
    const polyMesh& src,
    const polyMesh& tgt,
    const labelListList& srcToTargetAddr
) const
{
    Pout<< "Source size = " << src.nCells() << endl;
    Pout<< "Target size = " << tgt.nCells() << endl;

    word fName("addressing_" + src.name() + "_to_" + tgt.name());

    if (Pstream::parRun())
    {
        fName = fName +  "_proc" + Foam::name(Pstream::myProcNo());
    }

    OFstream os(src.time().path()/fName + ".obj");

    label vertI = 0;
    forAll(srcToTargetAddr, i)
    {
        const labelList& tgtAddress = srcToTargetAddr[i];
        forAll(tgtAddress, j)
        {
            label tgtI = tgtAddress[j];
            const vector& c0 = src.cellCentres()[i];

            const cell& c = tgt.cells()[tgtI];
            const pointField pts(c.points(tgt.faces(), tgt.points()));
            forAll(pts, j)
            {
                const point& p = pts[j];
                os  << "v " << p.x() << ' ' << p.y() << ' ' << p.z() << nl;
                vertI++;
                os  << "v " << c0.x() << ' ' << c0.y() << ' ' << c0.z()
                    << nl;
                vertI++;
                os  << "l " << vertI - 1 << ' ' << vertI << nl;
            }
        }
    }
}


Foam::labelList Foam::meshToMeshNew::maskCells
(
    const polyMesh& src,
    const polyMesh& tgt
) const
{
    boundBox intersectBb
    (
        max(src.bounds().min(), tgt.bounds().min()),
        min(src.bounds().max(), tgt.bounds().max())
    );

    intersectBb.inflate(0.01);

    const cellList& srcCells = src.cells();
    const faceList& srcFaces = src.faces();
    const pointField& srcPts = src.points();

    DynamicList<label> cells(src.size());
    forAll(srcCells, srcI)
    {
        boundBox cellBb(srcCells[srcI].points(srcFaces, srcPts), false);
        if (intersectBb.overlaps(cellBb))
        {
            cells.append(srcI);
        }
    }

    if (debug)
    {
        Pout<< "participating source mesh cells: " << cells.size() << endl;
    }

    return cells;
}


bool Foam::meshToMeshNew::findInitialSeeds
(
    const polyMesh& src,
    const polyMesh& tgt,
    const labelList& srcCellIDs,
    const boolList& mapFlag,
    const label startSeedI,
    label& srcSeedI,
    label& tgtSeedI
) const
{
    const cellList& srcCells = src.cells();
    const faceList& srcFaces = src.faces();
    const pointField& srcPts = src.points();

    for (label i = startSeedI; i < srcCellIDs.size(); i++)
    {
        label srcI = srcCellIDs[i];

        if (mapFlag[srcI])
        {
            const pointField
                pts(srcCells[srcI].points(srcFaces, srcPts).xfer());

            forAll(pts, ptI)
            {
                const point& pt = pts[ptI];
                label tgtI = tgt.cellTree().findInside(pt);

                if (tgtI != -1 && intersect(src, tgt, srcI, tgtI))
                {
                    srcSeedI = srcI;
                    tgtSeedI = tgtI;

                    return true;
                }
            }
        }
    }

    if (debug)
    {
        Pout<< "could not find starting seed" << endl;
    }

    return false;
}


void Foam::meshToMeshNew::appendToDirectSeeds
(
    const polyMesh& src,
    const polyMesh& tgt,
    boolList& mapFlag,
    labelList& srcTgtSeed,
    DynamicList<label>& srcSeeds,
    label& srcSeedI,
    label& tgtSeedI
) const
{
    const labelList& srcNbr = src.cellCells()[srcSeedI];
    const labelList& tgtNbr = tgt.cellCells()[tgtSeedI];

    const vectorField& srcCentre = src.cellCentres();

    forAll(srcNbr, i)
    {
        label srcI = srcNbr[i];

        if (mapFlag[srcI] && (srcTgtSeed[srcI] == -1))
        {
            // source cell srcI not yet mapped

            // identfy if target cell exists for source cell srcI
            bool found = false;
            forAll(tgtNbr, j)
            {
                label tgtI = tgtNbr[j];
                if
                (
                    tgt.pointInCell
                    (
                        srcCentre[srcI],
                        tgtI,
                        polyMesh::FACEPLANES
                    )
                )
                {
                    // new match - append to lists
                    found = true;

                    srcTgtSeed[srcI] = tgtI;
                    srcSeeds.append(srcI);

                    break;
                }
            }

            if (!found)
            {
                // no match available for source cell srcI
                mapFlag[srcI] = false;
            }
        }
    }

    if (srcSeeds.size())
    {
        srcSeedI = srcSeeds.remove();
        tgtSeedI = srcTgtSeed[srcSeedI];
    }
    else
    {
        srcSeedI = -1;
        tgtSeedI = -1;
    }
}


void Foam::meshToMeshNew::calcDirect
(
    const polyMesh& src,
    const polyMesh& tgt,
    const label srcSeedI,
    const label tgtSeedI
)
{
    // store a list of src cells already mapped
    boolList srcSeedFlag(src.nCells(), true);
    labelList srcTgtSeed(src.nCells(), -1);

    List<DynamicList<label> > srcToTgt(src.nCells());
    List<DynamicList<label> > tgtToSrc(tgt.nCells());

    DynamicList<label> srcSeeds;

    const scalarField& srcVc = src.cellVolumes();
    const scalarField& tgtVc = tgt.cellVolumes();

    label srcCellI = srcSeedI;
    label tgtCellI = tgtSeedI;

    do
    {
        // store src/tgt cell pair
        srcToTgt[srcCellI].append(tgtCellI);
        tgtToSrc[tgtCellI].append(srcCellI);

        // mark source cell srcSeedI as matched
        srcSeedFlag[srcCellI] = false;

        // accumulate intersection volume
        V_ += srcVc[srcCellI];

        // find new source seed cell
        appendToDirectSeeds
        (
            src,
            tgt,
            srcSeedFlag,
            srcTgtSeed,
            srcSeeds,
            srcCellI,
            tgtCellI
        );
    }
    while (srcCellI >= 0);

    // transfer addressing into persistent storage
    forAll(srcToTgtCellAddr_, i)
    {
        scalar v = srcVc[i];
        srcToTgtCellAddr_[i].transfer(srcToTgt[i]);
        srcToTgtCellWght_[i] = scalarList(srcToTgtCellAddr_[i].size(), v);
    }

    forAll(tgtToSrcCellAddr_, i)
    {
        scalar v = tgtVc[i];
        tgtToSrcCellAddr_[i].transfer(tgtToSrc[i]);
        tgtToSrcCellWght_[i] = scalarList(tgtToSrcCellAddr_[i].size(), v);
    }
}


void Foam::meshToMeshNew::normaliseWeights
(
    const word& descriptor,
    const scalarField& cellVolumes,
    const labelListList& addr,
    scalarListList& wght
) const
{
    const label nCell = returnReduce(wght.size(), sumOp<label>());

    if (nCell > 0)
    {
        scalar minW = GREAT;
        scalar maxW = -GREAT;

        forAll(wght, cellI)
        {
            scalarList& w = wght[cellI];
            scalar s = sum(w);
            scalar Vc = cellVolumes[cellI];

            forAll(w, i)
            {
                w[i] /= Vc;
            }

            minW = min(minW, s/Vc);
            maxW = max(maxW, s/Vc);
        }

        Info<< type() << ": " << descriptor << " weights min/max = "
            << returnReduce(minW, minOp<scalar>()) << ", "
            << returnReduce(maxW, maxOp<scalar>()) << endl;
    }
}


void Foam::meshToMeshNew::appendNbrTgtCells
(
    const label tgtCellI,
    const polyMesh& tgt,
    const DynamicList<label>& visitedTgtCells,
    DynamicList<label>& nbrTgtCellIDs
) const
{
    const labelList& nbrCells = tgt.cellCells()[tgtCellI];

    // filter out cells already visited from cell neighbours
    forAll(nbrCells, i)
    {
        label nbrCellI = nbrCells[i];

        if
        (
            (findIndex(visitedTgtCells, nbrCellI) == -1)
         && (findIndex(nbrTgtCellIDs, nbrCellI) == -1)
        )
        {
            nbrTgtCellIDs.append(nbrCellI);
        }
    }
}


void Foam::meshToMeshNew::setNextCells
(
    label& startSeedI,
    label& srcCellI,
    label& tgtCellI,
    const polyMesh& src,
    const polyMesh& tgt,
    const labelList& srcCellIDs,
    const boolList& mapFlag,
    const DynamicList<label>& visitedCells,
    labelList& seedCells
) const
{
    const labelList& srcNbrCells = src.cellCells()[srcCellI];

    // set possible seeds for later use by querying all src cell neighbours
    // with all visited target cells
    bool valuesSet = false;
    forAll(srcNbrCells, i)
    {
        label cellS = srcNbrCells[i];

        if (mapFlag[cellS] && seedCells[cellS] == -1)
        {
            forAll(visitedCells, j)
            {
                label cellT = visitedCells[j];

                if (intersect(src, tgt, cellS, cellT))
                {
                    seedCells[cellS] = cellT;

                    if (!valuesSet)
                    {
                        srcCellI = cellS;
                        tgtCellI = cellT;
                        valuesSet = true;
                    }
                }
            }
        }
    }

    // set next src and tgt cells if not set above
    if (valuesSet)
    {
        return;
    }
    else
    {
        // try to use existing seed
        bool foundNextSeed = false;
        for (label i = startSeedI; i < srcCellIDs.size(); i++)
        {
            label cellS = srcCellIDs[i];

            if (mapFlag[cellS])
            {
                if (!foundNextSeed)
                {
                    startSeedI = i;
                    foundNextSeed = true;
                }

                if (seedCells[cellS] != -1)
                {
                    srcCellI = cellS;
                    tgtCellI = seedCells[cellS];

                    return;
                }
            }
        }

        // perform new search to find match
        if (debug)
        {
            Pout<< "Advancing front stalled: searching for new "
                << "target cell" << endl;
        }

        bool restart =
            findInitialSeeds
            (
                src,
                tgt,
                srcCellIDs,
                mapFlag,
                startSeedI,
                srcCellI,
                tgtCellI
            );

        if (restart)
        {
            // successfully found new starting seed-pair
            return;
        }
    }

    // if we have got to here, there are no more src/tgt cell intersections
    srcCellI = -1;
    tgtCellI = -1;
}


bool Foam::meshToMeshNew::intersect
(
    const polyMesh& src,
    const polyMesh& tgt,
    const label srcCellI,
    const label tgtCellI
) const
{
    bool result = false;

    switch (method_)
    {
        case imMap:
        {
            result =
                tgt.pointInCell
                (
                    src.cellCentres()[srcCellI],
                    tgtCellI,
                    polyMesh::FACEPLANES
                );
            break;
        }
        case imCellVolumeWeight:
        {
            scalar threshold = tolerance_*src.cellVolumes()[srcCellI];

            tetOverlapVolume overlapEngine;

            treeBoundBox bbTgtCell
            (
                pointField
                (
                    tgt.points(),
                    tgt.cellPoints()[tgtCellI]
                )
            );

            result = overlapEngine.cellCellOverlapMinDecomp
            (
                src,
                srcCellI,
                tgt,
                tgtCellI,
                bbTgtCell,
                threshold
            );

            break;
        }
        default:
        {
            FatalErrorIn
            (
                "bool Foam::meshToMeshNew::intersect"
                "("
                    "const polyMesh&, "
                    "const polyMesh&, "
                    "const label, "
                    "const label"
                ") const"
            )
                << "Unknown interpolation method"
                << abort(FatalError);
        }
    }

    return result;
}


Foam::scalar Foam::meshToMeshNew::interVol
(
    const polyMesh& src,
    const polyMesh& tgt,
    const label srcCellI,
    const label tgtCellI
) const
{
    tetOverlapVolume overlapEngine;

    treeBoundBox bbTgtCell
    (
        pointField
        (
            tgt.points(),
            tgt.cellPoints()[tgtCellI]
        )
    );

    scalar vol = overlapEngine.cellCellOverlapVolumeMinDecomp
    (
        src,
        srcCellI,
        tgt,
        tgtCellI,
        bbTgtCell
    );

    return vol;
}


void Foam::meshToMeshNew::calcIndirect
(
    const polyMesh& src,
    const polyMesh& tgt,
    const label srcSeedI,
    const label tgtSeedI,
    const labelList& srcCellIDs,
    boolList& mapFlag,
    label& startSeedI
)
{
    label srcCellI = srcSeedI;
    label tgtCellI = tgtSeedI;

    List<DynamicList<label> > srcToTgtAddr(src.nCells());
    List<DynamicList<scalar> > srcToTgtWght(src.nCells());

    List<DynamicList<label> > tgtToSrcAddr(tgt.nCells());
    List<DynamicList<scalar> > tgtToSrcWght(tgt.nCells());

    // list of tgt cell neighbour cells
    DynamicList<label> nbrTgtCells(10);

    // list of tgt cells currently visited for srcCellI to avoid multiple hits
    DynamicList<label> visitedTgtCells(10);

    // list to keep track of tgt cells used to seed src cells
    labelList seedCells(src.nCells(), -1);
    seedCells[srcCellI] = tgtCellI;

    const scalarField& srcVol = src.cellVolumes();

    do
    {
        nbrTgtCells.clear();
        visitedTgtCells.clear();

        // append initial target cell and neighbours
        nbrTgtCells.append(tgtCellI);
        appendNbrTgtCells(tgtCellI, tgt, visitedTgtCells, nbrTgtCells);

        do
        {
            tgtCellI = nbrTgtCells.remove();
            visitedTgtCells.append(tgtCellI);

            scalar vol = interVol(src, tgt, srcCellI, tgtCellI);

            // accumulate addressing and weights for valid intersection
            if (vol/srcVol[srcCellI] > tolerance_)
            {
                // store src/tgt cell pair
                srcToTgtAddr[srcCellI].append(tgtCellI);
                srcToTgtWght[srcCellI].append(vol);

                tgtToSrcAddr[tgtCellI].append(srcCellI);
                tgtToSrcWght[tgtCellI].append(vol);

                appendNbrTgtCells(tgtCellI, tgt, visitedTgtCells, nbrTgtCells);

                // accumulate intersection volume
                V_ += vol;
            }
        }
        while (!nbrTgtCells.empty());

        mapFlag[srcCellI] = false;

        // find new source seed cell
        setNextCells
        (
            startSeedI,
            srcCellI,
            tgtCellI,
            src,
            tgt,
            srcCellIDs,
            mapFlag,
            visitedTgtCells,
            seedCells
        );
    }
    while (srcCellI != -1);

    // transfer addressing into persistent storage
    forAll(srcToTgtCellAddr_, i)
    {
        srcToTgtCellAddr_[i].transfer(srcToTgtAddr[i]);
        srcToTgtCellWght_[i].transfer(srcToTgtWght[i]);
    }

    forAll(tgtToSrcCellAddr_, i)
    {
        tgtToSrcCellAddr_[i].transfer(tgtToSrcAddr[i]);
        tgtToSrcCellWght_[i].transfer(tgtToSrcWght[i]);
    }
}


void Foam::meshToMeshNew::calcAddressing
(
    const polyMesh& src,
    const polyMesh& tgt
)
{
    srcToTgtCellAddr_.setSize(src.nCells());
    srcToTgtCellWght_.setSize(src.nCells());

    tgtToSrcCellAddr_.setSize(tgt.nCells());
    tgtToSrcCellWght_.setSize(tgt.nCells());

    if (!src.nCells() || !tgt.nCells())
    {
        if (debug)
        {
            Pout<< "mesh interpolation: cells not on processor: Source cells = "
                << src.nCells() << ", target cells = " << tgt.nCells()
                << endl;
        }
    }

    if (!src.nCells())
    {
        return;
    }
    else if (!tgt.nCells())
    {
        if (debug)
        {
            Pout<< "mesh interpolation: hhave " << src.nCells() << " source "
                << " cells but no target cells" << endl;
        }

        return;
    }

    // (potentially) participating source mesh cells
    const labelList srcCellIDs = maskCells(src, tgt);

    // list to keep track of whether src cell can be mapped
    boolList mapFlag(src.nCells(), false);
    UIndirectList<bool>(mapFlag, srcCellIDs) = true;

    // find initial point in tgt mesh
    label srcSeedI = -1;
    label tgtSeedI = -1;
    label startSeedI = 0;

    bool startWalk =
        findInitialSeeds
        (
            src,
            tgt,
            srcCellIDs,
            mapFlag,
            startSeedI,
            srcSeedI,
            tgtSeedI
        );

    if (!startWalk)
    {
        // if meshes are collocated, after inflating the source mesh bounding
        // box tgt mesh cells may be transferred, but may still not overlap
        // with the source mesh
        return;
    }

    switch (method_)
    {
        case imMap:
        {
            calcDirect(src, tgt, srcSeedI, tgtSeedI);
            break;
        }
        case imCellVolumeWeight:
        {
            calcIndirect
            (
                src,
                tgt,
                srcSeedI,
                tgtSeedI,
                srcCellIDs,
                mapFlag,
                startSeedI
            );
            break;
        }
        default:
        {
            FatalErrorIn
            (
                "void Foam::meshToMeshNew::calcAddressing"
                "("
                    "const polyMesh&, "
                    "const polyMesh&"
                ")"
            )
                << "Unknown interpolation method"
                << abort(FatalError);
        }
    }


    if (debug)
    {
        writeConnectivity(src, tgt, srcToTgtCellAddr_);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshToMeshNew::meshToMeshNew
(
    const polyMesh& src,
    const polyMesh& tgt,
    const interpolationMethod& method
)
:
    srcRegionName_(src.name()),
    tgtRegionName_(tgt.name()),
    srcToTgtCellAddr_(),
    tgtToSrcCellAddr_(),
    srcToTgtCellWght_(),
    tgtToSrcCellWght_(),
    method_(method),
    V_(0.0),
    singleMeshProc_(-1),
    srcMapPtr_(NULL),
    tgtMapPtr_(NULL)
{
    Info<< "Creating mesh-to-mesh addressing for " << src.name()
        << " and " << tgt.name() << " regions" << endl;

    singleMeshProc_ = calcDistribution(src, tgt);

    if (singleMeshProc_ == -1)
    {
        // create global indexing for src and tgt meshes
        globalIndex globalSrcCells(src.nCells());
        globalIndex globalTgtCells(tgt.nCells());

        // Create processor map of overlapping cells. This map gets
        // (possibly remote) cells from the tgt mesh such that they (together)
        // cover all of the src mesh
        autoPtr<mapDistribute> mapPtr = calcProcMap(src, tgt);
        const mapDistribute& map = mapPtr();

        pointField newTgtPoints;
        faceList newTgtFaces;
        labelList newTgtFaceOwners;
        labelList newTgtFaceNeighbours;
        labelList newTgtCellIDs;

        distributeAndMergeCells
        (
            map,
            tgt,
            globalTgtCells,
            newTgtPoints,
            newTgtFaces,
            newTgtFaceOwners,
            newTgtFaceNeighbours,
            newTgtCellIDs
        );


        // create a new target mesh
        polyMesh newTgt
        (
            IOobject
            (
                "newTgt::" + Foam::name(Pstream::myProcNo()),
                tgt.time().timeName(),
                tgt.time(),
                IOobject::NO_READ
            ),
            xferMove(newTgtPoints),
            xferMove(newTgtFaces),
            xferMove(newTgtFaceOwners),
            xferMove(newTgtFaceNeighbours),
            false                                   // no parallel comms
        );

        // create some dummy patch info
        List<polyPatch*> patches(1);
        patches[0] = new polyPatch
        (
            "defaultFaces",
            newTgt.nFaces() - newTgt.nInternalFaces(),
            newTgt.nInternalFaces(),
            0,
            newTgt.boundaryMesh(),
            word::null
        );

        newTgt.addPatches(patches);

        // force calculation of tet-base points used for point-in-cell
        (void)newTgt.tetBasePtIs();

        // force construction of cell tree
//        (void)newTgt.cellTree();

        if (debug)
        {
            Pout<< "Created newTgt mesh:" << nl
                << " old cells = " << tgt.nCells()
                << ", new cells = " << newTgt.nCells() << nl
                << " old faces = " << tgt.nFaces()
                << ", new faces = " << newTgt.nFaces() << endl;

            if (debug > 1)
            {
                Pout<< "Writing newTgt mesh: " << newTgt.name() << endl;
                newTgt.write();
            }
        }

        calcAddressing(src, newTgt);

        // per source cell the target cell address in newTgt mesh
        forAll(srcToTgtCellAddr_, i)
        {
            labelList& addressing = srcToTgtCellAddr_[i];
            forAll(addressing, addrI)
            {
                addressing[addrI] = newTgtCellIDs[addressing[addrI]];
            }
        }

        // convert target addresses in newTgtMesh into global cell numbering
        forAll(tgtToSrcCellAddr_, i)
        {
            labelList& addressing = tgtToSrcCellAddr_[i];
            forAll(addressing, addrI)
            {
                addressing[addrI] = globalSrcCells.toGlobal(addressing[addrI]);
            }
        }

        // set up as a reverse distribute
        mapDistribute::distribute
        (
            Pstream::nonBlocking,
            List<labelPair>(),
            tgt.nCells(),
            map.constructMap(),
            map.subMap(),
            tgtToSrcCellAddr_,
            ListPlusEqOp<label>(),
            labelList()
        );

        // set up as a reverse distribute
        mapDistribute::distribute
        (
            Pstream::nonBlocking,
            List<labelPair>(),
            tgt.nCells(),
            map.constructMap(),
            map.subMap(),
            tgtToSrcCellWght_,
            ListPlusEqOp<scalar>(),
            scalarList()
        );

        // weights normalisation
        normaliseWeights
        (
            "source",
            src.cellVolumes(),
            srcToTgtCellAddr_,
            srcToTgtCellWght_
        );

        normaliseWeights
        (
            "target",
            tgt.cellVolumes(),
            tgtToSrcCellAddr_,
            tgtToSrcCellWght_
        );

        // cache maps and reset addresses
        List<Map<label> > cMap;
        srcMapPtr_.reset
        (
            new mapDistribute(globalSrcCells, tgtToSrcCellAddr_, cMap)
        );
        tgtMapPtr_.reset
        (
            new mapDistribute(globalTgtCells, srcToTgtCellAddr_, cMap)
        );

        // collect volume intersection contributions
        reduce(V_, sumOp<scalar>());
    }
    else
    {
        calcAddressing(src, tgt);

        normaliseWeights
        (
            "source",
            src.cellVolumes(),
            srcToTgtCellAddr_,
            srcToTgtCellWght_
        );

        normaliseWeights
        (
            "target",
            tgt.cellVolumes(),
            tgtToSrcCellAddr_,
            tgtToSrcCellWght_
        );
    }

    Info<< "    Overlap volume: " << V_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::meshToMeshNew::~meshToMeshNew()
{}


// ************************************************************************* //
