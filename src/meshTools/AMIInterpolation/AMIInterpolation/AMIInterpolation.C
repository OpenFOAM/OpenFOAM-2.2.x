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

#include "AMIInterpolation.H"
#include "meshTools.H"
#include "mapDistribute.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::writeIntersectionOBJ
(
    const scalar area,
    const face& f1,
    const face& f2,
    const pointField& f1Points,
    const pointField& f2Points
) const
{
    static label count = 1;

    const pointField f1pts = f1.points(f1Points);
    const pointField f2pts = f2.points(f2Points);

    Pout<< "Face intersection area (" << count <<  "):" << nl
        << "    f1 face = " << f1 << nl
        << "    f1 pts  = " << f1pts << nl
        << "    f2 face = " << f2 << nl
        << "    f2 pts  = " << f2pts << nl
        << "    area    = " << area
        << endl;

    OFstream os("areas" + name(count) + ".obj");

    forAll(f1pts, i)
    {
        meshTools::writeOBJ(os, f1pts[i]);
    }
    os<< "l";
    forAll(f1pts, i)
    {
        os<< " " << i + 1;
    }
    os<< " 1" << endl;


    forAll(f2pts, i)
    {
        meshTools::writeOBJ(os, f2pts[i]);
    }
    os<< "l";
    forAll(f2pts, i)
    {
        os<< " " << f1pts.size() + i + 1;
    }
    os<< " " << f1pts.size() + 1 << endl;

    count++;

}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::checkPatches
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch
) const
{
    const scalar maxBoundsError = 0.05;

    // check bounds of source and target
    boundBox bbSrc(srcPatch.points(), srcPatch.meshPoints(), true);
    boundBox bbTgt(tgtPatch.points(), tgtPatch.meshPoints(), true);

    boundBox bbTgtInf(bbTgt);
    bbTgtInf.inflate(maxBoundsError);

    if (!bbTgtInf.contains(bbSrc))
    {
        WarningIn
        (
            "AMIInterpolation<SourcePatch, TargetPatch>::checkPatches"
            "("
                "const SourcePatch&, "
                "const TargetPatch&"
            ")"
        )   << "Source and target patch bounding boxes are not similar" << nl
            << "    source box span     : " << bbSrc.span() << nl
            << "    target box span     : " << bbTgt.span() << nl
            << "    source box          : " << bbSrc << nl
            << "    target box          : " << bbTgt << nl
            << "    inflated target box : " << bbTgtInf << endl;
    }
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::resetTree
(
    const TargetPatch& tgtPatch
)
{
    // Clear the old octree
    treePtr_.clear();


    treeBoundBox bb(tgtPatch.points());
    bb.inflate(0.01);

    if (!treePtr_.valid())
    {
        treePtr_.reset
        (
            new indexedOctree<treeType>
            (
                treeType(false, tgtPatch),
                bb,                         // overall search domain
                8,                          // maxLevel
                10,                         // leaf size
                3.0                         // duplicity
            )
        );
    }
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::projectPointsToSurface
(
    const searchableSurface& surf,
    pointField& pts
) const
{
    if (debug)
    {
        Info<< "AMI: projecting points to surface" << endl;
    }

    List<pointIndexHit> nearInfo;

    surf.findNearest(pts, scalarField(pts.size(), GREAT), nearInfo);

    label nMiss = 0;
    forAll(nearInfo, i)
    {
        const pointIndexHit& pi = nearInfo[i];

        if (pi.hit())
        {
            pts[i] = pi.hitPoint();
        }
        else
        {
            pts[i] = pts[i];
            nMiss++;
        }
    }

    if (nMiss > 0)
    {
        FatalErrorIn
        (
            "void Foam::projectPointsToSurface"
            "("
                "const searchableSurface&, "
                "pointField&"
            ") const"
        )
            << "Error projecting points to surface: "
            << nMiss << " faces could not be determined"
            << abort(FatalError);
    }
}


template<class SourcePatch, class TargetPatch>
Foam::label Foam::AMIInterpolation<SourcePatch, TargetPatch>::findTargetFace
(
    const label srcFaceI,
    const SourcePatch& srcPatch
) const
{
    label targetFaceI = -1;

    const pointField& srcPts = srcPatch.points();
    const face& srcFace = srcPatch[srcFaceI];
    const point srcPt = srcFace.centre(srcPts);
    const scalar srcFaceArea = srcMagSf_[srcFaceI];

    pointIndexHit sample = treePtr_->findNearest(srcPt, 10.0*srcFaceArea);


    if (sample.hit())
    {
        targetFaceI = sample.index();

        if (debug)
        {
            Pout<< "Source point = " << srcPt << ", Sample point = "
                << sample.hitPoint() << ", Sample index = " << sample.index()
                << endl;
        }
    }

    return targetFaceI;
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::appendNbrFaces
(
    const label faceI,
    const TargetPatch& patch,
    const DynamicList<label>& visitedFaces,
    DynamicList<label>& faceIDs
) const
{
    const labelList& nbrFaces = patch.faceFaces()[faceI];

    // filter out faces already visited from src face neighbours
    forAll(nbrFaces, i)
    {
        label nbrFaceI = nbrFaces[i];
        bool valid = true;
        forAll(visitedFaces, j)
        {
            if (nbrFaceI == visitedFaces[j])
            {
                valid = false;
                break;
            }
        }

        if (valid)
        {
            forAll(faceIDs, j)
            {
                if (nbrFaceI == faceIDs[j])
                {
                    valid = false;
                    break;
                }
            }
        }

        // prevent addition of face if it is not on the same plane-ish
        if (valid)
        {
            const vector& n1 = patch.faceNormals()[faceI];
            const vector& n2 = patch.faceNormals()[nbrFaceI];

            scalar cosI = n1 & n2;

            if (cosI > Foam::cos(degToRad(89.0)))
            {
                faceIDs.append(nbrFaceI);
            }
        }
    }
}


template<class SourcePatch, class TargetPatch>
bool Foam::AMIInterpolation<SourcePatch, TargetPatch>::processSourceFace
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const label srcFaceI,
    const label tgtStartFaceI,

    // list of tgt face neighbour faces
    DynamicList<label>& nbrFaces,
    // list of faces currently visited for srcFaceI to avoid multiple hits
    DynamicList<label>& visitedFaces,

    // temporary storage for addressing and weights
    List<DynamicList<label> >& srcAddr,
    List<DynamicList<scalar> >& srcWght,
    List<DynamicList<label> >& tgtAddr,
    List<DynamicList<scalar> >& tgtWght
)
{
    nbrFaces.clear();
    visitedFaces.clear();

    // append initial target face and neighbours
    nbrFaces.append(tgtStartFaceI);

    appendNbrFaces(tgtStartFaceI, tgtPatch, visitedFaces, nbrFaces);

    bool faceProcessed = false;

    scalar sumArea = 0.0;

    do
    {
        // process new target face
        label tgtFaceI = nbrFaces.remove();

        visitedFaces.append(tgtFaceI);
        scalar area = interArea(srcFaceI, tgtFaceI, srcPatch, tgtPatch);

        // store when intersection area > 0
        if (area/srcMagSf_[srcFaceI] > faceAreaIntersect::tolerance())
        {
            srcAddr[srcFaceI].append(tgtFaceI);
            srcWght[srcFaceI].append(area);

            tgtAddr[tgtFaceI].append(srcFaceI);
            tgtWght[tgtFaceI].append(area);

            appendNbrFaces(tgtFaceI, tgtPatch, visitedFaces, nbrFaces);

            faceProcessed = true;

            sumArea += area;
        }


    } while (nbrFaces.size() > 0);

    return faceProcessed;
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::setNextFaces
(
    label& startSeedI,
    label& srcFaceI,
    label& tgtFaceI,
    const SourcePatch& srcPatch0,
    const TargetPatch& tgtPatch0,
    const boolList& mapFlag,
    labelList& seedFaces,
    const DynamicList<label>& visitedFaces
) const
{
    const labelList& srcNbrFaces = srcPatch0.faceFaces()[srcFaceI];

    // set possible seeds for later use
    bool valuesSet = false;
    forAll(srcNbrFaces, i)
    {
        label faceS = srcNbrFaces[i];

        if (mapFlag[faceS] && seedFaces[faceS] == -1)
        {
            forAll(visitedFaces, j)
            {
                label faceT = visitedFaces[j];
                scalar area = interArea(faceS, faceT, srcPatch0, tgtPatch0);

                // Check that faces have enough overlap for robust walking
                if (area/srcMagSf_[srcFaceI] > faceAreaIntersect::tolerance())
                {
                    // TODO - throwing area away - re-use in next iteration?

                    seedFaces[faceS] = faceT;

                    if (!valuesSet)
                    {
                        srcFaceI = faceS;
                        tgtFaceI = faceT;
                        valuesSet = true;
                    }
                }
            }
        }
    }

    // set next src and tgt faces if not set above
    if (valuesSet)
    {
        return;
    }
    else
    {
        // try to use existing seed
        bool foundNextSeed = false;
        for (label faceI = startSeedI; faceI < mapFlag.size(); faceI++)
        {
            if (mapFlag[faceI])
            {
                if (!foundNextSeed)
                {
                    startSeedI = faceI;
                    foundNextSeed = true;
                }

                if (seedFaces[faceI] != -1)
                {
                    srcFaceI = faceI;
                    tgtFaceI = seedFaces[faceI];
                    return;
                }
            }
        }

        // perform new search to find match
        if (debug)
        {
            Pout<< "Advancing front stalled: searching for new "
                << "target face" << endl;
        }

        foundNextSeed = false;
        for (label faceI = startSeedI; faceI < mapFlag.size(); faceI++)
        {
            if (mapFlag[faceI])
            {
                if (!foundNextSeed)
                {
                    startSeedI = faceI + 1;
                    foundNextSeed = true;
                }

                srcFaceI = faceI;
                tgtFaceI = findTargetFace(srcFaceI, srcPatch0);

                if (tgtFaceI >= 0)
                {
                    return;
                }
            }
        }

        FatalErrorIn
        (
            "void Foam::AMIInterpolation<SourcePatch, TargetPatch>::"
            "setNextFaces"
            "("
                "label&, "
                "label&, "
                "label&, "
                "const SourcePatch&, "
                "const TargetPatch&, "
                "const boolList&, "
                "labelList&, "
                "const DynamicList<label>&"
            ") const"
        )  << "Unable to set source and target faces" << abort(FatalError);
    }
}


template<class SourcePatch, class TargetPatch>
Foam::scalar Foam::AMIInterpolation<SourcePatch, TargetPatch>::interArea
(
    const label srcFaceI,
    const label tgtFaceI,
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch
) const
{
    const pointField& srcPoints = srcPatch.points();
    const pointField& tgtPoints = tgtPatch.points();

    // references to candidate faces
    const face& src = srcPatch[srcFaceI];
    const face& tgt = tgtPatch[tgtFaceI];

    // quick reject if either face has zero area
    // Note: do not used stored face areas for target patch
    const scalar tgtMag = tgt.mag(tgtPoints);
    if ((srcMagSf_[srcFaceI] < ROOTVSMALL) || ( tgtMag < ROOTVSMALL))
    {
        return 0.0;
    }

    // create intersection object
    faceAreaIntersect inter(srcPoints, tgtPoints, reverseTarget_);


    // crude resultant norm
    vector n(-srcPatch.faceNormals()[srcFaceI]);
    if (reverseTarget_)
    {
        n -= tgtPatch.faceNormals()[tgtFaceI];
    }
    else
    {
        n += tgtPatch.faceNormals()[tgtFaceI];
    }
    n *= 0.5;

    scalar area = 0;
    if (mag(n) > ROOTVSMALL)
    {
        area = inter.calc(src, tgt, n, triMode_);
    }
    else
    {
        WarningIn
        (
            "void Foam::AMIInterpolation<SourcePatch, TargetPatch>::"
            "interArea"
            "("
                "const label, "
                "const label, "
                "const SourcePatch&, "
                "const TargetPatch&"
            ") const"
        )   << "Invalid normal for source face " << srcFaceI
            << " points " << UIndirectList<point>(srcPoints, src)
            << " target face " << tgtFaceI
            << " points " << UIndirectList<point>(tgtPoints, tgt)
            << endl;
    }

    return area;
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::
restartUncoveredSourceFace
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    List<DynamicList<label> >& srcAddr,
    List<DynamicList<scalar> >& srcWght,
    List<DynamicList<label> >& tgtAddr,
    List<DynamicList<scalar> >& tgtWght
)
{
    // Collect all src faces with a low weight
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelHashSet lowWeightFaces(100);

    forAll(srcWght, srcFaceI)
    {
        scalar s = sum(srcWght[srcFaceI]);
        scalar t = s/srcMagSf_[srcFaceI];

        if (t < 0.5)
        {
            lowWeightFaces.insert(srcFaceI);
        }
    }

    if (debug)
    {
        Pout<< "AMIInterpolation: restarting search on "
            << lowWeightFaces.size() << " faces since sum of weights < 0.5"
            << endl;
    }

    if (lowWeightFaces.size() > 0)
    {
        // Erase all the lowWeight source faces from the target
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        DynamicList<label> okSrcFaces(10);
        DynamicList<scalar> okSrcWeights(10);
        forAll(tgtAddr, tgtFaceI)
        {
            okSrcFaces.clear();
            okSrcWeights.clear();
            DynamicList<label>& srcFaces = tgtAddr[tgtFaceI];
            DynamicList<scalar>& srcWeights = tgtWght[tgtFaceI];
            forAll(srcFaces, i)
            {
                if (!lowWeightFaces.found(srcFaces[i]))
                {
                    okSrcFaces.append(srcFaces[i]);
                    okSrcWeights.append(srcWeights[i]);
                }
            }
            if (okSrcFaces.size() < srcFaces.size())
            {
                srcFaces.transfer(okSrcFaces);
                srcWeights.transfer(okSrcWeights);
            }
        }


        // Restart search from best hit
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // list of tgt face neighbour faces
        DynamicList<label> nbrFaces(10);

        // list of faces currently visited for srcFaceI to avoid multiple hits
        DynamicList<label> visitedFaces(10);

        forAllConstIter(labelHashSet, lowWeightFaces, iter)
        {
            label srcFaceI = iter.key();
            label tgtFaceI = findTargetFace(srcFaceI, srcPatch);
            if (tgtFaceI != -1)
            {
                //bool faceProcessed =
                processSourceFace
                (
                    srcPatch,
                    tgtPatch,
                    srcFaceI,
                    tgtFaceI,

                    nbrFaces,
                    visitedFaces,

                    srcAddr,
                    srcWght,
                    tgtAddr,
                    tgtWght
                );
                // ? Check faceProcessed to see if restarting has worked.
            }
        }
    }
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::calcAddressing
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    label srcFaceI,
    label tgtFaceI
)
{
    // Pre-size to handle early exit
    srcAddress_.setSize(srcPatch.size());
    srcWeights_.setSize(srcPatch.size());
    tgtAddress_.setSize(tgtPatch.size());
    tgtWeights_.setSize(tgtPatch.size());


    if (debug && (!srcPatch.size() || !tgtPatch.size()))
    {
        Pout<< "AMI: Patches not on processor: Source faces = "
            << srcPatch.size() << ", target faces = " << tgtPatch.size()
            << endl;
    }

    if (!srcPatch.size())
    {
        return;
    }
    else if (!tgtPatch.size())
    {
        WarningIn
        (
            "void Foam::AMIInterpolation<SourcePatch, TargetPatch>::"
            "calcAddressing"
            "("
                "const SourcePatch&, "
                "const TargetPatch&, "
                "label, "
                "label"
            ")"
        ) << "Have " << srcPatch.size() << " source faces but no target faces."
          << endl;

        return;
    }

    resetTree(tgtPatch);

    // temporary storage for addressing and weights
    List<DynamicList<label> > srcAddr(srcPatch.size());
    List<DynamicList<scalar> > srcWght(srcPatch.size());
    List<DynamicList<label> > tgtAddr(tgtPatch.size());
    List<DynamicList<scalar> > tgtWght(tgtPatch.size());


    // find initial face match using brute force/octree search
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if ((srcFaceI == -1) || (tgtFaceI == -1))
    {
        srcFaceI = 0;
        tgtFaceI = 0;
        bool foundFace = false;
        forAll(srcPatch, faceI)
        {
            tgtFaceI = findTargetFace(faceI, srcPatch);
            if (tgtFaceI >= 0)
            {
                srcFaceI = faceI;
                foundFace = true;
                break;
            }
        }

        if (!foundFace)
        {
            FatalErrorIn
            (
                "void Foam::AMIInterpolation<SourcePatch, TargetPatch>::"
                "calcAddressing"
                "("
                    "const SourcePatch&, "
                    "const TargetPatch&, "
                    "label, "
                    "label"
                ")"
            )   << "Unable to find initial target face" << abort(FatalError);
        }
    }

    if (debug)
    {
        Pout<< "AMI: initial target face = " << tgtFaceI << endl;
    }


    // construct weights and addressing
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    label nFacesRemaining = srcPatch.size();

    // list of tgt face neighbour faces
    DynamicList<label> nbrFaces(10);

    // list of faces currently visited for srcFaceI to avoid multiple hits
    DynamicList<label> visitedFaces(10);

    // list to keep track of tgt faces used to seed src faces
    labelList seedFaces(nFacesRemaining, -1);
    seedFaces[srcFaceI] = tgtFaceI;

    // list to keep track of whether src face can be mapped
    boolList mapFlag(nFacesRemaining, true);

    // reset starting seed
    label startSeedI = 0;

    DynamicList<label> nonOverlapFaces;
    do
    {
        // Do advancing front starting from srcFaceI,tgtFaceI
        bool faceProcessed = processSourceFace
        (
            srcPatch,
            tgtPatch,
            srcFaceI,
            tgtFaceI,

            nbrFaces,
            visitedFaces,

            srcAddr,
            srcWght,
            tgtAddr,
            tgtWght
        );

        mapFlag[srcFaceI] = false;

        nFacesRemaining--;

        if (!faceProcessed)
        {
            nonOverlapFaces.append(srcFaceI);
        }

        // choose new src face from current src face neighbour
        if (nFacesRemaining > 0)
        {
            setNextFaces
            (
                startSeedI,
                srcFaceI,
                tgtFaceI,
                srcPatch,
                tgtPatch,
                mapFlag,
                seedFaces,
                visitedFaces
            );
        }
    } while (nFacesRemaining > 0);

    if (nonOverlapFaces.size() != 0)
    {
        Pout<< "AMI: " << nonOverlapFaces.size()
            << " non-overlap faces identified"
            << endl;

        srcNonOverlap_.transfer(nonOverlapFaces);
    }


    // Check for badly covered faces
    if (restartUncoveredSourceFace_)
    {
        restartUncoveredSourceFace
        (
            srcPatch,
            tgtPatch,
            srcAddr,
            srcWght,
            tgtAddr,
            tgtWght
        );
    }


    // transfer data to persistent storage
    forAll(srcAddr, i)
    {
        srcAddress_[i].transfer(srcAddr[i]);
        srcWeights_[i].transfer(srcWght[i]);
    }
    forAll(tgtAddr, i)
    {
        tgtAddress_[i].transfer(tgtAddr[i]);
        tgtWeights_[i].transfer(tgtWght[i]);
    }
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::normaliseWeights
(
    const scalarField& patchAreas,
    const word& patchName,
    const labelListList& addr,
    scalarListList& wght,
    scalarField& wghtSum,
    const bool output
)
{
    // Normalise the weights
    wghtSum.setSize(wght.size());
    forAll(wght, faceI)
    {
        scalar s = sum(wght[faceI]);
        scalar t = s/patchAreas[faceI];

        forAll(addr[faceI], i)
        {
            wght[faceI][i] /= s;
        }

        wghtSum[faceI] = t;
    }


    if (output)
    {
        const label nFace = returnReduce(wght.size(), sumOp<label>());

        if (nFace)
        {
            Info<< "AMI: Patch " << patchName << " weights min/max/average = "
                << gMin(wghtSum) << ", "
                << gMax(wghtSum) << ", "
                << gAverage(wghtSum) << endl;
        }
    }
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::agglomerate
(
    const autoPtr<mapDistribute>& targetMapPtr,
    const scalarField& fineSrcMagSf,
    const labelListList& fineSrcAddress,
    const scalarListList& fineSrcWeights,

    const labelList& sourceRestrictAddressing,
    const labelList& targetRestrictAddressing,

    scalarField& srcMagSf,
    labelListList& srcAddress,
    scalarListList& srcWeights,
    scalarField& srcWeightsSum,
    autoPtr<mapDistribute>& tgtMap
)
{
    label sourceCoarseSize =
    (
        sourceRestrictAddressing.size()
      ? max(sourceRestrictAddressing)+1
      : 0
    );

    label targetCoarseSize =
    (
        targetRestrictAddressing.size()
      ? max(targetRestrictAddressing)+1
      : 0
    );

    // Agglomerate face areas
    {
        srcMagSf.setSize(sourceRestrictAddressing.size(), 0.0);
        forAll(sourceRestrictAddressing, faceI)
        {
            label coarseFaceI = sourceRestrictAddressing[faceI];
            srcMagSf[coarseFaceI] += fineSrcMagSf[faceI];
        }
    }


    // Agglomerate weights and indices
    if (targetMapPtr.valid())
    {
        const mapDistribute& map = targetMapPtr();

        // Get all restriction addressing.
        labelList allRestrict(targetRestrictAddressing);
        map.distribute(allRestrict);

        // So now we have agglomeration of the target side in
        // allRestrict:
        //  0..size-1 : local agglomeration (= targetRestrictAddressing)
        //  size..    : agglomeration data from other processors

        labelListList tgtSubMap(Pstream::nProcs());

        // Local subMap is just identity
        {
            tgtSubMap[Pstream::myProcNo()] = identity(targetCoarseSize);
        }

        forAll(map.subMap(), procI)
        {
            if (procI != Pstream::myProcNo())
            {
                // Combine entries that point to the same coarse element. All
                // the elements refer to local data so index into
                // targetRestrictAddressing or allRestrict (since the same
                // for local data).
                const labelList& elems = map.subMap()[procI];
                labelList& newSubMap = tgtSubMap[procI];
                newSubMap.setSize(elems.size());

                labelList oldToNew(targetCoarseSize, -1);
                label newI = 0;

                forAll(elems, i)
                {
                    label fineElem = elems[i];
                    label coarseElem = allRestrict[fineElem];
                    if (oldToNew[coarseElem] == -1)
                    {
                        oldToNew[coarseElem] = newI;
                        newSubMap[newI] = coarseElem;
                        newI++;
                    }
                }
                newSubMap.setSize(newI);
            }
        }

        // Reconstruct constructMap by combining entries. Note that order
        // of handing out indices should be the same as loop above to compact
        // the sending map

        labelListList tgtConstructMap(Pstream::nProcs());
        labelList tgtCompactMap;

        // Local constructMap is just identity
        {
            tgtConstructMap[Pstream::myProcNo()] =
                identity(targetCoarseSize);

            tgtCompactMap = targetRestrictAddressing;
        }
        tgtCompactMap.setSize(map.constructSize());
        label compactI = targetCoarseSize;

        // Compact data from other processors
        forAll(map.constructMap(), procI)
        {
            if (procI != Pstream::myProcNo())
            {
                // Combine entries that point to the same coarse element. All
                // elements now are remote data so we cannot use any local
                // data here - use allRestrict instead.
                const labelList& elems = map.constructMap()[procI];

                labelList& newConstructMap = tgtConstructMap[procI];
                newConstructMap.setSize(elems.size());

                if (elems.size())
                {
                    // Get the maximum target coarse size for this set of
                    // received data.
                    label remoteTargetCoarseSize = labelMin;
                    forAll(elems, i)
                    {
                        remoteTargetCoarseSize = max
                        (
                            remoteTargetCoarseSize,
                            allRestrict[elems[i]]
                        );
                    }
                    remoteTargetCoarseSize += 1;

                    // Combine locally data coming from procI
                    labelList oldToNew(remoteTargetCoarseSize, -1);
                    label newI = 0;

                    forAll(elems, i)
                    {
                        label fineElem = elems[i];
                        // fineElem now points to section from procI
                        label coarseElem = allRestrict[fineElem];
                        if (oldToNew[coarseElem] == -1)
                        {
                            oldToNew[coarseElem] = newI;
                            tgtCompactMap[fineElem] = compactI;
                            newConstructMap[newI] = compactI++;
                            newI++;
                        }
                        else
                        {
                            // Get compact index
                            label compactI = oldToNew[coarseElem];
                            tgtCompactMap[fineElem] = newConstructMap[compactI];
                        }
                    }
                    newConstructMap.setSize(newI);
                }
            }
        }


        srcAddress.setSize(sourceCoarseSize);
        srcWeights.setSize(sourceCoarseSize);

        forAll(fineSrcAddress, faceI)
        {
            // All the elements contributing to faceI. Are slots in
            // mapDistribute'd data.
            const labelList& elems = fineSrcAddress[faceI];
            const scalarList& weights = fineSrcWeights[faceI];
            const scalar fineArea = fineSrcMagSf[faceI];

            label coarseFaceI = sourceRestrictAddressing[faceI];

            labelList& newElems = srcAddress[coarseFaceI];
            scalarList& newWeights = srcWeights[coarseFaceI];

            forAll(elems, i)
            {
                label elemI = elems[i];
                label coarseElemI = tgtCompactMap[elemI];

                label index = findIndex(newElems, coarseElemI);
                if (index == -1)
                {
                    newElems.append(coarseElemI);
                    newWeights.append(fineArea*weights[i]);
                }
                else
                {
                    newWeights[index] += fineArea*weights[i];
                }
            }
        }

        tgtMap.reset
        (
            new mapDistribute
            (
                compactI,
                tgtSubMap.xfer(),
                tgtConstructMap.xfer()
            )
        );
    }
    else
    {
        srcAddress.setSize(sourceCoarseSize);
        srcWeights.setSize(sourceCoarseSize);

        forAll(fineSrcAddress, faceI)
        {
            // All the elements contributing to faceI. Are slots in
            // mapDistribute'd data.
            const labelList& elems = fineSrcAddress[faceI];
            const scalarList& weights = fineSrcWeights[faceI];
            const scalar fineArea = fineSrcMagSf[faceI];

            label coarseFaceI = sourceRestrictAddressing[faceI];

            labelList& newElems = srcAddress[coarseFaceI];
            scalarList& newWeights = srcWeights[coarseFaceI];

            forAll(elems, i)
            {
                label elemI = elems[i];
                label coarseElemI = targetRestrictAddressing[elemI];

                label index = findIndex(newElems, coarseElemI);
                if (index == -1)
                {
                    newElems.append(coarseElemI);
                    newWeights.append(fineArea*weights[i]);
                }
                else
                {
                    newWeights[index] += fineArea*weights[i];
                }
            }
        }
    }

    // weights normalisation
    normaliseWeights
    (
        srcMagSf,
        "source",
        srcAddress,
        srcWeights,
        srcWeightsSum,
        false
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::AMIInterpolation
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const faceAreaIntersect::triangulationMode& triMode,
    const bool reverseTarget,
    const bool restartUncoveredSourceFace
)
:
    reverseTarget_(reverseTarget),
    restartUncoveredSourceFace_(restartUncoveredSourceFace),
    singlePatchProc_(-999),
    srcAddress_(),
    srcWeights_(),
    srcWeightsSum_(),
    srcNonOverlap_(),
    tgtAddress_(),
    tgtWeights_(),
    tgtWeightsSum_(),
    treePtr_(NULL),
    triMode_(triMode),
    srcMapPtr_(NULL),
    tgtMapPtr_(NULL)
{
    label srcSize = returnReduce(srcPatch.size(), sumOp<label>());
    label tgtSize = returnReduce(tgtPatch.size(), sumOp<label>());

    Info<< "AMI: Creating addressing and weights between "
        << srcSize << " source faces and " << tgtSize << " target faces"
        << endl;

    update(srcPatch, tgtPatch);
}


template<class SourcePatch, class TargetPatch>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::AMIInterpolation
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const autoPtr<searchableSurface>& surfPtr,
    const faceAreaIntersect::triangulationMode& triMode,
    const bool reverseTarget,
    const bool restartUncoveredSourceFace
)
:
    reverseTarget_(reverseTarget),
    restartUncoveredSourceFace_(restartUncoveredSourceFace),
    singlePatchProc_(-999),
    srcAddress_(),
    srcWeights_(),
    srcWeightsSum_(),
    srcNonOverlap_(),
    tgtAddress_(),
    tgtWeights_(),
    tgtWeightsSum_(),
    treePtr_(NULL),
    triMode_(triMode),
    srcMapPtr_(NULL),
    tgtMapPtr_(NULL)
{
    label srcSize = returnReduce(srcPatch.size(), sumOp<label>());
    label tgtSize = returnReduce(tgtPatch.size(), sumOp<label>());

    Info<< "AMI: Creating addressing and weights between "
        << srcSize << " source faces and " << tgtSize << " target faces"
        << endl;

    if (surfPtr.valid())
    {
        // create new patches for source and target
        pointField srcPoints = srcPatch.points();
        SourcePatch srcPatch0
        (
            SubList<face>
            (
                srcPatch,
                srcPatch.size(),
                0
            ),
            srcPoints
        );

        if (debug)
        {
            OFstream os("amiSrcPoints.obj");
            forAll(srcPoints, i)
            {
                meshTools::writeOBJ(os, srcPoints[i]);
            }
        }

        pointField tgtPoints = tgtPatch.points();
        TargetPatch tgtPatch0
        (
            SubList<face>
            (
                tgtPatch,
                tgtPatch.size(),
                0
            ),
            tgtPoints
        );

        if (debug)
        {
            OFstream os("amiTgtPoints.obj");
            forAll(tgtPoints, i)
            {
                meshTools::writeOBJ(os, tgtPoints[i]);
            }
        }


        // map source and target patches onto projection surface
        projectPointsToSurface(surfPtr(), srcPoints);
        projectPointsToSurface(surfPtr(), tgtPoints);


        // calculate AMI interpolation
        update(srcPatch0, tgtPatch0);
    }
    else
    {
        update(srcPatch, tgtPatch);
    }
}


template<class SourcePatch, class TargetPatch>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::AMIInterpolation
(
    const AMIInterpolation<SourcePatch, TargetPatch>& fineAMI,
    const labelList& sourceRestrictAddressing,
    const labelList& targetRestrictAddressing
)
:
    reverseTarget_(fineAMI.reverseTarget_),
    restartUncoveredSourceFace_(fineAMI.restartUncoveredSourceFace_),
    singlePatchProc_(fineAMI.singlePatchProc_),
    srcAddress_(),
    srcWeights_(),
    srcWeightsSum_(),
    srcNonOverlap_(),
    tgtAddress_(),
    tgtWeights_(),
    tgtWeightsSum_(),
    treePtr_(NULL),
    triMode_(fineAMI.triMode_),
    srcMapPtr_(NULL),
    tgtMapPtr_(NULL)
{
    label sourceCoarseSize =
    (
        sourceRestrictAddressing.size()
      ? max(sourceRestrictAddressing)+1
      : 0
    );

    label neighbourCoarseSize =
    (
        targetRestrictAddressing.size()
      ? max(targetRestrictAddressing)+1
      : 0
    );

    if (debug & 2)
    {
        Pout<< "AMI: Creating addressing and weights as agglomeration of AMI :"
            << " source:" << fineAMI.srcAddress().size()
            << " target:" << fineAMI.tgtAddress().size()
            << " coarse source size:" << sourceCoarseSize
            << " neighbour source size:" << neighbourCoarseSize
            << endl;
    }

    if
    (
        fineAMI.srcAddress().size() != sourceRestrictAddressing.size()
     || fineAMI.tgtAddress().size() != targetRestrictAddressing.size()
    )
    {
        FatalErrorIn
        (
            "AMIInterpolation<SourcePatch, TargetPatch>::AMIInterpolation"
            "("
            "    const AMIInterpolation<SourcePatch, TargetPatch>&, "
            "    const labelList&, "
            "    const labelList&"
            ")"
        )   << "Size mismatch." << nl
            << "Source patch size:" << fineAMI.srcAddress().size() << nl
            << "Source agglomeration size:"
            << sourceRestrictAddressing.size() << nl
            << "Target patch size:" << fineAMI.tgtAddress().size() << nl
            << "Target agglomeration size:"
            << targetRestrictAddressing.size()
            << exit(FatalError);
    }


    // Agglomerate addresses and weights

    agglomerate
    (
        fineAMI.tgtMapPtr_,
        fineAMI.srcMagSf(),
        fineAMI.srcAddress(),
        fineAMI.srcWeights(),

        sourceRestrictAddressing,
        targetRestrictAddressing,

        srcMagSf_,
        srcAddress_,
        srcWeights_,
        srcWeightsSum_,
        tgtMapPtr_
    );

    //if (tgtMapPtr_.valid())
    //{
    //    Pout<< "tgtMap:" << endl;
    //    string oldPrefix = Pout.prefix();
    //    Pout.prefix() = oldPrefix + "  ";
    //    tgtMapPtr_().printLayout(Pout);
    //    Pout.prefix() = oldPrefix;
    //}

    agglomerate
    (
        fineAMI.srcMapPtr_,
        fineAMI.tgtMagSf(),
        fineAMI.tgtAddress(),
        fineAMI.tgtWeights(),

        targetRestrictAddressing,
        sourceRestrictAddressing,

        tgtMagSf_,
        tgtAddress_,
        tgtWeights_,
        tgtWeightsSum_,
        srcMapPtr_
    );

    //if (srcMapPtr_.valid())
    //{
    //    Pout<< "srcMap:" << endl;
    //    string oldPrefix = Pout.prefix();
    //    Pout.prefix() = oldPrefix + "  ";
    //    srcMapPtr_().printLayout(Pout);
    //    Pout.prefix() = oldPrefix;
    //}
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::~AMIInterpolation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::update
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch
)
{
    // Calculate face areas
    srcMagSf_.setSize(srcPatch.size());
    forAll(srcMagSf_, faceI)
    {
        srcMagSf_[faceI] = srcPatch[faceI].mag(srcPatch.points());
    }
    tgtMagSf_.setSize(tgtPatch.size());
    forAll(tgtMagSf_, faceI)
    {
        tgtMagSf_[faceI] = tgtPatch[faceI].mag(tgtPatch.points());
    }


    // Calculate if patches present on multiple processors
    singlePatchProc_ = calcDistribution(srcPatch, tgtPatch);

    if (singlePatchProc_ == -1)
    {
        // convert local addressing to global addressing
        globalIndex globalSrcFaces(srcPatch.size());
        globalIndex globalTgtFaces(tgtPatch.size());

        // Create processor map of overlapping faces. This map gets
        // (possibly remote) faces from the tgtPatch such that they (together)
        // cover all of the srcPatch
        autoPtr<mapDistribute> mapPtr = calcProcMap(srcPatch, tgtPatch);
        const mapDistribute& map = mapPtr();

        // create new target patch that fully encompasses source patch

        // faces and points
        faceList newTgtFaces;
        pointField newTgtPoints;
        // original faces from tgtPatch (in globalIndexing since might be
        // remote)
        labelList tgtFaceIDs;
        distributeAndMergePatches
        (
            map,
            tgtPatch,
            globalTgtFaces,
            newTgtFaces,
            newTgtPoints,
            tgtFaceIDs
        );

        TargetPatch
            newTgtPatch
            (
                SubList<face>
                (
                    newTgtFaces,
                    newTgtFaces.size()
                ),
                newTgtPoints
            );

        checkPatches(srcPatch, newTgtPatch);


        // calculate AMI interpolation
        calcAddressing(srcPatch, newTgtPatch);

        // Now
        // ~~~
        //  srcAddress_ :   per srcPatch face a list of the newTgtPatch (not
        //                  tgtPatch) faces it overlaps
        //  tgtAddress_ :   per newTgtPatch (not tgtPatch) face a list of the
        //                  srcPatch faces it overlaps


        // Rework newTgtPatch indices into globalIndices of tgtPatch
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


        forAll(srcAddress_, i)
        {
            labelList& addressing = srcAddress_[i];
            forAll(addressing, addrI)
            {
                addressing[addrI] = tgtFaceIDs[addressing[addrI]];
            }
        }

        forAll(tgtAddress_, i)
        {
            labelList& addressing = tgtAddress_[i];
            forAll(addressing, addrI)
            {
                addressing[addrI] = globalSrcFaces.toGlobal(addressing[addrI]);
            }
        }

        // send data back to originating procs. Note that contributions
        // from different processors get added (ListAppendEqOp)

        mapDistribute::distribute
        (
            Pstream::nonBlocking,
            List<labelPair>(),
            tgtPatch.size(),
            map.constructMap(),
            map.subMap(),
            tgtAddress_,
            ListAppendEqOp<label>(),
            labelList()
        );

        mapDistribute::distribute
        (
            Pstream::nonBlocking,
            List<labelPair>(),
            tgtPatch.size(),
            map.constructMap(),
            map.subMap(),
            tgtWeights_,
            ListAppendEqOp<scalar>(),
            scalarList()
        );

        // weights normalisation
        normaliseWeights
        (
            srcMagSf_,
            "source",
            srcAddress_,
            srcWeights_,
            srcWeightsSum_,
            true
        );
        normaliseWeights
        (
            tgtMagSf_,
            "target",
            tgtAddress_,
            tgtWeights_,
            tgtWeightsSum_,
            true
        );

        // cache maps and reset addresses
        List<Map<label> > cMap;
        srcMapPtr_.reset(new mapDistribute(globalSrcFaces, tgtAddress_, cMap));
        tgtMapPtr_.reset(new mapDistribute(globalTgtFaces, srcAddress_, cMap));


        if (debug)
        {
            writeFaceConnectivity(srcPatch, newTgtPatch, srcAddress_);
        }
    }
    else
    {
        checkPatches(srcPatch, tgtPatch);

        calcAddressing(srcPatch, tgtPatch);

        normaliseWeights
        (
            srcMagSf_,
            "source",
            srcAddress_,
            srcWeights_,
            srcWeightsSum_,
            true
        );
        normaliseWeights
        (
            tgtMagSf_,
            "target",
            tgtAddress_,
            tgtWeights_,
            tgtWeightsSum_,
            true
        );
    }

    if (debug)
    {
        Info<< "AMIInterpolation : Constructed addressing and weights" << nl
            << "    triMode        :" << triMode_ << nl
            << "    singlePatchProc:" << singlePatchProc_ << nl
            << "    srcMagSf       :" << gSum(srcMagSf_) << nl
            << "    tgtMagSf       :" << gSum(tgtMagSf_) << nl
            << endl;
    }
}


template<class SourcePatch, class TargetPatch>
template<class Type, class CombineOp>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToTarget
(
    const UList<Type>& fld,
    const CombineOp& cop,
    List<Type>& result
) const
{
    if (fld.size() != srcAddress_.size())
    {
        FatalErrorIn
        (
            "AMIInterpolation::interpolateToTarget"
            "("
                "const UList<Type>&, "
                "const CombineOp&, "
                "List<Type>&"
            ") const"
        )   << "Supplied field size is not equal to source patch size" << nl
            << "    source patch   = " << srcAddress_.size() << nl
            << "    target patch   = " << tgtAddress_.size() << nl
            << "    supplied field = " << fld.size()
            << abort(FatalError);
    }

    result.setSize(tgtAddress_.size());

    if (singlePatchProc_ == -1)
    {
        const mapDistribute& map = srcMapPtr_();

        List<Type> work(fld);
        map.distribute(work);

        forAll(result, faceI)
        {
            const labelList& faces = tgtAddress_[faceI];
            const scalarList& weights = tgtWeights_[faceI];

            forAll(faces, i)
            {
                cop(result[faceI], faceI, work[faces[i]], weights[i]);
            }
        }
    }
    else
    {
        forAll(result, faceI)
        {
            const labelList& faces = tgtAddress_[faceI];
            const scalarList& weights = tgtWeights_[faceI];

            forAll(faces, i)
            {
                cop(result[faceI], faceI, fld[faces[i]], weights[i]);
            }
        }
    }
}


template<class SourcePatch, class TargetPatch>
template<class Type, class CombineOp>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToSource
(
    const UList<Type>& fld,
    const CombineOp& cop,
    List<Type>& result
) const
{
    if (fld.size() != tgtAddress_.size())
    {
        FatalErrorIn
        (
            "AMIInterpolation::interpolateToSource"
            "("
                "const UList<Type>&, "
                "const CombineOp&, "
                "List<Type>&"
            ") const"
        )   << "Supplied field size is not equal to target patch size" << nl
            << "    source patch   = " << srcAddress_.size() << nl
            << "    target patch   = " << tgtAddress_.size() << nl
            << "    supplied field = " << fld.size()
            << abort(FatalError);
    }

    result.setSize(srcAddress_.size());

    if (singlePatchProc_ == -1)
    {
        const mapDistribute& map = tgtMapPtr_();

        List<Type> work(fld);
        map.distribute(work);

        forAll(result, faceI)
        {
            const labelList& faces = srcAddress_[faceI];
            const scalarList& weights = srcWeights_[faceI];

            forAll(faces, i)
            {
                cop(result[faceI], faceI, work[faces[i]], weights[i]);
            }
        }
    }
    else
    {
        forAll(result, faceI)
        {
            const labelList& faces = srcAddress_[faceI];
            const scalarList& weights = srcWeights_[faceI];

            forAll(faces, i)
            {
                cop(result[faceI], faceI, fld[faces[i]], weights[i]);
            }
        }
    }
}


template<class SourcePatch, class TargetPatch>
template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type> >
Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToSource
(
    const Field<Type>& fld,
    const CombineOp& cop
) const
{
    tmp<Field<Type> > tresult
    (
        new Field<Type>
        (
            srcAddress_.size(),
            pTraits<Type>::zero
        )
    );

    interpolateToSource
    (
        fld,
        multiplyWeightedOp<Type, CombineOp>(cop),
        tresult()
    );

    return tresult;
}


template<class SourcePatch, class TargetPatch>
template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type> >
Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToSource
(
    const tmp<Field<Type> >& tFld,
    const CombineOp& cop
) const
{
    return interpolateToSource(tFld(), cop);
}


template<class SourcePatch, class TargetPatch>
template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type> >
Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToTarget
(
    const Field<Type>& fld,
    const CombineOp& cop
) const
{
    tmp<Field<Type> > tresult
    (
        new Field<Type>
        (
            tgtAddress_.size(),
            pTraits<Type>::zero
        )
    );

    interpolateToTarget
    (
        fld,
        multiplyWeightedOp<Type, CombineOp>(cop),
        tresult()
    );

    return tresult;
}


template<class SourcePatch, class TargetPatch>
template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type> >
Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToTarget
(
    const tmp<Field<Type> >& tFld,
    const CombineOp& cop
) const
{
    return interpolateToTarget(tFld(), cop);
}


template<class SourcePatch, class TargetPatch>
template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToSource
(
    const Field<Type>& fld
) const
{
    return interpolateToSource(fld, plusEqOp<Type>());
}


template<class SourcePatch, class TargetPatch>
template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToSource
(
    const tmp<Field<Type> >& tFld
) const
{
    return interpolateToSource(tFld(), plusEqOp<Type>());
}


template<class SourcePatch, class TargetPatch>
template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToTarget
(
    const Field<Type>& fld
) const
{
    return interpolateToTarget(fld, plusEqOp<Type>());
}


template<class SourcePatch, class TargetPatch>
template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToTarget
(
    const tmp<Field<Type> >& tFld
) const
{
    return interpolateToTarget(tFld(), plusEqOp<Type>());
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::writeFaceConnectivity
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const labelListList& srcAddress
)
const
{
    OFstream os("faceConnectivity" + Foam::name(Pstream::myProcNo()) + ".obj");

    label ptI = 1;

    forAll(srcAddress, i)
    {
        const labelList& addr = srcAddress[i];
        const point& srcPt = srcPatch.faceCentres()[i];
        forAll(addr, j)
        {
            label tgtPtI = addr[j];
            const point& tgtPt = tgtPatch.faceCentres()[tgtPtI];

            meshTools::writeOBJ(os, srcPt);
            meshTools::writeOBJ(os, tgtPt);

            os  << "l " << ptI << " " << ptI + 1 << endl;

            ptI += 2;
        }
    }
}


template<class SourcePatch, class TargetPatch>
Foam::label Foam::AMIInterpolation<SourcePatch, TargetPatch>::srcPointFace
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const vector& n,
    const label tgtFaceI,
    point& tgtPoint
)
const
{
    const pointField& srcPoints = srcPatch.points();

    // source face addresses that intersect target face tgtFaceI
    const labelList& addr = tgtAddress_[tgtFaceI];

    forAll(addr, i)
    {
        label srcFaceI = addr[i];
        const face& f = srcPatch[tgtFaceI];

        pointHit ray = f.ray(tgtPoint, n, srcPoints);

        if (ray.hit())
        {
            tgtPoint = ray.rawPoint();

            return srcFaceI;
        }
    }

    // no hit registered - try with face normal instead of input normal
    forAll(addr, i)
    {
        label srcFaceI = addr[i];
        const face& f = srcPatch[tgtFaceI];

        vector nFace(-srcPatch.faceNormals()[srcFaceI]);
        nFace += tgtPatch.faceNormals()[tgtFaceI];

        pointHit ray = f.ray(tgtPoint, nFace, srcPoints);

        if (ray.hit())
        {
            tgtPoint = ray.rawPoint();

            return srcFaceI;
        }
    }

    return -1;
}


template<class SourcePatch, class TargetPatch>
Foam::label Foam::AMIInterpolation<SourcePatch, TargetPatch>::tgtPointFace
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const vector& n,
    const label srcFaceI,
    point& srcPoint
)
const
{
    const pointField& tgtPoints = tgtPatch.points();

    // target face addresses that intersect source face srcFaceI
    const labelList& addr = srcAddress_[srcFaceI];

    forAll(addr, i)
    {
        label tgtFaceI = addr[i];
        const face& f = tgtPatch[tgtFaceI];

        pointHit ray = f.ray(srcPoint, n, tgtPoints);

        if (ray.hit())
        {
            srcPoint = ray.rawPoint();

            return tgtFaceI;
        }
    }

    // no hit registered - try with face normal instead of input normal
    forAll(addr, i)
    {
        label tgtFaceI = addr[i];
        const face& f = tgtPatch[tgtFaceI];

        vector nFace(-srcPatch.faceNormals()[srcFaceI]);
        nFace += tgtPatch.faceNormals()[tgtFaceI];

        pointHit ray = f.ray(srcPoint, n, tgtPoints);

        if (ray.hit())
        {
            srcPoint = ray.rawPoint();

            return tgtFaceI;
        }
    }

    return -1;
}


// ************************************************************************* //
