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

#include "triSurfaceMesh.H"
#include "Random.H"
#include "addToRunTimeSelectionTable.H"
#include "EdgeMap.H"
#include "triSurfaceFields.H"
#include "Time.H"
#include "PackedBoolList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(triSurfaceMesh, 0);
addToRunTimeSelectionTable(searchableSurface, triSurfaceMesh, dict);

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//// Special version of Time::findInstance that does not check headerOk
//// to search for instances of raw files
//Foam::word Foam::triSurfaceMesh::findRawInstance
//(
//    const Time& runTime,
//    const fileName& dir,
//    const word& name
//)
//{
//    // Check current time first
//    if (isFile(runTime.path()/runTime.timeName()/dir/name))
//    {
//        return runTime.timeName();
//    }
//    instantList ts = runTime.times();
//    label instanceI;
//
//    for (instanceI = ts.size()-1; instanceI >= 0; --instanceI)
//    {
//        if (ts[instanceI].value() <= runTime.timeOutputValue())
//        {
//            break;
//        }
//    }
//
//    // continue searching from here
//    for (; instanceI >= 0; --instanceI)
//    {
//        if (isFile(runTime.path()/ts[instanceI].name()/dir/name))
//        {
//            return ts[instanceI].name();
//        }
//    }
//
//
//    // not in any of the time directories, try constant
//
//    // Note. This needs to be a hard-coded constant, rather than the
//    // constant function of the time, because the latter points to
//    // the case constant directory in parallel cases
//
//    if (isFile(runTime.path()/runTime.constant()/dir/name))
//    {
//        return runTime.constant();
//    }
//
//    FatalErrorIn
//    (
//        "searchableSurfaces::findRawInstance"
//        "(const Time&, const fileName&, const word&)"
//    )   << "Cannot find file \"" << name << "\" in directory "
//        << runTime.constant()/dir
//        << exit(FatalError);
//
//    return runTime.constant();
//}


//- Check file existence
const Foam::fileName& Foam::triSurfaceMesh::checkFile
(
    const fileName& fName,
    const fileName& objectName
)
{
    if (fName.empty())
    {
        FatalErrorIn
        (
            "triSurfaceMesh::checkFile(const fileName&, const fileName&)"
        )   << "Cannot find triSurfaceMesh starting from "
            << objectName << exit(FatalError);
    }
    return fName;
}


bool Foam::triSurfaceMesh::addFaceToEdge
(
    const edge& e,
    EdgeMap<label>& facesPerEdge
)
{
    EdgeMap<label>::iterator eFnd = facesPerEdge.find(e);
    if (eFnd != facesPerEdge.end())
    {
        if (eFnd() == 2)
        {
            return false;
        }
        eFnd()++;
    }
    else
    {
        facesPerEdge.insert(e, 1);
    }
    return true;
}


bool Foam::triSurfaceMesh::isSurfaceClosed() const
{
    // Construct pointFaces. Let's hope surface has compact point
    // numbering ...
    labelListList pointFaces;
    invertManyToMany(points()().size(), *this, pointFaces);

    // Loop over all faces surrounding point. Count edges emanating from point.
    // Every edge should be used by two faces exactly.
    // To prevent doing work twice per edge only look at edges to higher
    // point
    EdgeMap<label> facesPerEdge(100);
    forAll(pointFaces, pointI)
    {
        const labelList& pFaces = pointFaces[pointI];

        facesPerEdge.clear();
        forAll(pFaces, i)
        {
            const triSurface::FaceType& f = triSurface::operator[](pFaces[i]);
            label fp = findIndex(f, pointI);

            // Something weird: if I expand the code of addFaceToEdge in both
            // below instances it gives a segmentation violation on some
            // surfaces. Compiler (4.3.2) problem?


            // Forward edge
            label nextPointI = f[f.fcIndex(fp)];

            if (nextPointI > pointI)
            {
                bool okFace = addFaceToEdge
                (
                    edge(pointI, nextPointI),
                    facesPerEdge
                );

                if (!okFace)
                {
                    return false;
                }
            }
            // Reverse edge
            label prevPointI = f[f.rcIndex(fp)];

            if (prevPointI > pointI)
            {
                bool okFace = addFaceToEdge
                (
                    edge(pointI, prevPointI),
                    facesPerEdge
                );

                if (!okFace)
                {
                    return false;
                }
            }
        }

        // Check for any edges used only once.
        forAllConstIter(EdgeMap<label>, facesPerEdge, iter)
        {
            if (iter() != 2)
            {
                return false;
            }
        }
    }

    return true;
}


// Gets all intersections after initial one. Adds smallVec and starts tracking
// from there.
void Foam::triSurfaceMesh::getNextIntersections
(
    const indexedOctree<treeDataTriSurface>& octree,
    const point& start,
    const point& end,
    const vector& smallVec,
    DynamicList<pointIndexHit, 1, 1>& hits
)
{
    const vector dirVec(end-start);
    const scalar magSqrDirVec(magSqr(dirVec));

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
        pointIndexHit inter = octree.findLine(pt, end);

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
            hits.append(inter);
            // Restore perturbVec
            perturbVec = smallVec;
        }
    }
}


void Foam::triSurfaceMesh::calcBounds(boundBox& bb, label& nPoints) const
{
    // Unfortunately nPoints constructs meshPoints() so do compact version
    // ourselves

    const triSurface& s = static_cast<const triSurface&>(*this);

    PackedBoolList pointIsUsed(points()().size());

    nPoints = 0;
    bb = boundBox::invertedBox;

    forAll(s, faceI)
    {
        const triSurface::FaceType& f = s[faceI];

        forAll(f, fp)
        {
            label pointI = f[fp];
            if (pointIsUsed.set(pointI, 1u))
            {
                bb.min() = ::Foam::min(bb.min(), points()()[pointI]);
                bb.max() = ::Foam::max(bb.max(), points()()[pointI]);
                nPoints++;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::triSurfaceMesh::triSurfaceMesh(const IOobject& io, const triSurface& s)
:
    searchableSurface(io),
    objectRegistry
    (
        IOobject
        (
            io.name(),
            io.instance(),
            io.local(),
            io.db(),
            io.readOpt(),
            io.writeOpt(),
            false       // searchableSurface already registered under name
        )
    ),
    triSurface(s),
    tolerance_(indexedOctree<treeDataTriSurface>::perturbTol()),
    minQuality_(-1),
    maxTreeDepth_(10),
    surfaceClosed_(-1)
{
    bounds() = boundBox(points());
}


Foam::triSurfaceMesh::triSurfaceMesh(const IOobject& io)
:
    // Find instance for triSurfaceMesh
    searchableSurface(io),
    //(
    //    IOobject
    //    (
    //        io.name(),
    //        io.time().findInstance(io.local(), word::null),
    //        io.local(),
    //        io.db(),
    //        io.readOpt(),
    //        io.writeOpt(),
    //        io.registerObject()
    //    )
    //),
    // Reused found instance in objectRegistry
    objectRegistry
    (
        IOobject
        (
            io.name(),
            static_cast<const searchableSurface&>(*this).instance(),
            io.local(),
            io.db(),
            io.readOpt(),
            io.writeOpt(),
            false       // searchableSurface already registered under name
        )
    ),
    triSurface
    (
        checkFile
        (
            searchableSurface::filePath(),
            searchableSurface::objectPath()
        )
    ),
    tolerance_(indexedOctree<treeDataTriSurface>::perturbTol()),
    minQuality_(-1),
    maxTreeDepth_(10),
    surfaceClosed_(-1)
{
    bounds() = boundBox(points());
}


Foam::triSurfaceMesh::triSurfaceMesh
(
    const IOobject& io,
    const dictionary& dict
)
:
    searchableSurface(io),
    //(
    //    IOobject
    //    (
    //        io.name(),
    //        io.time().findInstance(io.local(), word::null),
    //        io.local(),
    //        io.db(),
    //        io.readOpt(),
    //        io.writeOpt(),
    //        io.registerObject()
    //    )
    //),
    // Reused found instance in objectRegistry
    objectRegistry
    (
        IOobject
        (
            io.name(),
            static_cast<const searchableSurface&>(*this).instance(),
            io.local(),
            io.db(),
            io.readOpt(),
            io.writeOpt(),
            false       // searchableSurface already registered under name
        )
    ),
    triSurface
    (
        checkFile
        (
            searchableSurface::filePath(),
            searchableSurface::objectPath()
        )
    ),
    tolerance_(indexedOctree<treeDataTriSurface>::perturbTol()),
    minQuality_(-1),
    maxTreeDepth_(10),
    surfaceClosed_(-1)
{
    scalar scaleFactor = 0;

    // allow rescaling of the surface points
    // eg, CAD geometries are often done in millimeters
    if (dict.readIfPresent("scale", scaleFactor) && scaleFactor > 0)
    {
        Info<< searchableSurface::name() << " : using scale " << scaleFactor
            << endl;
        triSurface::scalePoints(scaleFactor);
    }

    bounds() = boundBox(points());

    // Have optional non-standard search tolerance for gappy surfaces.
    if (dict.readIfPresent("tolerance", tolerance_) && tolerance_ > 0)
    {
        Info<< searchableSurface::name() << " : using intersection tolerance "
            << tolerance_ << endl;
    }

    // Have optional minimum quality for normal calculation
    if (dict.readIfPresent("minQuality", minQuality_) && minQuality_ > 0)
    {
        Info<< searchableSurface::name()
            << " : ignoring triangles with quality < "
            << minQuality_ << " for normals calculation." << endl;
    }

    // Have optional non-standard tree-depth to limit storage.
    if (dict.readIfPresent("maxTreeDepth", maxTreeDepth_) && maxTreeDepth_ > 0)
    {
        Info<< searchableSurface::name() << " : using maximum tree depth "
            << maxTreeDepth_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::triSurfaceMesh::~triSurfaceMesh()
{
    clearOut();
}


void Foam::triSurfaceMesh::clearOut()
{
    tree_.clear();
    edgeTree_.clear();
    triSurface::clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::triSurfaceMesh::coordinates() const
{
    tmp<pointField> tPts(new pointField(8));
    pointField& pt = tPts();

    // Use copy to calculate face centres so they don't get stored
    pt = PrimitivePatch<triSurface::FaceType, SubList, const pointField&>
    (
        SubList<triSurface::FaceType>(*this, triSurface::size()),
        triSurface::points()
    ).faceCentres();

    return tPts;
}


Foam::tmp<Foam::pointField> Foam::triSurfaceMesh::points() const
{
    return triSurface::points();
}


bool Foam::triSurfaceMesh::overlaps(const boundBox& bb) const
{
     const indexedOctree<treeDataTriSurface>& octree = tree();

     labelList indices = octree.findBox(treeBoundBox(bb));

     return !indices.empty();
}


void Foam::triSurfaceMesh::movePoints(const pointField& newPoints)
{
    tree_.clear();
    edgeTree_.clear();
    triSurface::movePoints(newPoints);
}


const Foam::indexedOctree<Foam::treeDataTriSurface>&
Foam::triSurfaceMesh::tree() const
{
    if (tree_.empty())
    {
        // Calculate bb without constructing local point numbering.
        treeBoundBox bb(vector::zero, vector::zero);

        if (size())
        {
            label nPoints;
            calcBounds(bb, nPoints);

            if (nPoints != points()().size())
            {
                WarningIn("triSurfaceMesh::tree() const")
                    << "Surface " << searchableSurface::name()
                    << " does not have compact point numbering."
                    << " Of " << points()().size() << " only " << nPoints
                    << " are used. This might give problems in some routines."
                    << endl;
            }


            // Random number generator. Bit dodgy since not exactly random ;-)
            Random rndGen(65431);

            // Slightly extended bb. Slightly off-centred just so on symmetric
            // geometry there are less face/edge aligned items.
            bb = bb.extend(rndGen, 1e-4);
            bb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
            bb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
        }

        scalar oldTol = indexedOctree<treeDataTriSurface>::perturbTol();
        indexedOctree<treeDataTriSurface>::perturbTol() = tolerance_;

        tree_.reset
        (
            new indexedOctree<treeDataTriSurface>
            (
                treeDataTriSurface(*this, tolerance_),
                bb,
                maxTreeDepth_,  // maxLevel
                10,             // leafsize
                3.0             // duplicity
            )
        );

        indexedOctree<treeDataTriSurface>::perturbTol() = oldTol;
    }

    return tree_();
}


const Foam::indexedOctree<Foam::treeDataEdge>&
Foam::triSurfaceMesh::edgeTree() const
{
    if (edgeTree_.empty())
    {
        // Boundary edges
        labelList bEdges
        (
            identity
            (
                nEdges()
               -nInternalEdges()
            )
          + nInternalEdges()
        );

        treeBoundBox bb(vector::zero, vector::zero);

        if (bEdges.size())
        {
            label nPoints;
            calcBounds(bb, nPoints);

            // Random number generator. Bit dodgy since not exactly random ;-)
            Random rndGen(65431);

            // Slightly extended bb. Slightly off-centred just so on symmetric
            // geometry there are less face/edge aligned items.

            bb = bb.extend(rndGen, 1e-4);
            bb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
            bb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
        }

        scalar oldTol = indexedOctree<treeDataTriSurface>::perturbTol();
        indexedOctree<treeDataEdge>::perturbTol() = tolerance_;

        edgeTree_.reset
        (
            new indexedOctree<treeDataEdge>
            (
                treeDataEdge
                (
                    false,          // cachebb
                    edges(),        // edges
                    localPoints(),  // points
                    bEdges          // selected edges
                ),
                bb,                 // bb
                maxTreeDepth_,      // maxLevel
                10,                 // leafsize
                3.0                 // duplicity
            )
        );

        indexedOctree<treeDataEdge>::perturbTol() = oldTol;
    }
    return edgeTree_();
}


Foam::scalar Foam::triSurfaceMesh::tolerance() const
{
    return tolerance_;
}


const Foam::wordList& Foam::triSurfaceMesh::regions() const
{
    if (regions_.empty())
    {
        regions_.setSize(patches().size());
        forAll(regions_, regionI)
        {
            regions_[regionI] = patches()[regionI].name();
        }
    }
    return regions_;
}


// Find out if surface is closed.
bool Foam::triSurfaceMesh::hasVolumeType() const
{
    if (surfaceClosed_ == -1)
    {
        if (isSurfaceClosed())
        {
            surfaceClosed_ = 1;
        }
        else
        {
            surfaceClosed_ = 0;
        }
    }

    return surfaceClosed_ == 1;
}


void Foam::triSurfaceMesh::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& info
) const
{
    const indexedOctree<treeDataTriSurface>& octree = tree();

    info.setSize(samples.size());

    scalar oldTol = indexedOctree<treeDataTriSurface>::perturbTol();
    indexedOctree<treeDataTriSurface>::perturbTol() = tolerance_;

    forAll(samples, i)
    {
        static_cast<pointIndexHit&>(info[i]) = octree.findNearest
        (
            samples[i],
            nearestDistSqr[i]
        );
    }

    indexedOctree<treeDataTriSurface>::perturbTol() = oldTol;
}


void Foam::triSurfaceMesh::findLine
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    const indexedOctree<treeDataTriSurface>& octree = tree();

    info.setSize(start.size());

    scalar oldTol = indexedOctree<treeDataTriSurface>::perturbTol();
    indexedOctree<treeDataTriSurface>::perturbTol() = tolerance_;

    forAll(start, i)
    {
        static_cast<pointIndexHit&>(info[i]) = octree.findLine
        (
            start[i],
            end[i]
        );
    }

    indexedOctree<treeDataTriSurface>::perturbTol() = oldTol;
}


void Foam::triSurfaceMesh::findLineAny
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    const indexedOctree<treeDataTriSurface>& octree = tree();

    info.setSize(start.size());

    scalar oldTol = indexedOctree<treeDataTriSurface>::perturbTol();
    indexedOctree<treeDataTriSurface>::perturbTol() = tolerance_;

    forAll(start, i)
    {
        static_cast<pointIndexHit&>(info[i]) = octree.findLineAny
        (
            start[i],
            end[i]
        );
    }

    indexedOctree<treeDataTriSurface>::perturbTol() = oldTol;
}


void Foam::triSurfaceMesh::findLineAll
(
    const pointField& start,
    const pointField& end,
    List<List<pointIndexHit> >& info
) const
{
    const indexedOctree<treeDataTriSurface>& octree = tree();

    info.setSize(start.size());

    scalar oldTol = indexedOctree<treeDataTriSurface>::perturbTol();
    indexedOctree<treeDataTriSurface>::perturbTol() = tolerance_;

    // Work array
    DynamicList<pointIndexHit, 1, 1> hits;

    // Tolerances:
    // To find all intersections we add a small vector to the last intersection
    // This is chosen such that
    // - it is significant (SMALL is smallest representative relative tolerance;
    //   we need something bigger since we're doing calculations)
    // - if the start-end vector is zero we still progress
    const vectorField dirVec(end-start);
    const vectorField smallVec
    (
        indexedOctree<treeDataTriSurface>::perturbTol()*dirVec
      + vector(ROOTVSMALL,ROOTVSMALL,ROOTVSMALL)
    );

    forAll(start, pointI)
    {
        // See if any intersection between pt and end
        pointIndexHit inter = octree.findLine(start[pointI], end[pointI]);

        if (inter.hit())
        {
            hits.clear();
            hits.append(inter);

            getNextIntersections
            (
                octree,
                start[pointI],
                end[pointI],
                smallVec[pointI],
                hits
            );

            info[pointI].transfer(hits);
        }
        else
        {
            info[pointI].clear();
        }
    }

    indexedOctree<treeDataTriSurface>::perturbTol() = oldTol;
}


void Foam::triSurfaceMesh::getRegion
(
    const List<pointIndexHit>& info,
    labelList& region
) const
{
    region.setSize(info.size());
    forAll(info, i)
    {
        if (info[i].hit())
        {
            region[i] = triSurface::operator[](info[i].index()).region();
        }
        else
        {
            region[i] = -1;
        }
    }
}


void Foam::triSurfaceMesh::getNormal
(
    const List<pointIndexHit>& info,
    vectorField& normal
) const
{
    const triSurface& s = static_cast<const triSurface&>(*this);

    normal.setSize(info.size());

    if (minQuality_ >= 0)
    {
        // Make sure we don't use triangles with low quality since
        // normal is not reliable.

        const labelListList& faceFaces = s.faceFaces();

        forAll(info, i)
        {
            if (info[i].hit())
            {
                label faceI = info[i].index();
                normal[i] = s[faceI].normal(points());

                scalar qual = s[faceI].tri(points()).quality();

                if (qual < minQuality_)
                {
                    // Search neighbouring triangles
                    const labelList& fFaces = faceFaces[faceI];

                    forAll(fFaces, j)
                    {
                        label nbrI = fFaces[j];
                        scalar nbrQual = s[nbrI].tri(points()).quality();
                        if (nbrQual > qual)
                        {
                            qual = nbrQual;
                            normal[i] = s[nbrI].normal(points());
                        }
                    }
                }

                normal[i] /= mag(normal[i]) + VSMALL;
            }
            else
            {
                // Set to what?
                normal[i] = vector::zero;
            }
        }
    }
    else
    {
        forAll(info, i)
        {
            if (info[i].hit())
            {
                label faceI = info[i].index();
                //- Cached:
                //normal[i] = faceNormals()[faceI];

                //- Uncached
                normal[i] = s[faceI].normal(points());
                normal[i] /= mag(normal[i]) + VSMALL;
            }
            else
            {
                // Set to what?
                normal[i] = vector::zero;
            }
        }
    }
}


void Foam::triSurfaceMesh::setField(const labelList& values)
{
    autoPtr<triSurfaceLabelField> fldPtr
    (
        new triSurfaceLabelField
        (
            IOobject
            (
                "values",
                objectRegistry::time().timeName(),  // instance
                "triSurface",                       // local
                *this,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            *this,
            dimless,
            labelField(values)
        )
    );

    // Store field on triMesh
    fldPtr.ptr()->store();
}


void Foam::triSurfaceMesh::getField
(
    const List<pointIndexHit>& info,
    labelList& values
) const
{
    if (foundObject<triSurfaceLabelField>("values"))
    {
        values.setSize(info.size());

        const triSurfaceLabelField& fld = lookupObject<triSurfaceLabelField>
        (
            "values"
        );

        forAll(info, i)
        {
            if (info[i].hit())
            {
                values[i] = fld[info[i].index()];
            }
        }
    }
}


void Foam::triSurfaceMesh::getVolumeType
(
    const pointField& points,
    List<volumeType>& volType
) const
{
    volType.setSize(points.size());

    scalar oldTol = indexedOctree<treeDataTriSurface>::perturbTol();
    indexedOctree<treeDataTriSurface>::perturbTol() = tolerance_;

    forAll(points, pointI)
    {
        const point& pt = points[pointI];

        // - use cached volume type per each tree node
        // - cheat conversion since same values
        volType[pointI] = static_cast<volumeType>(tree().getVolumeType(pt));
    }

    indexedOctree<treeDataTriSurface>::perturbTol() = oldTol;
}


//- Write using given format, version and compression
bool Foam::triSurfaceMesh::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    fileName fullPath(searchableSurface::objectPath());

    if (!mkDir(fullPath.path()))
    {
        return false;
    }

    triSurface::write(fullPath);

    if (!isFile(fullPath))
    {
        return false;
    }

    //return objectRegistry::writeObject(fmt, ver, cmp);
    return true;
}


// ************************************************************************* //
