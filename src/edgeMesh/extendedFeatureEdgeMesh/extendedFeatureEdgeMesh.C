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

#include "extendedFeatureEdgeMesh.H"
#include "surfaceFeatures.H"
#include "triSurface.H"
#include "Random.H"
#include "Time.H"
#include "meshTools.H"
#include "ListListOps.H"
#include "OFstream.H"
#include "IFstream.H"
#include "unitConversion.H"
#include "DynamicField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(extendedFeatureEdgeMesh, 0);

    template<>
    const char* Foam::NamedEnum
    <
        Foam::extendedFeatureEdgeMesh::pointStatus,
        4
    >::names[] =
    {
        "convex",
        "concave",
        "mixed",
        "nonFeature"
    };

    template<>
    const char* Foam::NamedEnum
    <
        Foam::extendedFeatureEdgeMesh::edgeStatus,
        6
    >::names[] =
    {
        "external",
        "internal",
        "flat",
        "open",
        "multiple",
        "none"
    };
}

const Foam::NamedEnum<Foam::extendedFeatureEdgeMesh::pointStatus, 4>
    Foam::extendedFeatureEdgeMesh::pointStatusNames_;

const Foam::NamedEnum<Foam::extendedFeatureEdgeMesh::edgeStatus, 6>
    Foam::extendedFeatureEdgeMesh::edgeStatusNames_;

Foam::scalar Foam::extendedFeatureEdgeMesh::cosNormalAngleTol_ =
    Foam::cos(degToRad(0.1));


Foam::label Foam::extendedFeatureEdgeMesh::convexStart_ = 0;


Foam::label Foam::extendedFeatureEdgeMesh::externalStart_ = 0;


Foam::label Foam::extendedFeatureEdgeMesh::nPointTypes = 4;


Foam::label Foam::extendedFeatureEdgeMesh::nEdgeTypes = 5;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::extendedFeatureEdgeMesh::extendedFeatureEdgeMesh(const IOobject& io)
:
    regIOobject(io),
    edgeMesh(pointField(0), edgeList(0)),
    concaveStart_(0),
    mixedStart_(0),
    nonFeatureStart_(0),
    internalStart_(0),
    flatStart_(0),
    openStart_(0),
    multipleStart_(0),
    normals_(0),
    edgeDirections_(0),
    edgeNormals_(0),
    featurePointNormals_(0),
    featurePointEdges_(0),
    regionEdges_(0),
    pointTree_(),
    edgeTree_(),
    edgeTreesByType_()
{
    if
    (
        io.readOpt() == IOobject::MUST_READ
     || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        if (readOpt() == IOobject::MUST_READ_IF_MODIFIED)
        {
            WarningIn
            (
                "extendedFeatureEdgeMesh::extendedFeatureEdgeMesh"
                "(const IOobject&)"
            )   << "Specified IOobject::MUST_READ_IF_MODIFIED but class"
                << " does not support automatic rereading."
                << endl;
        }

        Istream& is = readStream(typeName);

        is  >> *this
            >> concaveStart_
            >> mixedStart_
            >> nonFeatureStart_
            >> internalStart_
            >> flatStart_
            >> openStart_
            >> multipleStart_
            >> normals_
            >> edgeNormals_
            >> featurePointNormals_
            >> featurePointEdges_
            >> regionEdges_;

        close();

        {
            // Calculate edgeDirections

            const edgeList& eds(edges());

            const pointField& pts(points());

            edgeDirections_.setSize(eds.size());

            forAll(eds, eI)
            {
                edgeDirections_[eI] = eds[eI].vec(pts);
            }

            edgeDirections_ /= mag(edgeDirections_);
        }
    }

    if (debug)
    {
        Pout<< "extendedFeatureEdgeMesh::extendedFeatureEdgeMesh :"
            << " constructed from IOobject :"
            << " points:" << points().size()
            << " edges:" << edges().size()
            << endl;
    }
}


Foam::extendedFeatureEdgeMesh::extendedFeatureEdgeMesh
(
    const IOobject& io,
    const extendedFeatureEdgeMesh& fem
)
:
    regIOobject(io),
    edgeMesh(fem),
    concaveStart_(fem.concaveStart()),
    mixedStart_(fem.mixedStart()),
    nonFeatureStart_(fem.nonFeatureStart()),
    internalStart_(fem.internalStart()),
    flatStart_(fem.flatStart()),
    openStart_(fem.openStart()),
    multipleStart_(fem.multipleStart()),
    normals_(fem.normals()),
    edgeDirections_(fem.edgeDirections()),
    edgeNormals_(fem.edgeNormals()),
    featurePointNormals_(fem.featurePointNormals()),
    featurePointEdges_(fem.featurePointEdges()),
    regionEdges_(fem.regionEdges()),
    pointTree_(),
    edgeTree_(),
    edgeTreesByType_()
{}


Foam::extendedFeatureEdgeMesh::extendedFeatureEdgeMesh
(
    const IOobject& io,
    const Xfer<pointField>& pointLst,
    const Xfer<edgeList>& edgeLst
)
:
    regIOobject(io),
    edgeMesh(pointLst, edgeLst),
    concaveStart_(0),
    mixedStart_(0),
    nonFeatureStart_(0),
    internalStart_(0),
    flatStart_(0),
    openStart_(0),
    multipleStart_(0),
    normals_(0),
    edgeDirections_(0),
    edgeNormals_(0),
    featurePointNormals_(0),
    featurePointEdges_(0),
    regionEdges_(0),
    pointTree_(),
    edgeTree_(),
    edgeTreesByType_()
{}


Foam::extendedFeatureEdgeMesh::extendedFeatureEdgeMesh
(
    const surfaceFeatures& sFeat,
    const objectRegistry& obr,
    const fileName& sFeatFileName
)
:
    regIOobject
    (
        IOobject
        (
            sFeatFileName,
            obr.time().constant(),
            "extendedFeatureEdgeMesh",
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    edgeMesh(pointField(0), edgeList(0)),
    concaveStart_(-1),
    mixedStart_(-1),
    nonFeatureStart_(-1),
    internalStart_(-1),
    flatStart_(-1),
    openStart_(-1),
    multipleStart_(-1),
    normals_(0),
    edgeDirections_(0),
    edgeNormals_(0),
    featurePointNormals_(0),
    featurePointEdges_(0),
    regionEdges_(0),
    pointTree_(),
    edgeTree_(),
    edgeTreesByType_()
{
    // Extract and reorder the data from surfaceFeatures
    const triSurface& surf = sFeat.surface();
    const labelList& featureEdges = sFeat.featureEdges();
    const labelList& featurePoints = sFeat.featurePoints();

    // Get a labelList of all the featureEdges that are region edges
    const labelList regionFeatureEdges(identity(sFeat.nRegionEdges()));

    sortPointsAndEdges
    (
        surf,
        featureEdges,
        regionFeatureEdges,
        featurePoints
    );
}


Foam::extendedFeatureEdgeMesh::extendedFeatureEdgeMesh
(
    const IOobject& io,
    const PrimitivePatch<face, List, pointField, point>& surf,
    const labelList& featureEdges,
    const labelList& regionFeatureEdges,
    const labelList& featurePoints
)
:
    regIOobject(io),
    edgeMesh(pointField(0), edgeList(0)),
    concaveStart_(-1),
    mixedStart_(-1),
    nonFeatureStart_(-1),
    internalStart_(-1),
    flatStart_(-1),
    openStart_(-1),
    multipleStart_(-1),
    normals_(0),
    edgeDirections_(0),
    edgeNormals_(0),
    featurePointNormals_(0),
    featurePointEdges_(0),
    regionEdges_(0),
    pointTree_(),
    edgeTree_(),
    edgeTreesByType_()
{
    sortPointsAndEdges
    (
        surf,
        featureEdges,
        regionFeatureEdges,
        featurePoints
    );
}


Foam::extendedFeatureEdgeMesh::extendedFeatureEdgeMesh
(
    const IOobject& io,
    const pointField& pts,
    const edgeList& eds,
    label concaveStart,
    label mixedStart,
    label nonFeatureStart,
    label internalStart,
    label flatStart,
    label openStart,
    label multipleStart,
    const vectorField& normals,
    const vectorField& edgeDirections,
    const labelListList& edgeNormals,
    const labelListList& featurePointNormals,
    const labelListList& featurePointEdges,
    const labelList& regionEdges
)
:
    regIOobject(io),
    edgeMesh(pts, eds),
    concaveStart_(concaveStart),
    mixedStart_(mixedStart),
    nonFeatureStart_(nonFeatureStart),
    internalStart_(internalStart),
    flatStart_(flatStart),
    openStart_(openStart),
    multipleStart_(multipleStart),
    normals_(normals),
    edgeDirections_(edgeDirections),
    edgeNormals_(edgeNormals),
    featurePointNormals_(featurePointNormals),
    featurePointEdges_(featurePointEdges),
    regionEdges_(regionEdges),
    pointTree_(),
    edgeTree_(),
    edgeTreesByType_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::extendedFeatureEdgeMesh::~extendedFeatureEdgeMesh()
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::extendedFeatureEdgeMesh::pointStatus
Foam::extendedFeatureEdgeMesh::classifyFeaturePoint
(
    label ptI
) const
{
    labelList ptEds(pointEdges()[ptI]);

    label nPtEds = ptEds.size();
    label nExternal = 0;
    label nInternal = 0;

    if (nPtEds == 0)
    {
        // There are no edges attached to the point, this is a problem
        return NONFEATURE;
    }

    forAll(ptEds, i)
    {
        edgeStatus edStat = getEdgeStatus(ptEds[i]);

        if (edStat == EXTERNAL)
        {
            nExternal++;
        }
        else if (edStat == INTERNAL)
        {
            nInternal++;
        }
    }

    if (nExternal == nPtEds)
    {
        return CONVEX;
    }
    else if (nInternal == nPtEds)
    {
        return CONCAVE;
    }
    else
    {
        return MIXED;
    }
}


Foam::extendedFeatureEdgeMesh::edgeStatus
Foam::extendedFeatureEdgeMesh::classifyEdge
(
    const List<vector>& norms,
    const labelList& edNorms,
    const vector& fC0tofC1
) const
{
    label nEdNorms = edNorms.size();

    if (nEdNorms == 1)
    {
        return OPEN;
    }
    else if (nEdNorms == 2)
    {
        const vector n0(norms[edNorms[0]]);
        const vector n1(norms[edNorms[1]]);

        if ((n0 & n1) > cosNormalAngleTol_)
        {
            return FLAT;
        }
        else if ((fC0tofC1 & n0) > 0.0)
        {
            return INTERNAL;
        }
        else
        {
            return EXTERNAL;
        }
    }
    else if (nEdNorms > 2)
    {
        return MULTIPLE;
    }
    else
    {
        // There is a problem - the edge has no normals
        return NONE;
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::extendedFeatureEdgeMesh::nearestFeaturePoint
(
    const point& sample,
    scalar searchDistSqr,
    pointIndexHit& info
) const
{
    info = pointTree().findNearest
    (
        sample,
        searchDistSqr
    );
}


void Foam::extendedFeatureEdgeMesh::nearestFeatureEdge
(
    const point& sample,
    scalar searchDistSqr,
    pointIndexHit& info
) const
{
    info = edgeTree().findNearest
    (
        sample,
        searchDistSqr
    );
}


void Foam::extendedFeatureEdgeMesh::nearestFeatureEdge
(
    const pointField& samples,
    const scalarField& searchDistSqr,
    List<pointIndexHit>& info
) const
{
    info.setSize(samples.size());

    forAll(samples, i)
    {
        nearestFeatureEdge
        (
            samples[i],
            searchDistSqr[i],
            info[i]
        );
    }
}


void Foam::extendedFeatureEdgeMesh::nearestFeatureEdgeByType
(
    const point& sample,
    const scalarField& searchDistSqr,
    List<pointIndexHit>& info
) const
{
    const PtrList<indexedOctree<treeDataEdge> >& edgeTrees = edgeTreesByType();

    info.setSize(edgeTrees.size());

    labelList sliceStarts(edgeTrees.size());

    sliceStarts[0] = externalStart_;
    sliceStarts[1] = internalStart_;
    sliceStarts[2] = flatStart_;
    sliceStarts[3] = openStart_;
    sliceStarts[4] = multipleStart_;

    forAll(edgeTrees, i)
    {
        info[i] = edgeTrees[i].findNearest
        (
            sample,
            searchDistSqr[i]
        );

        // The index returned by the indexedOctree is local to the slice of
        // edges it was supplied with, return the index to the value in the
        // complete edge list

        info[i].setIndex(info[i].index() + sliceStarts[i]);
    }
}


void Foam::extendedFeatureEdgeMesh::allNearestFeaturePoints
(
    const point& sample,
    scalar searchRadiusSqr,
    List<pointIndexHit>& info
) const
{
    // Pick up all the feature points that intersect the search sphere
    labelList elems = pointTree().findSphere
    (
        sample,
        searchRadiusSqr
    );

    DynamicList<pointIndexHit> dynPointHit(elems.size());

    forAll(elems, elemI)
    {
        label index = elems[elemI];
        label ptI = pointTree().shapes().pointLabels()[index];
        const point& pt = points()[ptI];

        pointIndexHit nearHit(true, pt, index);

        dynPointHit.append(nearHit);
    }

    info.transfer(dynPointHit);
}


void Foam::extendedFeatureEdgeMesh::allNearestFeatureEdges
(
    const point& sample,
    const scalar searchRadiusSqr,
    List<pointIndexHit>& info
) const
{
    const PtrList<indexedOctree<treeDataEdge> >& edgeTrees = edgeTreesByType();

    info.setSize(edgeTrees.size());

    labelList sliceStarts(edgeTrees.size());

    sliceStarts[0] = externalStart_;
    sliceStarts[1] = internalStart_;
    sliceStarts[2] = flatStart_;
    sliceStarts[3] = openStart_;
    sliceStarts[4] = multipleStart_;

    DynamicList<pointIndexHit> dynEdgeHit(edgeTrees.size()*3);

    // Loop over all the feature edge types
    forAll(edgeTrees, i)
    {
        // Pick up all the edges that intersect the search sphere
        labelList elems = edgeTrees[i].findSphere
        (
            sample,
            searchRadiusSqr
        );

        forAll(elems, elemI)
        {
            label index = elems[elemI];
            label edgeI = edgeTrees[i].shapes().edgeLabels()[index];
            const edge& e = edges()[edgeI];

            pointHit hitPoint = e.line(points()).nearestDist(sample);

            label hitIndex = index + sliceStarts[i];

            pointIndexHit nearHit
            (
                hitPoint.hit(),
                hitPoint.rawPoint(),
                hitIndex
            );

            dynEdgeHit.append(nearHit);
        }
    }

    info.transfer(dynEdgeHit);
}


const Foam::indexedOctree<Foam::treeDataPoint>&
Foam::extendedFeatureEdgeMesh::pointTree() const
{
    if (pointTree_.empty())
    {
        Random rndGen(17301893);

        // Slightly extended bb. Slightly off-centred just so on symmetric
        // geometry there are less face/edge aligned items.
        treeBoundBox bb
        (
            treeBoundBox(points()).extend(rndGen, 1e-4)
        );

        bb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
        bb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

        const labelList featurePointLabels = identity(nonFeatureStart_);

        pointTree_.reset
        (
            new indexedOctree<treeDataPoint>
            (
                treeDataPoint
                (
                    points(),
                    featurePointLabels
                ),
                bb,     // bb
                8,      // maxLevel
                10,     // leafsize
                3.0     // duplicity
            )
        );
    }

    return pointTree_();
}


const Foam::indexedOctree<Foam::treeDataEdge>&
Foam::extendedFeatureEdgeMesh::edgeTree() const
{
    if (edgeTree_.empty())
    {
        Random rndGen(17301893);

        // Slightly extended bb. Slightly off-centred just so on symmetric
        // geometry there are less face/edge aligned items.
        treeBoundBox bb
        (
            treeBoundBox(points()).extend(rndGen, 1e-4)
        );

        bb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
        bb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

        labelList allEdges(identity(edges().size()));

        edgeTree_.reset
        (
            new indexedOctree<treeDataEdge>
            (
                treeDataEdge
                (
                    false,          // cachebb
                    edges(),        // edges
                    points(),       // points
                    allEdges        // selected edges
                ),
                bb,     // bb
                8,      // maxLevel
                10,     // leafsize
                3.0     // duplicity
            )
        );
    }

    return edgeTree_();
}


const Foam::PtrList<Foam::indexedOctree<Foam::treeDataEdge> >&
Foam::extendedFeatureEdgeMesh::edgeTreesByType() const
{
    if (edgeTreesByType_.size() == 0)
    {
        edgeTreesByType_.setSize(nEdgeTypes);

        Random rndGen(872141);

        // Slightly extended bb. Slightly off-centred just so on symmetric
        // geometry there are less face/edge aligned items.
        treeBoundBox bb
        (
            treeBoundBox(points()).extend(rndGen, 1e-4)
        );

        bb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
        bb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

        labelListList sliceEdges(nEdgeTypes);

        // External edges
        sliceEdges[0] =
            identity(internalStart_ - externalStart_) + externalStart_;

        // Internal edges
        sliceEdges[1] = identity(flatStart_ - internalStart_) + internalStart_;

        // Flat edges
        sliceEdges[2] = identity(openStart_ - flatStart_) + flatStart_;

        // Open edges
        sliceEdges[3] = identity(multipleStart_ - openStart_) + openStart_;

        // Multiple edges
        sliceEdges[4] =
            identity(edges().size() - multipleStart_) + multipleStart_;

        forAll(edgeTreesByType_, i)
        {
            edgeTreesByType_.set
            (
                i,
                new indexedOctree<treeDataEdge>
                (
                    treeDataEdge
                    (
                        false,          // cachebb
                        edges(),        // edges
                        points(),       // points
                        sliceEdges[i]   // selected edges
                    ),
                    bb,     // bb
                    8,      // maxLevel
                    10,     // leafsize
                    3.0     // duplicity
                )
            );
        }
    }

    return edgeTreesByType_;
}


void Foam::extendedFeatureEdgeMesh::add(const extendedFeatureEdgeMesh& fem)
{
    // Points
    // ~~~~~~

    // From current points into combined points
    labelList reversePointMap(points().size());
    labelList reverseFemPointMap(fem.points().size());

    label newPointI = 0;
    for (label i = 0; i < concaveStart(); i++)
    {
        reversePointMap[i] = newPointI++;
    }
    for (label i = 0; i < fem.concaveStart(); i++)
    {
        reverseFemPointMap[i] = newPointI++;
    }

    // Concave
    label newConcaveStart = newPointI;
    for (label i = concaveStart(); i < mixedStart(); i++)
    {
        reversePointMap[i] = newPointI++;
    }
    for (label i = fem.concaveStart(); i < fem.mixedStart(); i++)
    {
        reverseFemPointMap[i] = newPointI++;
    }

    // Mixed
    label newMixedStart = newPointI;
    for (label i = mixedStart(); i < nonFeatureStart(); i++)
    {
        reversePointMap[i] = newPointI++;
    }
    for (label i = fem.mixedStart(); i < fem.nonFeatureStart(); i++)
    {
        reverseFemPointMap[i] = newPointI++;
    }

    // Non-feature
    label newNonFeatureStart = newPointI;
    for (label i = nonFeatureStart(); i < points().size(); i++)
    {
        reversePointMap[i] = newPointI++;
    }
    for (label i = fem.nonFeatureStart(); i < fem.points().size(); i++)
    {
        reverseFemPointMap[i] = newPointI++;
    }

    pointField newPoints(newPointI);
    newPoints.rmap(points(), reversePointMap);
    newPoints.rmap(fem.points(), reverseFemPointMap);


    // Edges
    // ~~~~~

    // From current edges into combined edges
    labelList reverseEdgeMap(edges().size());
    labelList reverseFemEdgeMap(fem.edges().size());

    // External
    label newEdgeI = 0;
    for (label i = 0; i < internalStart(); i++)
    {
        reverseEdgeMap[i] = newEdgeI++;
    }
    for (label i = 0; i < fem.internalStart(); i++)
    {
        reverseFemEdgeMap[i] = newEdgeI++;
    }

    // Internal
    label newInternalStart = newEdgeI;
    for (label i = internalStart(); i < flatStart(); i++)
    {
        reverseEdgeMap[i] = newEdgeI++;
    }
    for (label i = fem.internalStart(); i < fem.flatStart(); i++)
    {
        reverseFemEdgeMap[i] = newEdgeI++;
    }

    // Flat
    label newFlatStart = newEdgeI;
    for (label i = flatStart(); i < openStart(); i++)
    {
        reverseEdgeMap[i] = newEdgeI++;
    }
    for (label i = fem.flatStart(); i < fem.openStart(); i++)
    {
        reverseFemEdgeMap[i] = newEdgeI++;
    }

    // Open
    label newOpenStart = newEdgeI;
    for (label i = openStart(); i < multipleStart(); i++)
    {
        reverseEdgeMap[i] = newEdgeI++;
    }
    for (label i = fem.openStart(); i < fem.multipleStart(); i++)
    {
        reverseFemEdgeMap[i] = newEdgeI++;
    }

    // Multiple
    label newMultipleStart = newEdgeI;
    for (label i = multipleStart(); i < edges().size(); i++)
    {
        reverseEdgeMap[i] = newEdgeI++;
    }
    for (label i = fem.multipleStart(); i < fem.edges().size(); i++)
    {
        reverseFemEdgeMap[i] = newEdgeI++;
    }

    edgeList newEdges(newEdgeI);
    forAll(edges(), i)
    {
        const edge& e = edges()[i];
        newEdges[reverseEdgeMap[i]] = edge
        (
            reversePointMap[e[0]],
            reversePointMap[e[1]]
        );
    }
    forAll(fem.edges(), i)
    {
        const edge& e = fem.edges()[i];
        newEdges[reverseFemEdgeMap[i]] = edge
        (
            reverseFemPointMap[e[0]],
            reverseFemPointMap[e[1]]
        );
    }

    pointField newEdgeDirections(newEdgeI);
    newEdgeDirections.rmap(edgeDirections(), reverseEdgeMap);
    newEdgeDirections.rmap(fem.edgeDirections(), reverseFemEdgeMap);




    // Normals
    // ~~~~~~~

    // Combine normals
    DynamicField<point> newNormals(normals().size()+fem.normals().size());
    newNormals.append(normals());
    newNormals.append(fem.normals());


    // Combine and re-index into newNormals
    labelListList newEdgeNormals(edgeNormals().size()+fem.edgeNormals().size());
    UIndirectList<labelList>(newEdgeNormals, reverseEdgeMap) =
        edgeNormals();
    UIndirectList<labelList>(newEdgeNormals, reverseFemEdgeMap) =
        fem.edgeNormals();
    forAll(reverseFemEdgeMap, i)
    {
        label mapI = reverseFemEdgeMap[i];
        labelList& en = newEdgeNormals[mapI];
        forAll(en, j)
        {
            en[j] += normals().size();
        }
    }


    // Combine and re-index into newFeaturePointNormals
    labelListList newFeaturePointNormals
    (
       featurePointNormals().size()
     + fem.featurePointNormals().size()
    );

    // Note: featurePointNormals only go up to nonFeatureStart
    UIndirectList<labelList>
    (
        newFeaturePointNormals,
        SubList<label>(reversePointMap, featurePointNormals().size())
    ) = featurePointNormals();
    UIndirectList<labelList>
    (
        newFeaturePointNormals,
        SubList<label>(reverseFemPointMap, fem.featurePointNormals().size())
    ) = fem.featurePointNormals();
    forAll(fem.featurePointNormals(), i)
    {
        label mapI = reverseFemPointMap[i];
        labelList& fn = newFeaturePointNormals[mapI];
        forAll(fn, j)
        {
            fn[j] += normals().size();
        }
    }


    // Combine regionEdges
    DynamicList<label> newRegionEdges
    (
        regionEdges().size()
      + fem.regionEdges().size()
    );
    forAll(regionEdges(), i)
    {
        newRegionEdges.append(reverseEdgeMap[regionEdges()[i]]);
    }
    forAll(fem.regionEdges(), i)
    {
        newRegionEdges.append(reverseFemEdgeMap[fem.regionEdges()[i]]);
    }


    // Assign
    // ~~~~~~

    // Transfer
    concaveStart_ = newConcaveStart;
    mixedStart_ = newMixedStart;
    nonFeatureStart_ = newNonFeatureStart;

    // Reset points and edges
    reset(xferMove(newPoints), newEdges.xfer());

    // Transfer
    internalStart_ = newInternalStart;
    flatStart_ = newFlatStart;
    openStart_ = newOpenStart;
    multipleStart_ = newMultipleStart;

    edgeDirections_.transfer(newEdgeDirections);

    normals_.transfer(newNormals);
    edgeNormals_.transfer(newEdgeNormals);
    featurePointNormals_.transfer(newFeaturePointNormals);

    regionEdges_.transfer(newRegionEdges);

    pointTree_.clear();
    edgeTree_.clear();
    edgeTreesByType_.clear();
}


void Foam::extendedFeatureEdgeMesh::flipNormals()
{
    // Points
    // ~~~~~~

    // From current points into new points
    labelList reversePointMap(identity(points().size()));

    // Flip convex and concave points

    label newPointI = 0;
    // Concave points become convex
    for (label i = concaveStart(); i < mixedStart(); i++)
    {
        reversePointMap[i] = newPointI++;
    }
    // Convex points become concave
    label newConcaveStart = newPointI;
    for (label i = 0; i < concaveStart(); i++)
    {
        reversePointMap[i] = newPointI++;
    }


    // Edges
    // ~~~~~~

    // From current edges into new edges
    labelList reverseEdgeMap(identity(edges().size()));

    // Flip external and internal edges

    label newEdgeI = 0;
    // Internal become external
    for (label i = internalStart(); i < flatStart(); i++)
    {
        reverseEdgeMap[i] = newEdgeI++;
    }
    // External become internal
    label newInternalStart = newEdgeI;
    for (label i = 0; i < internalStart(); i++)
    {
        reverseEdgeMap[i] = newEdgeI++;
    }


    pointField newPoints(points().size());
    newPoints.rmap(points(), reversePointMap);

    edgeList newEdges(edges().size());
    forAll(edges(), i)
    {
        const edge& e = edges()[i];
        newEdges[reverseEdgeMap[i]] = edge
        (
            reversePointMap[e[0]],
            reversePointMap[e[1]]
        );
    }


    // Normals are flipped
    // ~~~~~~~~~~~~~~~~~~~

    pointField newEdgeDirections(edges().size());
    newEdgeDirections.rmap(-1.0*edgeDirections(), reverseEdgeMap);

    pointField newNormals(-1.0*normals());

    labelListList newEdgeNormals(edgeNormals().size());
    UIndirectList<labelList>(newEdgeNormals, reverseEdgeMap) = edgeNormals();

    labelListList newFeaturePointNormals(featurePointNormals().size());

    // Note: featurePointNormals only go up to nonFeatureStart
    UIndirectList<labelList>
    (
        newFeaturePointNormals,
        SubList<label>(reversePointMap, featurePointNormals().size())
    ) = featurePointNormals();

    labelList newRegionEdges(regionEdges().size());
    forAll(regionEdges(), i)
    {
        newRegionEdges[i] = reverseEdgeMap[regionEdges()[i]];
    }

    // Transfer
    concaveStart_ = newConcaveStart;

    // Reset points and edges
    reset(xferMove(newPoints), newEdges.xfer());

    // Transfer
    internalStart_ = newInternalStart;

    edgeDirections_.transfer(newEdgeDirections);
    normals_.transfer(newNormals);
    edgeNormals_.transfer(newEdgeNormals);
    featurePointNormals_.transfer(newFeaturePointNormals);
    regionEdges_.transfer(newRegionEdges);

    pointTree_.clear();
    edgeTree_.clear();
    edgeTreesByType_.clear();
}
//XXXXX

void Foam::extendedFeatureEdgeMesh::writeObj
(
    const fileName& prefix
) const
{
    Info<< nl << "Writing extendedFeatureEdgeMesh components to " << prefix
        << endl;

    label verti = 0;

    edgeMesh::write(prefix + "_edgeMesh.obj");

    OFstream convexFtPtStr(prefix + "_convexFeaturePts.obj");
    Info<< "Writing convex feature points to " << convexFtPtStr.name() << endl;

    for(label i = 0; i < concaveStart_; i++)
    {
        meshTools::writeOBJ(convexFtPtStr, points()[i]);
    }

    OFstream concaveFtPtStr(prefix + "_concaveFeaturePts.obj");
    Info<< "Writing concave feature points to "
        << concaveFtPtStr.name() << endl;

    for(label i = concaveStart_; i < mixedStart_; i++)
    {
        meshTools::writeOBJ(concaveFtPtStr, points()[i]);
    }

    OFstream mixedFtPtStr(prefix + "_mixedFeaturePts.obj");
    Info<< "Writing mixed feature points to " << mixedFtPtStr.name() << endl;

    for(label i = mixedStart_; i < nonFeatureStart_; i++)
    {
        meshTools::writeOBJ(mixedFtPtStr, points()[i]);
    }

    OFstream mixedFtPtStructureStr(prefix + "_mixedFeaturePtsStructure.obj");
    Info<< "Writing mixed feature point structure to "
        << mixedFtPtStructureStr.name() << endl;

    verti = 0;
    for(label i = mixedStart_; i < nonFeatureStart_; i++)
    {
        const labelList& ptEds = pointEdges()[i];

        forAll(ptEds, j)
        {
            const edge& e = edges()[ptEds[j]];

            meshTools::writeOBJ(mixedFtPtStructureStr, points()[e[0]]); verti++;
            meshTools::writeOBJ(mixedFtPtStructureStr, points()[e[1]]); verti++;
            mixedFtPtStructureStr << "l " << verti-1 << ' ' << verti << endl;
        }
    }

    OFstream externalStr(prefix + "_externalEdges.obj");
    Info<< "Writing external edges to " << externalStr.name() << endl;

    verti = 0;
    for (label i = externalStart_; i < internalStart_; i++)
    {
        const edge& e = edges()[i];

        meshTools::writeOBJ(externalStr, points()[e[0]]); verti++;
        meshTools::writeOBJ(externalStr, points()[e[1]]); verti++;
        externalStr << "l " << verti-1 << ' ' << verti << endl;
    }

    OFstream internalStr(prefix + "_internalEdges.obj");
    Info<< "Writing internal edges to " << internalStr.name() << endl;

    verti = 0;
    for (label i = internalStart_; i < flatStart_; i++)
    {
        const edge& e = edges()[i];

        meshTools::writeOBJ(internalStr, points()[e[0]]); verti++;
        meshTools::writeOBJ(internalStr, points()[e[1]]); verti++;
        internalStr << "l " << verti-1 << ' ' << verti << endl;
    }

    OFstream flatStr(prefix + "_flatEdges.obj");
    Info<< "Writing flat edges to " << flatStr.name() << endl;

    verti = 0;
    for (label i = flatStart_; i < openStart_; i++)
    {
        const edge& e = edges()[i];

        meshTools::writeOBJ(flatStr, points()[e[0]]); verti++;
        meshTools::writeOBJ(flatStr, points()[e[1]]); verti++;
        flatStr << "l " << verti-1 << ' ' << verti << endl;
    }

    OFstream openStr(prefix + "_openEdges.obj");
    Info<< "Writing open edges to " << openStr.name() << endl;

    verti = 0;
    for (label i = openStart_; i < multipleStart_; i++)
    {
        const edge& e = edges()[i];

        meshTools::writeOBJ(openStr, points()[e[0]]); verti++;
        meshTools::writeOBJ(openStr, points()[e[1]]); verti++;
        openStr << "l " << verti-1 << ' ' << verti << endl;
    }

    OFstream multipleStr(prefix + "_multipleEdges.obj");
    Info<< "Writing multiple edges to " << multipleStr.name() << endl;

    verti = 0;
    for (label i = multipleStart_; i < edges().size(); i++)
    {
        const edge& e = edges()[i];

        meshTools::writeOBJ(multipleStr, points()[e[0]]); verti++;
        meshTools::writeOBJ(multipleStr, points()[e[1]]); verti++;
        multipleStr << "l " << verti-1 << ' ' << verti << endl;
    }

    OFstream regionStr(prefix + "_regionEdges.obj");
    Info<< "Writing region edges to " << regionStr.name() << endl;

    verti = 0;
    forAll(regionEdges_, i)
    {
        const edge& e = edges()[regionEdges_[i]];

        meshTools::writeOBJ(regionStr, points()[e[0]]); verti++;
        meshTools::writeOBJ(regionStr, points()[e[1]]); verti++;
        regionStr << "l " << verti-1 << ' ' << verti << endl;
    }
}


bool Foam::extendedFeatureEdgeMesh::writeData(Ostream& os) const
{
    os  << "// points" << nl
        << points() << nl
        << "// edges" << nl
        << edges() << nl
        << "// concaveStart mixedStart nonFeatureStart" << nl
        << concaveStart_ << token::SPACE
        << mixedStart_ << token::SPACE
        << nonFeatureStart_ << nl
        << "// internalStart flatStart openStart multipleStart" << nl
        << internalStart_ << token::SPACE
        << flatStart_ << token::SPACE
        << openStart_ << token::SPACE
        << multipleStart_ << nl
        << "// normals" << nl
        << normals_ << nl
        << "// edgeNormals" << nl
        << edgeNormals_ << nl
        << "// featurePointNormals" << nl
        << featurePointNormals_ << nl
        << "// featurePointEdges" << nl
        << featurePointEdges_ << nl
        << "// regionEdges" << nl
        << regionEdges_
        << endl;

    return os.good();
}


// ************************************************************************* //
