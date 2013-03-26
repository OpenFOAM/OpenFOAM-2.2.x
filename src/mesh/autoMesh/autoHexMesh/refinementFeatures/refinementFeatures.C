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

#include "refinementFeatures.H"
#include "Time.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::refinementFeatures::read
(
    const objectRegistry& io,
    const PtrList<dictionary>& featDicts
)
{
    forAll(featDicts, featI)
    {
        const dictionary& dict = featDicts[featI];

        fileName featFileName(dict.lookup("file"));

        {
            IOobject featObj
            (
                featFileName,                       // name
                io.time().constant(),               // instance
                "triSurface",                       // local
                io.time(),                          // registry
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            );

            autoPtr<edgeMesh> eMeshPtr = edgeMesh::New(featObj.filePath());

            set
            (
                featI,
                new featureEdgeMesh
                (
                    featObj,
                    eMeshPtr->points(),
                    eMeshPtr->edges()
                )
            );
        }

        const featureEdgeMesh& eMesh = operator[](featI);

        //eMesh.mergePoints(meshRefiner_.mergeDistance());

        if (dict.found("levels"))
        {
            List<Tuple2<scalar, label> > distLevels(dict["levels"]);

            if (dict.size() < 1)
            {
                FatalErrorIn
                (
                    "refinementFeatures::read"
                    "(const objectRegistry&"
                    ", const PtrList<dictionary>&)"
                )   << " : levels should be at least size 1" << endl
                    << "levels : "  << dict["levels"]
                    << exit(FatalError);
            }

            distances_[featI].setSize(distLevels.size());
            levels_[featI].setSize(distLevels.size());

            forAll(distLevels, j)
            {
                distances_[featI][j] = distLevels[j].first();
                levels_[featI][j] = distLevels[j].second();

                // Check in incremental order
                if (j > 0)
                {
                    if
                    (
                        (distances_[featI][j] <= distances_[featI][j-1])
                     || (levels_[featI][j] > levels_[featI][j-1])
                    )
                    {
                        FatalErrorIn
                        (
                            "refinementFeatures::read"
                            "(const objectRegistry&"
                            ", const PtrList<dictionary>&)"
                        )   << " : Refinement should be specified in order"
                            << " of increasing distance"
                            << " (and decreasing refinement level)." << endl
                            << "Distance:" << distances_[featI][j]
                            << " refinementLevel:" << levels_[featI][j]
                            << exit(FatalError);
                    }
                }
            }
        }
        else
        {
            // Look up 'level' for single level
            levels_[featI] = labelList(1, readLabel(dict.lookup("level")));
            distances_[featI] = scalarField(1, 0.0);
        }

        Info<< "Refinement level according to distance to "
            << featFileName << " (" << eMesh.points().size() << " points, "
            << eMesh.edges().size() << " edges)." << endl;
        forAll(levels_[featI], j)
        {
            Info<< "    level " << levels_[featI][j]
                << " for all cells within " << distances_[featI][j]
                << " meter." << endl;
        }
    }
}


void Foam::refinementFeatures::buildTrees
(
    const label featI,
    const labelList& featurePoints
)
{
    const featureEdgeMesh& eMesh = operator[](featI);
    const pointField& points = eMesh.points();
    const edgeList& edges = eMesh.edges();

    // Calculate bb of all points
    treeBoundBox bb(points);

    // Random number generator. Bit dodgy since not exactly random ;-)
    Random rndGen(65431);

    // Slightly extended bb. Slightly off-centred just so on symmetric
    // geometry there are less face/edge aligned items.
    bb = bb.extend(rndGen, 1e-4);
    bb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
    bb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

    edgeTrees_.set
    (
        featI,
        new indexedOctree<treeDataEdge>
        (
            treeDataEdge
            (
                false,                  // do not cache bb
                edges,
                points,
                identity(edges.size())
            ),
            bb,     // overall search domain
            8,      // maxLevel
            10,     // leafsize
            3.0     // duplicity
        )
    );

    pointTrees_.set
    (
        featI,
        new indexedOctree<treeDataPoint>
        (
            treeDataPoint(points, featurePoints),
            bb,     // overall search domain
            8,      // maxLevel
            10,     // leafsize
            3.0     // duplicity
        )
    );
}



// Find maximum level of a shell.
void Foam::refinementFeatures::findHigherLevel
(
    const pointField& pt,
    const label featI,
    labelList& maxLevel
) const
{
    const labelList& levels = levels_[featI];

    const scalarField& distances = distances_[featI];

    // Collect all those points that have a current maxLevel less than
    // (any of) the shell. Also collect the furthest distance allowable
    // to any shell with a higher level.

    pointField candidates(pt.size());
    labelList candidateMap(pt.size());
    scalarField candidateDistSqr(pt.size());
    label candidateI = 0;

    forAll(maxLevel, pointI)
    {
        forAllReverse(levels, levelI)
        {
            if (levels[levelI] > maxLevel[pointI])
            {
                candidates[candidateI] = pt[pointI];
                candidateMap[candidateI] = pointI;
                candidateDistSqr[candidateI] = sqr(distances[levelI]);
                candidateI++;
                break;
            }
        }
    }
    candidates.setSize(candidateI);
    candidateMap.setSize(candidateI);
    candidateDistSqr.setSize(candidateI);

    // Do the expensive nearest test only for the candidate points.
    const indexedOctree<treeDataEdge>& tree = edgeTrees_[featI];

    List<pointIndexHit> nearInfo(candidates.size());
    forAll(candidates, candidateI)
    {
        nearInfo[candidateI] = tree.findNearest
        (
            candidates[candidateI],
            candidateDistSqr[candidateI]
        );
    }

    // Update maxLevel
    forAll(nearInfo, candidateI)
    {
        if (nearInfo[candidateI].hit())
        {
            // Check which level it actually is in.
            label minDistI = findLower
            (
                distances,
                mag(nearInfo[candidateI].hitPoint()-candidates[candidateI])
            );

            label pointI = candidateMap[candidateI];

            // pt is inbetween shell[minDistI] and shell[minDistI+1]
            maxLevel[pointI] = levels[minDistI+1];
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinementFeatures::refinementFeatures
(
    const objectRegistry& io,
    const PtrList<dictionary>& featDicts
)
:
    PtrList<featureEdgeMesh>(featDicts.size()),
    distances_(featDicts.size()),
    levels_(featDicts.size()),
    edgeTrees_(featDicts.size()),
    pointTrees_(featDicts.size())
{
    // Read features
    read(io, featDicts);

    // Search engines
    forAll(*this, i)
    {
        const featureEdgeMesh& eMesh = operator[](i);
        const labelListList& pointEdges = eMesh.pointEdges();

        DynamicList<label> featurePoints;
        forAll(pointEdges, pointI)
        {
            if (pointEdges[pointI].size() > 2)
            {
                featurePoints.append(pointI);
            }
        }

        Info<< "Detected " << featurePoints.size()
            << " featurePoints out of " << pointEdges.size()
            << " on feature " << eMesh.name() << endl;

        buildTrees(i, featurePoints);
    }
}


Foam::refinementFeatures::refinementFeatures
(
    const objectRegistry& io,
    const PtrList<dictionary>& featDicts,
    const scalar minCos
)
:
    PtrList<featureEdgeMesh>(featDicts.size()),
    distances_(featDicts.size()),
    levels_(featDicts.size()),
    edgeTrees_(featDicts.size()),
    pointTrees_(featDicts.size())
{
    // Read features
    read(io, featDicts);

    // Search engines
    forAll(*this, i)
    {
        const featureEdgeMesh& eMesh = operator[](i);
        const pointField& points = eMesh.points();
        const edgeList& edges = eMesh.edges();
        const labelListList& pointEdges = eMesh.pointEdges();

        DynamicList<label> featurePoints;
        forAll(pointEdges, pointI)
        {
            const labelList& pEdges = pointEdges[pointI];
            if (pEdges.size() > 2)
            {
                featurePoints.append(pointI);
            }
            else if (pEdges.size() == 2)
            {
                // Check the angle
                const edge& e0 = edges[pEdges[0]];
                const edge& e1 = edges[pEdges[1]];

                const point& p = points[pointI];
                const point& p0 = points[e0.otherVertex(pointI)];
                const point& p1 = points[e1.otherVertex(pointI)];

                vector v0 = p-p0;
                scalar v0Mag = mag(v0);

                vector v1 = p1-p;
                scalar v1Mag = mag(v1);

                if
                (
                    v0Mag > SMALL
                 && v1Mag > SMALL
                 && ((v0/v0Mag & v1/v1Mag) < minCos)
                )
                {
                    featurePoints.append(pointI);
                }
            }
        }

        Info<< "Detected " << featurePoints.size()
            << " featurePoints out of " << points.size()
            << " on feature " << eMesh.name()
            << " when using feature cos " << minCos << endl;

        buildTrees(i, featurePoints);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::refinementFeatures::findNearestEdge
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    labelList& nearFeature,
    List<pointIndexHit>& nearInfo
) const
{
    nearFeature.setSize(samples.size());
    nearFeature = -1;
    nearInfo.setSize(samples.size());

    forAll(edgeTrees_, featI)
    {
        const indexedOctree<treeDataEdge>& tree = edgeTrees_[featI];

        if (tree.shapes().size() > 0)
        {
            forAll(samples, sampleI)
            {
                const point& sample = samples[sampleI];

                scalar distSqr;
                if (nearInfo[sampleI].hit())
                {
                    distSqr = magSqr(nearInfo[sampleI].hitPoint()-sample);
                }
                else
                {
                    distSqr = nearestDistSqr[sampleI];
                }

                pointIndexHit info = tree.findNearest(sample, distSqr);

                if (info.hit())
                {
                    nearInfo[sampleI] = info;
                    nearFeature[sampleI] = featI;
                }
            }
        }
    }
}


void Foam::refinementFeatures::findNearestPoint
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    labelList& nearFeature,
    labelList& nearIndex
) const
{
    nearFeature.setSize(samples.size());
    nearFeature = -1;
    nearIndex.setSize(samples.size());
    nearIndex = -1;

    forAll(pointTrees_, featI)
    {
        const indexedOctree<treeDataPoint>& tree = pointTrees_[featI];

        if (tree.shapes().pointLabels().size() > 0)
        {
            forAll(samples, sampleI)
            {
                const point& sample = samples[sampleI];

                scalar distSqr;
                if (nearFeature[sampleI] != -1)
                {
                    label nearFeatI = nearFeature[sampleI];
                    const indexedOctree<treeDataPoint>& nearTree =
                        pointTrees_[nearFeatI];
                    label featPointI =
                        nearTree.shapes().pointLabels()[nearIndex[sampleI]];
                    const point& featPt =
                        operator[](nearFeatI).points()[featPointI];
                    distSqr = magSqr(featPt-sample);
                }
                else
                {
                    distSqr = nearestDistSqr[sampleI];
                }

                pointIndexHit info = tree.findNearest(sample, distSqr);

                if (info.hit())
                {
                    nearFeature[sampleI] = featI;
                    nearIndex[sampleI] = info.index();
                }
            }
        }
    }
}


void Foam::refinementFeatures::findHigherLevel
(
    const pointField& pt,
    const labelList& ptLevel,
    labelList& maxLevel
) const
{
    // Maximum level of any shell. Start off with level of point.
    maxLevel = ptLevel;

    forAll(*this, featI)
    {
        findHigherLevel(pt, featI, maxLevel);
    }
}


Foam::scalar Foam::refinementFeatures::maxDistance() const
{
    scalar overallMax = -GREAT;
    forAll(distances_, featI)
    {
        overallMax = max(overallMax, max(distances_[featI]));
    }
    return overallMax;
}


// ************************************************************************* //
