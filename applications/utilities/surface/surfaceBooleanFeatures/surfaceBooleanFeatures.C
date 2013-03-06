/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    surfaceBooleanFeatures

Description

    Generates the extendedFeatureEdgeMesh for the interface between a boolean
    operation on two surfaces.  Assumes that the orientation of the surfaces is
    correct:

    + if the operation is union or intersection, that both surface's normals
      (n) have the same orientation with respect to a point, i.e. surfaces and b
      are orientated the same with respect to point x:

    @verbatim
       _______
      |       |--> n
      |    ___|___             x
      |a  |   |   |--> n
      |___|___|  b|
          |       |
          |_______|

    @endverbatim

    + if the operation is a subtraction, the surfaces should be oppositely
    oriented with respect to a point, i.e. for (a - b), then b's orientation
    should be such that x is "inside", and a's orientation such that x is
    "outside"

    @verbatim
       _______
      |       |--> n
      |    ___|___             x
      |a  |   |   |
      |___|___|  b|
          |  n <--|
          |_______|

    @endverbatim

    When the operation is peformed - for union, all of the edges generates where
    one surfaces cuts another are all "internal" for union, and "external" for
    intersection, b - a and a - b.  This has been assumed, formal (dis)proof is
    invited.

\*---------------------------------------------------------------------------*/

#include "triSurface.H"
#include "argList.H"
#include "Time.H"
#include "featureEdgeMesh.H"
#include "extendedFeatureEdgeMesh.H"
#include "triSurfaceSearch.H"
#include "OFstream.H"
#include "booleanSurface.H"
#include "edgeIntersections.H"
#include "meshTools.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Keep on shuffling surface points until no more degenerate intersections.
// Moves both surfaces and updates set of edge cuts.
bool intersectSurfaces
(
    triSurface& surf1,
    edgeIntersections& edgeCuts1,
    triSurface& surf2,
    edgeIntersections& edgeCuts2
)
{
    bool hasMoved1 = false;
    bool hasMoved2 = false;

    for (label iter = 0; iter < 10; iter++)
    {
        Info<< "Determining intersections of surf1 edges with surf2"
            << " faces" << endl;

        // Determine surface1 edge intersections. Allow surface to be moved.

        // Number of iterations needed to resolve degenerates
        label nIters1 = 0;
        {
            triSurfaceSearch querySurf2(surf2);

            scalarField surf1PointTol
            (
                1e-3*edgeIntersections::minEdgeLength(surf1)
            );

            // Determine raw intersections
            edgeCuts1 = edgeIntersections
            (
                surf1,
                querySurf2,
                surf1PointTol
            );

            // Shuffle a bit to resolve degenerate edge-face hits
            {
                pointField points1(surf1.points());

                nIters1 =
                    edgeCuts1.removeDegenerates
                    (
                        5,              // max iterations
                        surf1,
                        querySurf2,
                        surf1PointTol,
                        points1         // work array
                    );

                if (nIters1 != 0)
                {
                    // Update geometric quantities
                    surf1.movePoints(points1);
                    hasMoved1 = true;
                }
            }
        }

        Info<< "Determining intersections of surf2 edges with surf1"
            << " faces" << endl;

        label nIters2 = 0;
        {
            triSurfaceSearch querySurf1(surf1);

            scalarField surf2PointTol
            (
                1e-3*edgeIntersections::minEdgeLength(surf2)
            );

            // Determine raw intersections
            edgeCuts2 = edgeIntersections
            (
                surf2,
                querySurf1,
                surf2PointTol
            );

            // Shuffle a bit to resolve degenerate edge-face hits
            {
                pointField points2(surf2.points());

                nIters2 =
                    edgeCuts2.removeDegenerates
                    (
                        5,              // max iterations
                        surf2,
                        querySurf1,
                        surf2PointTol,
                        points2         // work array
                    );

                if (nIters2 != 0)
                {
                    // Update geometric quantities
                    surf2.movePoints(points2);
                    hasMoved2 = true;
                }
            }
        }

        if (nIters1 == 0 && nIters2 == 0)
        {
            Info<< "** Resolved all intersections to be proper edge-face pierce"
                << endl;
            break;
        }
    }

    if (hasMoved1)
    {
        fileName newFile("surf1.obj");
        Info<< "Surface 1 has been moved. Writing to " << newFile
            << endl;
        surf1.write(newFile);
    }

    if (hasMoved2)
    {
        fileName newFile("surf2.obj");
        Info<< "Surface 2 has been moved. Writing to " << newFile
            << endl;
        surf2.write(newFile);
    }

    return hasMoved1 || hasMoved2;
}


int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append("action");
    argList::validArgs.append("surface file");
    argList::validArgs.append("surface file");

    argList::addBoolOption
    (
        "perturb",
        "Perturb surface points to escape degenerate intersections"
    );

    argList::addBoolOption
    (
        "invertedSpace",
        "do the surfaces have inverted space orientation, "
        "i.e. a point at infinity is considered inside. "
        "This is only sensible for union and intersection."
    );

    #   include "setRootCase.H"
    #   include "createTime.H"

    word action(args.args()[1]);

    HashTable<booleanSurface::booleanOpType> validActions;
    validActions.insert("intersection", booleanSurface::INTERSECTION);
    validActions.insert("union", booleanSurface::UNION);
    validActions.insert("difference", booleanSurface::DIFFERENCE);

    if (!validActions.found(action))
    {
        FatalErrorIn(args.executable())
            << "Unsupported action " << action << endl
            << "Supported actions:" << validActions.toc() << exit(FatalError);
    }

    fileName surf1Name(args[2]);
    Info<< "Reading surface " << surf1Name << endl;
    triSurface surf1(surf1Name);

    Info<< surf1Name << " statistics:" << endl;
    surf1.writeStats(Info);
    Info<< endl;

    fileName surf2Name(args[3]);
    Info<< "Reading surface " << surf2Name << endl;
    triSurface surf2(surf2Name);

    Info<< surf2Name << " statistics:" << endl;
    surf2.writeStats(Info);
    Info<< endl;

    edgeIntersections edge1Cuts;
    edgeIntersections edge2Cuts;

    bool invertedSpace = args.optionFound("invertedSpace");

    if (invertedSpace && validActions[action] == booleanSurface::DIFFERENCE)
    {
        FatalErrorIn(args.executable())
            << "Inverted space only makes sense for union or intersection."
            << exit(FatalError);
    }

    if (args.optionFound("perturb"))
    {
        intersectSurfaces
        (
            surf1,
            edge1Cuts,
            surf2,
            edge2Cuts
        );
    }
    else
    {
        triSurfaceSearch querySurf2(surf2);

        Info<< "Determining intersections of surf1 edges with surf2 faces"
            << endl;

        edge1Cuts = edgeIntersections
        (
            surf1,
            querySurf2,
            1e-3*edgeIntersections::minEdgeLength(surf1)
        );

        triSurfaceSearch querySurf1(surf1);

        Info<< "Determining intersections of surf2 edges with surf1 faces"
            << endl;

        edge2Cuts = edgeIntersections
        (
            surf2,
            querySurf1,
            1e-3*edgeIntersections::minEdgeLength(surf2)
        );
    }

    // Determine intersection edges
    surfaceIntersection inter(surf1, edge1Cuts, surf2, edge2Cuts);

    fileName sFeatFileName =
        surf1Name.lessExt().name()
      + "_"
      + surf2Name.lessExt().name()
      + "_"
      + action;

    label nFeatEds = inter.cutEdges().size();

    vectorField normals(2*nFeatEds, vector::zero);
    vectorField edgeDirections(nFeatEds, vector::zero);
    labelListList edgeNormals(nFeatEds, labelList(2, label(-1)));

    triSurfaceSearch querySurf1(surf1);
    triSurfaceSearch querySurf2(surf2);

    OFstream normalFile(sFeatFileName + "_normals.obj");

    scalar scale = 0.05*min
    (
        querySurf1.tree().bb().mag(),
        querySurf2.tree().bb().mag()
    );

    forAll(inter.cutEdges(), i)
    {
        const edge& fE(inter.cutEdges()[i]);

        point fEC = fE.centre(inter.cutPoints());

        pointIndexHit nearest1 = querySurf1.tree().findNearest(fEC, sqr(GREAT));
        pointIndexHit nearest2 = querySurf2.tree().findNearest(fEC, sqr(GREAT));

        normals[2*i] = surf1.faceNormals()[nearest1.index()];
        normals[2*i + 1] = surf2.faceNormals()[nearest2.index()];

        edgeNormals[i][0] = 2*i;
        edgeNormals[i][1] = 2*i + 1;

        edgeDirections[i] = fE.vec(inter.cutPoints());

        {
            meshTools::writeOBJ(normalFile, inter.cutPoints()[fE.start()]);
            meshTools::writeOBJ(normalFile, inter.cutPoints()[fE.end()]);

            normalFile<< "l " << (5*i) + 1 << " " << (5*i) + 2<< endl;

            meshTools::writeOBJ(normalFile, fEC);
            meshTools::writeOBJ(normalFile, fEC + scale*normals[2*i]);
            meshTools::writeOBJ(normalFile, fEC + scale*normals[2*i + 1]);

            normalFile<< "l " << (5*i) + 3 << " " << (5*i) + 4 << endl;
            normalFile<< "l " << (5*i) + 3 << " " << (5*i) + 5 << endl;
        }
    }

    label internalStart = -1;

    if (validActions[action] == booleanSurface::UNION)
    {
        if (!invertedSpace)
        {
            // All edges are internal
            internalStart = 0;
        }
        else
        {
            // All edges are external
            internalStart = nFeatEds;
        }
    }
    else if (validActions[action] == booleanSurface::INTERSECTION)
    {
        if (!invertedSpace)
        {
            // All edges are external
            internalStart = nFeatEds;
        }
        else
        {
            // All edges are internal
            internalStart = 0;
        }
    }
    else if (validActions[action] == booleanSurface::DIFFERENCE)
    {
        // All edges are external
        internalStart = nFeatEds;
    }
    else
    {
        FatalErrorIn(args.executable())
            << "Unsupported booleanSurface:booleanOpType and space "
            << action << " " << invertedSpace
            << abort(FatalError);
    }

    // There are no feature points supported by surfaceIntersection
    // Flat, open or multiple edges are assumed to be impossible
    // Region edges are not explicitly supported by surfaceIntersection

    extendedFeatureEdgeMesh feMesh
    (
        IOobject
        (
            sFeatFileName + ".extendedFeatureEdgeMesh",
            runTime.constant(),
            "featureEdgeMesh",
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        inter.cutPoints(),
        inter.cutEdges(),
        0,                  // concaveStart,
        0,                  // mixedStart,
        0,                  // nonFeatureStart,
        internalStart,      // internalStart,
        nFeatEds,           // flatStart,
        nFeatEds,           // openStart,
        nFeatEds,           // multipleStart,
        normals,
        edgeDirections,
        edgeNormals,
        labelListList(0),   // featurePointNormals,
        labelListList(0),   // featurePointEdges,
        labelList(0)        // regionEdges
    );

    feMesh.write();

    feMesh.writeObj(sFeatFileName);

    {
        // Write a featureEdgeMesh for backwards compatibility
        featureEdgeMesh bfeMesh
        (
            IOobject
            (
                sFeatFileName + ".eMesh",   // name
                runTime.constant(),                         // instance
                "triSurface",
                runTime,                                    // registry
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            feMesh.points(),
            feMesh.edges()
        );

        Info<< nl << "Writing featureEdgeMesh to "
            << bfeMesh.objectPath() << endl;

        bfeMesh.regIOobject::write();
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
