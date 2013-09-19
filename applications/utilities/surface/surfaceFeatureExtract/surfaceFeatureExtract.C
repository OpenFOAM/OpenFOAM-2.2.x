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
    surfaceFeatureExtract

Description
    Extracts and writes surface features to file. All but the basic feature
    extraction is WIP.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "triSurface.H"
#include "surfaceFeatures.H"
#include "featureEdgeMesh.H"
#include "extendedFeatureEdgeMesh.H"
#include "treeBoundBox.H"
#include "meshTools.H"
#include "OFstream.H"
#include "triSurfaceMesh.H"
#include "vtkSurfaceWriter.H"
#include "triSurfaceFields.H"
#include "indexedOctree.H"
#include "treeDataEdge.H"
#include "unitConversion.H"
#include "plane.H"

#ifdef ENABLE_CURVATURE
#include "buildCGALPolyhedron.H"
#include "CGALPolyhedronRings.H"
#include <CGAL/Monge_via_jet_fitting.h>
#include <CGAL/Lapack/Linear_algebra_lapack.h>
#include <CGAL/property_map.h>
#endif

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef ENABLE_CURVATURE
scalarField calcCurvature(const triSurface& surf)
{
    scalarField k(surf.points().size(), 0);

    Polyhedron P;

    buildCGALPolyhedron convert(surf);
    P.delegate(convert);

    // Info<< "Created CGAL Polyhedron with " << label(P.size_of_vertices())
    //     << " vertices and " << label(P.size_of_facets())
    //     << " facets. " << endl;

    // The rest of this function adapted from
    //     CGAL-3.7/examples/Jet_fitting_3/Mesh_estimation.cpp

     //Vertex property map, with std::map
    typedef std::map<Vertex*, int> Vertex2int_map_type;
    typedef boost::associative_property_map< Vertex2int_map_type >
        Vertex_PM_type;
    typedef T_PolyhedralSurf_rings<Polyhedron, Vertex_PM_type > Poly_rings;

    typedef CGAL::Monge_via_jet_fitting<Kernel>         Monge_via_jet_fitting;
    typedef Monge_via_jet_fitting::Monge_form           Monge_form;

    std::vector<Point_3> in_points;  //container for data points

    // default parameter values and global variables
    unsigned int d_fitting = 2;
    unsigned int d_monge = 2;
    unsigned int min_nb_points = (d_fitting + 1)*(d_fitting + 2)/2;

    //initialize the tag of all vertices to -1
    Vertex_iterator vitb = P.vertices_begin();
    Vertex_iterator vite = P.vertices_end();

    Vertex2int_map_type vertex2props;
    Vertex_PM_type vpm(vertex2props);

    CGAL_For_all(vitb, vite)
    {
        put(vpm, &(*vitb), -1);
    }

    vite = P.vertices_end();

    label vertI = 0;

    for (vitb = P.vertices_begin(); vitb != vite; vitb++)
    {
        //initialize
        Vertex* v = &(*vitb);

        //gather points around the vertex using rings
        // From: gather_fitting_points(v, in_points, vpm);
        {
            std::vector<Vertex*> gathered;
            in_points.clear();

            Poly_rings::collect_enough_rings(v, min_nb_points, gathered, vpm);

            //store the gathered points
            std::vector<Vertex*>::iterator itb = gathered.begin();
            std::vector<Vertex*>::iterator ite = gathered.end();

            CGAL_For_all(itb, ite)
            {
                in_points.push_back((*itb)->point());
            }
        }

        //skip if the nb of points is to small
        if ( in_points.size() < min_nb_points )
        {
            std::cerr
                << "not enough pts for fitting this vertex"
                << in_points.size()
                << std::endl;

            continue;
        }

        // perform the fitting
        Monge_via_jet_fitting monge_fit;

        Monge_form monge_form = monge_fit
        (
            in_points.begin(),
            in_points.end(),
            d_fitting,
            d_monge
        );

//        std::cout<< monge_form;;
//        std::cout<< "condition number : "
//                 << monge_fit.condition_number() << nl << std::endl;

        // Use the maximum curvature to give smaller cell sizes later.
        k[vertI++] =
            max
            (
                mag(monge_form.principal_curvatures(0)),
                mag(monge_form.principal_curvatures(1))
            );
    }

    return k;
}
#endif


bool edgesConnected(const edge& e1, const edge& e2)
{
    if
    (
        e1.start() == e2.start()
     || e1.start() == e2.end()
     || e1.end() == e2.start()
     || e1.end() == e2.end()
    )
    {
        return true;
    }

    return false;
}


scalar calcProximityOfFeaturePoints
(
    const List<pointIndexHit>& hitList,
    const scalar defaultCellSize
)
{
    scalar minDist = defaultCellSize;

    for
    (
        label hI1 = 0;
        hI1 < hitList.size() - 1;
        ++hI1
    )
    {
        const pointIndexHit& pHit1 = hitList[hI1];

        if (pHit1.hit())
        {
            for
            (
                label hI2 = hI1 + 1;
                hI2 < hitList.size();
                ++hI2
            )
            {
                const pointIndexHit& pHit2 = hitList[hI2];

                if (pHit2.hit())
                {
                    scalar curDist = mag(pHit1.hitPoint() - pHit2.hitPoint());

                    minDist = min(curDist, minDist);
                }
            }
        }
    }

    return minDist;
}


scalar calcProximityOfFeatureEdges
(
    const extendedFeatureEdgeMesh& efem,
    const List<pointIndexHit>& hitList,
    const scalar defaultCellSize
)
{
    scalar minDist = defaultCellSize;

    for
    (
        label hI1 = 0;
        hI1 < hitList.size() - 1;
        ++hI1
    )
    {
        const pointIndexHit& pHit1 = hitList[hI1];

        if (pHit1.hit())
        {
            const edge& e1 = efem.edges()[pHit1.index()];

            for
            (
                label hI2 = hI1 + 1;
                hI2 < hitList.size();
                ++hI2
            )
            {
                const pointIndexHit& pHit2 = hitList[hI2];

                if (pHit2.hit())
                {
                    const edge& e2 = efem.edges()[pHit2.index()];

                    // Don't refine if the edges are connected to each other
                    if (!edgesConnected(e1, e2))
                    {
                        scalar curDist =
                            mag(pHit1.hitPoint() - pHit2.hitPoint());

                        minDist = min(curDist, minDist);
                    }
                }
            }
        }
    }

    return minDist;
}


void dumpBox(const treeBoundBox& bb, const fileName& fName)
{
    OFstream str(fName);

    Info<< "Dumping bounding box " << bb << " as lines to obj file "
        << str.name() << endl;


    pointField boxPoints(bb.points());

    forAll(boxPoints, i)
    {
        meshTools::writeOBJ(str, boxPoints[i]);
    }

    forAll(treeBoundBox::edges, i)
    {
        const edge& e = treeBoundBox::edges[i];

        str<< "l " << e[0]+1 <<  ' ' << e[1]+1 << nl;
    }
}


// Deletes all edges inside/outside bounding box from set.
void deleteBox
(
    const triSurface& surf,
    const treeBoundBox& bb,
    const bool removeInside,
    List<surfaceFeatures::edgeStatus>& edgeStat
)
{
    forAll(edgeStat, edgeI)
    {
        const point eMid = surf.edges()[edgeI].centre(surf.localPoints());

        if (removeInside ? bb.contains(eMid) : !bb.contains(eMid))
        {
            edgeStat[edgeI] = surfaceFeatures::NONE;
        }
    }
}


bool onLine(const point& p, const linePointRef& line)
{
    const point& a = line.start();
    const point& b = line.end();

    if
    (
        ( p.x() < min(a.x(), b.x()) || p.x() > max(a.x(), b.x()) )
     || ( p.y() < min(a.y(), b.y()) || p.y() > max(a.y(), b.y()) )
     || ( p.z() < min(a.z(), b.z()) || p.z() > max(a.z(), b.z()) )
    )
    {
        return false;
    }

    return true;
}


// Deletes all edges inside/outside bounding box from set.
void deleteEdges
(
    const triSurface& surf,
    const plane& cutPlane,
    List<surfaceFeatures::edgeStatus>& edgeStat
)
{
    const pointField& points = surf.points();
    const labelList& meshPoints = surf.meshPoints();

    forAll(edgeStat, edgeI)
    {
        const edge& e = surf.edges()[edgeI];
        const point& p0 = points[meshPoints[e.start()]];
        const point& p1 = points[meshPoints[e.end()]];
        const linePointRef line(p0, p1);

        // If edge does not intersect the plane, delete.
        scalar intersect = cutPlane.lineIntersect(line);

        point featPoint = intersect * (p1 - p0) + p0;

        if (!onLine(featPoint, line))
        {
            edgeStat[edgeI] = surfaceFeatures::NONE;
        }
    }
}


void drawHitProblem
(
    label fI,
    const triSurface& surf,
    const pointField& start,
    const pointField& faceCentres,
    const pointField& end,
    const List<pointIndexHit>& hitInfo
)
{
    Info<< nl << "# findLineAll did not hit its own face."
        << nl << "# fI " << fI
        << nl << "# start " << start[fI]
        << nl << "# f centre " << faceCentres[fI]
        << nl << "# end " << end[fI]
        << nl << "# hitInfo " << hitInfo
        << endl;

    meshTools::writeOBJ(Info, start[fI]);
    meshTools::writeOBJ(Info, faceCentres[fI]);
    meshTools::writeOBJ(Info, end[fI]);

    Info<< "l 1 2 3" << endl;

    meshTools::writeOBJ(Info, surf.points()[surf[fI][0]]);
    meshTools::writeOBJ(Info, surf.points()[surf[fI][1]]);
    meshTools::writeOBJ(Info, surf.points()[surf[fI][2]]);

    Info<< "f 4 5 6" << endl;

    forAll(hitInfo, hI)
    {
        label hFI = hitInfo[hI].index();

        meshTools::writeOBJ(Info, surf.points()[surf[hFI][0]]);
        meshTools::writeOBJ(Info, surf.points()[surf[hFI][1]]);
        meshTools::writeOBJ(Info, surf.points()[surf[hFI][2]]);

        Info<< "f "
            << 3*hI + 7 << " "
            << 3*hI + 8 << " "
            << 3*hI + 9
            << endl;
    }
}


// Unmark non-manifold edges if individual triangles are not features
void unmarkBaffles
(
    const triSurface& surf,
    const scalar includedAngle,
    List<surfaceFeatures::edgeStatus>& edgeStat
)
{
    scalar minCos = Foam::cos(degToRad(180.0 - includedAngle));

    const labelListList& edgeFaces = surf.edgeFaces();

    forAll(edgeFaces, edgeI)
    {
        const labelList& eFaces = edgeFaces[edgeI];

        if (eFaces.size() > 2)
        {
            label i0 = eFaces[0];
            //const labelledTri& f0 = surf[i0];
            const Foam::vector& n0 = surf.faceNormals()[i0];

            //Pout<< "edge:" << edgeI << " n0:" << n0 << endl;

            bool same = true;

            for (label i = 1; i < eFaces.size(); i++)
            {
                //const labelledTri& f = surf[i];
                const Foam::vector& n = surf.faceNormals()[eFaces[i]];

                //Pout<< "    mag(n&n0): " << mag(n&n0) << endl;

                if (mag(n&n0) < minCos)
                {
                    same = false;
                    break;
                }
            }

            if (same)
            {
                edgeStat[edgeI] = surfaceFeatures::NONE;
            }
        }
    }
}


//- Divide into multiple normal bins
//  - return REGION if != 2 normals
//  - return REGION if 2 normals that make feature angle
//  - otherwise return NONE and set normals,bins
surfaceFeatures::edgeStatus checkFlatRegionEdge
(
    const triSurface& surf,
    const scalar tol,
    const scalar includedAngle,
    const label edgeI
)
{
    const edge& e = surf.edges()[edgeI];
    const labelList& eFaces = surf.edgeFaces()[edgeI];

    // Bin according to normal

    DynamicList<Foam::vector> normals(2);
    DynamicList<labelList> bins(2);

    forAll(eFaces, eFaceI)
    {
        const Foam::vector& n = surf.faceNormals()[eFaces[eFaceI]];

        // Find the normal in normals
        label index = -1;
        forAll(normals, normalI)
        {
            if (mag(n&normals[normalI]) > (1-tol))
            {
                index = normalI;
                break;
            }
        }

        if (index != -1)
        {
            bins[index].append(eFaceI);
        }
        else if (normals.size() >= 2)
        {
            // Would be third normal. Mark as feature.
            //Pout<< "** at edge:" << surf.localPoints()[e[0]]
            //    << surf.localPoints()[e[1]]
            //    << " have normals:" << normals
            //    << " and " << n << endl;
            return surfaceFeatures::REGION;
        }
        else
        {
            normals.append(n);
            bins.append(labelList(1, eFaceI));
        }
    }


    // Check resulting number of bins
    if (bins.size() == 1)
    {
        // Note: should check here whether they are two sets of faces
        // that are planar or indeed 4 faces al coming together at an edge.
        //Pout<< "** at edge:"
        //    << surf.localPoints()[e[0]]
        //    << surf.localPoints()[e[1]]
        //    << " have single normal:" << normals[0]
        //    << endl;
        return surfaceFeatures::NONE;
    }
    else
    {
        // Two bins. Check if normals make an angle

        //Pout<< "** at edge:"
        //    << surf.localPoints()[e[0]]
        //    << surf.localPoints()[e[1]] << nl
        //    << "    normals:" << normals << nl
        //    << "    bins   :" << bins << nl
        //    << endl;

        if (includedAngle >= 0)
        {
            scalar minCos = Foam::cos(degToRad(180.0 - includedAngle));

            forAll(eFaces, i)
            {
                const Foam::vector& ni = surf.faceNormals()[eFaces[i]];
                for (label j=i+1; j<eFaces.size(); j++)
                {
                    const Foam::vector& nj = surf.faceNormals()[eFaces[j]];
                    if (mag(ni & nj) < minCos)
                    {
                        //Pout<< "have sharp feature between normal:" << ni
                        //    << " and " << nj << endl;

                        // Is feature. Keep as region or convert to
                        // feature angle? For now keep as region.
                        return surfaceFeatures::REGION;
                    }
                }
            }
        }


        // So now we have two normals bins but need to make sure both
        // bins have the same regions in it.

         // 1. store + or - region number depending
        //    on orientation of triangle in bins[0]
        const labelList& bin0 = bins[0];
        labelList regionAndNormal(bin0.size());
        forAll(bin0, i)
        {
            const labelledTri& t = surf.localFaces()[eFaces[bin0[i]]];
            int dir = t.edgeDirection(e);

            if (dir > 0)
            {
                regionAndNormal[i] = t.region()+1;
            }
            else if (dir == 0)
            {
                FatalErrorIn("problem.")
                    << exit(FatalError);
            }
            else
            {
                regionAndNormal[i] = -(t.region()+1);
            }
        }

        // 2. check against bin1
        const labelList& bin1 = bins[1];
        labelList regionAndNormal1(bin1.size());
        forAll(bin1, i)
        {
            const labelledTri& t = surf.localFaces()[eFaces[bin1[i]]];
            int dir = t.edgeDirection(e);

            label myRegionAndNormal;
            if (dir > 0)
            {
                myRegionAndNormal = t.region()+1;
            }
            else
            {
                myRegionAndNormal = -(t.region()+1);
            }

            regionAndNormal1[i] = myRegionAndNormal;

            label index = findIndex(regionAndNormal, -myRegionAndNormal);
            if (index == -1)
            {
                // Not found.
                //Pout<< "cannot find region " << myRegionAndNormal
                //    << " in regions " << regionAndNormal << endl;

                return surfaceFeatures::REGION;
            }
        }

        // Passed all checks, two normal bins with the same contents.
        //Pout<< "regionAndNormal:" << regionAndNormal << endl;
        //Pout<< "myRegionAndNormal:" << regionAndNormal1 << endl;

        return surfaceFeatures::NONE;
    }
}


void writeStats(const extendedFeatureEdgeMesh& fem, Ostream& os)
{
    os  << "    points : " << fem.points().size() << nl
        << "    of which" << nl
        << "        convex             : "
        << fem.concaveStart() << nl
        << "        concave            : "
        << (fem.mixedStart()-fem.concaveStart()) << nl
        << "        mixed              : "
        << (fem.nonFeatureStart()-fem.mixedStart()) << nl
        << "        non-feature        : "
        << (fem.points().size()-fem.nonFeatureStart()) << nl
        << "    edges  : " << fem.edges().size() << nl
        << "    of which" << nl
        << "        external edges     : "
        << fem.internalStart() << nl
        << "        internal edges     : "
        << (fem.flatStart()- fem.internalStart()) << nl
        << "        flat edges         : "
        << (fem.openStart()- fem.flatStart()) << nl
        << "        open edges         : "
        << (fem.multipleStart()- fem.openStart()) << nl
        << "        multiply connected : "
        << (fem.edges().size()- fem.multipleStart()) << nl;
}



int main(int argc, char *argv[])
{
    argList::addNote
    (
        "extract and write surface features to file"
    );
    argList::noParallel();

#   include "addDictOption.H"

#   include "setRootCase.H"
#   include "createTime.H"

    const word dictName("surfaceFeatureExtractDict");
#   include "setSystemRunTimeDictionaryIO.H"

    Info<< "Reading " << dictName << nl << endl;

    const IOdictionary dict(dictIO);

    forAllConstIter(dictionary, dict, iter)
    {
        if (!iter().isDict())
        {
            continue;
        }

        const dictionary& surfaceDict = iter().dict();

        if (!surfaceDict.found("extractionMethod"))
        {
            continue;
        }

        const word extractionMethod = surfaceDict.lookup("extractionMethod");

        const fileName surfFileName = iter().keyword();
        const fileName sFeatFileName = surfFileName.lessExt().name();

        Info<< "Surface            : " << surfFileName << nl << endl;

        const Switch writeVTK =
            surfaceDict.lookupOrDefault<Switch>("writeVTK", "off");
        const Switch writeObj =
            surfaceDict.lookupOrDefault<Switch>("writeObj", "off");

        const Switch curvature =
            surfaceDict.lookupOrDefault<Switch>("curvature", "off");
        const Switch featureProximity =
            surfaceDict.lookupOrDefault<Switch>("featureProximity", "off");
        const Switch closeness =
            surfaceDict.lookupOrDefault<Switch>("closeness", "off");


#ifndef ENABLE_CURVATURE
        if (curvature)
        {
            WarningIn(args.executable())
                << "Curvature calculation has been requested but "
                << args.executable() << " has not " << nl
                << "    been compiled with CGAL. "
                << "Skipping the curvature calculation." << endl;
        }
#else
        if (curvature && env("FOAM_SIGFPE"))
        {
            WarningIn(args.executable())
                << "Detected floating point exception trapping (FOAM_SIGFPE)."
                << " This might give" << nl
                << "    problems when calculating curvature on straight angles"
                << " (infinite curvature)" << nl
                << "    Switch it off in case of problems." << endl;
        }
#endif


        Info<< nl << "Feature line extraction is only valid on closed manifold "
            << "surfaces." << endl;

        // Read
        // ~~~~

        triSurface surf(runTime.constantPath()/"triSurface"/surfFileName);

        Info<< "Statistics:" << endl;
        surf.writeStats(Info);
        Info<< endl;

        faceList faces(surf.size());

        forAll(surf, fI)
        {
            faces[fI] = surf[fI].triFaceFace();
        }


        // Either construct features from surface & featureAngle or read set.
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        autoPtr<surfaceFeatures> set;

        scalar includedAngle = -1;

        if (extractionMethod == "extractFromFile")
        {
            const fileName featureEdgeFile =
                surfaceDict.subDict("extractFromFileCoeffs").lookup
                (
                    "featureEdgeFile"
                );

            edgeMesh eMesh(featureEdgeFile);

            // Sometimes duplicate edges are present. Remove them.
            eMesh.mergeEdges();

            Info<< nl << "Reading existing feature edges from file "
                << featureEdgeFile << endl;

            set.set(new surfaceFeatures(surf, eMesh.points(), eMesh.edges()));
        }
        else if (extractionMethod == "extractFromSurface")
        {
            includedAngle = readScalar
            (
                surfaceDict.subDict("extractFromSurfaceCoeffs").lookup
                (
                    "includedAngle"
                )
            );

            Info<< nl << "Constructing feature set from included angle "
                << includedAngle << endl;

            set.set(new surfaceFeatures(surf, includedAngle));
        }
        else
        {
            FatalErrorIn(args.executable())
                << "No initial feature set. Provide either one"
                << " of extractFromFile (to read existing set)" << nl
                << " or extractFromSurface (to construct new set from angle)"
                << exit(FatalError);
        }


        // Trim set
        // ~~~~~~~~

        if (surfaceDict.isDict("trimFeatures"))
        {
            dictionary trimDict = surfaceDict.subDict("trimFeatures");

            scalar minLen =
                trimDict.lookupOrAddDefault<scalar>("minLen", -GREAT);

            label minElem = trimDict.lookupOrAddDefault<label>("minElem", 0);

            // Trim away small groups of features
            if (minElem > 0 || minLen > 0)
            {
                Info<< "Removing features of length < "
                    << minLen << endl;
                Info<< "Removing features with number of edges < "
                    << minElem << endl;

                set().trimFeatures(minLen, minElem);
            }
        }


        // Subset
        // ~~~~~~

        // Convert to marked edges, points
        List<surfaceFeatures::edgeStatus> edgeStat(set().toStatus());

        if (surfaceDict.isDict("subsetFeatures"))
        {
            const dictionary& subsetDict = surfaceDict.subDict
            (
                "subsetFeatures"
            );

            if (subsetDict.found("insideBox"))
            {
                treeBoundBox bb(subsetDict.lookup("insideBox")());

                Info<< "Removing all edges outside bb " << bb << endl;
                dumpBox(bb, "subsetBox.obj");

                deleteBox(surf, bb, false, edgeStat);
            }
            else if (subsetDict.found("outsideBox"))
            {
                treeBoundBox bb(subsetDict.lookup("outsideBox")());

                Info<< "Removing all edges inside bb " << bb << endl;
                dumpBox(bb, "deleteBox.obj");

                deleteBox(surf, bb, true, edgeStat);
            }

            const Switch nonManifoldEdges =
                subsetDict.lookupOrDefault<Switch>("nonManifoldEdges", "yes");

            if (!nonManifoldEdges)
            {
                Info<< "Removing all non-manifold edges"
                    << " (edges with > 2 connected faces) unless they"
                    << " cross multiple regions" << endl;

                forAll(edgeStat, edgeI)
                {
                    const labelList& eFaces = surf.edgeFaces()[edgeI];

                    if
                    (
                        eFaces.size() > 2
                     && edgeStat[edgeI] == surfaceFeatures::REGION
                     && (eFaces.size() % 2) == 0
                    )
                    {
                        edgeStat[edgeI] = checkFlatRegionEdge
                        (
                            surf,
                            1e-5,   //tol,
                            includedAngle,
                            edgeI
                        );
                    }
                }
            }

            const Switch openEdges =
                subsetDict.lookupOrDefault<Switch>("openEdges", "yes");

            if (!openEdges)
            {
                Info<< "Removing all open edges"
                    << " (edges with 1 connected face)" << endl;

                forAll(edgeStat, edgeI)
                {
                    if (surf.edgeFaces()[edgeI].size() == 1)
                    {
                        edgeStat[edgeI] = surfaceFeatures::NONE;
                    }
                }
            }

            if (subsetDict.found("plane"))
            {
                plane cutPlane(subsetDict.lookup("plane")());

                deleteEdges(surf, cutPlane, edgeStat);

                Info<< "Only edges that intersect the plane with normal "
                    << cutPlane.normal()
                    << " and base point " << cutPlane.refPoint()
                    << " will be included as feature edges."<< endl;
            }
        }


        surfaceFeatures newSet(surf);
        newSet.setFromStatus(edgeStat);

        Info<< nl
            << "Initial feature set:" << nl
            << "    feature points : " << newSet.featurePoints().size() << nl
            << "    feature edges  : " << newSet.featureEdges().size() << nl
            << "    of which" << nl
            << "        region edges   : " << newSet.nRegionEdges() << nl
            << "        external edges : " << newSet.nExternalEdges() << nl
            << "        internal edges : " << newSet.nInternalEdges() << nl
            << endl;

        //if (writeObj)
        //{
        //    newSet.writeObj("final");
        //}

        // Extracting and writing a extendedFeatureEdgeMesh
        extendedFeatureEdgeMesh feMesh
        (
            newSet,
            runTime,
            sFeatFileName + ".extendedFeatureEdgeMesh"
        );


        if (surfaceDict.isDict("addFeatures"))
        {
            const word addFeName = surfaceDict.subDict("addFeatures")["name"];
            Info<< "Adding (without merging) features from " << addFeName
                << nl << endl;

            extendedFeatureEdgeMesh addFeMesh
            (
                IOobject
                (
                    addFeName,
                    runTime.time().constant(),
                    "extendedFeatureEdgeMesh",
                    runTime.time(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );
            Info<< "Read " << addFeMesh.name() << nl;
            writeStats(addFeMesh, Info);

            feMesh.add(addFeMesh);
        }


        Info<< nl
            << "Final feature set:" << nl;
        writeStats(feMesh, Info);

        Info<< nl << "Writing extendedFeatureEdgeMesh to "
            << feMesh.objectPath() << endl;

        mkDir(feMesh.path());

        if (writeObj)
        {
            feMesh.writeObj(feMesh.path()/surfFileName.lessExt().name());
        }

        feMesh.write();

        // Write a featureEdgeMesh for backwards compatibility
        featureEdgeMesh bfeMesh
        (
            IOobject
            (
                surfFileName.lessExt().name() + ".eMesh",   // name
                runTime.constant(),                         // instance
                "triSurface",
                runTime,                                    // registry
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            feMesh.points(),
            feMesh.edges()
        );

        Info<< nl << "Writing featureEdgeMesh to "
            << bfeMesh.objectPath() << endl;

        bfeMesh.regIOobject::write();

    // Find close features

    // // Dummy trim operation to mark features
    // labelList featureEdgeIndexing = newSet.trimFeatures(-GREAT, 0);

    // scalarField surfacePtFeatureIndex(surf.points().size(), -1);

    // forAll(newSet.featureEdges(), eI)
    // {
    //     const edge& e = surf.edges()[newSet.featureEdges()[eI]];

    //     surfacePtFeatureIndex[surf.meshPoints()[e.start()]] =
    //     featureEdgeIndexing[newSet.featureEdges()[eI]];

    //     surfacePtFeatureIndex[surf.meshPoints()[e.end()]] =
    //     featureEdgeIndexing[newSet.featureEdges()[eI]];
    // }

    // if (writeVTK)
    // {
    //     vtkSurfaceWriter().write
    //     (
    //         runTime.constant()/"triSurface",    // outputDir
    //         sFeatFileName,                      // surfaceName
    //         surf.points(),
    //         faces,
    //         "surfacePtFeatureIndex",            // fieldName
    //         surfacePtFeatureIndex,
    //         true,                               // isNodeValues
    //         true                                // verbose
    //     );
    // }

    // Random rndGen(343267);

    // treeBoundBox surfBB
    // (
    //     treeBoundBox(searchSurf.bounds()).extend(rndGen, 1e-4)
    // );

    // surfBB.min() -= Foam::point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
    // surfBB.max() += Foam::point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

    // indexedOctree<treeDataEdge> ftEdTree
    // (
    //     treeDataEdge
    //     (
    //         false,
    //         surf.edges(),
    //         surf.localPoints(),
    //         newSet.featureEdges()
    //     ),
    //     surfBB,
    //     8,      // maxLevel
    //     10,     // leafsize
    //     3.0     // duplicity
    // );

    // labelList nearPoints = ftEdTree.findBox
    // (
    //     treeBoundBox
    //     (
    //         sPt - featureSearchSpan*Foam::vector::one,
    //         sPt + featureSearchSpan*Foam::vector::one
    //     )
    // );

        if (closeness)
        {
            Info<< nl << "Extracting internal and external closeness of "
                << "surface." << endl;


            triSurfaceMesh searchSurf
            (
                IOobject
                (
                    sFeatFileName + ".closeness",
                    runTime.constant(),
                    "triSurface",
                    runTime,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                surf
            );


            // Internal and external closeness

            // Prepare start and end points for intersection tests

            const vectorField& normals = searchSurf.faceNormals();

            scalar span = searchSurf.bounds().mag();

            scalar externalAngleTolerance = 10;
            scalar externalToleranceCosAngle =
                Foam::cos
                (
                    degToRad(180 - externalAngleTolerance)
                );

            scalar internalAngleTolerance = 45;
            scalar internalToleranceCosAngle =
                Foam::cos
                (
                    degToRad(180 - internalAngleTolerance)
                );

            Info<< "externalToleranceCosAngle: " << externalToleranceCosAngle
                << nl
                << "internalToleranceCosAngle: " << internalToleranceCosAngle
                << endl;

            // Info<< "span " << span << endl;

            pointField start(searchSurf.faceCentres() - span*normals);
            pointField end(searchSurf.faceCentres() + span*normals);
            const pointField& faceCentres = searchSurf.faceCentres();

            List<List<pointIndexHit> > allHitInfo;

            // Find all intersections (in order)
            searchSurf.findLineAll(start, end, allHitInfo);

            scalarField internalCloseness(start.size(), GREAT);
            scalarField externalCloseness(start.size(), GREAT);

            forAll(allHitInfo, fI)
            {
                const List<pointIndexHit>& hitInfo = allHitInfo[fI];

                if (hitInfo.size() < 1)
                {
                    drawHitProblem(fI, surf, start, faceCentres, end, hitInfo);

                    // FatalErrorIn(args.executable())
                    //     << "findLineAll did not hit its own face."
                    //     << exit(FatalError);
                }
                else if (hitInfo.size() == 1)
                {
                    if (!hitInfo[0].hit())
                    {
                        // FatalErrorIn(args.executable())
                        //     << "findLineAll did not hit any face."
                        //     << exit(FatalError);
                    }
                    else if (hitInfo[0].index() != fI)
                    {
                        drawHitProblem
                        (
                            fI,
                            surf,
                            start,
                            faceCentres,
                            end,
                            hitInfo
                        );

                        // FatalErrorIn(args.executable())
                        //     << "findLineAll did not hit its own face."
                        //     << exit(FatalError);
                    }
                }
                else
                {
                    label ownHitI = -1;

                    forAll(hitInfo, hI)
                    {
                        // Find the hit on the triangle that launched the ray

                        if (hitInfo[hI].index() == fI)
                        {
                            ownHitI = hI;

                            break;
                        }
                    }

                    if (ownHitI < 0)
                    {
                        drawHitProblem
                        (
                            fI,
                            surf,
                            start,
                            faceCentres,
                            end,
                            hitInfo
                        );

                        // FatalErrorIn(args.executable())
                        //     << "findLineAll did not hit its own face."
                        //     << exit(FatalError);
                    }
                    else if (ownHitI == 0)
                    {
                        // There are no internal hits, the first hit is the
                        // closest external hit

                        if
                        (
                            (
                                normals[fI]
                              & normals[hitInfo[ownHitI + 1].index()]
                            )
                          < externalToleranceCosAngle
                        )
                        {
                            externalCloseness[fI] =
                                mag
                                (
                                    faceCentres[fI]
                                  - hitInfo[ownHitI + 1].hitPoint()
                                );
                        }
                    }
                    else if (ownHitI == hitInfo.size() - 1)
                    {
                        // There are no external hits, the last but one hit is
                        // the closest internal hit

                        if
                        (
                            (
                                normals[fI]
                              & normals[hitInfo[ownHitI - 1].index()]
                            )
                          < internalToleranceCosAngle
                        )
                        {
                            internalCloseness[fI] =
                                mag
                                (
                                    faceCentres[fI]
                                  - hitInfo[ownHitI - 1].hitPoint()
                                );
                        }
                    }
                    else
                    {
                        if
                        (
                            (
                                normals[fI]
                              & normals[hitInfo[ownHitI + 1].index()]
                            )
                          < externalToleranceCosAngle
                        )
                        {
                            externalCloseness[fI] =
                                mag
                                (
                                    faceCentres[fI]
                                  - hitInfo[ownHitI + 1].hitPoint()
                                );
                        }

                        if
                        (
                            (
                                normals[fI]
                              & normals[hitInfo[ownHitI - 1].index()]
                            )
                          < internalToleranceCosAngle
                        )
                        {
                            internalCloseness[fI] =
                                mag
                                (
                                    faceCentres[fI]
                                  - hitInfo[ownHitI - 1].hitPoint()
                                );
                        }
                    }
                }
            }

            triSurfaceScalarField internalClosenessField
            (
                IOobject
                (
                    sFeatFileName + ".internalCloseness",
                    runTime.constant(),
                    "triSurface",
                    runTime,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                surf,
                dimLength,
                internalCloseness
            );

            internalClosenessField.write();

            triSurfaceScalarField externalClosenessField
            (
                IOobject
                (
                    sFeatFileName + ".externalCloseness",
                    runTime.constant(),
                    "triSurface",
                    runTime,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                surf,
                dimLength,
                externalCloseness
            );

            externalClosenessField.write();

            if (writeVTK)
            {
                vtkSurfaceWriter().write
                (
                    runTime.constantPath()/"triSurface",// outputDir
                    sFeatFileName,                      // surfaceName
                    surf.points(),
                    faces,
                    "internalCloseness",                // fieldName
                    internalCloseness,
                    false,                              // isNodeValues
                    true                                // verbose
                );

                vtkSurfaceWriter().write
                (
                    runTime.constantPath()/"triSurface",// outputDir
                    sFeatFileName,                      // surfaceName
                    surf.points(),
                    faces,
                    "externalCloseness",                // fieldName
                    externalCloseness,
                    false,                              // isNodeValues
                    true                                // verbose
                );
            }
        }


#ifdef ENABLE_CURVATURE
        if (curvature)
        {
            Info<< nl << "Extracting curvature of surface at the points."
                << endl;

            scalarField k = calcCurvature(surf);

        // Modify the curvature values on feature edges and points to be zero.

        //    forAll(newSet.featureEdges(), fEI)
        //    {
        //        const edge& e = surf.edges()[newSet.featureEdges()[fEI]];
        //
        //        k[surf.meshPoints()[e.start()]] = 0.0;
        //        k[surf.meshPoints()[e.end()]] = 0.0;
        //    }

            triSurfacePointScalarField kField
            (
                IOobject
                (
                    sFeatFileName + ".curvature",
                    runTime.constant(),
                    "triSurface",
                    runTime,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                surf,
                dimLength,
                k
            );

            kField.write();

            if (writeVTK)
            {
                vtkSurfaceWriter().write
                (
                    runTime.constantPath()/"triSurface",// outputDir
                    sFeatFileName,                      // surfaceName
                    surf.points(),
                    faces,
                    "curvature",                        // fieldName
                    k,
                    true,                               // isNodeValues
                    true                                // verbose
                );
            }
        }
#endif


        if (featureProximity)
        {
            Info<< nl << "Extracting proximity of close feature points and "
                << "edges to the surface" << endl;

            const scalar searchDistance =
                readScalar(surfaceDict.lookup("maxFeatureProximity"));

            scalarField featureProximity(surf.size(), searchDistance);

            forAll(surf, fI)
            {
                const triPointRef& tri = surf[fI].tri(surf.points());
                const point& triCentre = tri.circumCentre();

                const scalar radiusSqr = min
                (
                    sqr(4*tri.circumRadius()),
                    sqr(searchDistance)
                );

                List<pointIndexHit> hitList;

                feMesh.allNearestFeatureEdges(triCentre, radiusSqr, hitList);

                featureProximity[fI] =
                    calcProximityOfFeatureEdges
                    (
                        feMesh,
                        hitList,
                        featureProximity[fI]
                    );

                feMesh.allNearestFeaturePoints(triCentre, radiusSqr, hitList);

                featureProximity[fI] =
                    calcProximityOfFeaturePoints
                    (
                        hitList,
                        featureProximity[fI]
                    );
            }

            triSurfaceScalarField featureProximityField
            (
                IOobject
                (
                    sFeatFileName + ".featureProximity",
                    runTime.constant(),
                    "triSurface",
                    runTime,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                surf,
                dimLength,
                featureProximity
            );

            featureProximityField.write();

            if (writeVTK)
            {
                vtkSurfaceWriter().write
                (
                    runTime.constantPath()/"triSurface",// outputDir
                    sFeatFileName,                      // surfaceName
                    surf.points(),
                    faces,
                    "featureProximity",                 // fieldName
                    featureProximity,
                    false,                              // isNodeValues
                    true                                // verbose
                );
            }
        }

        Info<< endl;
    }

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
