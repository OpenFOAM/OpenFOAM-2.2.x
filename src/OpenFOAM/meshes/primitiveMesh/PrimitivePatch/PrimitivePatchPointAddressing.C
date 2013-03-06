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

Description
     Point addressing on the patch: pointEdges and pointFaces.

\*---------------------------------------------------------------------------*/

#include "PrimitivePatch.H"
#include "SLList.H"
#include "ListOps.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
void
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
calcPointEdges() const
{
    if (debug)
    {
        Info<< "PrimitivePatch<Face, FaceList, PointField, PointType>::"
            << "calcPointEdges() : calculating pointEdges"
            << endl;
    }

    if (pointEdgesPtr_)
    {
        // it is considered an error to attempt to recalculate
        // if already allocated
        FatalErrorIn
        (
            "PrimitivePatch<Face, FaceList, PointField, PointType>::"
            "calcPointEdges()"
        )   << "pointEdges already calculated"
            << abort(FatalError);
    }

    pointEdgesPtr_ = new labelListList(meshPoints().size());

    labelListList& pe = *pointEdgesPtr_;

    invertManyToMany(pe.size(), edges(), pe);

    if (debug)
    {
        Info<< "PrimitivePatch<Face, FaceList, PointField, PointType>::"
            << "calcPointEdges() finished calculating pointEdges"
            << endl;
    }

    // Now order the edges of each point according to whether they share a
    // face

//    DynamicList<label> newEdgeList;

//    forAll(pe, pointI)
//    {
//        const labelList& pEdges = pe[pointI];

//        label edgeI = pEdges[0];

//        label prevFaceI = edgeFaces()[edgeI][0];

//        newEdgeList.clear();
//        newEdgeList.setCapacity(pEdges.size());

//        do
//        {
//            newEdgeList.append(edgeI);

//            // Cross edge to next face
//            const labelList& eFaces = edgeFaces()[edgeI];

//            if (eFaces.size() != 2)
//            {
//                break;
//            }

//            label faceI = eFaces[0];
//            if (faceI == prevFaceI)
//            {
//                faceI = eFaces[1];
//            }

//            // Cross face to next edge
//            const labelList& fEdges = faceEdges()[faceI];

//            forAll(fEdges, feI)
//            {
//                const label nextEdgeI = fEdges[feI];
//                const edge& nextEdge = edges()[nextEdgeI];

//                if
//                (
//                    nextEdgeI != edgeI
//                 && (nextEdge.start() == pointI || nextEdge.end() == pointI)
//                )
//                {
//                    edgeI = nextEdgeI;
//                    break;
//                }
//            }

//            prevFaceI = faceI;

//        } while (edgeI != pEdges[0]);

//        if (newEdgeList.size() == pEdges.size())
//        {
//            pe[pointI] = newEdgeList;
//        }
//    }

//    if (debug)
//    {
//        Info<< "PrimitivePatch<Face, FaceList, PointField, PointType>::"
//            << "calcPointEdges() finished ordering pointEdges"
//            << endl;
//    }
}


template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
void
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
calcPointFaces() const
{
    if (debug)
    {
        Info<< "PrimitivePatch<Face, FaceList, PointField, PointType>::"
            << "calcPointFaces() : calculating pointFaces"
            << endl;
    }

    if (pointFacesPtr_)
    {
        // it is considered an error to attempt to recalculate
        // if already allocated
        FatalErrorIn
        (
            "PrimitivePatch<Face, FaceList, PointField, PointType>::"
            "calcPointFaces()"
        )   << "pointFaces already calculated"
            << abort(FatalError);
    }

    const List<Face>& f = localFaces();

    // set up storage for pointFaces
    List<SLList<label> > pointFcs(meshPoints().size());

    forAll(f, faceI)
    {
        const Face& curPoints = f[faceI];

        forAll(curPoints, pointI)
        {
            pointFcs[curPoints[pointI]].append(faceI);
        }
    }

    // sort out the list
    pointFacesPtr_ = new labelListList(pointFcs.size());

    labelListList& pf = *pointFacesPtr_;

    forAll(pointFcs, pointI)
    {
        pf[pointI].setSize(pointFcs[pointI].size());

        label i = 0;
        forAllIter(SLList<label>, pointFcs[pointI], curFacesIter)
        {
            pf[pointI][i++] = curFacesIter();
        }
    }

    if (debug)
    {
        Info<< "PrimitivePatch<Face, FaceList, PointField, PointType>::"
            << "calcPointFaces() finished calculating pointFaces"
            << endl;
    }
}


// ************************************************************************* //
