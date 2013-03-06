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

Class
    Foam::polyMeshFilter

Description
    Filter the edges and faces of a polyMesh whilst satisfying the given mesh
    quality criteria.

    Works on a copy of the mesh.

SourceFiles
    polyMeshFilter.C

\*---------------------------------------------------------------------------*/

#ifndef polyMeshFilter_H
#define polyMeshFilter_H

#include "IOdictionary.H"
#include "Time.H"
#include "List.H"
#include "autoPtr.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;
class fvMesh;
class PackedBoolList;

/*---------------------------------------------------------------------------*\
                       Class polyMeshFilter Declaration
\*---------------------------------------------------------------------------*/

class polyMeshFilter
{
    // Private data

        //- Reference to the original mesh
        const fvMesh& mesh_;

        //- Copy of the original mesh to perform the filtering on
        autoPtr<fvMesh> newMeshPtr_;

        //- Dictionary containing the coefficient sub-dictionaries
        const IOdictionary dict_;

        //- After collapsing, check the mesh quality and redo the collapsing
        //  iteration if there are too many bad faces in the mesh
        Switch controlMeshQuality_;

        //- Coefficients for collapsing edges
        const dictionary& collapseEdgesCoeffDict_;

        //- Coefficients for collapsing faces
        dictionary collapseFacesCoeffDict_;

        //- Coefficients for controlling the mesh quality
        dictionary meshQualityCoeffDict_;

        //- Remove edges shorter than this length
        const scalar minLen_;

        //- Merge points that are only attached to two edges and have an angle
        //  between the edge greater than this value
        const scalar maxCos_;

        //- The amount that the local minimum edge length will be reduced by if
        //  the edge is part of a collapse string that generates poor quality
        //  faces
        const scalar edgeReductionFactor_;

        //- Maximum number of outer iterations
        const label maxIterations_;

        //- Maximum number of smoothing iterations of minEdgeLen_ and
        //  faceFilterFactor_
        const label maxSmoothIters_;

        //- Initialisation value of faceFilterFactor_
        const scalar initialFaceLengthFactor_;

        //- The amount that the local face size factor will be reduced by if
        //  the face is part of a collapse string that generates poor quality
        //  faces
        const scalar faceReductionFactor_;

        //- Maximum number of times a deleted point can be associated with the
        //  creation of a bad face it is forced to be kept.
        const label maxPointErrorCount_;

        //- The minimum edge length for each edge
        scalarField minEdgeLen_;

        //- The face filter factor for each face
        scalarField faceFilterFactor_;


    // Private Member Functions

        //- Increment pointErrorCount for points attached to a bad face
        void updatePointErrorCount
        (
            const PackedBoolList& isErrorPoint,
            const labelList& oldToNewMesh,
            labelList& pointErrorCount
        ) const;


        //- Given the new points that are part of bad faces, and a map from the
        //  old mesh points to the new mesh points, relax minEdgeLen_
        void checkMeshEdgesAndRelaxEdges
        (
            const polyMesh& newMesh,
            const labelList& oldToNewMesh,
            const PackedBoolList& isErrorPoint,
            const labelList& pointErrorCount
        );

        //- Given the new points that are part of bad faces, and a map from the
        //  old mesh points to the new mesh points, relax faceFilterFactor_
        void checkMeshFacesAndRelaxEdges
        (
            const polyMesh& newMesh,
            const labelList& oldToNewMesh,
            const PackedBoolList& isErrorPoint,
            const labelList& pointErrorCount
        );

        // Mark boundary points
        // boundaryPoint:
        // + -1 : point not on boundary
        // +  0 : point on a real boundary
        // + >0 : point on a processor patch with that ID
        // @todo Need to mark boundaryEdges as well, as an edge may have two
        //       boundary points but not itself lie on a boundary
        labelList findBoundaryPoints(const polyMesh& mesh) const;

        //- Print min/mean/max data for a field
        void printScalarFieldStats
        (
            const string desc,
            const scalarField& fld
        ) const;

        //- Update minEdgeLen_ for the new mesh based upon the movement of the
        //  old points to the new points
        void mapOldMeshEdgeFieldToNewMesh
        (
            const polyMesh& newMesh,
            const labelList& pointMap,
            scalarField& newMeshMinEdgeLen
        ) const;

        //- Update faceFilterFactor_ for the new mesh based upon the movement
        //  of the old faces to the new faces
        void mapOldMeshFaceFieldToNewMesh
        (
            const polyMesh& newMesh,
            const labelList& faceMap,
            scalarField& newMeshFaceFilterFactor
        ) const;

        //- Maintain a map of the original mesh points to the latest version of
        //  the filtered mesh.
        void updateOldToNewPointMap
        (
            const labelList& currToNew,
            labelList& origToCurrentPointMap
        ) const;

        //- Disallow default bitwise copy construct
        polyMeshFilter(const polyMeshFilter&);

        //- Disallow default bitwise assignment
        void operator=(const polyMeshFilter&);


public:

    //- Runtime type information
    ClassName("polyMeshFilter");


    // Constructors

        //- Construct from fvMesh
        explicit polyMeshFilter(const fvMesh& mesh);


    //- Destructor
    ~polyMeshFilter();


    // Member Functions

        // Access

            //- Return reference to the filtered mesh. Does not check if the
            //  mesh has actually been filtered.
            const autoPtr<fvMesh>& filteredMesh() const;


        // Edit

            //- Return a copy of an fvMesh
            static autoPtr<fvMesh> copyMesh(const fvMesh& mesh);

            //- Filter edges and faces
            label filter(const label nOriginalBadFaces);

            //- Filter edges only.
            label filterEdges(const label nOriginalBadFaces);

            //- Filter all faces that are in the face zone indirectPatchFaces
            label filterIndirectPatchFaces();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
