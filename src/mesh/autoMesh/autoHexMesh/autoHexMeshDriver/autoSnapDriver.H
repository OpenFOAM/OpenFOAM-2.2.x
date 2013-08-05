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

Class
    Foam::autoSnapDriver

Description
    All to do with snapping to surface

SourceFiles
    autoSnapDriver.C
    autoSnapDriverFeature.C

\*---------------------------------------------------------------------------*/

#ifndef autoSnapDriver_H
#define autoSnapDriver_H

#include "meshRefinement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class motionSmoother;
class snapParameters;
class pointConstraint;

/*---------------------------------------------------------------------------*\
                           Class autoSnapDriver Declaration
\*---------------------------------------------------------------------------*/

class autoSnapDriver
{
    // Private data

        //- Mesh+surface
        meshRefinement& meshRefiner_;

        //- From global surface region to master side patch
        const labelList globalToMasterPatch_;

        //- From global surface region to slave side patch
        const labelList globalToSlavePatch_;


    // Private Member Functions


        // Snapping

            //- Calculates (geometric) shared points
            static label getCollocatedPoints
            (
                const scalar tol,
                const pointField&,
                PackedBoolList&
            );

            //- Calculate displacement per patch point to smooth out patch.
            //  Quite complicated in determining which points to move where.
            pointField smoothPatchDisplacement
            (
                const motionSmoother&,
                const List<labelPair>&
            ) const;

            //- Check that face zones are synced
            void checkCoupledFaceZones() const;

            //- Per edge distance to patch
            static tmp<scalarField> edgePatchDist
            (
                const pointMesh&,
                const indirectPrimitivePatch&
            );

            //- Write displacement as .obj file.
            static void dumpMove
            (
                const fileName&,
                const pointField&,
                const pointField&
            );

            //- Check displacement is outwards pointing
            static bool outwardsDisplacement
            (
                const indirectPrimitivePatch&,
                const vectorField&
            );

            // Feature line snapping

                //- Is point on two feature edges that make a largish angle?
                bool isFeaturePoint
                (
                    const scalar featureCos,
                    const indirectPrimitivePatch& pp,
                    const PackedBoolList& isFeatureEdge,
                    const label pointI
                ) const;

                void smoothAndConstrain
                (
                    const indirectPrimitivePatch& pp,
                    const List<pointConstraint>& constraints,
                    vectorField& disp
                ) const;
                void smoothAndConstrain2
                (
                    const bool applyConstraints,
                    const indirectPrimitivePatch& pp,
                    const List<pointConstraint>& constraints,
                    vectorField& disp
                ) const;
                void calcNearest
                (
                    const label iter,
                    const indirectPrimitivePatch& pp,
                    vectorField& pointDisp,
                    vectorField& pointSurfaceNormal,
                    vectorField& pointRotation
                ) const;
                void calcNearestFace
                (
                    const label iter,
                    const indirectPrimitivePatch& pp,
                    vectorField& faceDisp,
                    vectorField& faceSurfaceNormal,
                    labelList& faceSurfaceRegion,
                    vectorField& faceRotation
                ) const;
                void calcNearestFacePointProperties
                (
                    const label iter,
                    const indirectPrimitivePatch& pp,

                    const vectorField& faceDisp,
                    const vectorField& faceSurfaceNormal,
                    const labelList& faceSurfaceRegion,

                    List<List<point> >& pointFaceSurfNormals,
                    List<List<point> >& pointFaceDisp,
                    List<List<point> >& pointFaceCentres,
                    List<labelList>&    pointFacePatchID
                ) const;
                void correctAttraction
                (
                    const DynamicList<point>& surfacePoints,
                    const DynamicList<label>& surfaceCounts,
                    const point& edgePt,
                    const vector& edgeNormal,   // normalised normal
                    const point& pt,
                    vector& edgeOffset  // offset from pt to point on edge
                ) const;


                //- For any reverse (so from feature back to mesh) attraction:
                //  add attraction if diagonal points on face attracted
                void stringFeatureEdges
                (
                    const label iter,
                    const scalar featureCos,

                    const indirectPrimitivePatch& pp,
                    const scalarField& snapDist,

                    const vectorField& rawPatchAttraction,
                    const List<pointConstraint>& rawPatchConstraints,

                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints
                ) const;

                //- Return hit if on multiple points
                pointIndexHit findMultiPatchPoint
                (
                    const point& pt,
                    const labelList& patchIDs,
                    const List<point>& faceCentres
                ) const;
                void binFeatureFace
                (
                    const label iter,
                    const scalar featureCos,

                    const indirectPrimitivePatch& pp,
                    const scalar snapDist,

                    const point& fc,
                    const vector& faceSurfaceNormal,
                    const vector& faceDisp,

                    DynamicList<point>& surfacePoints,
                    DynamicList<vector>& surfaceNormals,
                    DynamicList<label>& surfaceCounts
                ) const;
                void binFeatureFaces
                (
                    const label iter,
                    const scalar featureCos,

                    const indirectPrimitivePatch& pp,
                    const scalarField& snapDist,
                    const label pointI,

                    const List<List<point> >& pointFaceSurfNormals,
                    const List<List<point> >& pointFaceDisp,
                    const List<List<point> >& pointFaceCentres,

                    DynamicList<point>& surfacePoints,
                    DynamicList<vector>& surfaceNormals,
                    DynamicList<label>& surfaceCounts
                ) const;


                void featureAttractionUsingReconstruction
                (
                    const label iter,
                    const scalar featureCos,

                    const indirectPrimitivePatch& pp,
                    const scalarField& snapDist,
                    const label pointI,

                    const List<List<point> >& pointFaceSurfNormals,
                    const List<List<point> >& pointFaceDisp,
                    const List<List<point> >& pointFaceCentres,
                    const labelListList& pointFacePatchID,

                    vector& patchAttraction,
                    pointConstraint& patchConstraint
                ) const;

                void featureAttractionUsingReconstruction
                (
                    const label iter,
                    const scalar featureCos,
                    const indirectPrimitivePatch& pp,
                    const scalarField& snapDist,

                    const List<List<point> >& pointFaceSurfNormals,
                    const List<List<point> >& pointFaceDisp,
                    const List<List<point> >& pointFaceCentres,
                    const labelListList& pointFacePatchID,

                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints
                ) const;

                void determineFeatures
                (
                    const label iter,
                    const scalar featureCos,
                    const bool multiRegionFeatureSnap,

                    const indirectPrimitivePatch&,
                    const scalarField&,

                    const List<List<point> >& pointFaceSurfNormals,
                    const List<List<point> >& pointFaceDisp,
                    const List<List<point> >& pointFaceCentres,
                    const labelListList& pointFacePatchID,

                    List<labelList>& pointAttractor,
                    List<List<pointConstraint> >& pointConstraints,
                    // Feature-edge to pp point
                    List<List<DynamicList<point> > >& edgeAttractors,
                    List<List<DynamicList<pointConstraint> > >& edgeConstraints,
                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints
                ) const;

                //- Find point on nearest feature edge (within searchDist).
                //  Return point and feature
                //  and store feature-edge to mesh-point and vice versa
                pointIndexHit findNearFeatureEdge
                (
                    const indirectPrimitivePatch& pp,
                    const scalarField& snapDist,
                    const label pointI,
                    const point& estimatedPt,

                    label& featI,
                    List<List<DynamicList<point> > >&,
                    List<List<DynamicList<pointConstraint> > >&,
                    vectorField&,
                    List<pointConstraint>&
                ) const;

                //- Find nearest feature point (within searchDist).
                //  Return feature point
                //  and store feature-point to mesh-point and vice versa.
                //  If another mesh point already referring to this feature
                //  point and further away, reset that one to a near feature
                //  edge (using findNearFeatureEdge above)
                labelPair findNearFeaturePoint
                (
                    const indirectPrimitivePatch& pp,
                    const scalarField& snapDist,
                    const label pointI,
                    const point& estimatedPt,

                    // Feature-point to pp point
                    List<labelList>& pointAttractor,
                    List<List<pointConstraint> >& pointConstraints,
                    // Feature-edge to pp point
                    List<List<DynamicList<point> > >& edgeAttractors,
                    List<List<DynamicList<pointConstraint> > >& edgeConstraints,
                    // pp point to nearest feature
                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints
                ) const;

                void featureAttractionUsingFeatureEdges
                (
                    const label iter,
                    const scalar featureCos,
                    const bool multiRegionFeatureSnap,
                    const indirectPrimitivePatch& pp,
                    const scalarField& snapDist,

                    const List<List<point> >& pointFaceSurfNormals,
                    const List<List<point> >& pointFaceDisp,
                    const List<List<point> >& pointFaceCentres,
                    const labelListList& pointFacePatchID,

                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints
                ) const;
                void preventFaceSqueeze
                (
                    const label iter,
                    const scalar featureCos,
                    const indirectPrimitivePatch& pp,
                    const scalarField& snapDist,

                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints
                ) const;

                vectorField calcNearestSurfaceFeature
                (
                    const snapParameters& snapParams,
                    const label iter,
                    const scalar featureCos,
                    const scalar featureAttract,
                    const scalarField& snapDist,
                    const vectorField& nearestDisp,
                    motionSmoother& meshMover
                ) const;


        //- Disallow default bitwise copy construct
        autoSnapDriver(const autoSnapDriver&);

        //- Disallow default bitwise assignment
        void operator=(const autoSnapDriver&);


public:

    //- Runtime type information
    ClassName("autoSnapDriver");


    // Constructors

        //- Construct from components
        autoSnapDriver
        (
            meshRefinement& meshRefiner,
            const labelList& globalToMasterPatch,
            const labelList& globalToSlavePatch
        );


    // Member Functions

        // Snapping

            //- Merge baffles.
            autoPtr<mapPolyMesh> mergeZoneBaffles(const List<labelPair>&);

            //- Calculate edge length per patch point.
            scalarField calcSnapDistance
            (
                const snapParameters& snapParams,
                const indirectPrimitivePatch&
            ) const;

            //- Smooth the mesh (patch and internal) to increase visibility
            //  of surface points (on castellated mesh) w.r.t. surface.
            void preSmoothPatch
            (
                const snapParameters& snapParams,
                const label nInitErrors,
                const List<labelPair>& baffles,
                motionSmoother&
            ) const;

            //- Get points both on patch and facezone.
            labelList getZoneSurfacePoints
            (
                const indirectPrimitivePatch&,
                const word& zoneName
            ) const;

            //- Per patch point calculate point on nearest surface. Set as
            //  boundary conditions of motionSmoother displacement field. Return
            //  displacement of patch points.
            vectorField calcNearestSurface
            (
                const scalarField& snapDist,
                motionSmoother& meshMover
            ) const;

            //- Smooth the displacement field to the internal.
            void smoothDisplacement
            (
                const snapParameters& snapParams,
                motionSmoother&
            ) const;

            //- Do the hard work: move the mesh according to displacement,
            //  locally relax the displacement. Return true if ended up with
            //  correct mesh, false if not.
            bool scaleMesh
            (
                const snapParameters& snapParams,
                const label nInitErrors,
                const List<labelPair>& baffles,
                motionSmoother&
            );

            //- Repatch faces according to surface nearest the face centre
            autoPtr<mapPolyMesh> repatchToSurface
            (
                const snapParameters& snapParams,
                const labelList& adaptPatchIDs,
                const labelList& preserveFaces
            );

            void doSnap
            (
                const dictionary& snapDict,
                const dictionary& motionDict,
                const scalar featureCos,
                const snapParameters& snapParams
            );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
