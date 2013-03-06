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
    Foam::refinementFeatures

Description
    Encapsulates queries for features.

SourceFiles
    refinementFeatures.C

\*---------------------------------------------------------------------------*/

#ifndef refinementFeatures_H
#define refinementFeatures_H

#include "featureEdgeMesh.H"
#include "indexedOctree.H"
#include "treeDataEdge.H"
#include "treeDataPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class refinementFeatures Declaration
\*---------------------------------------------------------------------------*/

class refinementFeatures
:
    public PtrList<featureEdgeMesh>
{
private:

    // Private data

        //- Per shell the list of ranges
        List<scalarField> distances_;

        //- Per shell per distance the refinement level
        labelListList levels_;

        //- Edge
        PtrList<indexedOctree<treeDataEdge> > edgeTrees_;

        //- Features points
        PtrList<indexedOctree<treeDataPoint> > pointTrees_;


    // Private Member Functions

        //- Read set of feature edge meshes
        void read(const objectRegistry&, const PtrList<dictionary>&);

        //- Build edge tree and feature point tree
        void buildTrees(const label, const labelList&);

        //- Find shell level higher than ptLevel
        void findHigherLevel
        (
            const pointField& pt,
            const label featI,
            labelList& maxLevel
        ) const;

public:

    // Constructors

        //- Construct from description
        refinementFeatures
        (
            const objectRegistry& io,
            const PtrList<dictionary>& featDicts
        );

        //- Construct from description and do geometric analysis to determine
        //  feature points
        refinementFeatures
        (
            const objectRegistry& io,
            const PtrList<dictionary>& featDicts,
            const scalar minCos
        );


    // Member Functions

        // Access

            //- Per featureEdgeMesh the list of level
            const labelListList& levels() const
            {
                return levels_;
            }

            //- Per featureEdgeMesh the list of ranges
            const List<scalarField>& distances() const
            {
                return distances_;
            }

            const PtrList<indexedOctree<treeDataEdge> >& edgeTrees() const
            {
                return edgeTrees_;
            }

            const PtrList<indexedOctree<treeDataPoint> >& pointTrees() const
            {
                return pointTrees_;
            }


        // Query

            //- Highest distance of all features
            scalar maxDistance() const;

            //- Find nearest point on nearest feature edge
            void findNearestEdge
            (
                const pointField& samples,
                const scalarField& nearestDistSqr,
                labelList& nearFeature,
                List<pointIndexHit>& nearInfo
            ) const;

            //- Find nearest feature point. Is an index into feature points
            //  which itself is an index into the edgeMesh points.
            //  So the point index is
            //      pointTrees()[nearFeature].shapes().pointLabels()[nearIndex]
            //  Wip.
            void findNearestPoint
            (
                const pointField& samples,
                const scalarField& nearestDistSqr,
                labelList& nearFeature,
                labelList& nearIndex
            ) const;

            //- Find shell level higher than ptLevel
            void findHigherLevel
            (
                const pointField& pt,
                const labelList& ptLevel,
                labelList& maxLevel
            ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
