/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::treeDataTriSurface

Description
    Encapsulates data for (indexedOc)tree searches on triSurface.

SourceFiles
    treeDataTriSurface.C

\*---------------------------------------------------------------------------*/

#ifndef treeDataTriSurface_H
#define treeDataTriSurface_H

#include "triSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class treeBoundBox;
class treeDataTriSurface;
template<class Type> class indexedOctree;

/*---------------------------------------------------------------------------*\
                           Class treeDataTriSurface Declaration
\*---------------------------------------------------------------------------*/

class treeDataTriSurface
{
    // Private data

        //- Reference to triSurface
        const triSurface& surface_;

        //- Tolerance to use for intersection tests
        const scalar planarTol_;

    // Private Member Functions

        // //- fast triangle nearest point calculation. Returns point in E0, E1
        // //  coordinate system:  base + s*E0 + t*E1
        // static scalar nearestCoords
        // (
        //     const point& base,
        //     const point& E0,
        //     const point& E1,
        //     const scalar a,
        //     const scalar b,
        //     const scalar c,
        //     const point& P,
        //     scalar& s,
        //     scalar& t
        // );

public:

    // Declare name of the class and its debug switch
    ClassName("treeDataTriSurface");


    // Constructors

        //- Construct from triSurface and tolerance for intersection
        //  tests. Holds reference.
        treeDataTriSurface(const triSurface&, const scalar planarTol);


    // Member Functions

        // Access

            inline const triSurface& surface() const
            {
                return surface_;
            }

            inline label size() const
            {
                return surface_.size();
            }

            //- Get representative point cloud for all shapes inside
            //  (one point per shape)
            pointField shapePoints() const;


        // Search

            //- Get type (inside,outside,mixed,unknown) of point w.r.t. surface.
            //  Only makes sense for closed surfaces.
            label getVolumeType
            (
                const indexedOctree<treeDataTriSurface>&,
                const point&
            ) const;

            //- Does (bb of) shape at index overlap bb
            bool overlaps
            (
                const label index,
                const treeBoundBox& sampleBb
            ) const;

            //- Calculates nearest (to sample) point in shape.
            //  Returns actual point and distance (squared)
            void findNearest
            (
                const labelUList& indices,
                const point& sample,

                scalar& nearestDistSqr,
                label& nearestIndex,
                point& nearestPoint
            ) const;

            //- Calculates nearest (to line) point in shape.
            //  Returns point and distance (squared)
            void findNearest
            (
                const labelUList& indices,
                const linePointRef& ln,

                treeBoundBox& tightest,
                label& minIndex,
                point& linePoint,
                point& nearestPoint
            ) const;

            //- Calculate intersection of triangle with ray. Sets result
            //  accordingly
            bool intersects
            (
                const label index,
                const point& start,
                const point& end,
                point& result
            ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
