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
    Foam::layerParameters

Description
    Simple container to keep together layer specific information.

SourceFiles
    layerParameters.C

\*---------------------------------------------------------------------------*/

#ifndef layerParameters_H
#define layerParameters_H

#include "dictionary.H"
#include "scalarField.H"
#include "labelList.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class polyBoundaryMesh;
class refinementSurfaces;

/*---------------------------------------------------------------------------*\
                           Class layerParameters Declaration
\*---------------------------------------------------------------------------*/

class layerParameters
{
public:

    // Public data types

        //- Enumeration defining the layer specification:
        //  - first and total thickness specified
        //  - first and expansion ratio specified
        //  - final and total thickness specified
        //  - final and expansion ratio specified
        enum layerSpecification
        {
            ILLEGAL,
            FIRST_AND_TOTAL,
            FIRST_AND_EXPANSION,
            FINAL_AND_TOTAL,
            FINAL_AND_EXPANSION
        };


private:

    // Static data members

        //- Default angle for faces to be convcave
        static const scalar defaultConcaveAngle;


    // Private data

        // Per patch (not region!) information

            //- How many layers to add.
            labelList numLayers_;

            //- Are sizes relative to local cell size
            Switch relativeSizes_;

            //- How thickness is specified.
            layerSpecification layerSpec_;

                scalarField firstLayerThickness_;

                scalarField finalLayerThickness_;

                scalarField thickness_;

                scalarField expansionRatio_;

            //- Minimum total thickness
            scalarField minThickness_;


        scalar featureAngle_;

        scalar slipFeatureAngle_;

        scalar concaveAngle_;

        label nGrow_;

        label nSmoothSurfaceNormals_;

        label nSmoothNormals_;

        label nSmoothThickness_;

        scalar maxFaceThicknessRatio_;

        scalar layerTerminationCos_;

        scalar maxThicknessToMedialRatio_;

        scalar minMedianAxisAngleCos_;

        label nBufferCellsNoExtrude_;

        label nSnap_;

        label nLayerIter_;

        label nRelaxedIter_;

        Switch additionalReporting_;

    // Private Member Functions

        //- Calculate expansion ratio from overall size v.s. thickness of
        //  first layer.
        scalar layerExpansionRatio
        (
            const label n,
            const scalar totalOverFirst
        ) const;

        //- Disallow default bitwise copy construct
        layerParameters(const layerParameters&);

        //- Disallow default bitwise assignment
        void operator=(const layerParameters&);


public:

    // Constructors

        //- Construct from dictionary
        layerParameters(const dictionary& dict, const polyBoundaryMesh&);


    // Member Functions

        // Per patch information

            //- How many layers to add.
            //  -1 : no specification. Assume 0 layers but allow sliding
            //       to make layers
            //   0 : specified to have 0 layers. No sliding allowed.
            //  >0 : number of layers
            const labelList& numLayers() const
            {
                return numLayers_;
            }

            //- Are size parameters relative to inner cell size or
            //  absolute distances.
            bool relativeSizes() const
            {
                return relativeSizes_;
            }

            // Expansion factor for layer mesh
            const scalarField& expansionRatio() const
            {
                return expansionRatio_;
            }

            //- Wanted thickness of the layer furthest away
            //  from the wall (i.e. nearest the original mesh).
            //  If relativeSize() this number is relative to undistorted
            //  size of the cell outside layer.
            const scalarField& finalLayerThickness() const
            {
                return finalLayerThickness_;
            }

            //- Wanted thickness of the layer nearest to the wall.
            //  If relativeSize() this number is relative to undistorted
            //  size of the cell outside layer.
            const scalarField& firstLayerThickness() const
            {
                return firstLayerThickness_;
            }

            //- Wanted overall thickness of all layers.
            //  If relativeSize() this number is relative to undistorted
            //  size of the cell outside layer.
            const scalarField& thickness() const
            {
                return thickness_;
            }

            //- Minimum overall thickness of cell layer. If for any reason layer
            //  cannot be above minThickness do not add layer.
            //  If relativeSize() this number is relative to undistorted
            //  size of the cell outside layer.
            const scalarField& minThickness() const
            {
                return minThickness_;
            }


        scalar featureAngle() const
        {
            return featureAngle_;
        }

        //- At non-patched sides allow mesh to slip if extrusion
        //  direction makes angle larger than slipFeatureAngle.
        scalar slipFeatureAngle() const
        {
            return slipFeatureAngle_;
        }

        scalar concaveAngle() const
        {
            return concaveAngle_;
        }

        //- If points get not extruded do nGrow layers of connected faces
        //  that are not grown. Is used to not do layers at all close to
        //  features.
        label nGrow() const
        {
            return nGrow_;
        }

        //- Number of smoothing iterations of surface normals
        label nSmoothSurfaceNormals() const
        {
            return nSmoothSurfaceNormals_;
        }

        //- Number of smoothing iterations of interior mesh movement
        //  direction
        label nSmoothNormals() const
        {
            return nSmoothNormals_;
        }

        //- Stop layer growth on highly warped cells
        scalar maxFaceThicknessRatio() const
        {
            return maxFaceThicknessRatio_;
        }

        scalar layerTerminationCos() const
        {
            return layerTerminationCos_;
        }

        //- Smooth layer thickness over surface patches
        label nSmoothThickness() const
        {
            return nSmoothThickness_;
        }

        //- Reduce layer growth where ratio thickness to medial
        //  distance is large
        scalar maxThicknessToMedialRatio() const
        {
            return maxThicknessToMedialRatio_;
        }

        //- Angle used to pick up medial axis points
        scalar minMedianAxisAngleCos() const
        {
            return minMedianAxisAngleCos_;
        }

        //- Create buffer region for new layer terminations
        label nBufferCellsNoExtrude() const
        {
            return nBufferCellsNoExtrude_;
        }

        label nSnap() const
        {
            return nSnap_;
        }

        const Switch& additionalReporting() const
        {
            return additionalReporting_;
        }


        // Overall

            //- Number of overall layer addition iterations
            label nLayerIter() const
            {
                return nLayerIter_;
            }

            //- Number of iterations after which relaxed motion rules
            //  are to be used.
            label nRelaxedIter() const
            {
                return nRelaxedIter_;
            }


        // Helper

            //- Determine overall thickness. Uses two of the four parameters
            //  according to the layerSpecification
            scalar layerThickness
            (
                const label nLayers,
                const scalar firstLayerThickess,
                const scalar finalLayerThickess,
                const scalar totalThickness,
                const scalar expansionRatio
            ) const;

            //- Determine expansion ratio. Uses two of the four parameters
            //  according to the layerSpecification
            scalar layerExpansionRatio
            (
                const label nLayers,
                const scalar firstLayerThickess,
                const scalar finalLayerThickess,
                const scalar totalThickness,
                const scalar expansionRatio
            ) const;

            //- Determine first layer (near-wall) thickness. Uses two of the
            //  four parameters according to the layerSpecification
            scalar firstLayerThickness
            (
                const label nLayers,
                const scalar firstLayerThickess,
                const scalar finalLayerThickess,
                const scalar totalThickness,
                const scalar expansionRatio
            ) const;

            //- Determine ratio of final layer thickness to
            //  overall layer thickness
            scalar finalLayerThicknessRatio
            (
                const label nLayers,
                const scalar expansionRatio
            ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
