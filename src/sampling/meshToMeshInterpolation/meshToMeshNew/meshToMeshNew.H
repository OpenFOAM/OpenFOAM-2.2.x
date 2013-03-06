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
    Foam::meshToMeshNew

Description
    Class to calculate the cell-addressing between two overlapping meshes

SourceFiles
    meshToMeshNew.C
    meshToMeshNewTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef meshToMeshNew_H
#define meshToMeshNew_H

#include "polyMesh.H"
#include "boundBox.H"
#include "mapDistribute.H"
#include "volFieldsFwd.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class meshToMeshNew Declaration
\*---------------------------------------------------------------------------*/

class meshToMeshNew
{
public:

    // Public data types

        //- Enumeration specifying required accuracy
        enum interpolationMethod
        {
            imMap,
            imCellVolumeWeight
        };

        static const NamedEnum<interpolationMethod, 2>
            interpolationMethodNames_;

private:

    // Private data

        //- Name of source mesh region
        const word srcRegionName_;

        //- Name of target mesh region
        const word tgtRegionName_;

        //- Source to target cell addressing
        labelListList srcToTgtCellAddr_;

        //- Target to source cell addressing
        labelListList tgtToSrcCellAddr_;

        //- Source to target cell interplation weights
        scalarListList srcToTgtCellWght_;

        //- Target to source cell interpolation weights
        scalarListList tgtToSrcCellWght_;

        //- Interpolation method
        interpolationMethod method_;

        //- Cell total volume in overlap region [m3]
        scalar V_;

        //- Index of processor that holds all of both sides. -1 in all other
        //  cases
        label singleMeshProc_;

        //- Source map pointer - parallel running only
        autoPtr<mapDistribute> srcMapPtr_;

        //- Target map pointer - parallel running only
        autoPtr<mapDistribute> tgtMapPtr_;

        //- Tolerance used in volume overlap calculations
        static scalar tolerance_;


    // Private Member Functions

        //- Helper function to add a constant offset to a list
        template<class Type>
        void add(UList<Type>& fld, const label offset) const;

        //- Write the connectivity (debugging)
        void writeConnectivity
        (
            const polyMesh& src,
            const polyMesh& tgt,
            const labelListList& srcToTargetAddr
        ) const;

        //- Return src cell IDs for the overlap region
        labelList maskCells(const polyMesh& src, const polyMesh& tgt) const;

        //- Find indices of overlapping cells in src and tgt meshes - returns
        //  true if found a matching pair
        bool findInitialSeeds
        (
            const polyMesh& src,
            const polyMesh& tgt,
            const labelList& srcCellIDs,
            const boolList& mapFlag,
            const label startSeedI,
            label& srcSeedI,
            label& tgtSeedI
        ) const;


        // Direct (one-to-one) mapping

            //- Append to list of src mesh seed indices
            void appendToDirectSeeds
            (
                const polyMesh& src,
                const polyMesh& tgt,
                boolList& mapFlag,
                labelList& srcTgtSeed,
                DynamicList<label>& srcSeeds,
                label& srcSeedI,
                label& tgtSeedI
            ) const;

            //- Main driver routine for direct mapping
            void calcDirect
            (
                const polyMesh& src,
                const polyMesh& tgt,
                const label srcSeedI,
                const label tgtSeedI
            );


        // Indirect (non-conformal) mapping

            //- Normalise the interpolation weights
            void normaliseWeights
            (
                const word& descriptor,
                const scalarField& cellVolumes,
                const labelListList& addr,
                scalarListList& wght
            ) const;

            //- Append target cell neihgbour cells to cellIDs list
            void appendNbrTgtCells
            (
                const label tgtCellI,
                const polyMesh& tgt,
                const DynamicList<label>& visitedTgtCells,
                DynamicList<label>& nbrTgtCellIDs
            ) const;

            //- Set the next cells in the advancing front algorithm
            void setNextCells
            (
                label& startSeedI,
                label& srcCellI,
                label& tgtCellI,
                const polyMesh& src,
                const polyMesh& tgt,
                const labelList& srcCellIDs,
                const boolList& mapFlag,
                const DynamicList<label>& visitedCells,
                labelList& seedCells
            ) const;

            //- Return the true if cells intersect
            bool intersect
            (
                const polyMesh& src,
                const polyMesh& tgt,
                const label srcCellI,
                const label tgtCellI
            ) const;

            //- Return the intersection volume between two cells
            scalar interVol
            (
                const polyMesh& src,
                const polyMesh& tgt,
                const label srcCellI,
                const label tgtCellI
            ) const;

            //- Main driver routine for indirect mapping
            void calcIndirect
            (
                const polyMesh& src,
                const polyMesh& tgt,
                const label srcSeedI,
                const label tgtSeedI,
                const labelList& srcCellIDs,
                boolList& mapFlag,
                label& startSeedI
            );


        //- Calculate the addressing between overalping regions of src and tgt
        //  meshes - main driver function
        void calcAddressing(const polyMesh& src, const polyMesh& tgt);


        // Parallel operations

            //- Determine whether the meshes are split across multiple pocessors
            label calcDistribution
            (
                const polyMesh& src,
                const polyMesh& tgt
            ) const;

            //- Determine which processor bounding-boxes overlap
            label calcOverlappingProcs
            (
                const List<boundBox>& procBb,
                const boundBox& bb,
                boolList& overlaps
            ) const;

            //- Calculate the mapping between processors
            autoPtr<mapDistribute> calcProcMap
            (
                const polyMesh& src,
                const polyMesh& tgt
            ) const;

            //- Distribute mesh info from 'my' processor to others
            void distributeCells
            (
                const mapDistribute& map,
                const polyMesh& tgtMesh,
                const globalIndex& globalI,
                List<pointField>& points,
                List<label>& nInternalFaces,
                List<faceList>& faces,
                List<labelList>& faceOwner,
                List<labelList>& faceNeighbour,
                List<labelList>& cellIDs,
                List<labelList>& nbrProcIDs,
                List<labelList>& procLocalFaceIDs
            ) const;

            //- Collect pieces of tgt mesh from other procssors and restructure
            void distributeAndMergeCells
            (
                const mapDistribute& map,
                const polyMesh& tgt,
                const globalIndex& globalI,
                pointField& tgtPoints,
                faceList& tgtFaces,
                labelList& tgtFaceOwners,
                labelList& tgtFaceNeighbours,
                labelList& tgtCellIDs
            ) const;


        //- Disallow default bitwise copy construct
        meshToMeshNew(const meshToMeshNew&);

        //- Disallow default bitwise assignment
        void operator=(const meshToMeshNew&);


public:

    //- Run-time type information
    TypeName("meshToMeshNew");


    //- Construct from source and target meshes
    meshToMeshNew
    (
        const polyMesh& src,
        const polyMesh& tgt,
        const interpolationMethod& method
    );


    //- Destructor
    virtual ~meshToMeshNew();


    // Member Functions

        // Access

            //- Return const access to the source to target cell addressing
            inline const labelListList& srcToTgtCellAddr() const;

            //- Return const access to the target to source cell addressing
            inline const labelListList& tgtToSrcCellAddr() const;

            //- Return const access to the source to target cell weights
            inline const scalarListList& srcToTgtCellWght() const;

            //- Return const access to the target to source cell weights
            inline const scalarListList& tgtToSrcCellWght() const;

            //- Return const access to the overlap volume
            inline scalar V() const;


        // Evaluation

            // Source-to-target field mapping

                //- Map field from src to tgt mesh with defined operation
                //  Values passed in via 'result' are used to initialise the
                //  return value
                template<class Type, class CombineOp>
                void mapSrcToTgt
                (
                    const UList<Type>& srcFld,
                    const CombineOp& cop,
                    List<Type>& result
                ) const;

                //- Return the src field mapped to the tgt mesh with a defined
                //  operation.  Initial values of the result are set to zero
                template<class Type, class CombineOp>
                tmp<Field<Type> > mapSrcToTgt
                (
                    const Field<Type>& srcFld,
                    const CombineOp& cop
                ) const;

                //- Convenience function to map a tmp field to the tgt mesh
                //  with a defined operation
                template<class Type, class CombineOp>
                tmp<Field<Type> > mapSrcToTgt
                (
                    const tmp<Field<Type> >& tsrcFld,
                    const CombineOp& cop
                ) const;

                //- Convenience function to map a field to the tgt mesh with a
                //  default operation (plusEqOp)
                template<class Type>
                tmp<Field<Type> > mapSrcToTgt
                (
                    const Field<Type>& srcFld
                ) const;

                //- Convenience function to map a tmp field to the tgt mesh
                //  with a default operation (plusEqOp)
                template<class Type>
                tmp<Field<Type> > mapSrcToTgt
                (
                    const tmp<Field<Type> >& tsrcFld
                ) const;


            // Target-to-source field mapping

                //- Map field from tgt to src mesh with defined operation
                //  Values passed in via 'result' are used to initialise the
                //  return value
                template<class Type, class CombineOp>
                void mapTgtToSrc
                (
                    const UList<Type>& tgtFld,
                    const CombineOp& cop,
                    List<Type>& result
                ) const;

                //- Return the tgt field mapped to the src mesh with a defined
                //  operation.  Initial values of the result are set to zero
                template<class Type, class CombineOp>
                tmp<Field<Type> > mapTgtToSrc
                (
                    const Field<Type>& tgtFld,
                    const CombineOp& cop
                ) const;

                //- Convenience function to map a tmp field to the src mesh
                //  with a defined operation
                template<class Type, class CombineOp>
                tmp<Field<Type> > mapTgtToSrc
                (
                    const tmp<Field<Type> >& ttgtFld,
                    const CombineOp& cop
                ) const;

                //- Convenience function to map a field to the src mesh with a
                //  default operation (plusEqOp)
                template<class Type>
                tmp<Field<Type> > mapTgtToSrc
                (
                    const Field<Type>& tgtFld
                ) const;

                //- Convenience function to map a tmp field to the src mesh
                //  with a default operation (plusEqOp)
                template<class Type>
                tmp<Field<Type> > mapTgtToSrc
                (
                    const tmp<Field<Type> >& ttgtFld
                ) const;


            // Volume field mapping

                //- Interpolare a field with a defined operation.  Values
                //  passed in via 'result' are used to initialise the return
                //  value.  Optionally interpolate patch values
                template<class Type, class CombineOp>
                void interpolate
                (
                    const GeometricField<Type, fvPatchField, volMesh>& field,
                    const CombineOp& cop,
                    GeometricField<Type, fvPatchField, volMesh>& result,
                    const bool interpPatches = true
                ) const;

                //- Interpolare a field with a defined operation.  The initial
                //  values of the result are set to zero.  Optionally
                //  interpolate patch values
                template<class Type, class CombineOp>
                tmp<GeometricField<Type, fvPatchField, volMesh> > interpolate
                (
                    const GeometricField<Type, fvPatchField, volMesh>& field,
                    const CombineOp& cop,
                    const bool interpPatches = true
                ) const;

                //- Interpolare a tmp field with a defined operation.  The
                //  initial values of the result are set to zero.  Optionally
                //  interpolate patch values
                template<class Type, class CombineOp>
                tmp<GeometricField<Type, fvPatchField, volMesh> > interpolate
                (
                    const tmp<GeometricField<Type, fvPatchField, volMesh> >&
                        tfield,
                    const CombineOp& cop,
                    const bool interpPatches = true
                ) const;

                //- Convenience function to map a field with a default
                //  operation (plusEqOp).  Optionally interpolate patch values
                template<class Type>
                tmp<GeometricField<Type, fvPatchField, volMesh> > interpolate
                (
                    const GeometricField<Type, fvPatchField, volMesh>& field,
                    const bool interpPatches = true
                ) const;

                //- Convenience function to map a tmp field with a default
                //  operation (plusEqOp).  Optionally interpolate patch values
                template<class Type>
                tmp<GeometricField<Type, fvPatchField, volMesh> > interpolate
                (
                    const tmp<GeometricField<Type, fvPatchField, volMesh> >&
                        tfield,
                    const bool interpPatches = true
                ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "meshToMeshNewI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "meshToMeshNewTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
