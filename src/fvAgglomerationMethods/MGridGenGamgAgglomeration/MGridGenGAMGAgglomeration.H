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
    Foam::MGridGenGAMGAgglomeration

Description
    Agglomerate using the MGridGen algorithm.

SourceFiles
    MGridGenGAMGAgglomeration.C
    MGridGenGAMGAgglomerate.C

\*---------------------------------------------------------------------------*/

#ifndef MGridGenGAMGAgglomeration_H
#define MGridGenGAMGAgglomeration_H

#include "fvMesh.H"
#include "GAMGAgglomeration.H"

extern "C"
{
#   include "mgridgen.h"
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class MGridGenGAMGAgglomeration Declaration
\*---------------------------------------------------------------------------*/

class MGridGenGAMGAgglomeration
:
    public GAMGAgglomeration
{
    // Private data

        const fvMesh& fvMesh_;


    // Private Member Functions

        //- Construct the CSR format addressing
        void makeCompactCellFaceAddressingAndFaceWeights
        (
            const lduAddressing& fineAddressing,
            List<idxtype>& cellCells,
            List<idxtype>& cellCellOffsets,
            const vectorField& Si,
            List<scalar>& faceWeights
        );

        //- Calculate and return agglomeration
        tmp<labelField> agglomerate
        (
            label& nCoarseCells,
            const label minSize,
            const label maxSize,
            const lduAddressing& fineAddressing,
            const scalarField& V,
            const vectorField& Sf,
            const scalarField& Sb
        );


        //- Disallow default bitwise copy construct
        MGridGenGAMGAgglomeration(const MGridGenGAMGAgglomeration&);

        //- Disallow default bitwise assignment
        void operator=(const MGridGenGAMGAgglomeration&);


public:

    //- Runtime type information
    TypeName("MGridGen");


    // Constructors

        //- Construct given mesh and controls
        MGridGenGAMGAgglomeration
        (
            const lduMesh& mesh,
            const dictionary& controlDict
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
