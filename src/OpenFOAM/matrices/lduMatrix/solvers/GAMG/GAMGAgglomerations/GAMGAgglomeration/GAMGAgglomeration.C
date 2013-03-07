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

\*---------------------------------------------------------------------------*/

#include "GAMGAgglomeration.H"
#include "lduMesh.H"
#include "lduMatrix.H"
#include "Time.H"
#include "dlLibraryTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GAMGAgglomeration, 0);
    defineRunTimeSelectionTable(GAMGAgglomeration, lduMesh);
    defineRunTimeSelectionTable(GAMGAgglomeration, lduMatrix);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::GAMGAgglomeration::compactLevels(const label nCreatedLevels)
{
    nCells_.setSize(nCreatedLevels);
    restrictAddressing_.setSize(nCreatedLevels);
    meshLevels_.setSize(nCreatedLevels);
    interfaceLevels_.setSize(nCreatedLevels + 1);
}


bool Foam::GAMGAgglomeration::continueAgglomerating
(
    const label nCoarseCells
) const
{
    // Check the need for further agglomeration on all processors
    bool contAgg = nCoarseCells >= nCellsInCoarsestLevel_;
    reduce(contAgg, andOp<bool>());
    return contAgg;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GAMGAgglomeration::GAMGAgglomeration
(
    const lduMesh& mesh,
    const dictionary& controlDict
)
:
    MeshObject<lduMesh, Foam::GeometricMeshObject, GAMGAgglomeration>(mesh),

    maxLevels_(50),

    nCellsInCoarsestLevel_
    (
        readLabel(controlDict.lookup("nCellsInCoarsestLevel"))
    ),

    nCells_(maxLevels_),
    restrictAddressing_(maxLevels_),
    faceRestrictAddressing_(maxLevels_),

    meshLevels_(maxLevels_),
    interfaceLevels_(maxLevels_ + 1)
{}


const Foam::GAMGAgglomeration& Foam::GAMGAgglomeration::New
(
    const lduMesh& mesh,
    const dictionary& controlDict
)
{
    if
    (
        !mesh.thisDb().foundObject<GAMGAgglomeration>
        (
            GAMGAgglomeration::typeName
        )
    )
    {
        const word agglomeratorType(controlDict.lookup("agglomerator"));

        const_cast<Time&>(mesh.thisDb().time()).libs().open
        (
            controlDict,
            "geometricGAMGAgglomerationLibs",
            lduMeshConstructorTablePtr_
        );

        lduMeshConstructorTable::iterator cstrIter =
            lduMeshConstructorTablePtr_->find(agglomeratorType);

        if (cstrIter == lduMeshConstructorTablePtr_->end())
        {
            FatalErrorIn
            (
                "GAMGAgglomeration::New"
                "(const lduMesh& mesh, const dictionary& controlDict)"
            )   << "Unknown GAMGAgglomeration type "
                << agglomeratorType << ".\n"
                << "Valid algebraic GAMGAgglomeration types are :"
                << lduMatrixConstructorTablePtr_->sortedToc() << endl
                << "Valid algebraic GAMGAgglomeration types are :"
                << lduMeshConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }

        return store(cstrIter()(mesh, controlDict).ptr());
    }
    else
    {
        return mesh.thisDb().lookupObject<GAMGAgglomeration>
        (
            GAMGAgglomeration::typeName
        );
    }
}


const Foam::GAMGAgglomeration& Foam::GAMGAgglomeration::New
(
    const lduMatrix& matrix,
    const dictionary& controlDict
)
{
    const lduMesh& mesh = matrix.mesh();

    if
    (
        !mesh.thisDb().foundObject<GAMGAgglomeration>
        (
            GAMGAgglomeration::typeName
        )
    )
    {
        const word agglomeratorType(controlDict.lookup("agglomerator"));

        const_cast<Time&>(mesh.thisDb().time()).libs().open
        (
            controlDict,
            "algebraicGAMGAgglomerationLibs",
            lduMatrixConstructorTablePtr_
        );

        if
        (
            !lduMatrixConstructorTablePtr_
         || !lduMatrixConstructorTablePtr_->found(agglomeratorType)
        )
        {
            return New(mesh, controlDict);
        }
        else
        {
            lduMatrixConstructorTable::iterator cstrIter =
                lduMatrixConstructorTablePtr_->find(agglomeratorType);

            return store(cstrIter()(matrix, controlDict).ptr());
        }
    }
    else
    {
        return mesh.thisDb().lookupObject<GAMGAgglomeration>
        (
            GAMGAgglomeration::typeName
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::GAMGAgglomeration::~GAMGAgglomeration()
{
    // Clear the interface storage by hand.
    // It is a list of ptrs not a PtrList for consistency of the interface
    for (label leveli=1; leveli<interfaceLevels_.size(); leveli++)
    {
        lduInterfacePtrsList& curLevel = interfaceLevels_[leveli];

        forAll(curLevel, i)
        {
            if (curLevel.set(i))
            {
                delete curLevel(i);
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::lduMesh& Foam::GAMGAgglomeration::meshLevel
(
    const label i
) const
{
    if (i == 0)
    {
        return mesh_;
    }
    else
    {
        return meshLevels_[i - 1];
    }
}


const Foam::lduInterfacePtrsList& Foam::GAMGAgglomeration::interfaceLevel
(
    const label i
) const
{
    return interfaceLevels_[i];
}


// ************************************************************************* //
