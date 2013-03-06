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

\*---------------------------------------------------------------------------*/

#include "pairGAMGAgglomeration.H"
#include "lduInterfacePtrsList.H"
#include "GAMGInterface.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::pairGAMGAgglomeration::combineLevels(const label curLevel)
{
    label prevLevel = curLevel - 1;

    // Set the previous level nCells to the current
    nCells_[prevLevel] = nCells_[curLevel];

    // Map the restrictAddressing from the coarser level into the previous
    // finer level

    const labelList& curResAddr = restrictAddressing_[curLevel];
    labelList& prevResAddr = restrictAddressing_[prevLevel];

    const labelList& curFaceResAddr = faceRestrictAddressing_[curLevel];
    labelList& prevFaceResAddr = faceRestrictAddressing_[prevLevel];

    forAll(prevFaceResAddr, i)
    {
        if (prevFaceResAddr[i] >= 0)
        {
            prevFaceResAddr[i] = curFaceResAddr[prevFaceResAddr[i]];
        }
        else
        {
            prevFaceResAddr[i] = -curResAddr[-prevFaceResAddr[i] - 1] - 1;
        }
    }

    // Delete the restrictAddressing for the coarser level
    faceRestrictAddressing_.set(curLevel, NULL);

    forAll(prevResAddr, i)
    {
        prevResAddr[i] = curResAddr[prevResAddr[i]];
    }

    // Delete the restrictAddressing for the coarser level
    restrictAddressing_.set(curLevel, NULL);


    // Delete the matrix addressing and coefficients from the previous level
    // and replace with the corresponding entried from the coarser level
    meshLevels_.set(prevLevel, meshLevels_.set(curLevel, NULL));

    // Same for the lduInterfaceFields taking care to delete the sub-entries
    // held on List<T*>
    const lduInterfacePtrsList& curInterLevel = interfaceLevels_[curLevel+1];
    lduInterfacePtrsList& prevInterLevel = interfaceLevels_[prevLevel+1];

    forAll(prevInterLevel, inti)
    {
        if (prevInterLevel.set(inti))
        {
            refCast<GAMGInterface>(const_cast<lduInterface&>
            (
                prevInterLevel[inti]
            )).combine(refCast<const GAMGInterface>(curInterLevel[inti]));

            delete curInterLevel(inti);
        }
    }

    interfaceLevels_.set(curLevel+1, NULL);
}


// ************************************************************************* //
