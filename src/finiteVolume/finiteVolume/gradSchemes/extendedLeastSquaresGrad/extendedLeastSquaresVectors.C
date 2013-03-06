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

\*---------------------------------------------------------------------------*/

#include "extendedLeastSquaresVectors.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(extendedLeastSquaresVectors, 0);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::extendedLeastSquaresVectors::extendedLeastSquaresVectors
(
    const fvMesh& mesh,
    const scalar minDet
)
:
    MeshObject<fvMesh, extendedLeastSquaresVectors>(mesh),
    minDet_(minDet),
    pVectorsPtr_(NULL),
    nVectorsPtr_(NULL),
    additionalCellsPtr_(NULL),
    additionalVectorsPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::extendedLeastSquaresVectors::~extendedLeastSquaresVectors()
{
    deleteDemandDrivenData(pVectorsPtr_);
    deleteDemandDrivenData(nVectorsPtr_);

    deleteDemandDrivenData(additionalCellsPtr_);
    deleteDemandDrivenData(additionalVectorsPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::extendedLeastSquaresVectors::makeLeastSquaresVectors() const
{
    if (debug)
    {
        Info<< "extendedLeastSquaresVectors::makeLeastSquaresVectors() :"
            << "Constructing least square gradient vectors"
            << endl;
    }

    const fvMesh& mesh = mesh_;

    pVectorsPtr_ = new surfaceVectorField
    (
        IOobject
        (
            "LeastSquaresP",
            mesh_.pointsInstance(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedVector("zero", dimless/dimLength, vector::zero)
    );
    surfaceVectorField& lsP = *pVectorsPtr_;

    nVectorsPtr_ = new surfaceVectorField
    (
        IOobject
        (
            "LeastSquaresN",
            mesh_.pointsInstance(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedVector("zero", dimless/dimLength, vector::zero)
    );
    surfaceVectorField& lsN = *nVectorsPtr_;

    // Set local references to mesh data
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();


    // Determine number of dimensions and (for 2D) missing dimension
    label nDims = 0;
    label twoD = -1;
    for (direction dir = 0; dir < vector::nComponents; dir++)
    {
        if (mesh.geometricD()[dir] == 1)
        {
            nDims++;
        }
        else
        {
            twoD = dir;
        }
    }

    if (nDims == 1)
    {
        FatalErrorIn
        (
            "extendedLeastSquaresVectors::makeLeastSquaresVectors() const"
        )   << "Found a mesh with only one geometric dimension : "
            << mesh.geometricD()
            << exit(FatalError);
    }
    else if (nDims == 2)
    {
        Info<< "extendedLeastSquares : detected " << nDims
            << " valid directions. Missing direction " << twoD << endl;
    }


    const volVectorField& C = mesh.C();

    // Set up temporary storage for the dd tensor (before inversion)
    symmTensorField dd(mesh_.nCells(), symmTensor::zero);

    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        vector d = C[nei] - C[own];

        const symmTensor wdd(1.0/magSqr(d)*sqr(d));

        dd[own] += wdd;
        dd[nei] += wdd;
    }

    // Visit the boundaries. Coupled boundaries are taken into account
    // in the construction of d vectors.
    surfaceVectorField::GeometricBoundaryField& blsP = lsP.boundaryField();

    forAll(blsP, patchi)
    {
        const fvPatch& p = blsP[patchi].patch();
        const labelUList& faceCells = p.faceCells();

        // Build the d-vectors
        vectorField pd(p.delta());

        forAll(pd, patchFaceI)
        {
            dd[faceCells[patchFaceI]] +=
                (1.0/magSqr(pd[patchFaceI]))*sqr(pd[patchFaceI]);
        }
    }


    // Check for missing dimensions
    // Add the missing eigenvector (such that it does not
    // affect the determinant)
    if (nDims == 2)
    {
        forAll(dd, cellI)
        {
            if (twoD == 0)
            {
                dd[cellI].xx() = 1;
            }
            else if (twoD == 1)
            {
                dd[cellI].yy() = 1;
            }
            else
            {
                dd[cellI].zz() = 1;
            }
        }
    }
    scalarField detdd(det(dd));

    Info<< "max(detdd) = " << max(detdd) << nl
        << "min(detdd) = " << min(detdd) << nl
        << "average(detdd) = " << average(detdd) << endl;

    label nAdaptedCells = 0;
    label nAddCells = 0;
    label maxNaddCells = 4*detdd.size();
    additionalCellsPtr_ = new List<labelPair>(maxNaddCells);
    List<labelPair>& additionalCells_ = *additionalCellsPtr_;

    forAll(detdd, i)
    {
        label count = 0;

        label oldNAddCells = nAddCells;

        while (++count < 100 && detdd[i] < minDet_)
        {
            if (nAddCells == maxNaddCells)
            {
                FatalErrorIn
                (
                    "extendedLeastSquaresVectors::"
                    "makeLeastSquaresVectors() const"
                )   << "nAddCells exceeds maxNaddCells ("
                    << maxNaddCells << ")"
                    << exit(FatalError);
            }

            labelList pointLabels = mesh.cells()[i].labels(mesh.faces());

            scalar maxDetddij = 0.0;

            label addCell = -1;

            forAll(pointLabels, j)
            {
                forAll(mesh.pointCells()[pointLabels[j]], k)
                {
                    label cellj = mesh.pointCells()[pointLabels[j]][k];

                    if (cellj != i)
                    {
                        vector dCij = (mesh.C()[cellj] - mesh.C()[i]);

                        symmTensor ddij =
                            dd[i] + (1.0/magSqr(dCij))*sqr(dCij);

                        // Add the missing eigenvector (such that it does not
                        // affect the determinant)
                        if (nDims == 2)
                        {
                            if (twoD == 0)
                            {
                                ddij.xx() = 1;
                            }
                            else if (twoD == 1)
                            {
                                ddij.yy() = 1;
                            }
                            else
                            {
                                ddij.zz() = 1;
                            }
                        }

                        scalar detddij = det(ddij);

                        if (detddij > maxDetddij)
                        {
                            addCell = cellj;
                            maxDetddij = detddij;
                        }
                    }
                }
            }

            if (addCell != -1)
            {
                additionalCells_[nAddCells][0] = i;
                additionalCells_[nAddCells++][1] = addCell;
                vector dCij = mesh.C()[addCell] - mesh.C()[i];
                dd[i] += (1.0/magSqr(dCij))*sqr(dCij);

                // Add the missing eigenvector (such that it does not
                // affect the determinant)
                if (nDims == 2)
                {
                    if (twoD == 0)
                    {
                        dd[i].xx() = 1;
                    }
                    else if (twoD == 1)
                    {
                        dd[i].yy() = 1;
                    }
                    else
                    {
                        dd[i].zz() = 1;
                    }
                }

                detdd[i] = det(dd[i]);
            }
        }

        if (oldNAddCells < nAddCells)
        {
            nAdaptedCells++;
        }
    }

    additionalCells_.setSize(nAddCells);

    reduce(nAddCells, sumOp<label>());
    reduce(nAdaptedCells, sumOp<label>());
    if (nAddCells)
    {
        Info<< "max(detdd) = " << max(detdd) << nl
            << "min(detdd) = " << min(detdd) << nl
            << "average(detdd) = " << average(detdd) << nl
            << "nAdapted/nCells = "
            << scalar(nAdaptedCells)/mesh.globalData().nTotalCells() << nl
            << "nAddCells/nCells = "
            << scalar(nAddCells)/mesh.globalData().nTotalCells()
            << endl;
    }

    // Invert the dd tensor
    const symmTensorField invDd(inv(dd));


    // Revisit all faces and calculate the lsP and lsN vectors
    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        vector d = C[nei] - C[own];

        lsP[facei] = (1.0/magSqr(d))*(invDd[own] & d);
        lsN[facei] = ((-1.0)/magSqr(d))*(invDd[nei] & d);
    }

    forAll(blsP, patchI)
    {
        fvsPatchVectorField& patchLsP = blsP[patchI];

        const fvPatch& p = patchLsP.patch();
        const labelUList& faceCells = p.faceCells();

        // Build the d-vectors
        vectorField pd(p.delta());

        forAll(p, patchFaceI)
        {
            patchLsP[patchFaceI] =
                (1.0/magSqr(pd[patchFaceI]))
               *(invDd[faceCells[patchFaceI]] & pd[patchFaceI]);
        }
    }


    additionalVectorsPtr_ = new vectorField(additionalCells_.size());
    vectorField& additionalVectors_ = *additionalVectorsPtr_;

    forAll(additionalCells_, i)
    {
        vector dCij =
            mesh.C()[additionalCells_[i][1]] - mesh.C()[additionalCells_[i][0]];

        additionalVectors_[i] =
            (1.0/magSqr(dCij))*(invDd[additionalCells_[i][0]] & dCij);
    }

    if (debug)
    {
        Info<< "extendedLeastSquaresVectors::makeLeastSquaresVectors() :"
            << "Finished constructing least square gradient vectors"
            << endl;
    }
}


const Foam::surfaceVectorField&
Foam::extendedLeastSquaresVectors::pVectors() const
{
    if (!pVectorsPtr_)
    {
        makeLeastSquaresVectors();
    }

    return *pVectorsPtr_;
}


const Foam::surfaceVectorField&
Foam::extendedLeastSquaresVectors::nVectors() const
{
    if (!nVectorsPtr_)
    {
        makeLeastSquaresVectors();
    }

    return *nVectorsPtr_;
}


const Foam::List<Foam::labelPair>&
Foam::extendedLeastSquaresVectors::additionalCells() const
{
    if (!additionalCellsPtr_)
    {
        makeLeastSquaresVectors();
    }

    return *additionalCellsPtr_;
}


const Foam::vectorField&
Foam::extendedLeastSquaresVectors::additionalVectors() const
{
    if (!additionalVectorsPtr_)
    {
        makeLeastSquaresVectors();
    }

    return *additionalVectorsPtr_;
}


bool Foam::extendedLeastSquaresVectors::movePoints()
{
    deleteDemandDrivenData(pVectorsPtr_);
    deleteDemandDrivenData(nVectorsPtr_);

    deleteDemandDrivenData(additionalCellsPtr_);
    deleteDemandDrivenData(additionalVectorsPtr_);

    return true;
}


// ************************************************************************* //
