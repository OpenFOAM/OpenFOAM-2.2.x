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

#include "refinementParameters.H"
#include "unitConversion.H"
#include "polyMesh.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::refinementParameters::refinementParameters
(
    const dictionary& dict,
    const label dummy
)
:
    maxGlobalCells_(readLabel(dict.lookup("cellLimit"))),
    maxLocalCells_(readLabel(dict.lookup("procCellLimit"))),
    minRefineCells_(readLabel(dict.lookup("minimumRefine"))),
    curvature_(readScalar(dict.lookup("curvature"))),
    nBufferLayers_(readLabel(dict.lookup("nBufferLayers"))),
    keepPoints_(dict.lookup("keepPoints")),
    allowFreeStandingZoneFaces_(dict.lookup("allowFreeStandingZoneFaces")),
    maxLoadUnbalance_(dict.lookupOrDefault<scalar>("maxLoadUnbalance",0))
{}


Foam::refinementParameters::refinementParameters(const dictionary& dict)
:
    maxGlobalCells_(readLabel(dict.lookup("maxGlobalCells"))),
    maxLocalCells_(readLabel(dict.lookup("maxLocalCells"))),
    minRefineCells_(readLabel(dict.lookup("minRefinementCells"))),
    nBufferLayers_(readLabel(dict.lookup("nCellsBetweenLevels"))),
    keepPoints_(pointField(1, dict.lookup("locationInMesh"))),
    allowFreeStandingZoneFaces_(dict.lookup("allowFreeStandingZoneFaces")),
    maxLoadUnbalance_(dict.lookupOrDefault<scalar>("maxLoadUnbalance",0))
{
    scalar featAngle(readScalar(dict.lookup("resolveFeatureAngle")));

    if (featAngle < 0 || featAngle > 180)
    {
        curvature_ = -GREAT;
    }
    else
    {
        curvature_ = Foam::cos(degToRad(featAngle));
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::refinementParameters::findCells(const polyMesh& mesh)
const
{
    // Force calculation of tet-diag decomposition (for use in findCell)
    (void)mesh.tetBasePtIs();

    // Global calculation engine
    globalIndex globalCells(mesh.nCells());

    // Cell label per point
    labelList cellLabels(keepPoints_.size());

    forAll(keepPoints_, i)
    {
        const point& keepPoint = keepPoints_[i];

        label localCellI = mesh.findCell(keepPoint);

        label globalCellI = -1;

        if (localCellI != -1)
        {
            globalCellI = globalCells.toGlobal(localCellI);
        }

        reduce(globalCellI, maxOp<label>());

        if (globalCellI == -1)
        {
            FatalErrorIn
            (
                "refinementParameters::findCells(const polyMesh&) const"
            )   << "Point " << keepPoint
                << " is not inside the mesh or on a face or edge." << nl
                << "Bounding box of the mesh:" << mesh.bounds()
                << exit(FatalError);
        }


        label procI = globalCells.whichProcID(globalCellI);
        label procCellI = globalCells.toLocal(procI, globalCellI);

        Info<< "Found point " << keepPoint << " in cell " << procCellI
            << " on processor " << procI << endl;


        if (globalCells.isLocal(globalCellI))
        {
            cellLabels[i] = localCellI;
        }
        else
        {
            cellLabels[i] = -1;
        }
    }
    return cellLabels;
}


// ************************************************************************* //
