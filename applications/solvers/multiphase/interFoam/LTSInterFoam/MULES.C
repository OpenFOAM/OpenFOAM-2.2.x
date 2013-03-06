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

#include "MULES.H"
#include "upwind.H"
#include "uncorrectedSnGrad.H"
#include "gaussConvectionScheme.H"
#include "gaussLaplacianScheme.H"
#include "uncorrectedSnGrad.H"
#include "surfaceInterpolate.H"
#include "fvcSurfaceIntegrate.H"
#include "slicedSurfaceFields.H"
#include "syncTools.H"

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::MULES::explicitLTSSolve
(
    volScalarField& psi,
    const surfaceScalarField& phi,
    surfaceScalarField& phiPsi,
    const scalar psiMax,
    const scalar psiMin
)
{
    explicitLTSSolve
    (
        geometricOneField(),
        psi,
        phi,
        phiPsi,
        zeroField(), zeroField(),
        psiMax, psiMin
    );
}


void Foam::MULES::implicitSolve
(
    volScalarField& psi,
    const surfaceScalarField& phi,
    surfaceScalarField& phiPsi,
    const scalar psiMax,
    const scalar psiMin
)
{
    implicitSolve
    (
        geometricOneField(),
        psi,
        phi,
        phiPsi,
        zeroField(), zeroField(),
        psiMax, psiMin
    );
}


// ************************************************************************* //
