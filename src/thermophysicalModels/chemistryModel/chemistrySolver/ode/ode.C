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

#include "ode.H"
#include "chemistryModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ChemistryModel>
Foam::ode<ChemistryModel>::ode
(
    const fvMesh& mesh
)
:
    chemistrySolver<ChemistryModel>(mesh),
    coeffsDict_(this->subDict("odeCoeffs")),
    solverName_(coeffsDict_.lookup("solver")),
    odeSolver_(ODESolver::New(solverName_, *this)),
    eps_(readScalar(coeffsDict_.lookup("eps")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ChemistryModel>
Foam::ode<ChemistryModel>::~ode()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ChemistryModel>
Foam::scalar Foam::ode<ChemistryModel>::solve
(
    scalarField& c,
    const scalar T,
    const scalar p,
    const scalar t0,
    const scalar dt
) const
{
    label nSpecie = this->nSpecie();
    scalarField c1(this->nEqns(), 0.0);

    // copy the concentration, T and P to the total solve-vector
    for (label i = 0; i < nSpecie; i++)
    {
        c1[i] = c[i];
    }
    c1[nSpecie] = T;
    c1[nSpecie+1] = p;

    scalar dtEst = dt;

    odeSolver_->solve
    (
        *this,
        t0,
        t0 + dt,
        c1,
        eps_,
        dtEst
    );

    forAll(c, i)
    {
        c[i] = max(0.0, c1[i]);
    }

    return dtEst;
}


// ************************************************************************* //
