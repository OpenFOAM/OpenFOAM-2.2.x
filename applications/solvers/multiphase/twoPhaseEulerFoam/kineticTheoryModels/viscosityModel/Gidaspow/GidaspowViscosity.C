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

#include "GidaspowViscosity.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
    defineTypeNameAndDebug(GidaspowViscosity, 0);
    addToRunTimeSelectionTable(viscosityModel, GidaspowViscosity, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::GidaspowViscosity::GidaspowViscosity
(
    const dictionary& dict
)
:
    viscosityModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::GidaspowViscosity::~GidaspowViscosity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::GidaspowViscosity::mu1
(
    const volScalarField& alpha1,
    const volScalarField& Theta,
    const volScalarField& g0,
    const dimensionedScalar& rho1,
    const dimensionedScalar& da,
    const dimensionedScalar& e
) const
{
    const scalar sqrtPi = sqrt(constant::mathematical::pi);

    return rho1*da*sqrt(Theta)*
    (
        (4.0/5.0)*sqr(alpha1)*g0*(1.0 + e)/sqrtPi
      + (1.0/15.0)*sqrtPi*g0*(1.0 + e)*sqr(alpha1)
      + (1.0/6.0)*sqrtPi*alpha1
      + (10.0/96.0)*sqrtPi/((1.0 + e)*g0)
    );
}


// ************************************************************************* //
