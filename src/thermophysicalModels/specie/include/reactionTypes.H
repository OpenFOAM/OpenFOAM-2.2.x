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

Typedefs
    Foam::reactionTypes

Description
    Type definitions for reactions

\*---------------------------------------------------------------------------*/

#ifndef reactionTypes_H
#define reactionTypes_H

#include "thermoPhysicsTypes.H"
#include "Reaction.H"

#include "icoPolynomial.H"
#include "hPolynomialThermo.H"
#include "polynomialTransport.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef Reaction<constGasThermoPhysics> constGasReaction;

    typedef Reaction<gasThermoPhysics> gasReaction;

    typedef Reaction<constIncompressibleGasThermoPhysics>
        constIncompressibleGasReaction;

    typedef Reaction<incompressibleGasThermoPhysics> incompressibleGasReaction;

    typedef Reaction<icoPoly8ThermoPhysics> icoPoly8Reaction;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
