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
    Foam::thermoPhysicsTypes

Description
    Type definitions for thermo-physics models

\*---------------------------------------------------------------------------*/

#ifndef thermoPhysicsTypes_H
#define thermoPhysicsTypes_H

#include "specie.H"
#include "perfectGas.H"
#include "incompressiblePerfectGas.H"
#include "hConstThermo.H"
#include "janafThermo.H"
#include "sensibleEnthalpy.H"
#include "thermo.H"
#include "sutherlandTransport.H"
#include "constTransport.H"

#include "icoPolynomial.H"
#include "hPolynomialThermo.H"
#include "polynomialTransport.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef
    constTransport
    <
        species::thermo
        <
            hConstThermo
            <
                perfectGas<specie>
            >,
            sensibleEnthalpy
        >
    > constGasThermoPhysics;

    typedef
    sutherlandTransport
    <
        species::thermo
        <
            janafThermo
            <
                perfectGas<specie>
            >,
            sensibleEnthalpy
        >
    > gasThermoPhysics;

    typedef
    constTransport
    <
        species::thermo
        <
            hConstThermo
            <
                incompressiblePerfectGas<specie>
            >,
            sensibleEnthalpy
        >
    > constIncompressibleGasThermoPhysics;

    typedef
    sutherlandTransport
    <
        species::thermo
        <
            janafThermo
            <
                incompressiblePerfectGas<specie>
            >,
            sensibleEnthalpy
        >
    > incompressibleGasThermoPhysics;

    typedef
    polynomialTransport
    <
        species::thermo
        <
            hPolynomialThermo
            <
                icoPolynomial<specie, 8>,
                8
            >,
            sensibleEnthalpy
        >,
        8
    > icoPoly8ThermoPhysics;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
