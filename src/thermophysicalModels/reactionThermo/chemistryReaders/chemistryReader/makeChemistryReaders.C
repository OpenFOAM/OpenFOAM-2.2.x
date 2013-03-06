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

#include "makeReactionThermo.H"
#include "thermoPhysicsTypes.H"
#include "solidThermoPhysicsTypes.H"

#include "chemistryReader.H"
#include "foamChemistryReader.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeChemistryReader(constGasThermoPhysics);
makeChemistryReader(gasThermoPhysics);
makeChemistryReader(constIncompressibleGasThermoPhysics);
makeChemistryReader(incompressibleGasThermoPhysics);
makeChemistryReader(icoPoly8ThermoPhysics);
makeChemistryReader(hConstSolidThermoPhysics);
makeChemistryReader(hExponentialSolidThermoPhysics);

makeChemistryReaderType(foamChemistryReader, constGasThermoPhysics);
makeChemistryReaderType(foamChemistryReader, gasThermoPhysics);
makeChemistryReaderType
(
    foamChemistryReader,
    constIncompressibleGasThermoPhysics
);
makeChemistryReaderType(foamChemistryReader, incompressibleGasThermoPhysics);
makeChemistryReaderType(foamChemistryReader, icoPoly8ThermoPhysics);
makeChemistryReaderType(foamChemistryReader, hConstSolidThermoPhysics);
makeChemistryReaderType(foamChemistryReader, hExponentialSolidThermoPhysics);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
