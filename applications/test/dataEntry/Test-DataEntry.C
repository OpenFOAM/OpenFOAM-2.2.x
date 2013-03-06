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

Application
    testDataEntry

Description
    Tests DataEntry

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "DataEntry.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    IOdictionary dataEntryProperties
    (
        IOobject
        (
            "dataEntryProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    autoPtr<DataEntry<scalar> > dataEntry
    (
        DataEntry<scalar>::New
        (
            "dataEntry",
            dataEntryProperties
        )
    );

    scalar x0 = readScalar(dataEntryProperties.lookup("x0"));
    scalar x1 = readScalar(dataEntryProperties.lookup("x1"));

    Info<< "Data entry type: " << dataEntry().type() << nl << endl;

    Info<< "Inputs" << nl
        << "    x0 = " << x0 << nl
        << "    x1 = " << x1 << nl
        << endl;

    Info<< "Interpolation" << nl
        << "    f(x0) = " << dataEntry().value(x0) << nl
        << "    f(x1) = " << dataEntry().value(x1) << nl
        << endl;

    Info<< "Integration" << nl
        << "    int(f(x)) lim(x0->x1) = " << dataEntry().integrate(x0, x1) << nl
        << endl;

    return 0;
}


// ************************************************************************* //
