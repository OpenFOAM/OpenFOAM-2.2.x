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
    surfaceSmooth

Description
    Example of a simple laplacian smoother.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "OFstream.H"
#include "boundBox.H"

#include "MeshedSurfaces.H"

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validOptions.clear();
    argList::validArgs.append("surfaceFile");
    argList::validArgs.append("underrelax factor (0..1)");
    argList::validArgs.append("iterations");
    argList::validArgs.append("output surfaceFile");
    argList args(argc, argv);

    const fileName surfFileName = args[1];
    const scalar relax = args.argRead<scalar>(2);
    const label  iters = args.argRead<label>(3);
    const fileName outFileName = args[4];

    if (relax <= 0 || relax > 1)
    {
        FatalErrorIn(args.executable()) << "Illegal relaxation factor "
            << relax << endl
            << "0: no change   1: move vertices to average of neighbours"
            << exit(FatalError);
    }

    Info<< "Relax:" << relax << nl
        << "Iters:" << iters << nl
        << "Reading surface from " << surfFileName << " ..." << endl;

    meshedSurface surf1(surfFileName);

    Info<< "Faces        : " << surf1.size() << nl
        << "Vertices     : " << surf1.nPoints() << nl
        << "Bounding Box : " << boundBox(surf1.localPoints()) << endl;

    pointField newPoints(surf1.localPoints());

    const labelListList& pointEdges = surf1.pointEdges();

    for (label iter = 0; iter < iters; iter++)
    {
        forAll(pointEdges, vertI)
        {
            vector avgPos(vector::zero);

            const labelList& myEdges = pointEdges[vertI];

            forAll(myEdges, myEdgeI)
            {
                const edge& e = surf1.edges()[myEdges[myEdgeI]];

                label otherVertI = e.otherVertex(vertI);

                avgPos += surf1.localPoints()[otherVertI];
            }
            avgPos /= myEdges.size();

            newPoints[vertI] = (1-relax)*newPoints[vertI] + relax*avgPos;
        }
    }

    Info<< "Writing surface to " << outFileName << " ..." << endl;

    meshedSurface
    (
        xferMove(newPoints),
        xferCopy(surf1.localFaces()),
        xferCopy(surf1.surfZones())
    ).write(outFileName);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
