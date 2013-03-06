/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

Global
    execFlowFunctionObjects

Description
    Execute the set of functionObjects specified in the selected dictionary
    (which defaults to system/controlDict) for the selected set of times.
    Alternative dictionaries should be placed in the system/ folder.

    The flow (p-U) and optionally turbulence fields are available for the
    function objects to operate on allowing forces and other related properties
    to be calculated in addition to cutting planes etc.

\*---------------------------------------------------------------------------*/

#include "calc.H"

#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "ReadFields.H"

#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"

#include "incompressible/RAS/RASModel/RASModel.H"
#include "incompressible/LES/LESModel/LESModel.H"

#include "fluidThermo.H"
#include "compressible/RAS/RASModel/RASModel.H"
#include "compressible/LES/LESModel/LESModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    void execFlowFunctionObjects(const argList& args, const Time& runTime)
    {
        if (args.optionFound("dict"))
        {
            IOdictionary dict
            (
                IOobject
                (
                    args["dict"],
                    runTime,
                    IOobject::MUST_READ_IF_MODIFIED
                )
            );

            functionObjectList fol(runTime, dict);
            fol.start();
            fol.execute(true);  // override outputControl - force writing
        }
        else
        {
            functionObjectList fol(runTime);
            fol.start();
            fol.execute(true);  // override outputControl - force writing
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::calc(const argList& args, const Time& runTime, const fvMesh& mesh)
{
    if (args.optionFound("noFlow"))
    {
        Info<< "    Operating in no-flow mode; no models will be loaded."
            << " All vol, surface and point fields will be loaded." << endl;

        // Read objects in time directory
        IOobjectList objects(mesh, runTime.timeName());

        // Read vol fields.

        PtrList<volScalarField> vsFlds;
        ReadFields(mesh, objects, vsFlds);

        PtrList<volVectorField> vvFlds;
        ReadFields(mesh, objects, vvFlds);

        PtrList<volSphericalTensorField> vstFlds;
        ReadFields(mesh, objects, vstFlds);

        PtrList<volSymmTensorField> vsymtFlds;
        ReadFields(mesh, objects, vsymtFlds);

        PtrList<volTensorField> vtFlds;
        ReadFields(mesh, objects, vtFlds);

        // Read surface fields.

        PtrList<surfaceScalarField> ssFlds;
        ReadFields(mesh, objects, ssFlds);

        PtrList<surfaceVectorField> svFlds;
        ReadFields(mesh, objects, svFlds);

        PtrList<surfaceSphericalTensorField> sstFlds;
        ReadFields(mesh, objects, sstFlds);

        PtrList<surfaceSymmTensorField> ssymtFlds;
        ReadFields(mesh, objects, ssymtFlds);

        PtrList<surfaceTensorField> stFlds;
        ReadFields(mesh, objects, stFlds);

        // Read point fields.
        const pointMesh& pMesh = pointMesh::New(mesh);

        PtrList<pointScalarField> psFlds;
        ReadFields(pMesh, objects, psFlds);

        PtrList<pointVectorField> pvFlds;
        ReadFields(pMesh, objects, pvFlds);

        PtrList<pointSphericalTensorField> pstFlds;
        ReadFields(pMesh, objects, pstFlds);

        PtrList<pointSymmTensorField> psymtFlds;
        ReadFields(pMesh, objects, psymtFlds);

        PtrList<pointTensorField> ptFlds;
        ReadFields(pMesh, objects, ptFlds);

        execFlowFunctionObjects(args, runTime);
    }
    else
    {
        Info<< "    Reading phi" << endl;
        surfaceScalarField phi
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ
            ),
            mesh
        );

        Info<< "    Reading U" << endl;
        volVectorField U
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ
            ),
            mesh
        );

        Info<< "    Reading p" << endl;
        volScalarField p
        (
            IOobject
            (
                "p",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ
            ),
            mesh
        );

        if (phi.dimensions() == dimensionSet(0, 3, -1, 0, 0))
        {
            IOobject RASPropertiesHeader
            (
                "RASProperties",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            );

            IOobject LESPropertiesHeader
            (
                "LESProperties",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            );

            if (RASPropertiesHeader.headerOk())
            {
                IOdictionary RASProperties(RASPropertiesHeader);

                singlePhaseTransportModel laminarTransport(U, phi);

                autoPtr<incompressible::RASModel> RASModel
                (
                    incompressible::RASModel::New
                    (
                        U,
                        phi,
                        laminarTransport
                    )
                );
                execFlowFunctionObjects(args, runTime);
            }
            else if (LESPropertiesHeader.headerOk())
            {
                IOdictionary LESProperties(LESPropertiesHeader);

                singlePhaseTransportModel laminarTransport(U, phi);

                autoPtr<incompressible::LESModel> sgsModel
                (
                    incompressible::LESModel::New(U, phi, laminarTransport)
                );

                execFlowFunctionObjects(args, runTime);
            }
            else
            {
                IOdictionary transportProperties
                (
                    IOobject
                    (
                        "transportProperties",
                        runTime.constant(),
                        mesh,
                        IOobject::MUST_READ_IF_MODIFIED,
                        IOobject::NO_WRITE
                    )
                );

                dimensionedScalar nu(transportProperties.lookup("nu"));

                execFlowFunctionObjects(args, runTime);
            }
        }
        else if (phi.dimensions() == dimensionSet(1, 0, -1, 0, 0))
        {
            autoPtr<fluidThermo> thermo(fluidThermo::New(mesh));

            volScalarField rho
            (
                IOobject
                (
                    "rho",
                    runTime.timeName(),
                    mesh
                ),
                thermo->rho()
            );

            IOobject RASPropertiesHeader
            (
                "RASProperties",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            );

            IOobject LESPropertiesHeader
            (
                "LESProperties",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            );

            if (RASPropertiesHeader.headerOk())
            {
                IOdictionary RASProperties(RASPropertiesHeader);

                autoPtr<compressible::RASModel> RASModel
                (
                    compressible::RASModel::New
                    (
                        rho,
                        U,
                        phi,
                        thermo()
                    )
                );

                execFlowFunctionObjects(args, runTime);
            }
            else if (LESPropertiesHeader.headerOk())
            {
                IOdictionary LESProperties(LESPropertiesHeader);

                autoPtr<compressible::LESModel> sgsModel
                (
                    compressible::LESModel::New(rho, U, phi, thermo())
                );

                execFlowFunctionObjects(args, runTime);
            }
            else
            {
                IOdictionary transportProperties
                (
                    IOobject
                    (
                        "transportProperties",
                        runTime.constant(),
                        mesh,
                        IOobject::MUST_READ_IF_MODIFIED,
                        IOobject::NO_WRITE
                    )
                );

                dimensionedScalar mu(transportProperties.lookup("mu"));

                execFlowFunctionObjects(args, runTime);
            }
        }
        else
        {
            FatalErrorIn(args.executable())
                << "Incorrect dimensions of phi: " << phi.dimensions()
                << nl << exit(FatalError);
        }
    }
}


// ************************************************************************* //
