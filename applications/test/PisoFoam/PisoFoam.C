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

Application
    pisoFoam

Description
    Transient solver for incompressible flow.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"

#include "LduMatrix.H"
#include "diagTensorField.H"

typedef LduMatrix<vector, scalar, scalar> lduVectorMatrix;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "readPISOControls.H"
        #include "CourantNo.H"

        // Pressure-velocity PISO corrector
        {
            // Momentum predictor

            fvVectorMatrix UEqn
            (
                fvm::ddt(U)
              + fvm::div(phi, U)
              + turbulence->divDevReff(U)
            );

            //UEqn.relax();

            fvVectorMatrix UEqnp(UEqn == -fvc::grad(p));

            lduVectorMatrix U3Eqnp(mesh);
            U3Eqnp.diag() = UEqnp.diag();
            U3Eqnp.upper() = UEqnp.upper();
            U3Eqnp.lower() = UEqnp.lower();
            U3Eqnp.source() = UEqnp.source();

            UEqnp.addBoundaryDiag(U3Eqnp.diag(), 0);
            UEqnp.addBoundarySource(U3Eqnp.source(), false);

            U3Eqnp.interfaces() = U.boundaryField().interfaces();
            U3Eqnp.interfacesUpper() = UEqnp.boundaryCoeffs().component(0);
            U3Eqnp.interfacesLower() = UEqnp.internalCoeffs().component(0);

            autoPtr<lduVectorMatrix::solver> U3EqnpSolver =
            lduVectorMatrix::solver::New
            (
                U.name(),
                U3Eqnp,
                dictionary
                (
                    IStringStream
                    (
                        "{"
                        "    /*solver          SmoothSolver;*/"
                        "    smoother        GaussSeidel;"
                        "    solver           PBiCCCG;"
                        "    preconditioner   none;"
                        "    tolerance        (1e-7 1e-7 1);"
                        "    relTol           (0 0 0);"
                        "}"
                    )()
                )
            );

            //for (int i=0; i<3; i++)
            {
                U3EqnpSolver->solve(U).print(Info);
                U.correctBoundaryConditions();
            }
            //solve(UEqnp);

            // --- PISO loop

            for (int corr=0; corr<nCorr; corr++)
            {
                volScalarField rAU(1.0/UEqn.A());

                volVectorField HbyA("HbyA", U);
                HbyA = rAU*UEqn.H();
                surfaceScalarField phiHbyA
                (
                    "phiHbyA",
                    (fvc::interpolate(HbyA) & mesh.Sf())
                  + fvc::ddtPhiCorr(rAU, U, phi)
                );

                adjustPhi(phiHbyA, U, p);

                // Non-orthogonal pressure corrector loop
                for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
                {
                    // Pressure corrector

                    fvScalarMatrix pEqn
                    (
                        fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                    );

                    pEqn.setReference(pRefCell, pRefValue);

                    if
                    (
                        corr == nCorr-1
                     && nonOrth == nNonOrthCorr
                    )
                    {
                        pEqn.solve(mesh.solver("pFinal"));
                    }
                    else
                    {
                        pEqn.solve();
                    }

                    if (nonOrth == nNonOrthCorr)
                    {
                        phi = phiHbyA - pEqn.flux();
                    }
                }

                #include "continuityErrs.H"

                U = HbyA - rAU*fvc::grad(p);
                U.correctBoundaryConditions();
            }
        }

        turbulence->correct();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
