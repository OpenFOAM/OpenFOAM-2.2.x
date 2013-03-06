rho = thermo.rho();

volScalarField rAU(1.0/UEqn.A());
volVectorField HbyA("HbyA", U);
HbyA = rAU*UEqn.H();

surfaceScalarField phid
(
    "phid",
    fvc::interpolate(psi)
   *(
        (fvc::interpolate(HbyA) & mesh.Sf())
      - fvc::meshPhi(rho, U)
    )
);

volScalarField Dp("Dp", rho*rAU);

for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
{
    fvScalarMatrix pEqn
    (
        fvm::ddt(psi, p)
      + fvm::div(phid, p)
      - fvm::laplacian(Dp, p)
    );

    pEqn.solve();

    phi = pEqn.flux();
}

#include "compressibleContinuityErrs.H"

U = HbyA - rAU*fvc::grad(p);
U.correctBoundaryConditions();
