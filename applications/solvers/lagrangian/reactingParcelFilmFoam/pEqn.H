rho = thermo.rho();

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("Dp", fvc::interpolate(rho*rAU));
volVectorField HbyA("HbyA", U);
HbyA = rAU*UEqn.H();

surfaceScalarField phig(-rhorAUf*ghf*fvc::snGrad(rho)*mesh.magSf());

surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::interpolate(rho)
   *(
        (fvc::interpolate(HbyA) & mesh.Sf())
      + fvc::ddtPhiCorr(rAU, rho, U, phi)
    )
  + phig
);

fvOptions.relativeFlux(fvc::interpolate(rho), phiHbyA);

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix p_rghEqn
    (
        fvc::ddt(psi, rho)*gh
      + fvc::div(phiHbyA)
      + fvm::ddt(psi, p_rgh)
      - fvm::laplacian(rhorAUf, p_rgh)
     ==
        parcels.Srho()
      + surfaceFilm.Srho()
      + fvOptions(psi, p_rgh, rho.name())
    );

    fvOptions.constrain(p_rghEqn);

    p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA + p_rghEqn.flux();
        U = HbyA + rAU*fvc::reconstruct((p_rghEqn.flux() + phig)/rhorAUf);
        U.correctBoundaryConditions();
        fvOptions.correct(U);
    }
}

p = p_rgh + rho*gh;

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

K = 0.5*magSqr(U);

if (thermo.dpdt())
{
    dpdt = fvc::ddt(p);
}
