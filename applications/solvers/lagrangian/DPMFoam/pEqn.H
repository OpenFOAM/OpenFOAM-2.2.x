{
    volVectorField HbyA("HbyA", Uc);
    HbyA = rAUc*UcEqn.H();

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
           (fvc::interpolate(HbyA) & mesh.Sf())
         + fvc::ddtPhiCorr(rAUc, alphac, Uc, phic)
         + phicForces
        )
    );

    // Update phi BCs before pEqn
    forAll(phic.boundaryField(), patchi)
    {
        phic.boundaryField()[patchi] =
            mesh.Sf().boundaryField()[patchi] & Uc.boundaryField()[patchi];
    }

    // Non-orthogonal pressure corrector loop
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(alphacf*rAUcf, p)
         ==
            fvc::ddt(alphac) + fvc::div(alphacf*phiHbyA)
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phic = phiHbyA - pEqn.flux()/alphacf;

            p.relax();

            Uc = HbyA
              + rAUc*fvc::reconstruct((phicForces - pEqn.flux()/alphacf)/rAUcf);
            Uc.correctBoundaryConditions();
        }
    }
}

#include "continuityErrs.H"
