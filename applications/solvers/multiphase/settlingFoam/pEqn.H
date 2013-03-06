{
    volScalarField rAU("rAU", 1.0/UEqn.A());

    surfaceScalarField rAUf("Dp", fvc::interpolate(rho*rAU));

    volVectorField HbyA("HbyA", U);
    HbyA = rAU*UEqn.H();

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::interpolate(rho)
       *(
           (fvc::interpolate(HbyA) & mesh.Sf())
         + fvc::ddtPhiCorr(rAU, rho, U, phi)
        )
    );
    phi = phiHbyA;

    surfaceScalarField phig
    (
        - ghf*fvc::snGrad(rho)*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh) == fvc::ddt(rho) + fvc::div(phiHbyA)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - p_rghEqn.flux();

            U = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
            U.correctBoundaryConditions();
        }
    }

    p == p_rgh + rho*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }

    #include "rhoEqn.H"
    #include "compressibleContinuityErrs.H"
}
