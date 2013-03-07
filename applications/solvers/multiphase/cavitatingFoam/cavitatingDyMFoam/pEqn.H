{
    if (pimple.nCorrPIMPLE() == 1)
    {
        p =
        (
            rho
          - gamma2*rhol0
          - ((gamma*psiv + gamma2*psil) - psi)*pSat
        )/psi;
    }

    surfaceScalarField rhof("rhof", fvc::interpolate(rho));

    volScalarField rAU(1.0/UEqn.A());
    surfaceScalarField rAUf("Dp", rhof*fvc::interpolate(rAU));

    volVectorField HbyA("HbyA", U);
    HbyA = rAU*UEqn.H();

    phiv = (fvc::interpolate(HbyA) & mesh.Sf())
         + fvc::ddtPhiCorr(rAU, rho, U, phivAbs);
    fvc::makeRelative(phiv, U);

    surfaceScalarField phiGradp(rAUf*mesh.magSf()*fvc::snGrad(p));

    phiv -= phiGradp/rhof;

    volScalarField rho0(rho - psi*p);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvc::ddt(rho)
          + psi*correction(fvm::ddt(p))
          + fvc::div(phiv, rho)
          + fvc::div(phiGradp)
          - fvm::laplacian(rAUf, p)
        );

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phiv += (phiGradp + pEqn.flux())/rhof;
            phivAbs = fvc::absolute(phiv, U);
        }
    }

    Info<< "Predicted p max-min : " << max(p).value()
        << " " << min(p).value() << endl;

    rho == max(rho0 + psi*p, rhoMin);

    #include "gammaPsi.H"

    p =
    (
        rho
      - gamma2*rhol0
      - ((gamma*psiv + gamma2*psil) - psi)*pSat
    )/psi;

    p.correctBoundaryConditions();

    Info<< "Phase-change corrected p max-min : " << max(p).value()
        << " " << min(p).value() << endl;

    // Correct velocity

    U = HbyA - rAU*fvc::grad(p);

    // Remove the swirl component of velocity for "wedge" cases
    if (pimple.dict().found("removeSwirl"))
    {
        label swirlCmpt(readLabel(pimple.dict().lookup("removeSwirl")));

        Info<< "Removing swirl component-" << swirlCmpt << " of U" << endl;
        U.field().replace(swirlCmpt, 0.0);
    }

    U.correctBoundaryConditions();

    Info<< "max(U) " << max(mag(U)).value() << endl;
}
