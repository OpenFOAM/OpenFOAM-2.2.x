{
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("Dp", fvc::interpolate(rAU));

    volVectorField HbyA("HbyA", U);
    HbyA = rAU*UEqn.H();

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (fvc::interpolate(HbyA) & mesh.Sf())
      + fvc::ddtPhiCorr(rAU, rho, U, phi)
    );
    phi = phiHbyA;

    surfaceScalarField phig
    (
        (
            fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    tmp<fvScalarMatrix> p_rghEqnComp1;
    tmp<fvScalarMatrix> p_rghEqnComp2;

    //if (transonic)
    //{
    //}
    //else
    {
        surfaceScalarField phid1("phid1", fvc::interpolate(psi1)*phi);
        surfaceScalarField phid2("phid2", fvc::interpolate(psi2)*phi);

        p_rghEqnComp1 =
            fvc::ddt(rho1) + psi1*correction(fvm::ddt(p_rgh))
          + fvc::div(phid1, p_rgh)
          - fvc::Sp(fvc::div(phid1), p_rgh);

        p_rghEqnComp2 =
            fvc::ddt(rho2) + psi2*correction(fvm::ddt(p_rgh))
          + fvc::div(phid2, p_rgh)
          - fvc::Sp(fvc::div(phid2), p_rgh);
    }

    // Thermodynamic density needs to be updated by psi*d(p) after the
    // pressure solution - done in 2 parts. Part 1:
    //twoPhaseProperties.rho() -= psi*p_rgh;

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqnIncomp
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rAUf, p_rgh)
        );

        solve
        (
            (
                (max(alpha1, scalar(0))/rho1)*p_rghEqnComp1()
              + (max(alpha2, scalar(0))/rho2)*p_rghEqnComp2()
            )
          + p_rghEqnIncomp,
            mesh.solver(p_rgh.select(pimple.finalInnerIter()))
        );

        if (pimple.finalNonOrthogonalIter())
        {
            // Second part of thermodynamic density update
            //twoPhaseProperties.rho() += psi*p_rgh;

            dgdt =
            (
                pos(alpha2)*(p_rghEqnComp2 & p_rgh)/rho2
              - pos(alpha1)*(p_rghEqnComp1 & p_rgh)/rho1
            );

            phi = phiHbyA + p_rghEqnIncomp.flux();

            U = HbyA
              + rAU*fvc::reconstruct((phig + p_rghEqnIncomp.flux())/rAUf);
            U.correctBoundaryConditions();
        }
    }

    p = max(p_rgh + (alpha1*rho1 + alpha2*rho2)*gh, pMin);

    twoPhaseProperties.correct();

    K = 0.5*magSqr(U);

    Info<< "max(U) " << max(mag(U)).value() << endl;
    Info<< "min(p_rgh) " << min(p_rgh).value() << endl;
}
