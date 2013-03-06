{
    rho1 = thermo1.rho();
    rho2 = thermo2.rho();

    surfaceScalarField alpha1f(fvc::interpolate(alpha1));
    surfaceScalarField alpha2f(scalar(1) - alpha1f);

    rAU1 = 1.0/U1Eqn.A();
    volScalarField rAU2(1.0/U2Eqn.A());

    surfaceScalarField rAlphaAU1f(fvc::interpolate(alpha1*rAU1));
    surfaceScalarField rAlphaAU2f(fvc::interpolate(alpha2*rAU2));

    volVectorField HbyA1("HbyA" + phase1Name, U1);
    HbyA1 = rAU1*U1Eqn.H();

    volVectorField HbyA2("HbyA" + phase2Name, U2);
    HbyA2 = rAU2*U2Eqn.H();

    mrfZones.absoluteFlux(phi1.oldTime());
    mrfZones.absoluteFlux(phi1);
    mrfZones.absoluteFlux(phi2.oldTime());
    mrfZones.absoluteFlux(phi2);

    surfaceScalarField ppDrag("ppDrag", 0.0*phi1);

    if (g0.value() > 0.0)
    {
        ppDrag -= ppMagf*fvc::snGrad(alpha1)*mesh.magSf();
    }

    if (kineticTheory.on())
    {
        ppDrag -=
            fvc::interpolate(1.0/rho1)*rAlphaAU1f
           *fvc::snGrad(kineticTheory.pa())*mesh.magSf();
    }

    surfaceScalarField phiHbyA1
    (
        "phiHbyA" + phase1Name,
        (fvc::interpolate(HbyA1) & mesh.Sf())
      + fvc::ddtPhiCorr(rAU1, alpha1, U1, phi1)
    );

    surfaceScalarField phiHbyA2
    (
        "phiHbyA" + phase2Name,
        (fvc::interpolate(HbyA2) & mesh.Sf())
      + fvc::ddtPhiCorr(rAU2, alpha2, U2, phi2)
    );

    phi = alpha1f*phiHbyA1 + alpha2f*phiHbyA2;
    mrfZones.relativeFlux(phi);

    phiHbyA1 +=
    (
        fvc::interpolate((1.0/rho1)*rAU1*dragCoeff)*phi2
      + ppDrag
      + rAlphaAU1f*(g & mesh.Sf())
    );
    mrfZones.relativeFlux(phiHbyA1);

    phiHbyA2 +=
    (
        fvc::interpolate((1.0/rho2)*rAU2*dragCoeff)*phi1
      + rAlphaAU2f*(g & mesh.Sf())
    );
    mrfZones.relativeFlux(phiHbyA2);

    mrfZones.relativeFlux(phi1.oldTime());
    mrfZones.relativeFlux(phi1);
    mrfZones.relativeFlux(phi2.oldTime());
    mrfZones.relativeFlux(phi2);

    surfaceScalarField phiHbyA("phiHbyA", alpha1f*phiHbyA1 + alpha2f*phiHbyA2);

    HbyA1 += (1.0/rho1)*rAU1*dragCoeff*U2;
    HbyA2 += (1.0/rho2)*rAU2*dragCoeff*U1;

    surfaceScalarField Dp
    (
        "Dp",
        mag
        (
            alpha1f*rAlphaAU1f/fvc::interpolate(rho1)
          + alpha2f*rAlphaAU2f/fvc::interpolate(rho2)
        )
    );

    tmp<fvScalarMatrix> pEqnComp1;
    tmp<fvScalarMatrix> pEqnComp2;

    //if (transonic)
    //{
    //}
    //else
    {
        surfaceScalarField phid1
        (
            "phid" + phase1Name,
            fvc::interpolate(psi1)*phi1
        );
        surfaceScalarField phid2
        (
            "phid" + phase2Name,
            fvc::interpolate(psi2)*phi2
        );

        pEqnComp1 =
            fvc::ddt(rho1) + psi1*correction(fvm::ddt(p))
          + fvc::div(phid1, p)
          - fvc::Sp(fvc::div(phid1), p);

        pEqnComp2 =
            fvc::ddt(rho2) + psi2*correction(fvm::ddt(p))
          + fvc::div(phid2, p)
          - fvc::Sp(fvc::div(phid2), p);
    }

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqnIncomp
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(Dp, p)
        );

        solve
        (
            (
                (alpha1/rho1)*pEqnComp1()
              + (alpha2/rho2)*pEqnComp2()
            )
          + pEqnIncomp,
            mesh.solver(p.select(pimple.finalInnerIter()))
        );

        if (pimple.finalNonOrthogonalIter())
        {
            surfaceScalarField mSfGradp(pEqnIncomp.flux()/Dp);

            phi1 = phiHbyA1 + rAlphaAU1f*mSfGradp/fvc::interpolate(rho1);
            phi2 = phiHbyA2 + rAlphaAU2f*mSfGradp/fvc::interpolate(rho2);
            phi = alpha1f*phi1 + alpha2f*phi2;

            dgdt =
            (
                pos(alpha2)*(pEqnComp2 & p)/rho2
              - pos(alpha1)*(pEqnComp1 & p)/rho1
            );

            p.relax();
            mSfGradp = pEqnIncomp.flux()/Dp;

            U1 = HbyA1
               + fvc::reconstruct
                 (
                     ppDrag
                   + rAlphaAU1f
                    *(
                         (g & mesh.Sf())
                       + mSfGradp/fvc::interpolate(rho1)
                     )
                 );
            U1.correctBoundaryConditions();

            U2 = HbyA2
               + fvc::reconstruct
                 (
                     rAlphaAU2f
                    *(
                        (g & mesh.Sf())
                      + mSfGradp/fvc::interpolate(rho2)
                    )
                 );
            U2.correctBoundaryConditions();

            U = alpha1*U1 + alpha2*U2;
        }
    }

    p = max(p, pMin);

    thermo1.correct();
    thermo2.correct();
    rho1 = thermo1.rho();
    rho2 = thermo2.rho();

    K1 = 0.5*magSqr(U1);
    K2 = 0.5*magSqr(U2);

    if (thermo1.dpdt())
    {
        dalpha1pdt = fvc::ddt(alpha1, p);
    }

    if (thermo2.dpdt())
    {
        dalpha2pdt = fvc::ddt(alpha2, p);
    }
}
