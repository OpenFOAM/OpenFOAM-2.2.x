{
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

    if (pimple.transonic())
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
            fvc::ddt(rho1)
          + fvc::div(phi1, rho1) - fvc::Sp(fvc::div(phi1), rho1)
          + correction
            (
                psi1*fvm::ddt(p)
              + fvm::div(phid1, p) - fvm::Sp(fvc::div(phid1), p)
            );
        deleteDemandDrivenData(pEqnComp1().faceFluxCorrectionPtr());
        pEqnComp1().relax();

        pEqnComp2 =
            fvc::ddt(rho2)
          + fvc::div(phi2, rho2) - fvc::Sp(fvc::div(phi2), rho2)
          + correction
            (
                psi2*fvm::ddt(p)
              + fvm::div(phid2, p) - fvm::Sp(fvc::div(phid2), p)
            );
        deleteDemandDrivenData(pEqnComp2().faceFluxCorrectionPtr());
        pEqnComp2().relax();
    }
    else
    {
        pEqnComp1 =
            fvc::ddt(rho1) + psi1*correction(fvm::ddt(p))
          + fvc::div(phi1, rho1) - fvc::Sp(fvc::div(phi1), rho1);

        pEqnComp2 =
            fvc::ddt(rho2) + psi2*correction(fvm::ddt(p))
          + fvc::div(phi2, rho2) - fvc::Sp(fvc::div(phi2), rho2);
    }

    // Cache p prior to solve for density update
    volScalarField p_0(p);

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

            mrfZones.absoluteFlux(phi1);
            phi1 == (mesh.Sf() & fvc::interpolate(U1));
            mrfZones.relativeFlux(phi1);
            phi1 = phiHbyA1 + rAlphaAU1f*mSfGradp/fvc::interpolate(rho1);

            mrfZones.absoluteFlux(phi2);
            phi2 == (mesh.Sf() & fvc::interpolate(U2));
            mrfZones.relativeFlux(phi2);
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

    // Update densities from change in p
    rho1 += psi1*(p - p_0);
    rho2 += psi2*(p - p_0);

    K1 = 0.5*magSqr(U1);
    K2 = 0.5*magSqr(U2);

    if (thermo1.dpdt() || thermo2.dpdt())
    {
        dpdt = fvc::ddt(p);
    }
}
