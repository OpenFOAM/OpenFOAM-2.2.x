{
    surfaceScalarField alpha1f(fvc::interpolate(alpha1));
    surfaceScalarField alpha2f(scalar(1) - alpha1f);

    volScalarField rAU1(1.0/U1Eqn.A());
    volScalarField rAU2(1.0/U2Eqn.A());

    rAU1f = 1.0/fvc::interpolate(U1Eqn.A());
    surfaceScalarField rAU2f(1.0/fvc::interpolate(U2Eqn.A()));

    volVectorField HbyA1("HbyA1", U1);
    HbyA1 = rAU1*U1Eqn.H();

    volVectorField HbyA2("HbyA2", U2);
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
        ppDrag -= rAU1f*fvc::snGrad(kineticTheory.pa()/rho1)*mesh.magSf();
    }

    surfaceScalarField phiHbyA1
    (
        "phiHbyA1",
        (fvc::interpolate(HbyA1) & mesh.Sf())
      + fvc::ddtPhiCorr(rAU1, U1, phi1)
    );

    surfaceScalarField phiHbyA2
    (
        "phiHbyA2",
        (fvc::interpolate(HbyA2) & mesh.Sf())
      + fvc::ddtPhiCorr(rAU2, U2, phi2)
    );

    phi = alpha1f*phiHbyA1 + alpha2f*phiHbyA2;
    mrfZones.relativeFlux(phi);

    phiHbyA1 +=
    (
        fvc::interpolate(alpha2/rho1*K*rAU1)*phi2
      + ppDrag
      + rAU1f*(g & mesh.Sf())
    );
    mrfZones.relativeFlux(phiHbyA1);

    phiHbyA2 +=
    (
        fvc::interpolate(alpha1/rho2*K*rAU2)*phi1
      + rAU2f*(g & mesh.Sf())
    );
    mrfZones.relativeFlux(phiHbyA2);

    mrfZones.relativeFlux(phi1.oldTime());
    mrfZones.relativeFlux(phi1);
    mrfZones.relativeFlux(phi2.oldTime());
    mrfZones.relativeFlux(phi2);

    surfaceScalarField phiHbyA("phiHbyA", alpha1f*phiHbyA1 + alpha2f*phiHbyA2);

    HbyA1 += alpha2*(1.0/rho1)*rAU1*K*U2;
    HbyA2 += alpha1*(1.0/rho2)*rAU2*K*U1;

    surfaceScalarField Dp
    (
        "Dp",
        alpha1f*rAU1f/rho1 + alpha2f*rAU2f/rho2
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(Dp, p) == fvc::div(phiHbyA)
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            surfaceScalarField SfGradp(pEqn.flux()/Dp);

            phi1 = phiHbyA1 - rAU1f*SfGradp/rho1;
            phi2 = phiHbyA2 - rAU2f*SfGradp/rho2;
            phi = alpha1f*phi1 + alpha2f*phi2;

            p.relax();
            SfGradp = pEqn.flux()/Dp;

            U1 = HbyA1
               + fvc::reconstruct
                 (
                     ppDrag
                   + rAU1f*((g & mesh.Sf()) - SfGradp/rho1)
                 );
            U1.correctBoundaryConditions();

            U2 = HbyA2
               + fvc::reconstruct
                 (
                     rAU2f*((g & mesh.Sf()) - SfGradp/rho2)
                 );
            U2.correctBoundaryConditions();

            U = alpha1*U1 + alpha2*U2;
        }
    }
}

#include "continuityErrs.H"
