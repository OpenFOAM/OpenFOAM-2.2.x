fvVectorMatrix U1Eqn(U1, U1.dimensions()*dimVol/dimTime);
fvVectorMatrix U2Eqn(U2, U2.dimensions()*dimVol/dimTime);

{
    {
        volTensorField gradU1T(T(fvc::grad(U1)));

        if (kineticTheory.on())
        {
            kineticTheory.solve(gradU1T);
            nuEff1 = kineticTheory.mu1()/rho1;
        }
        else // If not using kinetic theory is using Ct model
        {
            nuEff1 = sqr(Ct)*nut2 + nu1;
        }

        volTensorField Rc1
        (
            "Rc1",
            (((2.0/3.0)*I)*nuEff1)*tr(gradU1T) - nuEff1*gradU1T
        );

        if (kineticTheory.on())
        {
            Rc1 -= ((kineticTheory.lambda()/rho1)*tr(gradU1T))*tensor(I);
        }

        surfaceScalarField phiR1
        (
            -fvc::interpolate(nuEff1)*mesh.magSf()*fvc::snGrad(alpha1)
            /fvc::interpolate(alpha1 + scalar(0.001))
        );

        U1Eqn =
        (
            (scalar(1) + Cvm*rho2*alpha2/rho1)*
            (
                fvm::ddt(U1)
              + fvm::div(phi1, U1, "div(phi1,U1)")
              - fvm::Sp(fvc::div(phi1), U1)
            )

          - fvm::laplacian(nuEff1, U1)
          + fvc::div(Rc1)

          + fvm::div(phiR1, U1, "div(phi1,U1)")
          - fvm::Sp(fvc::div(phiR1), U1)
          + (fvc::grad(alpha1)/(fvc::average(alpha1) + scalar(0.001)) & Rc1)
         ==
        //  g                          // Buoyancy term transfered to p-equation
          - fvm::Sp(alpha2/rho1*K, U1)
        //+ alpha2/rho1*K*U2           // Explicit drag transfered to p-equation
          - alpha2/rho1*(liftCoeff - Cvm*rho2*DDtU2)
        );
        mrfZones.addCoriolis(scalar(1) + Cvm*rho2*alpha2/rho1, U1Eqn);
        U1Eqn.relax();
    }

    {
        volTensorField gradU2T(T(fvc::grad(U2)));
        volTensorField Rc2
        (
            "Rc2",
            (((2.0/3.0)*I)*nuEff2)*tr(gradU2T) - nuEff2*gradU2T
        );

        surfaceScalarField phiR2
        (
            -fvc::interpolate(nuEff2)*mesh.magSf()*fvc::snGrad(alpha2)
            /fvc::interpolate(alpha2 + scalar(0.001))
        );

        U2Eqn =
        (
            (scalar(1) + Cvm*rho2*alpha1/rho2)*
            (
                fvm::ddt(U2)
              + fvm::div(phi2, U2, "div(phi2,U2)")
              - fvm::Sp(fvc::div(phi2), U2)
            )

          - fvm::laplacian(nuEff2, U2)
          + fvc::div(Rc2)

          + fvm::div(phiR2, U2, "div(phi2,U2)")
          - fvm::Sp(fvc::div(phiR2), U2)

          + (fvc::grad(alpha2)/(fvc::average(alpha2) + scalar(0.001)) & Rc2)
         ==
        //  g                          // Buoyancy term transfered to p-equation
          - fvm::Sp(alpha1/rho2*K, U2)
        //+ alpha1/rho2*K*U1           // Explicit drag transfered to p-equation
          + alpha1/rho2*(liftCoeff + Cvm*rho2*DDtU1)
        );
        mrfZones.addCoriolis(scalar(1) + Cvm*rho2*alpha1/rho2, U2Eqn);
        U2Eqn.relax();
    }
}
