mrfZones.correctBoundaryVelocity(U1);
mrfZones.correctBoundaryVelocity(U2);
mrfZones.correctBoundaryVelocity(U);

fvVectorMatrix U1Eqn(U1, U1.dimensions()*dimVol/dimTime);
fvVectorMatrix U2Eqn(U2, U2.dimensions()*dimVol/dimTime);

{
    {
        volTensorField gradU1T(fvc::grad(U1)().T());

        if (kineticTheory.on())
        {
            kineticTheory.solve(gradU1T);
            nuEff1 = kineticTheory.mu1()/rho1;
        }
        else // If not using kinetic theory is using Ct model
        {
            nuEff1 = sqr(Ct)*nut2 + thermo1.mu()/rho1;
        }

        volTensorField Rc1
        (
            "Rc",
            (((2.0/3.0)*I)*nuEff1)*tr(gradU1T) - nuEff1*gradU1T
        );

        if (kineticTheory.on())
        {
            Rc1 -= ((kineticTheory.lambda()/rho1)*tr(gradU1T))*tensor(I);
        }

        U1Eqn =
        (
            fvm::ddt(alpha1, U1)
          + fvm::div(alphaPhi1, U1)

            // Compressibity correction
          - fvm::Sp(fvc::ddt(alpha1) + fvc::div(alphaPhi1), U1)

          + Cvm*rho2*alpha1*alpha2/rho1*
            (
                fvm::ddt(U1)
              + fvm::div(phi1, U1)
              - fvm::Sp(fvc::div(phi1), U1)
            )

          - fvm::laplacian(alpha1*nuEff1, U1)
          + fvc::div(alpha1*Rc1)
         ==
          - fvm::Sp(dragCoeff/rho1, U1)
          - alpha1*alpha2/rho1*(liftForce - Cvm*rho2*DDtU2)
        );
        mrfZones.addCoriolis(alpha1*(1 + Cvm*rho2*alpha2/rho1), U1Eqn);
        U1Eqn.relax();
    }

    {
        volTensorField gradU2T(fvc::grad(U2)().T());
        volTensorField Rc2
        (
            "Rc",
            (((2.0/3.0)*I)*nuEff2)*tr(gradU2T) - nuEff2*gradU2T
        );

        U2Eqn =
        (
            fvm::ddt(alpha2, U2)
          + fvm::div(alphaPhi2, U2)

            // Compressibity correction
          - fvm::Sp(fvc::ddt(alpha2) + fvc::div(alphaPhi2), U2)

          + Cvm*rho2*alpha1*alpha2/rho2*
            (
                fvm::ddt(U2)
              + fvm::div(phi2, U2)
              - fvm::Sp(fvc::div(phi2), U2)
            )

          - fvm::laplacian(alpha2*nuEff2, U2)
          + fvc::div(alpha2*Rc2)
         ==
          - fvm::Sp(dragCoeff/rho2, U2)
          + alpha1*alpha2/rho2*(liftForce + Cvm*rho2*DDtU1)
        );
        mrfZones.addCoriolis(alpha2*(1 + Cvm*rho2*alpha1/rho2), U2Eqn);
        U2Eqn.relax();
    }
}
