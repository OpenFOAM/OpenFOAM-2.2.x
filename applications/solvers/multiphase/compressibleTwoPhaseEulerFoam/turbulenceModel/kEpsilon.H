if (turbulence)
{
    if (mesh.changing())
    {
        y.correct();
    }

    tmp<volTensorField> tgradU2(fvc::grad(U2));
    volScalarField G(2*nut2*(tgradU2() && dev(symm(tgradU2()))));
    tgradU2.clear();

    #include "wallFunctions.H"

    // Dissipation equation
    fvScalarMatrix epsEqn
    (
        fvm::ddt(alpha2, epsilon)
      + fvm::div(alphaPhi2, epsilon)

        // Compressibity correction
      - fvm::Sp(fvc::ddt(alpha2) + fvc::div(alphaPhi2), epsilon)

      - fvm::laplacian
        (
            alpha1Eps*nuEff2, epsilon,
            "laplacian(DepsilonEff,epsilon)"
        )
      ==
         C1*alpha2*G*epsilon/k
       - fvm::Sp(C2*alpha2*epsilon/k, epsilon)
    );

    #include "wallDissipation.H"



    epsEqn.relax();
    epsEqn.solve();

    epsilon.max(dimensionedScalar("zero", epsilon.dimensions(), 1.0e-15));


    // Turbulent kinetic energy equation
    fvScalarMatrix kEqn
    (
        fvm::ddt(alpha2, k)
      + fvm::div(alphaPhi2, k)

        // Compressibity correction
      - fvm::Sp(fvc::ddt(alpha2) + fvc::div(alphaPhi2), k)

      - fvm::laplacian
        (
            alpha1k*nuEff2, k,
            "laplacian(DkEff,k)"
        )
      ==
        alpha2*G
      - fvm::Sp(alpha2*epsilon/k, k)
    );
    kEqn.relax();
    kEqn.solve();

    k.max(dimensionedScalar("zero", k.dimensions(), 1.0e-8));

    //- Re-calculate turbulence viscosity
    nut2 = Cmu*sqr(k)/epsilon;

    #include "wallViscosity.H"
}

nuEff2 = nut2 + thermo2.mu()/rho2;
