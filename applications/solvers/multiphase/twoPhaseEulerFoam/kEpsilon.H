if (turbulence)
{
    if (mesh.changing())
    {
        y.correct();
    }

    tmp<volTensorField> tgradU2 = fvc::grad(U2);
    volScalarField G(2*nut2*(tgradU2() && dev(symm(tgradU2()))));
    tgradU2.clear();

    #include "wallFunctions.H"

    // Dissipation equation
    fvScalarMatrix epsEqn
    (
        fvm::ddt(epsilon)
      + fvm::div(phi2, epsilon)
      - fvm::Sp(fvc::div(phi2), epsilon)
      - fvm::laplacian
        (
            alpha1Eps*nuEff2, epsilon,
            "laplacian(DepsilonEff,epsilon)"
        )
      ==
         C1*G*epsilon/k
       - fvm::Sp(C2*epsilon/k, epsilon)
    );

    #include "wallDissipation.H"

    epsEqn.relax();
    epsEqn.solve();

    epsilon.max(dimensionedScalar("zero", epsilon.dimensions(), 1.0e-15));


    // Turbulent kinetic energy equation
    fvScalarMatrix kEqn
    (
        fvm::ddt(k)
      + fvm::div(phi2, k)
      - fvm::Sp(fvc::div(phi2), k)
      - fvm::laplacian
        (
            alpha1k*nuEff2, k,
            "laplacian(DkEff,k)"
        )
      ==
        G
      - fvm::Sp(epsilon/k, k)
    );
    kEqn.relax();
    kEqn.solve();

    k.max(dimensionedScalar("zero", k.dimensions(), 1.0e-8));

    //- Re-calculate turbulence viscosity
    nut2 = Cmu*sqr(k)/epsilon;

    #include "wallViscosity.H"
}

nuEff2 = nut2 + nu2;
