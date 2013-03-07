{
    kappat = turbulence->nut()/Prt;
    kappat.correctBoundaryConditions();

    volScalarField kappaEff("kappaEff", turbulence->nu()/Pr + kappat);

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(kappaEff, T)
     ==
        radiation->ST(rhoCpRef, T)
    );

    TEqn.relax();
    TEqn.solve();

    radiation->correct();

    rhok = 1.0 - beta*(T - TRef);
}
