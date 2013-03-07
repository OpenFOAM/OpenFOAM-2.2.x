{
    solve
    (
        fvm::ddt(rho, T)
      + fvm::div(rhoPhi, T)
      - fvm::laplacian(twoPhaseProperties.alphaEff(turbulence->mut()), T)
      + (
            p*fvc::div(phi)
          + fvc::ddt(rho, K) + fvc::div(rhoPhi, K)
        )
       *(
           alpha1/twoPhaseProperties.thermo1().Cv()
         + alpha2/twoPhaseProperties.thermo2().Cv()
        )
    );

    twoPhaseProperties.correct();
}
