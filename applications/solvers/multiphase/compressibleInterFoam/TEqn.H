{
    volScalarField kByCv
    (
        "kByCv",
        (alpha1*k1/Cv1 + alpha2*k2/Cv2)
      + (alpha1*rho1 + alpha2*rho2)*turbulence->nut()
    );

    solve
    (
        fvm::ddt(rho, T)
      + fvm::div(rhoPhi, T)
      - fvm::laplacian(kByCv, T)
      + p*fvc::div(phi)*(alpha1/Cv1 + alpha2/Cv2)
    );

    // Update compressibilities
    psi1 = eos1->psi(p, T);
    psi2 = eos2->psi(p, T);
}
