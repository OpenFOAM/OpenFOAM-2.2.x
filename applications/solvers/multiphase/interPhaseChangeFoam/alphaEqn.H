{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phir("phir", phic*interface.nHatf());

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        surfaceScalarField phiAlpha
        (
            fvc::flux
            (
                phi,
                alpha1,
                alphaScheme
            )
          + fvc::flux
            (
                -fvc::flux(-phir, alpha2, alpharScheme),
                alpha1,
                alpharScheme
            )
        );

        Pair<tmp<volScalarField> > vDotAlphal =
            twoPhaseProperties->vDotAlphal();
        const volScalarField& vDotcAlphal = vDotAlphal[0]();
        const volScalarField& vDotvAlphal = vDotAlphal[1]();

        //MULES::explicitSolve
        //(
        //    geometricOneField(),
        //    alpha1,
        //    phi,
        //    phiAlpha,
        //    (vDotvAlphal - vDotcAlphal)(),
        //    (divU*alpha1 + vDotcAlphal)(),
        //    1,
        //    0
        //);

        MULES::implicitSolve
        (
            geometricOneField(),
            alpha1,
            phi,
            phiAlpha,
            (divU + vDotvAlphal - vDotcAlphal)(),
            vDotcAlphal,
            1,
            0
        );

        alpha2 = 1.0 - alpha1;
        rhoPhi = phiAlpha*(rho1 - rho2) + phi*rho2;
    }

    Info<< "Liquid phase volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(alpha1) = " << min(alpha1).value()
        << "  Max(alpha1) = " << max(alpha1).value()
        << endl;
}
