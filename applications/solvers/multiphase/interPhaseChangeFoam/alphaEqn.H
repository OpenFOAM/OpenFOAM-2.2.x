{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phir("phir", phic*interface.nHatf());

    for (int gCorr=0; gCorr<nAlphaCorr; gCorr++)
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
                -fvc::flux(-phir, scalar(1) - alpha1, alpharScheme),
                alpha1,
                alpharScheme
            )
        );

        Pair<tmp<volScalarField> > vDotAlphal =
            twoPhaseProperties->vDotAlphal();
        const volScalarField& vDotcAlphal = vDotAlphal[0]();
        const volScalarField& vDotvAlphal = vDotAlphal[1]();

        volScalarField Sp
        (
            IOobject
            (
                "Sp",
                runTime.timeName(),
                mesh
            ),
            vDotvAlphal - vDotcAlphal
        );

        volScalarField Su
        (
            IOobject
            (
                "Su",
                runTime.timeName(),
                mesh
            ),
            // Divergence term is handled explicitly to be
            // consistent with the explicit transport solution
            divU*alpha1
          + vDotcAlphal
        );

        //MULES::explicitSolve
        //(
        //    geometricOneField(),
        //    alpha1,
        //    phi,
        //    phiAlpha,
        //    Sp,
        //    Su,
        //    1,
        //    0
        //);

        MULES::implicitSolve
        (
            geometricOneField(),
            alpha1,
            phi,
            phiAlpha,
            Sp,
            Su,
            1,
            0
        );

        rhoPhi +=
            (runTime.deltaT()/totalDeltaT)
           *(phiAlpha*(rho1 - rho2) + phi*rho2);
    }

    Info<< "Liquid phase volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(alpha1) = " << min(alpha1).value()
        << "  Max(alpha1) = " << max(alpha1).value()
        << endl;
}
