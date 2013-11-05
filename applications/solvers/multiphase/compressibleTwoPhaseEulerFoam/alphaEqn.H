{
    word alphaScheme("div(phi," + alpha1.name() + ')');
    word alpharScheme("div(phir," + alpha1.name() + ')');

    surfaceScalarField phic("phic", phi);
    surfaceScalarField phir("phir", phi1 - phi2);

    if (g0.value() > 0.0)
    {
        surfaceScalarField alpha1f(fvc::interpolate(alpha1));
        surfaceScalarField phipp(ppMagf*fvc::snGrad(alpha1)*mesh.magSf());
        phir += phipp;
        phic += alpha1f*phipp;
    }

    for (int acorr=0; acorr<nAlphaCorr; acorr++)
    {
        volScalarField::DimensionedInternalField Sp
        (
            IOobject
            (
                "Sp",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("Sp", dgdt.dimensions(), 0.0)
        );

        volScalarField::DimensionedInternalField Su
        (
            IOobject
            (
                "Su",
                runTime.timeName(),
                mesh
            ),
            // Divergence term is handled explicitly to be
            // consistent with the explicit transport solution
            fvc::div(phi)*min(alpha1, scalar(1))
        );

        forAll(dgdt, celli)
        {
            if (dgdt[celli] > 0.0 && alpha1[celli] > 0.0)
            {
                Sp[celli] -= dgdt[celli]*alpha1[celli];
                Su[celli] += dgdt[celli]*alpha1[celli];
            }
            else if (dgdt[celli] < 0.0 && alpha1[celli] < 1.0)
            {
                Sp[celli] += dgdt[celli]*(1.0 - alpha1[celli]);
            }
        }

        dimensionedScalar totalDeltaT = runTime.deltaT();
        if (nAlphaSubCycles > 1)
        {
            alphaPhi1 = dimensionedScalar("0", alphaPhi1.dimensions(), 0);
        }

        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            surfaceScalarField alphaPhic1
            (
                fvc::flux
                (
                    phic,
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

            // Ensure that the flux at inflow BCs is preserved
            forAll(alphaPhic1.boundaryField(), patchi)
            {
                fvsPatchScalarField& alphaPhic1p =
                    alphaPhic1.boundaryField()[patchi];

                if (!alphaPhic1p.coupled())
                {
                    const scalarField& phi1p = phi1.boundaryField()[patchi];
                    const scalarField& alpha1p = alpha1.boundaryField()[patchi];

                    forAll(alphaPhic1p, facei)
                    {
                        if (phi1p[facei] < 0)
                        {
                            alphaPhic1p[facei] = alpha1p[facei]*phi1p[facei];
                        }
                    }
                }
            }

            MULES::explicitSolve
            (
                geometricOneField(),
                alpha1,
                phi,
                alphaPhic1,
                Sp,
                Su,
                (g0.value() > 0 ? alphaMax : 1),
                0
            );

            if (nAlphaSubCycles > 1)
            {
                alphaPhi1 += (runTime.deltaT()/totalDeltaT)*alphaPhic1;
            }
            else
            {
                alphaPhi1 = alphaPhic1;
            }

            /*
            // Legacy semi-implicit and potentially unbounded form
            fvScalarMatrix alpha1Eqn
            (
                fvm::ddt(alpha1)
              + fvm::div(phic, alpha1, alphaScheme)
              + fvm::div
                (
                    -fvc::flux(-phir, alpha2, alpharScheme),
                    alpha1,
                    alpharScheme
                )
             ==
               fvm::Sp(Sp, alpha1) + Su
            );

            alpha1Eqn.relax();
            alpha1Eqn.solve();

            if (nAlphaSubCycles > 1)
            {
                alphaPhi1 += (runTime.deltaT()/totalDeltaT)*alpha1Eqn.flux();
            }
            else
            {
                alphaPhi1 = alpha1Eqn.flux();
            }
            */
        }

        if (g0.value() > 0.0)
        {
            surfaceScalarField alpha1f(fvc::interpolate(alpha1));

            ppMagf =
                fvc::interpolate((1.0/rho1)*rAU1)
               *g0*min(exp(preAlphaExp*(alpha1f - alphaMax)), expMax);

            fvScalarMatrix alpha1Eqn
            (
                fvm::ddt(alpha1) - fvc::ddt(alpha1)
              - fvm::laplacian
                (
                    alpha1f*ppMagf,
                    alpha1,
                    "laplacian(alpha1PpMag,alpha1)"
                )
            );

            alpha1Eqn.relax();
            alpha1Eqn.solve();

            #include "packingLimiter.H"

            alphaPhi1 += alpha1Eqn.flux();
        }

        alphaPhi2 = phi - alphaPhi1;
        alpha2 = scalar(1) - alpha1;

        Info<< "Dispersed phase volume fraction = "
            << alpha1.weightedAverage(mesh.V()).value()
            << "  Min(alpha1) = " << min(alpha1).value()
            << "  Max(alpha1) = " << max(alpha1).value()
            << endl;
    }
}

rho = alpha1*rho1 + alpha2*rho2;
