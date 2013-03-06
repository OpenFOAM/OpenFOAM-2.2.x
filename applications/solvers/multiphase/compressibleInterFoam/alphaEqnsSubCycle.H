{
    label nAlphaCorr(readLabel(pimple.dict().lookup("nAlphaCorr")));

    label nAlphaSubCycles(readLabel(pimple.dict().lookup("nAlphaSubCycles")));

    surfaceScalarField phic(mag(phi/mesh.magSf()));
    phic = min(interface.cAlpha()*phic, max(phic));

    volScalarField divU(fvc::div(fvc::absolute(phi, U)));

    if (nAlphaSubCycles > 1)
    {
        dimensionedScalar totalDeltaT = runTime.deltaT();
        surfaceScalarField rhoPhiSum(0.0*rhoPhi);

        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            #include "alphaEqns.H"
            rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
        }

        rhoPhi = rhoPhiSum;
    }
    else
    {
        #include "alphaEqns.H"
    }
}
