const dictionary& pimpleDict = pimple.dict();
label nAlphaCorr(readLabel(pimpleDict.lookup("nAlphaCorr")));
label nAlphaSubCycles(readLabel(pimpleDict.lookup("nAlphaSubCycles")));

if (nAlphaSubCycles > 1)
{
    surfaceScalarField rhoPhiSum(0.0*rhoPhi);
    dimensionedScalar totalDeltaT = runTime.deltaT();

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

{
    volScalarField rhoNew(alpha1*rho1 + alpha2*rho2 + alpha3*rho3);

    //solve(fvm::ddt(rho) + fvc::div(rhoPhi));
    //Info<< "density error = "
    //    << max((mag(rho - rhoNew)/mag(rhoNew))().internalField()) << endl;

    rho == rhoNew;
}
