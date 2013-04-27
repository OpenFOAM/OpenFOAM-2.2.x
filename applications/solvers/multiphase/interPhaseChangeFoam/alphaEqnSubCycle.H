surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("0", dimMass/dimTime, 0)
);

{
    const dictionary& pimpleDict = pimple.dict();

    label nAlphaCorr(readLabel(pimpleDict.lookup("nAlphaCorr")));

    label nAlphaSubCycles(readLabel(pimpleDict.lookup("nAlphaSubCycles")));

    surfaceScalarField phic(mag(phi/mesh.magSf()));
    phic = min(interface.cAlpha()*phic, max(phic));

    volScalarField divU(fvc::div(phi));

    dimensionedScalar totalDeltaT = runTime.deltaT();

    if (nAlphaSubCycles > 1)
    {
        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            #include "alphaEqn.H"
        }
    }
    else
    {
        #include "alphaEqn.H"
    }

    rho == alpha1*rho1 + alpha2*rho2;
}
