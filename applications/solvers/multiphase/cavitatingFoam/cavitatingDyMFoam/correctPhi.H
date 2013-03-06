{
    wordList pcorrTypes
    (
        p.boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
    );

    forAll (p.boundaryField(), i)
    {
        if (p.boundaryField()[i].fixesValue())
        {
            pcorrTypes[i] = fixedValueFvPatchScalarField::typeName;
        }
    }

    volScalarField pcorr
    (
        IOobject
        (
            "pcorr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("pcorr", p.dimensions(), 0.0),
        pcorrTypes
    );

    surfaceScalarField rhof(fvc::interpolate(rho, "div(phiv,rho)"));
    dimensionedScalar rAUf("(1|A(U))", dimTime, 1.0);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pcorrEqn
        (
            fvm::laplacian(rAUf, pcorr) == fvc::ddt(rho) + fvc::div(phiv*rhof)
        );

        pcorrEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phiv -= pcorrEqn.flux()/rhof;
        }
    }
}
