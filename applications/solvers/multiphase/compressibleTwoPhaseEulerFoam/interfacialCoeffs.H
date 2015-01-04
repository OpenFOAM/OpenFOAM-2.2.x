volScalarField dragCoeff
(
    IOobject
    (
        "dragCoeff",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("dragCoeff", dimensionSet(1, -3, -1, 0, 0), 0)
);

volVectorField liftForce
(
    IOobject
    (
        "liftForce",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedVector("liftForce", dimensionSet(1, -2, -2, 0, 0), vector::zero)
);

volScalarField heatTransferCoeff
(
    IOobject
    (
        "heatTransferCoeff",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("heatTransferCoeff", dimensionSet(1, -1, -3, -1, 0), 0)
);

{
    volVectorField Ur(U1 - U2);
    volScalarField magUr(mag(Ur) + residualSlip);

    if (dispersedPhase == phase1Name)
    {
        volScalarField alphaCoeff(max(alpha1, residualPhaseFraction));

        dragCoeff = alphaCoeff*drag1->K(magUr);
        heatTransferCoeff = alphaCoeff*heatTransfer1->K(magUr);
    }
    else if (dispersedPhase == phase2Name)
    {
        volScalarField alphaCoeff(max(alpha2, residualPhaseFraction));

        dragCoeff = alphaCoeff*drag2->K(magUr);
        heatTransferCoeff = alphaCoeff*heatTransfer2->K(magUr);
    }
    else if (dispersedPhase == "both")
    {
        // Blend with the phase-fraction of the continuous phase
        volScalarField alphaCoeff(max(alpha1*alpha2, residualPhaseFraction));

        dragCoeff = alphaCoeff*(drag1->K(magUr) + drag2->K(magUr));

        heatTransferCoeff = alphaCoeff*
        (
            heatTransfer1->K(magUr)
          + heatTransfer2->K(magUr)
        );
    }
    else
    {
        FatalErrorIn(args.executable())
            << "dispersedPhase: " << dispersedPhase << " is incorrect"
            << exit(FatalError);
    }

    liftForce = Cl*(alpha1*rho1 + alpha2*rho2)*(Ur ^ fvc::curl(U));

    // Remove lift, drag and phase heat-transfer at fixed-flux boundaries
    forAll(phi1.boundaryField(), patchi)
    {
        if (isA<fixedValueFvsPatchScalarField>(phi1.boundaryField()[patchi]))
        {
            dragCoeff.boundaryField()[patchi] = 0.0;
            heatTransferCoeff.boundaryField()[patchi] = 0.0;
            liftForce.boundaryField()[patchi] = vector::zero;
        }
    }
}
