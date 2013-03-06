    volVectorField Ur(U1 - U2);
    volScalarField magUr(mag(Ur) + residualSlip);

    volScalarField Ka(drag1->K(magUr));
    volScalarField K(Ka);

    if (dragPhase == "2")
    {
        volScalarField Kb(drag2->K(magUr));
        K = Kb;
    }
    else if (dragPhase == "blended")
    {
        volScalarField Kb(drag2->K(magUr));
        K = (alpha2*Ka + alpha1*Kb);
    }

    volVectorField liftCoeff
    (
        Cl*(alpha2*rho2 + alpha1*rho1)*(Ur ^ fvc::curl(U))
    );

    // Remove lift and drag at fixed-flux boundaries
    forAll(phi1.boundaryField(), patchi)
    {
        if (isA<fixedValueFvsPatchScalarField>(phi1.boundaryField()[patchi]))
        {
            K.boundaryField()[patchi] = 0.0;
            liftCoeff.boundaryField()[patchi] = vector::zero;
        }
    }
