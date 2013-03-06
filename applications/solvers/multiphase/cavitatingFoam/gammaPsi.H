{
    gamma = max(min((rho - rholSat)/(rhovSat - rholSat), scalar(1)), scalar(0));

    Info<< "max-min gamma: " << max(gamma).value()
        << " " << min(gamma).value() << endl;

    psiModel->correct();

    //Info<< "min a: " << 1.0/sqrt(max(psi)).value() << endl;
}
