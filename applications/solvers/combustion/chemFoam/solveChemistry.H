    dtChem = chemistry.solve
    (
        runTime.value() - runTime.deltaT().value(),
        runTime.deltaT().value()
    );
    scalar Sh = chemistry.Sh()()[0]/rho[0];
    integratedHeat += Sh*runTime.deltaT().value();
