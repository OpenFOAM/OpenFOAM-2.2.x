Info<< "\nConstructing coal cloud" << endl;
coalCloud coalParcels
(
    "coalCloud1",
    rho,
    U,
    g,
    slgThermo
);

Info<< "\nConstructing limestone cloud" << endl;
fluidThermoCloud limestoneParcels
(
    "limestoneCloud1",
    rho,
    U,
    g,
    slgThermo
);
