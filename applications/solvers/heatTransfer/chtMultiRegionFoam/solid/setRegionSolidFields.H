    fvMesh& mesh = solidRegions[i];
    solidThermo& thermo = thermos[i];
    const radiation::radiationModel& radiation = radiations[i];

    tmp<volScalarField> trho = thermo.rho();
    const volScalarField& rho = trho();

    tmp<volScalarField> tcp = thermo.Cp();
    const volScalarField& cp = tcp();

    tmp<volScalarField> talpha = thermo.alpha();
    const volScalarField& alpha = talpha();
    tmp<volScalarField> tkappa = thermo.kappa();
    const volScalarField& kappa = tkappa();

    volScalarField& h = thermo.he();

    const volScalarField& betav = betavSolid[i];

    fv::IOoptionList& fvOptions = solidHeatSources[i];
