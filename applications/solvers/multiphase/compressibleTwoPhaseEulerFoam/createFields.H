    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    word phase1Name
    (
        transportProperties.found("phases")
      ? wordList(transportProperties.lookup("phases"))[0]
      : "1"
    );

    word phase2Name
    (
        transportProperties.found("phases")
      ? wordList(transportProperties.lookup("phases"))[1]
      : "2"
    );

    autoPtr<phaseModel> phase1 = phaseModel::New
    (
        mesh,
        transportProperties,
        phase1Name
    );

    autoPtr<phaseModel> phase2 = phaseModel::New
    (
        mesh,
        transportProperties,
        phase2Name
    );

    volScalarField& alpha1 = phase1();
    volScalarField& alpha2 = phase2();
    alpha2 = scalar(1) - alpha1;

    volVectorField& U1 = phase1->U();
    surfaceScalarField& phi1 = phase1->phi();

    volVectorField& U2 = phase2->U();
    surfaceScalarField& phi2 = phase2->phi();

    dimensionedScalar pMin
    (
        "pMin",
        dimPressure,
        transportProperties.lookup("pMin")
    );

    rhoThermo& thermo1 = phase1->thermo();
    rhoThermo& thermo2 = phase2->thermo();

    volScalarField& p = thermo1.p();

    volScalarField& rho1 = thermo1.rho();
    const volScalarField& psi1 = thermo1.psi();

    volScalarField& rho2 = thermo2.rho();
    const volScalarField& psi2 = thermo2.psi();

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1*U1 + alpha2*U2
    );

    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(alpha1)*phi1 + fvc::interpolate(alpha2)*phi2
    );

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1*rho1 + alpha2*rho2
    );

    Info<< "Calculating field DDtU1 and DDtU2\n" << endl;

    volVectorField DDtU1
    (
        fvc::ddt(U1)
      + fvc::div(phi1, U1)
      - fvc::div(phi1)*U1
    );

    volVectorField DDtU2
    (
        fvc::ddt(U2)
      + fvc::div(phi2, U2)
      - fvc::div(phi2)*U2
    );


    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());

    dimensionedScalar Cvm
    (
        "Cvm",
        dimless,
        transportProperties.lookup("Cvm")
    );

    dimensionedScalar Cl
    (
        "Cl",
        dimless,
        transportProperties.lookup("Cl")
    );

    dimensionedScalar Ct
    (
        "Ct",
        dimless,
        transportProperties.lookup("Ct")
    );

    #include "createRASTurbulence.H"

    IOdictionary interfacialProperties
    (
        IOobject
        (
            "interfacialProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    autoPtr<dragModel> drag1 = dragModel::New
    (
        interfacialProperties,
        alpha1,
        phase1,
        phase2
    );

    autoPtr<dragModel> drag2 = dragModel::New
    (
        interfacialProperties,
        alpha2,
        phase2,
        phase1
    );

    autoPtr<heatTransferModel> heatTransfer1 = heatTransferModel::New
    (
        interfacialProperties,
        alpha1,
        phase1,
        phase2
    );

    autoPtr<heatTransferModel> heatTransfer2 = heatTransferModel::New
    (
        interfacialProperties,
        alpha2,
        phase2,
        phase1
    );

    word dispersedPhase(interfacialProperties.lookup("dispersedPhase"));

    if
    (
        !(
            dispersedPhase == phase1Name
         || dispersedPhase == phase2Name
         || dispersedPhase == "both"
        )
    )
    {
        FatalErrorIn(args.executable())
            << "invalid dispersedPhase " << dispersedPhase
            << exit(FatalError);
    }

    Info << "dispersedPhase is " << dispersedPhase << endl;

    scalar residualPhaseFraction
    (
        readScalar
        (
            interfacialProperties.lookup("residualPhaseFraction")
        )
    );

    dimensionedScalar residualSlip
    (
        "residualSlip",
        dimVelocity,
        interfacialProperties.lookup("residualSlip")
    );

    kineticTheoryModel kineticTheory
    (
        phase1,
        U2,
        alpha1,
        drag1
    );

    volScalarField rAU1
    (
        IOobject
        (
            "rAU" + phase1Name,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, 0, 1, 0, 0), 0.0)
    );

    surfaceScalarField ppMagf
    (
        IOobject
        (
            "ppMagf",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.0)
    );


    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);


    volScalarField dgdt
    (
        pos(alpha2)*fvc::div(phi)/max(alpha2, scalar(0.0001))
    );


    Info<< "Creating field dpdt\n" << endl;
    volScalarField dpdt
    (
        IOobject
        (
            "dpdt",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
    );


    Info<< "Creating field kinetic energy K\n" << endl;
    volScalarField K1("K" + phase1Name, 0.5*magSqr(U1));
    volScalarField K2("K" + phase2Name, 0.5*magSqr(U2));

    surfaceScalarField alphaPhi1
    (
        IOobject
        (
            "alphaPhi" + phase1Name,
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(0, 3, -1, 0, 0), 0)
    );

    surfaceScalarField alphaPhi2
    (
        IOobject
        (
            "alphaPhi" + phase2Name,
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(0, 3, -1, 0, 0), 0)
    );
