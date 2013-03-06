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

    autoPtr<phaseModel> phase1 = phaseModel::New
    (
        mesh,
        transportProperties,
        "1"
    );

    autoPtr<phaseModel> phase2 = phaseModel::New
    (
        mesh,
        transportProperties,
        "2"
    );

    volVectorField& U1 = phase1->U();
    surfaceScalarField& phi1 = phase1->phi();
    const dimensionedScalar& rho1 = phase1->rho();
    const dimensionedScalar& nu1 = phase1->nu();

    volVectorField& U2 = phase2->U();
    surfaceScalarField& phi2 = phase2->phi();
    const dimensionedScalar& rho2 = phase2->rho();
    const dimensionedScalar& nu2 = phase2->nu();

    Info<< "Reading field alpha1\n" << endl;
    volScalarField alpha1
    (
        IOobject
        (
            "alpha1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField alpha2
    (
        IOobject
        (
            "alpha2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        scalar(1) - alpha1
        //,alpha1.boundaryField().types()
    );

    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

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
            mesh
        ),
        alpha1*rho1 + alpha2*rho2
    );

    #include "createRASTurbulence.H"

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

    word dragPhase("blended");
    if (interfacialProperties.found("dragPhase"))
    {
        dragPhase = word(interfacialProperties.lookup("dragPhase"));

        bool validDrag =
            dragPhase == "1" || dragPhase == "2" || dragPhase == "blended";

        if (!validDrag)
        {
            FatalErrorIn(args.executable())
                << "invalid dragPhase " << dragPhase
                << exit(FatalError);
        }
    }

    dimensionedScalar residualSlip
    (
        dimensionedScalar::lookupOrDefault
        (
            "residualSlip",
            interfacialProperties,
            0,
            dimVelocity
        )
    );

    Info << "dragPhase is " << dragPhase << endl;
    kineticTheoryModel kineticTheory
    (
        phase1,
        U2,
        alpha1,
        drag1
    );

    surfaceScalarField rAU1f
    (
        IOobject
        (
            "rAU1f",
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
