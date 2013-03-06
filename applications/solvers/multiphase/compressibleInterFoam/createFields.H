    Info<< "Reading field p_rgh\n" << endl;
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "createPhi.H"

    Info<< "Reading field T\n" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        p_rgh
    );


    Info<< "Reading transportProperties\n" << endl;
    twoPhaseMixture twoPhaseProperties(U, phi);

    volScalarField& alpha1(twoPhaseProperties.alpha1());

    Info<< "Calculating phase-fraction alpha" << twoPhaseProperties.phase2Name()
        << nl << endl;
    volScalarField alpha2
    (
        "alpha" + twoPhaseProperties.phase2Name(),
        scalar(1) - alpha1
    );

    dimensionedScalar k1
    (
        "k",
        dimensionSet(1, 1, -3, -1, 0),
        twoPhaseProperties.subDict
        (
            twoPhaseProperties.phase1Name()
        ).lookup("k")
    );

    dimensionedScalar k2
    (
        "k",
        dimensionSet(1, 1, -3, -1, 0),
        twoPhaseProperties.subDict
        (
            twoPhaseProperties.phase2Name()
        ).lookup("k")
    );

    dimensionedScalar Cv1
    (
        "Cv",
        dimensionSet(0, 2, -2, -1, 0),
        twoPhaseProperties.subDict
        (
            twoPhaseProperties.phase1Name()
        ).lookup("Cv")
    );

    dimensionedScalar Cv2
    (
        "Cv",
        dimensionSet(0, 2, -2, -1, 0),
        twoPhaseProperties.subDict
        (
            twoPhaseProperties.phase2Name()
        ).lookup("Cv")
    );

    autoPtr<phaseEquationOfState> eos1
    (
        phaseEquationOfState::New
        (
            twoPhaseProperties.subDict
            (
                twoPhaseProperties.phase1Name()
            )
        )
    );

    autoPtr<phaseEquationOfState> eos2
    (
        phaseEquationOfState::New
        (
            twoPhaseProperties.subDict
            (
                twoPhaseProperties.phase2Name()
            )
        )
    );

    volScalarField psi1
    (
        IOobject
        (
            "psi1",
            runTime.timeName(),
            mesh
        ),
        eos1->psi(p, T)
    );
    psi1.oldTime();

    volScalarField psi2
    (
        IOobject
        (
            "psi2",
            runTime.timeName(),
            mesh
        ),
        eos2->psi(p, T)
    );
    psi2.oldTime();

    dimensionedScalar pMin(twoPhaseProperties.lookup("pMin"));

    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("ghf", g & mesh.Cf());

    volScalarField rho1("rho1", eos1->rho(p, T));
    volScalarField rho2("rho2", eos2->rho(p, T));

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        alpha1*rho1 + alpha2*rho2
    );


    // Mass flux
    // Initialisation does not matter because rhoPhi is reset after the
    // alpha1 solution before it is used in the U equation.
    surfaceScalarField rhoPhi
    (
        IOobject
        (
            "rho*phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(rho)*phi
    );

    volScalarField dgdt
    (
        pos(alpha2)*fvc::div(phi)/max(alpha2, scalar(0.0001))
    );

    // Construct interface from alpha1 distribution
    interfaceProperties interface(alpha1, U, twoPhaseProperties);

    // Construct incompressible turbulence model
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, twoPhaseProperties)
    );
