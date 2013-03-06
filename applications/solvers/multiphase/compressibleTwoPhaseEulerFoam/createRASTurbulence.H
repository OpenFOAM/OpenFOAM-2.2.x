    IOdictionary RASProperties
    (
        IOobject
        (
            "RASProperties",
            runTime.constant(),
            mesh,
            //IOobject::MUST_READ_IF_MODIFIED,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );


    Switch turbulence
    (
        RASProperties.lookup("turbulence")
    );

    dictionary kEpsilonDict
    (
        RASProperties.subDictPtr("kEpsilonCoeffs")
    );

    dimensionedScalar Cmu
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Cmu",
            kEpsilonDict,
            0.09
        )
    );

    dimensionedScalar C1
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "C1",
            kEpsilonDict,
            1.44
        )
    );

    dimensionedScalar C2
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "C2",
            kEpsilonDict,
            1.92
        )
    );

    dimensionedScalar alpha1k
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "alpha1k",
            kEpsilonDict,
            1.0
        )
    );

    dimensionedScalar alpha1Eps
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "alpha1Eps",
            kEpsilonDict,
            0.76923
        )
    );

    dimensionedScalar Prt
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Prt",
            kEpsilonDict,
            1.0
        )
    );

    dictionary wallFunctionDict
    (
        RASProperties.subDictPtr("wallFunctionCoeffs")
    );

    dimensionedScalar kappa
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "kappa",
            wallFunctionDict,
            0.41
        )
    );

    dimensionedScalar E
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "E",
            wallFunctionDict,
            9.8
        )
    );

    if (RASProperties.lookupOrDefault("printCoeffs", false))
    {
        Info<< "kEpsilonCoeffs" << kEpsilonDict << nl
            << "wallFunctionCoeffs" << wallFunctionDict << endl;
    }


    nearWallDist y(mesh);


    Info<< "Reading field k\n" << endl;
    volScalarField k
    (
        IOobject
        (
            "k",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field epsilon\n" << endl;
    volScalarField epsilon
    (
        IOobject
        (
            "epsilon",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Calculating field nut2\n" << endl;
    volScalarField nut2
    (
        IOobject
        (
            "nut" + phase2Name,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Cmu*sqr(k)/epsilon
    );

    Info<< "Calculating field nuEff1\n" << endl;
    volScalarField nuEff1
    (
        IOobject
        (
            "nuEff" + phase1Name,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sqr(Ct)*nut2 + thermo1.mu()/rho1
    );

    Info<< "Calculating field nuEff2\n" << endl;
    volScalarField nuEff2
    (
        IOobject
        (
            "nuEff" + phase2Name,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        nut2 + thermo2.mu()/rho2
    );
