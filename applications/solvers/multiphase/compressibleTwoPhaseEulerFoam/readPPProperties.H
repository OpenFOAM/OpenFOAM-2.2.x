    IOdictionary ppProperties
    (
        IOobject
        (
            "ppProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    scalar preAlphaExp
    (
        readScalar(ppProperties.lookup("preAlphaExp"))
    );

    scalar alphaMax
    (
        readScalar(ppProperties.lookup("alphaMax"))
    );

    scalar expMax
    (
        readScalar(ppProperties.lookup("expMax"))
    );

    dimensionedScalar g0
    (
        ppProperties.lookup("g0")
    );

    Switch packingLimiter
    (
        ppProperties.lookup("packingLimiter")
    );
