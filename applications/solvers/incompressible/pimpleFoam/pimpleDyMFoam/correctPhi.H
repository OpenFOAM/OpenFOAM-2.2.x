{
    if (mesh.changing())
    {
        forAll(U.boundaryField(), patchI)
        {
            if (U.boundaryField()[patchI].fixesValue())
            {
                U.boundaryField()[patchI].initEvaluate();
            }
        }

        forAll(U.boundaryField(), patchI)
        {
            if (U.boundaryField()[patchI].fixesValue())
            {
                U.boundaryField()[patchI].evaluate();

                phi.boundaryField()[patchI] =
                    U.boundaryField()[patchI]
                  & mesh.Sf().boundaryField()[patchI];
            }
        }
    }

    wordList pcorrTypes
    (
        p.boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
    );

    forAll(p.boundaryField(), patchI)
    {
        if (p.boundaryField()[patchI].fixesValue())
        {
            pcorrTypes[patchI] = fixedValueFvPatchScalarField::typeName;
        }
    }

    volScalarField pcorr
    (
        IOobject
        (
            "pcorr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("pcorr", p.dimensions(), 0.0),
        pcorrTypes
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pcorrEqn
        (
            fvm::laplacian(rAU, pcorr) == fvc::div(phi)
        );

        pcorrEqn.setReference(pRefCell, pRefValue);
        pcorrEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi -= pcorrEqn.flux();
        }
    }
}

phi.oldTime() = phi;

#include "continuityErrs.H"
