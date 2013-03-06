{
    labelList cellBoundaryFaceCount(epsilon.size(), 0);

    scalar Cmu25 = ::pow(Cmu.value(), 0.25);
    scalar Cmu75 = ::pow(Cmu.value(), 0.75);
    scalar kappa_ = kappa.value();

    const fvPatchList& patches = mesh.boundary();

    //- Initialise the near-wall P field to zero
    forAll(patches, patchi)
    {
        const fvPatch& currPatch = patches[patchi];

        if (isA<wallFvPatch>(currPatch))
        {
            forAll(currPatch, facei)
            {
                label faceCelli = currPatch.faceCells()[facei];

                epsilon[faceCelli] = 0.0;
                G[faceCelli] = 0.0;
            }
        }
    }

    //- Accumulate the wall face contributions to epsilon and G
    //  Increment cellBoundaryFaceCount for each face for averaging
    forAll(patches, patchi)
    {
        const fvPatch& currPatch = patches[patchi];
        const scalarField& mu2_ = thermo2.mu().boundaryField()[patchi];
        const scalarField& rho2_ = rho2.boundaryField()[patchi];

        if (isA<wallFvPatch>(currPatch))
        {
            const scalarField& nut2w = nut2.boundaryField()[patchi];

            scalarField magFaceGradU(mag(U2.boundaryField()[patchi].snGrad()));

            forAll(currPatch, facei)
            {
                label faceCelli = currPatch.faceCells()[facei];

                // For corner cells (with two boundary or more faces),
                // epsilon and G in the near-wall cell are calculated
                // as an average

                cellBoundaryFaceCount[faceCelli]++;

                epsilon[faceCelli] +=
                     Cmu75*::pow(k[faceCelli], 1.5)
                    /(kappa_*y[patchi][facei]);

                G[faceCelli] +=
                    (nut2w[facei] + mu2_[facei]/rho2_[facei])
                   *magFaceGradU[facei]
                   *Cmu25*::sqrt(k[faceCelli])
                   /(kappa_*y[patchi][facei]);
            }
        }
    }


    // perform the averaging

    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isA<wallFvPatch>(curPatch))
        {
            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];

                epsilon[faceCelli] /= cellBoundaryFaceCount[faceCelli];
                G[faceCelli] /= cellBoundaryFaceCount[faceCelli];
                cellBoundaryFaceCount[faceCelli] = 1;
            }
        }
    }
}
