{
    scalar Cmu25 = ::pow(Cmu.value(), 0.25);
    scalar kappa_ = kappa.value();
    scalar E_ = E.value();

    const fvPatchList& patches = mesh.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& currPatch = patches[patchi];
        const scalarField& mu2_ = thermo2.mu().boundaryField()[patchi];
        const scalarField& rho2_ = rho2.boundaryField()[patchi];

        if (isA<wallFvPatch>(currPatch))
        {
            scalarField& nutw = nut2.boundaryField()[patchi];

            forAll(currPatch, facei)
            {
                label faceCelli = currPatch.faceCells()[facei];

                // calculate yPlus
                scalar yPlus =
                    Cmu25*y[patchi][facei]*::sqrt(k[faceCelli])
                   /(mu2_[facei]/rho2_[facei]);

                if (yPlus > 11.6)
                {
                    nutw[facei] =
                        (mu2_[facei]/rho2_[facei])
                       *(yPlus*kappa_/::log(E_*yPlus) - 1);
                }
                else
                {
                    nutw[facei] = 0.0;
                }
            }
        }
    }
}
