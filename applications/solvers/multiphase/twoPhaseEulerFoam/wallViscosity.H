{
    scalar Cmu25 = ::pow(Cmu.value(), 0.25);
    scalar kappa_ = kappa.value();
    scalar E_ = E.value();
    scalar nu2_ = nu2.value();

    const fvPatchList& patches = mesh.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& currPatch = patches[patchi];

        if (isA<wallFvPatch>(currPatch))
        {
            scalarField& nutw = nut2.boundaryField()[patchi];

            forAll(currPatch, facei)
            {
                label faceCelli = currPatch.faceCells()[facei];

                // calculate yPlus
                scalar yPlus =
                    Cmu25*y[patchi][facei]*::sqrt(k[faceCelli])/nu2_;

                if (yPlus > 11.6)
                {
                    nutw[facei] = nu2_*(yPlus*kappa_/::log(E_*yPlus) -1);
                }
                else
                {
                    nutw[facei] = 0.0;
                }
            }
        }
    }
}
