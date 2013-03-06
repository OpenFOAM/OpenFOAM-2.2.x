    #include "readTimeControls.H"

    label nAlphaCorr(readLabel(pimple.dict().lookup("nAlphaCorr")));

    label nAlphaSubCycles(readLabel(pimple.dict().lookup("nAlphaSubCycles")));

    if (nAlphaSubCycles > 1 && pimple.nCorrPIMPLE() != 1)
    {
        FatalErrorIn(args.executable())
            << "Sub-cycling alpha is only allowed for PISO, "
               "i.e. when the number of outer-correctors = 1"
            << exit(FatalError);
    }

    bool correctPhi =
        pimple.dict().lookupOrDefault<Switch>("correctPhi", true);

    bool checkMeshCourantNo =
        pimple.dict().lookupOrDefault<Switch>("checkMeshCourantNo", false);
