    #include "readTimeControls.H"

    const dictionary& pimpleDict = pimple.dict();

    const bool correctPhi =
        pimpleDict.lookupOrDefault("correctPhi", false);

    const bool checkMeshCourantNo =
        pimpleDict.lookupOrDefault("checkMeshCourantNo", false);

    const bool ddtPhiCorr =
        pimpleDict.lookupOrDefault("ddtPhiCorr", true);
