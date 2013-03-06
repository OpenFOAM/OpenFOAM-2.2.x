    #include "readTimeControls.H"

    int nAlphaCorr(readInt(pimple.dict().lookup("nAlphaCorr")));
    int nAlphaSubCycles(readInt(pimple.dict().lookup("nAlphaSubCycles")));
    Switch correctAlpha(pimple.dict().lookup("correctAlpha"));
