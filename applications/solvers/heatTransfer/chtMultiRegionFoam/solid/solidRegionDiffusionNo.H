    scalar DiNum = -GREAT;

    forAll(solidRegions, i)
    {
        #include "setRegionSolidFields.H"

        DiNum = max
        (
            solidRegionDiffNo
            (
                solidRegions[i],
                runTime,
                rho*cp,
                kappa
            ),
            DiNum
        );
    }
