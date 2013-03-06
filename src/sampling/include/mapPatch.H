    if (planeMapping && targetPatchNamePresent && cut.cut())
    {
        //Map variables U, phi & k
        //Info<< typeInfo("isoLESmodel") << endl;
        //Info<< "SGS model" << sgsModel.type() << endl;

        scalarField planek(cutCells.size());
        scalarField planeNuTilda(cutCells.size());
        vectorField planeU(cutCells.size());

        forAll(cutCells, cCellsI)
        {
            if (sgsModel.type() == "SpalartAllmaras")
            {
                planeNuTilda[cCellsI] = sgsModel.nuTilda()()[cutCells[cCellsI]];
            }
            else
            {
                planek[cCellsI] = sgsModel.k()()[cutCells[cCellsI]];
            }
            planeU[cCellsI] = U[cutCells[cCellsI]];
        }

        if (sgsModel.type() == "SpalartAllmaras")
        {
            Info<< "Mapping NuTilda." << endl;
            sgsModel.nuTilda()().boundaryField()[targetPatchNumber] ==
            interPatch.faceInterpolate(planeNuTilda)();
        }
        else
        {
            Info<< "Mapping k." << endl;
            sgsModel.k()().boundaryField()[targetPatchNumber] ==
            interPatch.faceInterpolate(planek)();
        }

        U.boundaryField()[targetPatchNumber] ==
            interPatch.faceInterpolate(planeU)();

/*         sgsModel.k()().boundaryField()[targetPatchNumber] ==
         toPatchInter.pointToFaceInterpolate
         (
             interPatch.pointInterpolate
             (
                 planePatchInter.faceToPointInterpolate
                 (
                     planek
                 )
             )
         );

         U.boundaryField()[targetPatchNumber] ==
         toPatchInter.pointToFaceInterpolate
         (
             interPatch.pointInterpolate
             (
                 planePatchInter.faceToPointInterpolate
                 (
                     planeU
                 )
             )
         );

*/

        scalar Q = sum
        (
            mesh.Sf().boundaryField()[targetPatchNumber]
            & U.boundaryField()[targetPatchNumber]
        );

        scalar Qbar = sum
        (
            mesh.Sf().boundaryField()[targetPatchNumber] &
            Ubar.value()
        );

        U.boundaryField()[targetPatchNumber] ==
            U.boundaryField()[targetPatchNumber] * (Qbar/Q);

        Info<< "Mass flux correction factor: " << (Qbar/Q) << endl;

        phi.boundaryField()[targetPatchNumber] ==
        (
            mesh.Sf().boundaryField()[targetPatchNumber] &
            U.boundaryField()[targetPatchNumber]
        );
    }


