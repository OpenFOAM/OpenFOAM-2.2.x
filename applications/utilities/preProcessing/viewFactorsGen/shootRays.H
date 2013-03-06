// All rays expressed as start face (local) index end end face (global)
// Pre-size by assuming a certain percentage is visible.

// Maximum lenght for dynamicList
const label maxDynListLength = 10000;

for (label procI = 0; procI < Pstream::nProcs(); procI++)
{
    // Shoot rays from me to procI. Note that even if processor has
    // 0 faces we still need to call findLine to keep calls synced.

    DynamicField<point> start(coarseMesh.nFaces());
    DynamicField<point> end(start.size());
    DynamicList<label> startIndex(start.size());
    DynamicList<label> endIndex(start.size());


    const pointField& myFc = remoteCoarseCf[Pstream::myProcNo()];
    const vectorField& myArea = remoteCoarseSf[Pstream::myProcNo()];

    const pointField& remoteArea = remoteCoarseSf[procI];
    const pointField& remoteFc = remoteCoarseCf[procI];

    label i = 0;
    label j = 0;
    do
    {
        for (; i < myFc.size(); i++)
        {
            const point& fc = myFc[i];
            const vector& fA = myArea[i];

            for (; j < remoteFc.size(); j++)//
            {
                if (procI != Pstream::myProcNo() || i != j)
                {
                    const point& remFc = remoteFc[j];
                    const vector& remA = remoteArea[j];
                    const vector& d = remFc - fc;

                    if (((d & fA) < 0.) && ((d & remA) > 0))
                    {
                        start.append(fc + 0.0001*d);
                        startIndex.append(i);
                        end.append(fc + 0.9999*d);
                        label globalI = globalNumbering.toGlobal(procI, j);
                        endIndex.append(globalI);
                        if (startIndex.size() > maxDynListLength)
                        {
                            break;
                        }
                    }
                }
            }
            if (startIndex.size() > maxDynListLength)
            {
                break;
            }

            if (j == remoteFc.size())
            {
                j = 0;
            }
        }

        List<pointIndexHit> hitInfo(startIndex.size());
        surfacesMesh.findLine(start, end, hitInfo);

        forAll (hitInfo, rayI)
        {
            if (!hitInfo[rayI].hit())
            {
                rayStartFace.append(startIndex[rayI]);
                rayEndFace.append(endIndex[rayI]);
            }
        }

        start.clear();
        startIndex.clear();
        end.clear();
        endIndex.clear();

    }while (returnReduce(i < myFc.size(), orOp<bool>()));

}
