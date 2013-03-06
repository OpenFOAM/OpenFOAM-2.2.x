    if (packingLimiter)
    {
        // Calculating exceeding volume fractions
        volScalarField alpha1Ex(max(alpha1 - alphaMax, scalar(0)));

        // Finding neighbouring cells of the whole domain
        labelListList neighbour = mesh.cellCells();
        scalarField cellVolumes(mesh.cellVolumes());

        forAll (alpha1Ex, celli)
        {
            // Finding the labels of the neighbouring cells
            labelList neighbourCell = neighbour[celli];

            // Initializing neighbouring cells contribution
            scalar neighboursEx = 0.0;

            forAll (neighbourCell, cellj)
            {
                labelList neighboursNeighbour = neighbour[neighbourCell[cellj]];
                scalar neighboursNeighbourCellVolumes = 0.0;

                forAll (neighboursNeighbour, cellk)
                {
                    neighboursNeighbourCellVolumes +=
                        cellVolumes[neighboursNeighbour[cellk]];
                }

                neighboursEx +=
                    alpha1Ex[neighbourCell[cellj]]*cellVolumes[celli]
                   /neighboursNeighbourCellVolumes;
            }

            alpha1[celli] += neighboursEx - alpha1Ex[celli];
        }
    }
