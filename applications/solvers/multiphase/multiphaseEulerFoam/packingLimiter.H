    if (packingLimiter)
    {
        // Calculating exceeding volume fractions
        volScalarField alphaEx(max(alpha1 - alphaMax, scalar(0)));

        // Finding neighbouring cells of the whole domain
        labelListList neighbour = mesh.cellCells();
        scalarField cellVolumes(mesh.cellVolumes());

        forAll (alphaEx, celli)
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
                    alphaEx[neighbourCell[cellj]]*cellVolumes[celli]
                   /neighboursNeighbourCellVolumes;
            }

            alpha1[celli] += neighboursEx - alphaEx[celli];
        }
    }
