/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "domainDecomposition.H"
#include "decompositionMethod.H"
#include "cpuTime.H"
#include "cellSet.H"
#include "regionSplit.H"
#include "Tuple2.H"
#include "faceSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::domainDecomposition::distributeCells()
{
    Info<< "\nCalculating distribution of cells" << endl;

    cpuTime decompositionTime;


    // See if any faces need to have owner and neighbour on same processor
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelHashSet sameProcFaces;

    if (decompositionDict_.found("preservePatches"))
    {
        wordList pNames(decompositionDict_.lookup("preservePatches"));

        Info<< nl
            << "Keeping owner of faces in patches " << pNames
            << " on same processor. This only makes sense for cyclics." << endl;

        const polyBoundaryMesh& patches = boundaryMesh();

        forAll(pNames, i)
        {
            const label patchI = patches.findPatchID(pNames[i]);

            if (patchI == -1)
            {
                FatalErrorIn("domainDecomposition::distributeCells()")
                    << "Unknown preservePatch " << pNames[i]
                    << endl << "Valid patches are " << patches.names()
                    << exit(FatalError);
            }

            const polyPatch& pp = patches[patchI];

            forAll(pp, i)
            {
                sameProcFaces.insert(pp.start() + i);
            }
        }
    }
    if (decompositionDict_.found("preserveFaceZones"))
    {
        wordList zNames(decompositionDict_.lookup("preserveFaceZones"));

        Info<< nl
            << "Keeping owner and neighbour of faces in zones " << zNames
            << " on same processor" << endl;

        const faceZoneMesh& fZones = faceZones();

        forAll(zNames, i)
        {
            label zoneI = fZones.findZoneID(zNames[i]);

            if (zoneI == -1)
            {
                FatalErrorIn("domainDecomposition::distributeCells()")
                    << "Unknown preserveFaceZone " << zNames[i]
                    << endl << "Valid faceZones are " << fZones.names()
                    << exit(FatalError);
            }

            const faceZone& fz = fZones[zoneI];

            forAll(fz, i)
            {
                sameProcFaces.insert(fz[i]);
            }
        }
    }


    // Specified processor for owner and neighbour of faces
    Map<label> specifiedProcessorFaces;
    List<Tuple2<word, label> > zNameAndProcs;

    if (decompositionDict_.found("singleProcessorFaceSets"))
    {
        decompositionDict_.lookup("singleProcessorFaceSets") >> zNameAndProcs;

        label nCells = 0;

        Info<< endl;

        forAll(zNameAndProcs, i)
        {
            Info<< "Keeping all cells connected to faceSet "
                << zNameAndProcs[i].first()
                << " on processor " << zNameAndProcs[i].second() << endl;

            // Read faceSet
            faceSet fz(*this, zNameAndProcs[i].first());
            nCells += fz.size();
        }


        // Size
        specifiedProcessorFaces.resize(2*nCells);


        // Fill
        forAll(zNameAndProcs, i)
        {
            faceSet fz(*this, zNameAndProcs[i].first());

            label procI = zNameAndProcs[i].second();

            forAllConstIter(faceSet, fz, iter)
            {
                label faceI = iter.key();

                specifiedProcessorFaces.insert(faceI, procI);
            }
        }
    }


    // Construct decomposition method and either do decomposition on
    // cell centres or on agglomeration


    autoPtr<decompositionMethod> decomposePtr = decompositionMethod::New
    (
        decompositionDict_
    );


    if (sameProcFaces.empty() && specifiedProcessorFaces.empty())
    {
        if (decompositionDict_.found("weightField"))
        {
            word weightName = decompositionDict_.lookup("weightField");

            volScalarField weights
            (
                IOobject
                (
                    weightName,
                    time().timeName(),
                    *this,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                *this
            );

            cellToProc_ = decomposePtr().decompose
            (
                *this,
                cellCentres(),
                weights.internalField()
            );
        }
        else
        {
            cellToProc_ = decomposePtr().decompose(*this, cellCentres());
        }

    }
    else
    {
        Info<< "Constrained decomposition:" << endl
            << "    faces with same processor owner and neighbour : "
            << sameProcFaces.size() << endl
            << "    faces all on same processor                   : "
            << specifiedProcessorFaces.size() << endl << endl;

        // Faces where owner and neighbour are not 'connected' (= all except
        // sameProcFaces)
        boolList blockedFace(nFaces(), true);

        forAllConstIter(labelHashSet, sameProcFaces, iter)
        {
            blockedFace[iter.key()] = false;
        }


        // For specifiedProcessorFaces add all point connected faces
        {
            forAllConstIter(Map<label>, specifiedProcessorFaces, iter)
            {
                const face& f = faces()[iter.key()];
                forAll(f, fp)
                {
                    const labelList& pFaces = pointFaces()[f[fp]];
                    forAll(pFaces, i)
                    {
                        blockedFace[pFaces[i]] = false;
                    }
                }
            }
        }


        // Connect coupled boundary faces
        const polyBoundaryMesh& patches =  boundaryMesh();

        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            if (pp.coupled())
            {
                forAll(pp, i)
                {
                    blockedFace[pp.start()+i] = false;
                }
            }
        }

        // Determine global regions, separated by blockedFaces
        regionSplit globalRegion(*this, blockedFace);


        // Determine region cell centres
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // This just takes the first cell in the region. Otherwise the problem
        // is with cyclics - if we'd average the region centre might be
        // somewhere in the middle of the domain which might not be anywhere
        // near any of the cells.

        pointField regionCentres(globalRegion.nRegions(), point::max);

        forAll(globalRegion, cellI)
        {
            label regionI = globalRegion[cellI];

            if (regionCentres[regionI] == point::max)
            {
                regionCentres[regionI] = cellCentres()[cellI];
            }
        }

        // Do decomposition on agglomeration
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        scalarField regionWeights(globalRegion.nRegions(), 0);

        if (decompositionDict_.found("weightField"))
        {
            word weightName = decompositionDict_.lookup("weightField");

            volScalarField weights
            (
                IOobject
                (
                    weightName,
                    time().timeName(),
                    *this,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                *this
            );

            forAll(globalRegion, cellI)
            {
                label regionI = globalRegion[cellI];

                regionWeights[regionI] += weights.internalField()[cellI];
            }
        }
        else
        {
            forAll(globalRegion, cellI)
            {
                label regionI = globalRegion[cellI];

                regionWeights[regionI] += 1.0;
            }
        }

        cellToProc_ = decomposePtr().decompose
        (
            *this,
            globalRegion,
            regionCentres,
            regionWeights
        );


        // For specifiedProcessorFaces rework the cellToProc to enforce
        // all on one processor since we can't guarantee that the input
        // to regionSplit was a single region.
        // E.g. faceSet 'a' with the cells split into two regions
        // by a notch formed by two walls
        //
        //          \   /
        //           \ /
        //    ---a----+-----a-----
        //
        //
        // Note that reworking the cellToProc might make the decomposition
        // unbalanced.
        if (specifiedProcessorFaces.size())
        {
            forAll(zNameAndProcs, i)
            {
                faceSet fz(*this, zNameAndProcs[i].first());

                if (fz.size())
                {
                    label procI = zNameAndProcs[i].second();
                    if (procI == -1)
                    {
                        // If no processor specified use the one from the
                        // 0th element
                        procI = cellToProc_[faceOwner()[fz[0]]];
                    }

                    forAllConstIter(faceSet, fz, iter)
                    {
                        label faceI = iter.key();

                        cellToProc_[faceOwner()[faceI]] = procI;
                        if (isInternalFace(faceI))
                        {
                            cellToProc_[faceNeighbour()[faceI]] = procI;
                        }
                    }
                }
            }
        }
    }

    Info<< "\nFinished decomposition in "
        << decompositionTime.elapsedCpuTime()
        << " s" << endl;
}


// ************************************************************************* //
