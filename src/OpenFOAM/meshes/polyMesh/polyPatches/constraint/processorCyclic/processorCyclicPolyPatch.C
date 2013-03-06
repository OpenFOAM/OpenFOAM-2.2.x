/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "processorCyclicPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "SubField.H"
#include "cyclicPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(processorCyclicPolyPatch, 0);
    addToRunTimeSelectionTable(polyPatch, processorCyclicPolyPatch, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::processorCyclicPolyPatch::processorCyclicPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const int myProcNo,
    const int neighbProcNo,
    const word& referPatchName,
    const transformType transform
)
:
    processorPolyPatch
    (
        name,
        size,
        start,
        index,
        bm,
        myProcNo,
        neighbProcNo,
        transform
    ),
    tag_
    (
        Pstream::nProcs()*max(myProcNo, neighbProcNo)
      + min(myProcNo, neighbProcNo)
    ),
    referPatchName_(referPatchName),
    referPatchID_(-1)
{
    if (debug)
    {
        Pout<< "processorCyclicPolyPatch " << name << " uses tag " << tag_
            << endl;
    }
}


Foam::processorCyclicPolyPatch::processorCyclicPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    processorPolyPatch(name, dict, index, bm, patchType),
    tag_
    (
        Pstream::nProcs()*max(myProcNo(), neighbProcNo())
      + min(myProcNo(), neighbProcNo())
    ),
    referPatchName_(dict.lookup("referPatch")),
    referPatchID_(-1)
{
    if (debug)
    {
        Pout<< "processorCyclicPolyPatch " << name << " uses tag " << tag_
            << endl;
    }
}


Foam::processorCyclicPolyPatch::processorCyclicPolyPatch
(
    const processorCyclicPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    processorPolyPatch(pp, bm),
    tag_(pp.tag_),
    referPatchName_(pp.referPatchName()),
    referPatchID_(-1)
{}


Foam::processorCyclicPolyPatch::processorCyclicPolyPatch
(
    const processorCyclicPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart,
    const word& referPatchName
)
:
    processorPolyPatch(pp, bm, index, newSize, newStart),
    tag_(pp.tag_),
    referPatchName_(referPatchName),
    referPatchID_(-1)
{}


Foam::processorCyclicPolyPatch::processorCyclicPolyPatch
(
    const processorCyclicPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    processorPolyPatch(pp, bm, index, mapAddressing, newStart),
    tag_(pp.tag_),
    referPatchName_(pp.referPatchName()),
    referPatchID_(-1)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::processorCyclicPolyPatch::~processorCyclicPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::processorCyclicPolyPatch::initGeometry(PstreamBuffers& pBufs)
{
    // Send over processorPolyPatch data
    processorPolyPatch::initGeometry(pBufs);
}


void Foam::processorCyclicPolyPatch::calcGeometry(PstreamBuffers& pBufs)
{
    // Receive and initialise processorPolyPatch data
    processorPolyPatch::calcGeometry(pBufs);

    if (Pstream::parRun())
    {

        // Where do we store the calculated transformation?
        // - on the processor patch?
        // - on the underlying cyclic patch?
        // - or do we not auto-calculate the transformation but
        //   have option of reading it.

        // Update underlying cyclic halves. Need to do both since only one
        // half might be present as a processorCyclic.
        coupledPolyPatch& pp = const_cast<coupledPolyPatch&>(referPatch());
        pp.calcGeometry
        (
            *this,
            faceCentres(),
            faceAreas(),
            faceCellCentres(),
            neighbFaceCentres(),
            neighbFaceAreas(),
            neighbFaceCellCentres()
        );

        if (isA<cyclicPolyPatch>(pp))
        {
            const cyclicPolyPatch& cpp = refCast<const cyclicPolyPatch>(pp);
            const_cast<cyclicPolyPatch&>(cpp.neighbPatch()).calcGeometry
            (
                *this,
                neighbFaceCentres(),
                neighbFaceAreas(),
                neighbFaceCellCentres(),
                faceCentres(),
                faceAreas(),
                faceCellCentres()
            );
        }
    }
}


void Foam::processorCyclicPolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    // Recalculate geometry
    initGeometry(pBufs);
}


void Foam::processorCyclicPolyPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField&
)
{
    calcGeometry(pBufs);
}


void Foam::processorCyclicPolyPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    processorPolyPatch::initUpdateMesh(pBufs);
}


void Foam::processorCyclicPolyPatch::updateMesh(PstreamBuffers& pBufs)
{
     referPatchID_ = -1;
     processorPolyPatch::updateMesh(pBufs);
}


void Foam::processorCyclicPolyPatch::initOrder
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp
) const
{
    // For now use the same algorithm as processorPolyPatch
    processorPolyPatch::initOrder(pBufs, pp);
}


// Return new ordering. Ordering is -faceMap: for every face index
// the new face -rotation:for every new face the clockwise shift
// of the original face. Return false if nothing changes (faceMap
// is identity, rotation is 0)
bool Foam::processorCyclicPolyPatch::order
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    // For now use the same algorithm as processorPolyPatch
    return processorPolyPatch::order(pBufs, pp, faceMap, rotation);
}


void Foam::processorCyclicPolyPatch::write(Ostream& os) const
{
    processorPolyPatch::write(os);
    os.writeKeyword("referPatch") << referPatchName_
        << token::END_STATEMENT << nl;
}


// ************************************************************************* //
