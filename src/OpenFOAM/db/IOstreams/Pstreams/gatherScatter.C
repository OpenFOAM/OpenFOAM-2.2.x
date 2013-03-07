/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Description
    Gather data from all processors onto single processor according to some
    communication schedule (usually linear-to-master or tree-to-master).
    The gathered data will be a single value constructed from the values
    on individual processors using a user-specified operator.

\*---------------------------------------------------------------------------*/

#include "UOPstream.H"
#include "OPstream.H"
#include "UIPstream.H"
#include "IPstream.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, class BinaryOp>
void Pstream::gather
(
    const List<UPstream::commsStruct>& comms,
    T& Value,
    const BinaryOp& bop,
    const int tag
)
{
    if (UPstream::parRun())
    {
        // Get my communication order
        const commsStruct& myComm = comms[UPstream::myProcNo()];

        // Receive from my downstairs neighbours
        forAll(myComm.below(), belowI)
        {
            T value;

            if (contiguous<T>())
            {
                UIPstream::read
                (
                    UPstream::scheduled,
                    myComm.below()[belowI],
                    reinterpret_cast<char*>(&value),
                    sizeof(T),
                    tag
                );
            }
            else
            {
                IPstream fromBelow
                (
                    UPstream::scheduled,
                    myComm.below()[belowI],
                    0,
                    tag
                );
                fromBelow >> value;
            }

            Value = bop(Value, value);
        }

        // Send up Value
        if (myComm.above() != -1)
        {
            if (contiguous<T>())
            {
                UOPstream::write
                (
                    UPstream::scheduled,
                    myComm.above(),
                    reinterpret_cast<const char*>(&Value),
                    sizeof(T),
                    tag
                );
            }
            else
            {
                OPstream toAbove(UPstream::scheduled, myComm.above(), 0, tag);
                toAbove << Value;
            }
        }
    }
}


template<class T, class BinaryOp>
void Pstream::gather(T& Value, const BinaryOp& bop, const int tag)
{
    if (UPstream::nProcs() < UPstream::nProcsSimpleSum)
    {
        gather(UPstream::linearCommunication(), Value, bop, tag);
    }
    else
    {
        gather(UPstream::treeCommunication(), Value, bop, tag);
    }
}


template<class T>
void Pstream::scatter
(
    const List<UPstream::commsStruct>& comms,
    T& Value,
    const int tag
)
{
    if (UPstream::parRun())
    {
        // Get my communication order
        const commsStruct& myComm = comms[UPstream::myProcNo()];

        // Reveive from up
        if (myComm.above() != -1)
        {
            if (contiguous<T>())
            {
                UIPstream::read
                (
                    UPstream::scheduled,
                    myComm.above(),
                    reinterpret_cast<char*>(&Value),
                    sizeof(T),
                    tag
                );
            }
            else
            {
                IPstream fromAbove(UPstream::scheduled, myComm.above(), 0, tag);
                fromAbove >> Value;
            }
        }

        // Send to my downstairs neighbours
        forAll(myComm.below(), belowI)
        {
            if (contiguous<T>())
            {
                UOPstream::write
                (
                    UPstream::scheduled,
                    myComm.below()[belowI],
                    reinterpret_cast<const char*>(&Value),
                    sizeof(T),
                    tag
                );
            }
            else
            {
                OPstream toBelow
                (
                    UPstream::scheduled,
                    myComm.below()[belowI],
                    0,
                    tag
                );
                toBelow << Value;
            }
        }
    }
}


template<class T>
void Pstream::scatter(T& Value, const int tag)
{
    if (UPstream::nProcs() < UPstream::nProcsSimpleSum)
    {
        scatter(UPstream::linearCommunication(), Value, tag);
    }
    else
    {
        scatter(UPstream::treeCommunication(), Value, tag);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
