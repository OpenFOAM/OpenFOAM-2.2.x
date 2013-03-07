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

InClass
    Foam

Description
    Combination-Reduction operation for a parallel run.  The
    information from all nodes is collected on the master node,
    combined using the given combination function and the result is
    broadcast to all nodes


\*---------------------------------------------------------------------------*/

#ifndef PstreamCombineReduceOps_H
#define PstreamCombineReduceOps_H

#include "UPstream.H"
#include "Pstream.H"
#include "ops.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T, class CombineOp>
void combineReduce
(
    const List<UPstream::commsStruct>& comms,
    T& Value,
    const CombineOp& cop,
    const int tag
)
{
    Pstream::combineGather(comms, Value, cop, tag);
    Pstream::combineScatter(comms, Value, tag);
}


template<class T, class CombineOp>
void combineReduce
(
    T& Value,
    const CombineOp& cop,
    const int tag = Pstream::msgType()
)
{
    if (UPstream::nProcs() < UPstream::nProcsSimpleSum)
    {
        Pstream::combineGather
        (
            UPstream::linearCommunication(),
            Value,
            cop,
            tag
        );
        Pstream::combineScatter(UPstream::linearCommunication(), Value, tag);
    }
    else
    {
        Pstream::combineGather(UPstream::treeCommunication(), Value, cop, tag);
        Pstream::combineScatter(UPstream::treeCommunication(), Value, tag);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
