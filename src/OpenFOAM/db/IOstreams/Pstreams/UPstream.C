/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "UPstream.H"
#include "debug.H"
#include "dictionary.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(UPstream, 0);

    template<>
    const char* Foam::NamedEnum
    <
        Foam::UPstream::commsTypes,
        3
    >::names[] =
    {
        "blocking",
        "scheduled",
        "nonBlocking"
    };
}


const Foam::NamedEnum<Foam::UPstream::commsTypes, 3>
    Foam::UPstream::commsTypeNames;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::UPstream::setParRun()
{
    parRun_ = true;

    Pout.prefix() = '[' +  name(myProcNo()) + "] ";
    Perr.prefix() = '[' +  name(myProcNo()) + "] ";
}


void Foam::UPstream::calcLinearComm(const label nProcs)
{
    linearCommunication_.setSize(nProcs);

    // Master
    labelList belowIDs(nProcs - 1);
    forAll(belowIDs, i)
    {
        belowIDs[i] = i + 1;
    }

    linearCommunication_[0] = commsStruct
    (
        nProcs,
        0,
        -1,
        belowIDs,
        labelList(0)
    );

    // Slaves. Have no below processors, only communicate up to master
    for (label procID = 1; procID < nProcs; procID++)
    {
        linearCommunication_[procID] = commsStruct
        (
            nProcs,
            procID,
            0,
            labelList(0),
            labelList(0)
        );
    }
}


// Append my children (and my children children etc.) to allReceives.
void Foam::UPstream::collectReceives
(
    const label procID,
    const List<DynamicList<label> >& receives,
    DynamicList<label>& allReceives
)
{
    const DynamicList<label>& myChildren = receives[procID];

    forAll(myChildren, childI)
    {
        allReceives.append(myChildren[childI]);
        collectReceives(myChildren[childI], receives, allReceives);
    }
}


// Tree like schedule. For 8 procs:
// (level 0)
//      0 receives from 1
//      2 receives from 3
//      4 receives from 5
//      6 receives from 7
// (level 1)
//      0 receives from 2
//      4 receives from 6
// (level 2)
//      0 receives from 4
//
// The sends/receives for all levels are collected per processor (one send per
// processor; multiple receives possible) creating a table:
//
// So per processor:
// proc     receives from   sends to
// ----     -------------   --------
//  0       1,2,4           -
//  1       -               0
//  2       3               0
//  3       -               2
//  4       5               0
//  5       -               4
//  6       7               4
//  7       -               6
void Foam::UPstream::calcTreeComm(label nProcs)
{
    label nLevels = 1;
    while ((1 << nLevels) < nProcs)
    {
        nLevels++;
    }

    List<DynamicList<label> > receives(nProcs);
    labelList sends(nProcs, -1);

    // Info<< "Using " << nLevels << " communication levels" << endl;

    label offset = 2;
    label childOffset = offset/2;

    for (label level = 0; level < nLevels; level++)
    {
        label receiveID = 0;
        while (receiveID < nProcs)
        {
            // Determine processor that sends and we receive from
            label sendID = receiveID + childOffset;

            if (sendID < nProcs)
            {
                receives[receiveID].append(sendID);
                sends[sendID] = receiveID;
            }

            receiveID += offset;
        }

        offset <<= 1;
        childOffset <<= 1;
    }

    // For all processors find the processors it receives data from
    // (and the processors they receive data from etc.)
    List<DynamicList<label> > allReceives(nProcs);
    for (label procID = 0; procID < nProcs; procID++)
    {
        collectReceives(procID, receives, allReceives[procID]);
    }


    treeCommunication_.setSize(nProcs);

    for (label procID = 0; procID < nProcs; procID++)
    {
        treeCommunication_[procID] = commsStruct
        (
            nProcs,
            procID,
            sends[procID],
            receives[procID].shrink(),
            allReceives[procID].shrink()
        );
    }
}


// Callback from UPstream::init() : initialize linear and tree communication
// schedules now that nProcs is known.
void Foam::UPstream::initCommunicationSchedule()
{
    calcLinearComm(nProcs());
    calcTreeComm(nProcs());
}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Initialise my process number to 0 (the master)
int Foam::UPstream::myProcNo_(0);

// By default this is not a parallel run
bool Foam::UPstream::parRun_(false);

// List of process IDs
Foam::List<int> Foam::UPstream::procIDs_(label(1), 0);

// Standard transfer message type
int Foam::UPstream::msgType_(1);

// Linear communication schedule
Foam::List<Foam::UPstream::commsStruct> Foam::UPstream::linearCommunication_(0);

// Multi level communication schedule
Foam::List<Foam::UPstream::commsStruct> Foam::UPstream::treeCommunication_(0);

// Should compact transfer be used in which floats replace doubles
// reducing the bandwidth requirement at the expense of some loss
// in accuracy
bool Foam::UPstream::floatTransfer
(
    debug::optimisationSwitch("floatTransfer", 0)
);
registerOptSwitchWithName
(
    Foam::UPstream::floatTransfer,
    floatTransfer,
    "floatTransfer"
);

// Number of processors at which the reduce algorithm changes from linear to
// tree
int Foam::UPstream::nProcsSimpleSum
(
    debug::optimisationSwitch("nProcsSimpleSum", 16)
);
registerOptSwitchWithName
(
    Foam::UPstream::nProcsSimpleSum,
    nProcsSimpleSum,
    "nProcsSimpleSum"
);

// Default commsType
Foam::UPstream::commsTypes Foam::UPstream::defaultCommsType
(
    commsTypeNames.read(debug::optimisationSwitches().lookup("commsType"))
);
// Register re-reader
class addcommsTypeToOpt
:
    public ::Foam::simpleRegIOobject
{
public:
    addcommsTypeToOpt(const char* name)
    :
        ::Foam::simpleRegIOobject(Foam::debug::addOptimisationObject, name)
    {}
    virtual ~addcommsTypeToOpt()
    {}
    virtual void readData(Foam::Istream& is)
    {
        Foam::UPstream::defaultCommsType = Foam::UPstream::commsTypeNames.read
        (
            is
        );
    }
    virtual void writeData(Foam::Ostream& os) const
    {
        os << Foam::UPstream::commsTypeNames[Foam::UPstream::defaultCommsType];
    }
};
addcommsTypeToOpt addcommsTypeToOpt_("commsType");


// Number of polling cycles in processor updates
int Foam::UPstream::nPollProcInterfaces
(
    debug::optimisationSwitch("nPollProcInterfaces", 0)
);
registerOptSwitchWithName
(
    Foam::UPstream::nPollProcInterfaces,
    nPollProcInterfaces,
    "nPollProcInterfaces"
);

// ************************************************************************* //
