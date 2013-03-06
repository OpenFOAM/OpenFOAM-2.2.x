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

Class
    Foam::UPstream

Description
    Inter-processor communications stream

SourceFiles
    UPstream.C
    UPstreamsPrint.C
    UPstreamCommsStruct.C
    gatherScatter.C
    combineGatherScatter.C
    gatherScatterList.C

\*---------------------------------------------------------------------------*/

#ifndef UPstream_H
#define UPstream_H

#include "labelList.H"
#include "DynamicList.H"
#include "HashTable.H"
#include "string.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class UPstream Declaration
\*---------------------------------------------------------------------------*/

class UPstream
{

public:

    //- Types of communications
    enum commsTypes
    {
        blocking,
        scheduled,
        nonBlocking
    };

    static const NamedEnum<commsTypes, 3> commsTypeNames;

    // Public classes

        //- Structure for communicating between processors
        class commsStruct
        {
            // Private data

                //- procID of above processor
                label above_;

                //- procIDs of processors directly below me
                labelList below_;

                //- procIDs of all processors below (so not just directly below)
                labelList allBelow_;

                //- procIDs of all processors not below. (inverse set of
                //  allBelow_ and minus myProcNo)
                labelList allNotBelow_;


        public:

            // Constructors

                //- Construct null
                commsStruct();

                //- Construct from components
                commsStruct
                (
                    const label,
                    const labelList&,
                    const labelList&,
                    const labelList&
                );

                //- Construct from components; construct allNotBelow_
                commsStruct
                (
                    const label nProcs,
                    const label myProcID,
                    const label,
                    const labelList&,
                    const labelList&
                );


            // Member Functions

                // Access

                    label above() const
                    {
                        return above_;
                    }

                    const labelList& below() const
                    {
                        return below_;
                    }

                    const labelList& allBelow() const
                    {
                        return allBelow_;
                    }

                    const labelList& allNotBelow() const
                    {
                        return allNotBelow_;
                    }


            // Member operators

                bool operator==(const commsStruct&) const;

                bool operator!=(const commsStruct&) const;


             // Ostream Operator

                friend Ostream& operator<<(Ostream&, const commsStruct&);
        };


        //- combineReduce operator for lists. Used for counting.
        class listEq
        {

        public:

            template<class T>
            void operator()(T& x, const T& y) const
            {
                forAll(y, i)
                {
                    if (y[i].size())
                    {
                        x[i] = y[i];
                    }
                }
            }
        };


private:

    // Private data

        static int myProcNo_;
        static bool parRun_;

        static List<int> procIDs_;
        static int msgType_;

        static List<commsStruct> linearCommunication_;
        static List<commsStruct> treeCommunication_;


    // Private Member Functions

        //- Set data for parallel running
        static void setParRun();

        //- Calculate linear communication schedule
        static void calcLinearComm(const label nProcs);

        //- Calculate tree communication schedule
        static void calcTreeComm(const label nProcs);

        //- Helper function for tree communication schedule determination
        //  Collects all processorIDs below a processor
        static void collectReceives
        (
            const label procID,
            const List<DynamicList<label> >& receives,
            DynamicList<label>& allReceives
        );

        //- Initialize all communication schedules. Callback from
        //  UPstream::init()
        static void initCommunicationSchedule();


protected:

    // Protected data

        //- Communications type of this stream
        commsTypes commsType_;

public:

    // Declare name of the class and its debug switch
    ClassName("UPstream");


    // Static data

        //- Should compact transfer be used in which floats replace doubles
        //  reducing the bandwidth requirement at the expense of some loss
        //  in accuracy
        static bool floatTransfer;

        //- Number of processors at which the sum algorithm changes from linear
        //  to tree
        static int nProcsSimpleSum;

        //- Default commsType
        static commsTypes defaultCommsType;

        //- Number of polling cycles in processor updates
        static int nPollProcInterfaces;

    // Constructors

        //- Construct given optional buffer size
        UPstream(const commsTypes commsType)
        :
            commsType_(commsType)
        {}


    // Member functions

        //- Add the valid option this type of communications library
        //  adds/requires on the command line
        static void addValidParOptions(HashTable<string>& validParOptions);

        //- Initialisation function called from main
        //  Spawns slave processes and initialises inter-communication
        static bool init(int& argc, char**& argv);

        // Non-blocking comms

            //- Get number of outstanding requests
            static label nRequests();

            //- Truncate number of outstanding requests
            static void resetRequests(const label sz);

            //- Wait until all requests (from start onwards) have finished.
            static void waitRequests(const label start = 0);

            //- Wait until request i has finished.
            static void waitRequest(const label i);

            //- Non-blocking comms: has request i finished?
            static bool finishedRequest(const label i);


        //- Is this a parallel run?
        static bool& parRun()
        {
            return parRun_;
        }

        //- Number of processes in parallel run
        static label nProcs()
        {
            return procIDs_.size();
        }

        //- Am I the master process
        static bool master()
        {
            return myProcNo_ == 0;
        }

        //- Process index of the master
        static int masterNo()
        {
            return 0;
        }

        //- Number of this process (starting from masterNo() = 0)
        static int myProcNo()
        {
            return myProcNo_;
        }

        //- Process IDs
        static const List<int>& procIDs()
        {
            return procIDs_;
        }

        //- Process ID of given process index
        static int procID(int procNo)
        {
            return procIDs_[procNo];
        }

        //- Process index of first slave
        static int firstSlave()
        {
            return 1;
        }

        //- Process index of last slave
        static int lastSlave()
        {
            return nProcs() - 1;
        }

        //- Communication schedule for linear all-to-master (proc 0)
        static const List<commsStruct>& linearCommunication()
        {
            return linearCommunication_;
        }

        //- Communication schedule for tree all-to-master (proc 0)
        static const List<commsStruct>& treeCommunication()
        {
            return treeCommunication_;
        }

        //- Message tag of standard messages
        static int& msgType()
        {
            return msgType_;
        }


            //- Get the communications type of the stream
            commsTypes commsType() const
            {
                return commsType_;
            }

            //- Set the communications type of the stream
            commsTypes commsType(const commsTypes ct)
            {
                commsTypes oldCommsType = commsType_;
                commsType_ = ct;
                return oldCommsType;
            }


        //- Exit program
        static void exit(int errnum = 1);

        //- Abort program
        static void abort();


};


Ostream& operator<<(Ostream&, const UPstream::commsStruct&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
