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

Class
    Foam::globalIndex

Description
    Calculates a unique integer (label so might not have enough room - 2G max)
    for processor + local index. E.g.

    globalIndex globalFaces(mesh.nFaces());
    label globalFaceI = globalFaces.toGlobal(faceI);


SourceFiles
    globalIndex.C

\*---------------------------------------------------------------------------*/

#ifndef globalIndex_H
#define globalIndex_H

#include "Pstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class globalIndex;

Istream& operator>>(Istream& is, globalIndex& gi);
Ostream& operator<<(Ostream& os, const globalIndex& gi);


/*---------------------------------------------------------------------------*\
                           Class globalIndex Declaration
\*---------------------------------------------------------------------------*/

class globalIndex
{
    // Private data

        //- Start of procI. Size is nProcs()+1. (so like CompactListList)
        labelList offsets_;


public:

    // Constructors

        //- Construct from local max size
        globalIndex
        (
            const label localSize,
            const int tag = Pstream::msgType(),
            const bool parallel = true              // use parallel comms
        );

        //- Construct from Istream
        globalIndex(Istream& is);


    // Member Functions

        // Queries relating to my processor

            //- my local size
            inline label localSize() const;

            //- From local to global
            inline label toGlobal(const label i) const;

            //- Is on local processor
            inline bool isLocal(const label i) const;

            //- From global to local on current processor.
            //  FatalError if not on local processor.
            inline label toLocal(const label i) const;


        // Global queries

            //- Global sum of localSizes
            inline label size() const;

            //- Size of procI data
            inline label localSize(const label procI) const;

            //- From local to global on procI
            inline label toGlobal(const label procI, const label i) const;

            //- Is on processor procI
            inline bool isLocal(const label procI, const label i) const;

            //- From global to local on procI
            inline label toLocal(const label procI, const label i) const;

            //- Which processor does global come from? Binary search.
            inline label whichProcID(const label i) const;

            //- Start of procI data
            inline label offset(const label procI) const;



    // IOstream Operators

        friend Istream& operator>>(Istream& is, globalIndex& gi);
        friend Ostream& operator<<(Ostream& os, const globalIndex& gi);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "globalIndexI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
