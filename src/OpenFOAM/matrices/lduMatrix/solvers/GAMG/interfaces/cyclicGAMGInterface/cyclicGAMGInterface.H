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
    Foam::cyclicGAMGInterface

Description
    GAMG agglomerated cyclic interface.

SourceFiles
    cyclicGAMGInterface.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicGAMGInterface_H
#define cyclicGAMGInterface_H

#include "GAMGInterface.H"
#include "cyclicLduInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class cyclicGAMGInterface Declaration
\*---------------------------------------------------------------------------*/

class cyclicGAMGInterface
:
    public GAMGInterface,
    virtual public cyclicLduInterface
{
    // Private data

        //- Reference tor the cyclicLduInterface from which this is
        //  agglomerated
        const cyclicLduInterface& fineCyclicInterface_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        cyclicGAMGInterface(const cyclicGAMGInterface&);

        //- Disallow default bitwise assignment
        void operator=(const cyclicGAMGInterface&);


public:

    //- Runtime type information
    TypeName("cyclic");


    // Constructors

        //- Construct from fine level interface,
        //  local and neighbour restrict addressing
        cyclicGAMGInterface
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            const lduInterface& fineInterface,
            const labelField& restrictAddressing,
            const labelField& neighbourRestrictAddressing,
            const label fineLevelIndex
        );


    //- Destructor
    virtual ~cyclicGAMGInterface();


    // Member Functions

        // Interface transfer functions

            //- Transfer and return internal field adjacent to the interface
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const;


        //- Cyclic interface functions

            //- Return neigbour processor number
            virtual label neighbPatchID() const
            {
                return fineCyclicInterface_.neighbPatchID();
            }

            virtual bool owner() const
            {
                return fineCyclicInterface_.owner();
            }

            virtual const cyclicGAMGInterface& neighbPatch() const
            {
                return dynamic_cast<const cyclicGAMGInterface&>
                (
                    coarseInterfaces_[neighbPatchID()]
                );
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return fineCyclicInterface_.forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return fineCyclicInterface_.reverseT();
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
