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
    Foam::processorGAMGInterface

Description
    GAMG agglomerated processor interface.

SourceFiles
    processorGAMGInterface.C
    processorGAMGInterfaceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef processorGAMGInterface_H
#define processorGAMGInterface_H

#include "GAMGInterface.H"
#include "processorLduInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class processorGAMGInterface Declaration
\*---------------------------------------------------------------------------*/

class processorGAMGInterface
:
    public GAMGInterface,
    public processorLduInterface
{
    // Private data

        //- Reference tor the processorLduInterface from which this is
        //  agglomerated
        const processorLduInterface& fineProcInterface_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        processorGAMGInterface(const processorGAMGInterface&);

        //- Disallow default bitwise assignment
        void operator=(const processorGAMGInterface&);


public:

    //- Runtime type information
    TypeName("processor");

    // Constructors

        //- Construct from fine-level interface,
        //  local and neighbour restrict addressing
        processorGAMGInterface
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            const lduInterface& fineInterface,
            const labelField& restrictAddressing,
            const labelField& neighbourRestrictAddressing,
            const label fineLevelIndex
        );


    //- Destructor
    virtual ~processorGAMGInterface();


    // Member Functions

        // Interface transfer functions

            //- Initialise neighbour field transfer
            virtual void initInternalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const;

            //- Transfer and return internal field adjacent to the interface
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const;


        //- Processor interface functions

            //- Return processor number
            virtual int myProcNo() const
            {
                return fineProcInterface_.myProcNo();
            }

            //- Return neigbour processor number
            virtual int neighbProcNo() const
            {
                return fineProcInterface_.neighbProcNo();
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return fineProcInterface_.forwardT();
            }

            //- Return message tag used for sending
            virtual int tag() const
            {
                return fineProcInterface_.tag();
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
