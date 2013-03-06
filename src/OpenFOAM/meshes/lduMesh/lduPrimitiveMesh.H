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
    Foam::lduPrimitiveMesh

Description
    Simplest contrete lduMesh which stores the addressing needed by lduMatrix.

\*---------------------------------------------------------------------------*/

#ifndef lduPrimitiveMesh_H
#define lduPrimitiveMesh_H

#include "lduMesh.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class lduPrimitiveMesh Declaration
\*---------------------------------------------------------------------------*/

class lduPrimitiveMesh
:
    public lduMesh,
    public lduAddressing
{
    // Private data

        //- Lower addressing
        labelList lowerAddr_;

        //- Upper addressing
        labelList upperAddr_;

        //- Patch to internal addressing
        labelListList patchAddr_;

        //- List of pointers for each patch
        //  with only those pointing to interfaces being set
        lduInterfacePtrsList interfaces_;

        //- Patch field evaluation schedule.
        //  Note this does not need to be held as a copy because it is invariant
        const lduSchedule& patchSchedule_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        lduPrimitiveMesh(const lduPrimitiveMesh&);

        //- Disallow default bitwise assignment
        void operator=(const lduPrimitiveMesh&);


public:

    // Constructors

        //- Construct from components as copies
        lduPrimitiveMesh
        (
            const label nCells,
            const labelUList& l,
            const labelUList& u,
            const labelListList& pa,
            lduInterfacePtrsList interfaces,
            const lduSchedule& ps
        )
        :
            lduAddressing(nCells),
            lowerAddr_(l),
            upperAddr_(u),
            patchAddr_(pa),
            interfaces_(interfaces),
            patchSchedule_(ps)
        {}


        //- Construct from components and  re-use storage as specified.
        lduPrimitiveMesh
        (
            const label nCells,
            labelList& l,
            labelList& u,
            labelListList& pa,
            lduInterfacePtrsList interfaces,
            const lduSchedule& ps,
            bool reUse
        )
        :
            lduAddressing(nCells),
            lowerAddr_(l, reUse),
            upperAddr_(u, reUse),
            patchAddr_(pa, reUse),
            interfaces_(interfaces, reUse),
            patchSchedule_(ps)
        {}


    //- Destructor
    virtual ~lduPrimitiveMesh()
    {}


    // Member Functions

        // Access

            //- Return ldu addressing
            virtual const lduAddressing& lduAddr() const
            {
                return *this;
            }

            //- Return a list of pointers for each patch
            //  with only those pointing to interfaces being set
            virtual lduInterfacePtrsList interfaces() const
            {
                return interfaces_;
            }

            //- Return Lower addressing
            virtual const labelUList& lowerAddr() const
            {
                return lowerAddr_;
            }

            //- Return Upper addressing
            virtual const labelUList& upperAddr() const
            {
                return upperAddr_;
            }

            //- Return patch addressing
            virtual const labelUList& patchAddr(const label i) const
            {
                return patchAddr_[i];
            }

            //- Return patch evaluation schedule
            virtual const lduSchedule& patchSchedule() const
            {
                return patchSchedule_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
