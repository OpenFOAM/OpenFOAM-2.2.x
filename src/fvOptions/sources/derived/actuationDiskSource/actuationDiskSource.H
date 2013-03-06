/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::fv::actuationDiskSource

Description
    Actuation disk source

    Constant values for momentum source for actuation disk

        T = 2*rho*A*sqr(Uo)*a*(1-a)
        U1 = (1 -a)Uo

    where:
        A: disk area
        Uo: upstream velocity
        a: 1 - Cp/Ct
        U1: velocity at the disk

    Sources described by:

        actuationDiskSourceCoeffs
        {
            fieldNames      (U);        // names of fields to apply source
            diskDir         (-1 0 0);   // disk direction
            Cp              0.1;        // power coefficient
            Ct              0.5;        // thrust coefficient
            diskArea        5.0;        // disk area
            upstreamPoint   (0 0 0);    // upstream point
        }


SourceFiles
    actuationDiskSource.C
    actuationDiskSourceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef actuationDiskSource_H
#define actuationDiskSource_H

#include "fvOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                     Class actuationDiskSource Declaration
\*---------------------------------------------------------------------------*/

class actuationDiskSource
:
    public option
{

protected:

    // Protected data

        //- Disk area normal
        vector diskDir_;

        //- Power coefficient
        scalar Cp_;

        //- Thrust coefficient
        scalar Ct_;

        //- Disk area
        scalar diskArea_;

        //- Upstream point sample
        point upstreamPoint_;

        //- Upstream cell ID
        label upstreamCellId_;


private:

    // Private Member Functions

        //- Check data
        void checkData() const;

        //- Add resistance to the UEqn
        template<class RhoFieldType>
        void addActuationDiskAxialInertialResistance
        (
            vectorField& Usource,
            const labelList& cells,
            const scalarField& V,
            const RhoFieldType& rho,
            const vectorField& U
        ) const;

        //- Disallow default bitwise copy construct
        actuationDiskSource(const actuationDiskSource&);

        //- Disallow default bitwise assignment
        void operator=(const actuationDiskSource&);


public:

    //- Runtime type information
    TypeName("actuationDiskSource");


    // Constructors

        //- Construct from components
        actuationDiskSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~actuationDiskSource()
    {}


    // Member Functions

        // Access

            //- Return Cp
            scalar Cp() const
            {
                return Cp_;
            }

            //- Return Ct
            scalar Ct() const
            {
                return Ct_;
            }

            //- Normal disk direction
            const vector& diskDir() const
            {
                return diskDir_;
            }

            //- Disk area
            scalar diskArea() const
            {
                return diskArea_;
            }


        // Public Functions

            //- Source term to fvMatrix<vector>
            virtual void addSup(fvMatrix<vector>& eqn, const label fieldI);


        // I-O

            //- Write data
            virtual void writeData(Ostream&) const;

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "actuationDiskSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
