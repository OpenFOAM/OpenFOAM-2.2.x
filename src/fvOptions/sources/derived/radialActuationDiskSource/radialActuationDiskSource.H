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
    Foam::fv::radialActuationDiskSource

Description
    Actuation disk zone definition.
    Constant values for momentum source for actuation disk

        T = 2*rho*A*sqr(Uo)*a*(1-a)
        U1 = (1 -a)Uo

    where:
        A: disk area
        Uo: upstream velocity
        a: 1 - Cp/Ct
        U1: velocity at the disk


    The thrust is distributed by a radial function:

        thrust(r) = T*(C0 + C1*r^2 + C2*r^4)

    Sources described by:

        actuationDiskSourceCoeffs
        {
            fieldName       U;          // name of field to apply source
            diskDir         (-1 0 0);   // disk direction
            Cp              0.1;        // power coefficient
            Ct              0.5;        // thrust coefficient
            diskArea        5.0;        // disk area
            coeffs          (0.1 0.5 0.01); // radial distribution coefficients
            upstreamPoint   (0 0 0);    // upstream point
        }


SourceFiles
    radialActuationDiskSource.C
    radialActuationDiskSourceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef radialActuationDiskSource_H
#define radialActuationDiskSource_H

#include "actuationDiskSource.H"
#include "FixedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                  Class radialActuationDiskSource Declaration
\*---------------------------------------------------------------------------*/

class radialActuationDiskSource
:
    public actuationDiskSource
{
    // Private data

        //- Coeffcients for the radial distribution
        FixedList<scalar, 3> radialCoeffs_;


    // Private Member Functions

        //- Add resistance to the UEqn
        template<class RhoFieldType>
        void addRadialActuationDiskAxialInertialResistance
        (
            vectorField& Usource,
            const labelList& cells,
            const scalarField& V,
            const RhoFieldType& rho,
            const vectorField& U
        ) const;

        //- Disallow default bitwise copy construct
        radialActuationDiskSource(const radialActuationDiskSource&);

        //- Disallow default bitwise assignment
        void operator=(const radialActuationDiskSource&);


public:

    //- Runtime type information
    TypeName("radialActuationDiskSource");


    // Constructors

        //- Construct from components
        radialActuationDiskSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~radialActuationDiskSource()
    {}


    // Member Functions

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
#   include "radialActuationDiskSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
