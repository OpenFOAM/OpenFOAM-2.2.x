/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::DPMDragModels::GidaspowSchillerNaumann

Description
    Gidaspow Schiller Naumann drag model for solid spheres

SourceFiles
    GidaspowSchillerNaumann.C

\*---------------------------------------------------------------------------*/

#ifndef GidaspowSchillerNaumann_H
#define GidaspowSchillerNaumann_H

#include "ParticleForce.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace DPMDragModels
{

/*---------------------------------------------------------------------------*\
                       Class GidaspowSchillerNaumann Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class GidaspowSchillerNaumann
:
    public ParticleForce<CloudType>
{
    // Private Data

        const volScalarField& alphac_;


    // Private Member Functions

        //- Drag coefficient multiplied by Reynolds number
        scalar CdRe(const scalar Re) const;


public:

    //- Runtime type information
    TypeName("GidaspowSchillerNaumannDrag");


    // Constructors

        //- Construct from mesh
        GidaspowSchillerNaumann
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        GidaspowSchillerNaumann(const GidaspowSchillerNaumann<CloudType>& df);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType> > clone() const
        {
            return autoPtr<ParticleForce<CloudType> >
            (
                new GidaspowSchillerNaumann<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~GidaspowSchillerNaumann();


    // Member Functions

        // Evaluation

            //- Calculate the coupled force
            virtual forceSuSp calcCoupled
            (
                const typename CloudType::parcelType& p,
                const scalar dt,
                const scalar mass,
                const scalar Re,
                const scalar muc
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace DPMDragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "GidaspowSchillerNaumann.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
