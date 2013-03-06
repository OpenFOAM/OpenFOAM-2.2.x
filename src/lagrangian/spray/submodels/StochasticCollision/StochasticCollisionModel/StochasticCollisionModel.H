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
    Foam::StochasticCollisionModel

Description
    Templated stochastic collision model class

SourceFiles
    StochasticCollisionModel.C
    StochasticCollisionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef StochasticCollisionModel_H
#define StochasticCollisionModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "SubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class StochasticCollisionModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class StochasticCollisionModel
:
    public SubModelBase<CloudType>
{
public:

    //- Runtime type information
    TypeName("collisionModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        StochasticCollisionModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        StochasticCollisionModel(CloudType& owner);

        //- Construct from dictionary
        StochasticCollisionModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Construct copy
        StochasticCollisionModel(const StochasticCollisionModel<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<StochasticCollisionModel<CloudType> > clone() const
        {
            return autoPtr<StochasticCollisionModel<CloudType> >
            (
                new StochasticCollisionModel<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~StochasticCollisionModel();


    //- Selector
    static autoPtr<StochasticCollisionModel<CloudType> > New
    (
        const dictionary& dict,
        CloudType& owner
    );


    // Member Functions

        virtual bool update
        (
            const scalar dt,
            cachedRandom& rndGen,
            vector& pos1,
            scalar& m1,
            scalar& d1,
            scalar& N1,
            vector& U1,
            scalar& rho1,
            scalar& T1,
            scalarField& Y1,
            const scalar sigma1,
            const label celli,
            const scalar voli,
            vector& pos2,
            scalar& m2,
            scalar& d2,
            scalar& N2,
            vector& U2,
            scalar& rho2,
            scalar& T2,
            scalarField& Y2,
            const scalar sigma2,
            const label cellj,
            const scalar volj
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeStochasticCollisionModel(CloudType)                                         \
                                                                              \
    typedef CloudType::sprayCloudType sprayCloudType;                         \
    defineNamedTemplateTypeNameAndDebug                                       \
    (                                                                         \
        StochasticCollisionModel<sprayCloudType>,                             \
        0                                                                     \
    );                                                                        \
    defineTemplateRunTimeSelectionTable                                       \
    (                                                                         \
        StochasticCollisionModel<sprayCloudType>,                             \
        dictionary                                                            \
    );


#define makeStochasticCollisionModelType(SS, CloudType)                       \
                                                                              \
    typedef CloudType::sprayCloudType sprayCloudType;                         \
    defineNamedTemplateTypeNameAndDebug(SS<sprayCloudType>, 0);               \
                                                                              \
    StochasticCollisionModel<sprayCloudType>::                                \
        adddictionaryConstructorToTable<SS<sprayCloudType> >                  \
            add##SS##CloudType##sprayCloudType##ConstructorToTable_;



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "StochasticCollisionModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
