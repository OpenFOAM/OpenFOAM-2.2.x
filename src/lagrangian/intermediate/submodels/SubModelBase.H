/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::SubModelBase

Description
    Base class for cloud sub-models

SourceFiles
    SubModelBase.C

\*---------------------------------------------------------------------------*/

#ifndef SubModelBase_H
#define SubModelBase_H

#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
template<class CloudType>
class SubModelBase;

/*---------------------------------------------------------------------------*\
                       Class SubModelBase Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class SubModelBase
{
protected:

    // Protected Data

        //- Reference to the cloud
        CloudType& owner_;

        //- Reference to the cloud dictionary
        const dictionary dict_;

        //- Name of the sub-model base class
        const word baseName_;

        //- Type of the sub-model
        const word modelType_;

        //- Name of the sub-model
        const word modelName_;

        //- Coefficients dictionary
        const dictionary coeffDict_;


    // Protected Member Functions

        //- Flag to indicate whether data is/was read in-line
        bool inLine() const;


public:

    // Constructors

        //- Construct null from owner cloud
        SubModelBase(CloudType& owner);

        //- Construct from owner cloud, dictionary, and model type name
        SubModelBase
        (
            CloudType& owner,
            const dictionary& dict,
            const word& baseName,
            const word& modelType,
            const word& dictExt = "Coeffs"
        );

        //- Construct from owner cloud, dictionary, and model type name
        SubModelBase
        (
            const word& modelName,
            CloudType& owner,
            const dictionary& dict,
            const word& baseName,
            const word& modelType
        );

        //- Construct as copy
        SubModelBase(const SubModelBase<CloudType>& smb);


    //- Destructor
    virtual ~SubModelBase();

    //- Type of cloud this model was instantiated for
    typedef CloudType cloudType;


    // Member Functions

        // Access

            //- Return const access to the owner cloud
            const CloudType& owner() const;

            //- Return const access to the cloud dictionary
            const dictionary& dict() const;

            //- Return const access to the sub-model type
            const word& modelType() const;

            //- Return const access to the base name of the sub-model
            const word& baseName() const;

            //- Return const access to the name of the sub-model
            const word& modelName() const;

            //- Return const access to the coefficients dictionary
            const dictionary& coeffDict() const;

            //- Return const access to the properties dictionary
            const IOdictionary& properties() const;

            //- Returns true if defaultCoeffs is true and outputs on printMsg
            bool defaultCoeffs(const bool printMsg) const;

            //- Return the model 'active' status - default active = true
            virtual bool active() const;

            //- Cache dependant sub-model fields
            virtual void cacheFields(const bool store);

            //- Flag to indicate when to write a property
            bool outputTime() const;


        // Edit

            //- Return non-const access to the owner cloud for manipulation
            CloudType& owner();

            // Base properties

                //- Retrieve generic property from the base model
                template<class Type>
                Type getBaseProperty
                (
                    const word& entryName,
                    const Type& defaultValue = pTraits<Type>::zero
                ) const;

                //- Retrieve generic property from the base model
                template<class Type>
                void getBaseProperty(const word& entryName, Type& value) const;

                //- Add generic property to the base model
                template<class Type>
                void setBaseProperty(const word& entryName, const Type& value);


            // Model properties

                //- Retrieve generic property from the sub-model
                template<class Type>
                Type getModelProperty
                (
                    const word& entryName,
                    const Type& defaultValue = pTraits<Type>::zero
                ) const;

                //- Retrieve generic property from the sub-model
                template<class Type>
                void getModelProperty(const word& entryName, Type& value) const;

                //- Add generic property to the sub-model
                template<class Type>
                void setModelProperty(const word& entryName, const Type& value);


        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "SubModelBase.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
