/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "fieldAverageItem.H"
#include "volFields.H"
#include "OFstream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::fieldAverage::addMeanField
(
    const label fieldI,
    wordList& meanFieldList
) const
{
    if (faItems_[fieldI].mean())
    {
        typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

        const word& fieldName = faItems_[fieldI].fieldName();

        word meanFieldName = fieldName + EXT_MEAN;
        if
        (
            (faItems_[fieldI].window() > 0)
         && (faItems_[fieldI].windowName() != "")
        )
        {
            meanFieldName = meanFieldName + "_" + faItems_[fieldI].windowName();
        }

        Info<< "Reading/calculating field " << meanFieldName << nl << endl;

        if (obr_.foundObject<fieldType>(meanFieldName))
        {
            meanFieldList[fieldI] = meanFieldName;
        }
        else if (obr_.found(meanFieldName))
        {
            Info<< "Cannot allocate average field " << meanFieldName
                << " since an object with that name already exists."
                << " Disabling averaging." << nl << endl;
            meanFieldList[fieldI] = word::null;
        }
        else
        {
            const fieldType& baseField =
                obr_.lookupObject<fieldType>(fieldName);

            // Store on registry
            obr_.store
            (
                new fieldType
                (
                    IOobject
                    (
                        meanFieldName,
                        obr_.time().timeName(obr_.time().startTime().value()),
                        obr_,
                        IOobject::READ_IF_PRESENT,
                        IOobject::NO_WRITE
                    ),
                    1*baseField
                )
            );

            meanFieldList[fieldI] = meanFieldName;
        }
    }
}


template<class Type1, class Type2>
void Foam::fieldAverage::addPrime2MeanField
(
    const label fieldI,
    const wordList& meanFieldList,
    wordList& prime2MeanFieldList
) const
{
    if (faItems_[fieldI].mean() && meanFieldList[fieldI].size())
    {
        typedef GeometricField<Type1, fvPatchField, volMesh> fieldType1;
        typedef GeometricField<Type2, fvPatchField, volMesh> fieldType2;

        const word& fieldName = faItems_[fieldI].fieldName();

        word meanFieldName = fieldName + EXT_PRIME2MEAN;
        if
        (
            (faItems_[fieldI].window() > 0)
         && (faItems_[fieldI].windowName() != "")
        )
        {
            meanFieldName = meanFieldName + "_" + faItems_[fieldI].windowName();
        }

        Info<< "Reading/calculating field " << meanFieldName << nl << endl;

        if (obr_.foundObject<fieldType2>(meanFieldName))
        {
            prime2MeanFieldList[fieldI] = meanFieldName;
        }
        else if (obr_.found(meanFieldName))
        {
            Info<< "Cannot allocate average field " << meanFieldName
                << " since an object with that name already exists."
                << " Disabling averaging." << nl << endl;
            prime2MeanFieldList[fieldI] = word::null;
        }
        else
        {
            const fieldType1& baseField =
                obr_.lookupObject<fieldType1>(fieldName);
            const fieldType1& meanField =
                obr_.lookupObject<fieldType1>(meanFieldList[fieldI]);

            obr_.store
            (
                new fieldType2
                (
                    IOobject
                    (
                        meanFieldName,
                        obr_.time().timeName(obr_.time().startTime().value()),
                        obr_,
                        IOobject::READ_IF_PRESENT,
                        IOobject::NO_WRITE
                    ),
                    sqr(baseField) - sqr(meanField)
                )
            );

            prime2MeanFieldList[fieldI] = meanFieldName;
        }
    }
}


template<class Type>
void Foam::fieldAverage::calculateMeanFields(const wordList& meanFieldList)
const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    scalar dt = obr_.time().deltaTValue();

    forAll(faItems_, i)
    {
        if (faItems_[i].mean() && meanFieldList[i].size())
        {
            const word& fieldName = faItems_[i].fieldName();
            const fieldType& baseField =
                obr_.lookupObject<fieldType>(fieldName);
            fieldType& meanField = const_cast<fieldType&>
            (
                obr_.lookupObject<fieldType>(meanFieldList[i])
            );

            scalar Dt = totalTime_[i];
            if (faItems_[i].iterBase())
            {
                dt = 1.0;
                Dt = scalar(totalIter_[i]);
            }

            scalar alpha = (Dt - dt)/Dt;
            scalar beta = dt/Dt;
            if (faItems_[i].window() > 0)
            {
                const scalar w = faItems_[i].window();

                if (Dt - dt >= w)
                {
                    alpha = (w - dt)/w;
                    beta = dt/w;
                }
            }

            meanField = alpha*meanField + beta*baseField;
        }
    }
}


template<class Type1, class Type2>
void Foam::fieldAverage::calculatePrime2MeanFields
(
    const wordList& meanFieldList,
    const wordList& prime2MeanFieldList
) const
{
    typedef GeometricField<Type1, fvPatchField, volMesh> fieldType1;
    typedef GeometricField<Type2, fvPatchField, volMesh> fieldType2;

    scalar dt = obr_.time().deltaTValue();

    forAll(faItems_, i)
    {
        if
        (
            faItems_[i].prime2Mean()
         && meanFieldList[i].size()
         && prime2MeanFieldList[i].size()
        )
        {
            const word& fieldName = faItems_[i].fieldName();
            const fieldType1& baseField =
                obr_.lookupObject<fieldType1>(fieldName);
            const fieldType1& meanField =
                obr_.lookupObject<fieldType1>(meanFieldList[i]);
            fieldType2& prime2MeanField = const_cast<fieldType2&>
            (
                obr_.lookupObject<fieldType2>(prime2MeanFieldList[i])
            );

            scalar Dt = totalTime_[i];
            if (faItems_[i].iterBase())
            {
                dt = 1.0;
                Dt = scalar(totalIter_[i]);
            }

            scalar alpha = (Dt - dt)/Dt;
            scalar beta = dt/Dt;
            if (faItems_[i].window() > 0)
            {
                const scalar w = faItems_[i].window();

                if (Dt - dt >= w)
                {
                    alpha = (w - dt)/w;
                    beta = dt/w;
                }
            }

            prime2MeanField =
                alpha*prime2MeanField
              + beta*sqr(baseField)
              - sqr(meanField);
        }
    }
}


template<class Type1, class Type2>
void Foam::fieldAverage::addMeanSqrToPrime2Mean
(
    const wordList& meanFieldList,
    const wordList& prime2MeanFieldList
) const
{
    typedef GeometricField<Type1, fvPatchField, volMesh> fieldType1;
    typedef GeometricField<Type2, fvPatchField, volMesh> fieldType2;

    forAll(faItems_, i)
    {
        if
        (
            faItems_[i].prime2Mean()
         && meanFieldList[i].size()
         && prime2MeanFieldList[i].size()
        )
        {
            const fieldType1& meanField =
                obr_.lookupObject<fieldType1>(meanFieldList[i]);
            fieldType2& prime2MeanField = const_cast<fieldType2&>
            (
                obr_.lookupObject<fieldType2>(prime2MeanFieldList[i])
            );

            prime2MeanField += sqr(meanField);
        }
    }
}


template<class Type>
void Foam::fieldAverage::writeFieldList(const wordList& fieldList) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    forAll(fieldList, i)
    {
        if (fieldList[i].size())
        {
            const fieldType& f = obr_.lookupObject<fieldType>(fieldList[i]);
            f.write();
        }
    }
}


// ************************************************************************* //
