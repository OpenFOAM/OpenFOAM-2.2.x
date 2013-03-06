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

\*---------------------------------------------------------------------------*/

#include "steadyParticleTracksTemplates.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
bool fieldOk(const IOobjectList& cloudObjs, const word& name)
{
    IOobjectList objects(cloudObjs.lookupClass(IOField<Type>::typeName));

    return (objects.lookup(name) != NULL);
}


template<class Type>
tmp<Field<Type> > readParticleField
(
    const word& name,
    const IOobjectList cloudObjs
)
{
    IOobjectList objects(cloudObjs.lookupClass(IOField<Type>::typeName));

    const IOobject* obj = objects.lookup(name);
    if (obj != NULL)
    {
        IOField<Type> newField(*obj);
        return tmp<Field<Type> >(new Field<Type>(newField.xfer()));
    }

    Info<< "error: cloud field name " << name << " not found" << endl;

    return Field<Type>::null();
}


template<class Type>
PtrList<List<Type> > readFields
(
    PtrList<List<Type> >& values,
    const List<word>& fields,
    const IOobjectList& cloudObjs
)
{
    IOobjectList objects(cloudObjs.lookupClass(IOField<Type>::typeName));

    label fieldI = 0;
    forAllConstIter(IOobjectList, objects, iter)
    {
        const IOobject& obj = *iter();
        forAll(fields, j)
        {
            if (obj.name() == fields[j])
            {
                Info<< "        reading field " << obj.name() << endl;
                IOField<Type> newField(obj);
                values.set(fieldI++, new List<Type>(newField.xfer()));
                break;
            }
        }
    }

    return values;
}


template<class Type>
void writeVTK(OFstream& os, const Type& value)
{
    os  << value.component(0);
    for (label i=1; i<pTraits<Type>::nComponents; i++)
    {
        os  << ' ' << value.component(i);
    }
}


template<class Type>
void writeVTKFields
(
    OFstream& os,
    const PtrList<List<Type> >& values,
    const List<SortableList<scalar> >& agePerTrack,
    const List<word>& fieldNames
)
{
    label step = max(floor(8/pTraits<Type>::nComponents), 1);

    forAll(values, fieldI)
    {
        Info<< "        writing field " << fieldNames[fieldI] << endl;
        os  << nl << fieldNames[fieldI] << ' ' << pTraits<Type>::nComponents
            << ' ' << values[fieldI].size() << " float" << nl;
        label offset = 0;
        forAll(agePerTrack, trackI)
        {
            const List<label> ids = agePerTrack[trackI].indices() + offset;

            List<Type> data(UIndirectList<Type>(values[fieldI], ids));
            label nData = data.size() - 1;
            forAll(data, i)
            {
                writeVTK<Type>(os, data[i]);
                if (((i + 1) % step == 0) || (i == nData))
                {
                    os  << nl;
                }
                else
                {
                    os  << ' ';
                }
            }
            offset += ids.size();
        }
    }
}


template<class Type>
void processFields
(
    OFstream& os,
    const List<SortableList<scalar> >& agePerTrack,
    const List<word>& userFieldNames,
    const IOobjectList& cloudObjs
)
{
    IOobjectList objects(cloudObjs.lookupClass(IOField<Type>::typeName));

    if (objects.size())
    {
        DynamicList<word> fieldNames(objects.size());
        forAll(userFieldNames, i)
        {
            IOobject* obj = objects.lookup(userFieldNames[i]);
            if (obj != NULL)
            {
                fieldNames.append(obj->name());
            }
        }
        fieldNames.shrink();

        PtrList<List<Type> > values(fieldNames.size());
        readFields<Type>(values, fieldNames, cloudObjs);

        writeVTKFields<Type>
        (
            os,
            values,
            agePerTrack,
            fieldNames.xfer()
        );
    }
}

} // End namespace Foam

// ************************************************************************* //
