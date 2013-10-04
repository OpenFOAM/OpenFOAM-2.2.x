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

#include "nearWallFields.H"
#include "mappedFieldFvPatchFields.H"
//#include "interpolationCellPoint.H"
#include "cachedInterpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::nearWallFields::createFields
(
    PtrList<GeometricField<Type, fvPatchField, volMesh> >& sflds
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> vfType;

    HashTable<const vfType*> flds(obr_.lookupClass<vfType>());

    forAllConstIter(typename HashTable<const vfType*>, flds, iter)
    {
        const vfType& fld = *iter();

        if (fieldMap_.found(fld.name()))
        {
            const word& sampleFldName = fieldMap_[fld.name()];

            if (obr_.found(sampleFldName))
            {
                Info<< "    a field " << sampleFldName
                    << " already exists on the mesh."
                    << endl;
            }
            else
            {
                label sz = sflds.size();
                sflds.setSize(sz+1);

                IOobject io(fld);
                io.readOpt() = IOobject::NO_READ;
                io.writeOpt() = IOobject::NO_WRITE;
                io.rename(sampleFldName);

                sflds.set(sz, new vfType(io, fld));
                vfType& sampleFld = sflds[sz];

                // Reset the bcs to be mapped
                forAllConstIter(labelHashSet, patchSet_, iter)
                {
                    label patchI = iter.key();

                    sampleFld.boundaryField().set
                    (
                        patchI,
                        new mappedFieldFvPatchField<Type>
                        (
                            sampleFld.mesh().boundary()[patchI],
                            sampleFld.dimensionedInternalField(),

                            sampleFld.mesh().name(),
                            mappedPatchBase::NEARESTCELL,
                            word::null,     // samplePatch
                            -distance_,

                            sampleFld.name(),       // fieldName
                            false,                  // setAverage
                            pTraits<Type>::zero,    // average
                            cachedInterpolationCellPoint<Type>::typeName
                        )
                    );
                }

                Info<< "    created " << sampleFld.name() << " to sample "
                    << fld.name() << endl;
            }
        }
    }
}


template<class Type>
void Foam::nearWallFields::sampleFields
(
    PtrList<GeometricField<Type, fvPatchField, volMesh> >& sflds
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> vfType;

    forAll(sflds, i)
    {
        const word& fldName = reverseFieldMap_[sflds[i].name()];
        const vfType& fld = obr_.lookupObject<vfType>(fldName);

        // Take over internal and boundary values
        sflds[i] == fld;
        // Evaluate to update the mapped
        sflds[i].correctBoundaryConditions();
    }
}


// ************************************************************************* //
