/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "fvMeshTools.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class GeoField>
void Foam::fvMeshTools::addPatchFields
(
    fvMesh& mesh,
    const dictionary& patchFieldDict,
    const word& defaultPatchFieldType,
    const typename GeoField::value_type& defaultPatchValue
)
{
    HashTable<const GeoField*> flds
    (
        mesh.objectRegistry::lookupClass<GeoField>()
    );

    forAllConstIter(typename HashTable<const GeoField*>, flds, iter)
    {
        const GeoField& fld = *iter();

        typename GeoField::GeometricBoundaryField& bfld =
            const_cast<typename GeoField::GeometricBoundaryField&>
            (
                fld.boundaryField()
            );

        label sz = bfld.size();
        bfld.setSize(sz+1);

        if (patchFieldDict.found(fld.name()))
        {
            bfld.set
            (
                sz,
                GeoField::PatchFieldType::New
                (
                    mesh.boundary()[sz],
                    fld.dimensionedInternalField(),
                    patchFieldDict.subDict(fld.name())
                )
            );
        }
        else
        {
            bfld.set
            (
                sz,
                GeoField::PatchFieldType::New
                (
                    defaultPatchFieldType,
                    mesh.boundary()[sz],
                    fld.dimensionedInternalField()
                )
            );
            bfld[sz] == defaultPatchValue;
        }
    }
}


template<class GeoField>
void Foam::fvMeshTools::setPatchFields
(
    fvMesh& mesh,
    const label patchI,
    const dictionary& patchFieldDict
)
{
    HashTable<const GeoField*> flds
    (
        mesh.objectRegistry::lookupClass<GeoField>()
    );

    forAllConstIter(typename HashTable<const GeoField*>, flds, iter)
    {
        const GeoField& fld = *iter();

        typename GeoField::GeometricBoundaryField& bfld =
            const_cast<typename GeoField::GeometricBoundaryField&>
            (
                fld.boundaryField()
            );

        if (patchFieldDict.found(fld.name()))
        {
            bfld.set
            (
                patchI,
                GeoField::PatchFieldType::New
                (
                    mesh.boundary()[patchI],
                    fld.dimensionedInternalField(),
                    patchFieldDict.subDict(fld.name())
                )
            );
        }
    }
}




template<class GeoField>
void Foam::fvMeshTools::setPatchFields
(
    fvMesh& mesh,
    const label patchI,
    const typename GeoField::value_type& value
)
{
    HashTable<const GeoField*> flds
    (
        mesh.objectRegistry::lookupClass<GeoField>()
    );

    forAllConstIter(typename HashTable<const GeoField*>, flds, iter)
    {
        const GeoField& fld = *iter();

        typename GeoField::GeometricBoundaryField& bfld =
            const_cast<typename GeoField::GeometricBoundaryField&>
            (
                fld.boundaryField()
            );

        bfld[patchI] == value;
    }
}


// Remove last patch field
template<class GeoField>
void Foam::fvMeshTools::trimPatchFields(fvMesh& mesh, const label nPatches)
{
    HashTable<const GeoField*> flds
    (
        mesh.objectRegistry::lookupClass<GeoField>()
    );

    forAllConstIter(typename HashTable<const GeoField*>, flds, iter)
    {
        const GeoField& fld = *iter();

        const_cast<typename GeoField::GeometricBoundaryField&>
        (
            fld.boundaryField()
        ).setSize(nPatches);
    }
}


// Reorder patch field
template<class GeoField>
void Foam::fvMeshTools::reorderPatchFields
(
    fvMesh& mesh,
    const labelList& oldToNew
)
{
    HashTable<const GeoField*> flds
    (
        mesh.objectRegistry::lookupClass<GeoField>()
    );

    forAllConstIter(typename HashTable<const GeoField*>, flds, iter)
    {
        const GeoField& fld = *iter();

        typename GeoField::GeometricBoundaryField& bfld =
            const_cast<typename GeoField::GeometricBoundaryField&>
            (
                fld.boundaryField()
            );
        bfld.reorder(oldToNew);
    }
}


// ************************************************************************* //
