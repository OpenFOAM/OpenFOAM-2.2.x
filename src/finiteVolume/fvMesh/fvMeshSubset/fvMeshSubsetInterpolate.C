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

\*---------------------------------------------------------------------------*/

#include "fvMeshSubset.H"
#include "emptyFvsPatchField.H"
#include "emptyPointPatchField.H"
#include "emptyFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> > fvMeshSubset::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const fvMesh& sMesh,
    const labelList& patchMap,
    const labelList& cellMap,
    const labelList& faceMap
)
{
    // Create and map the internal-field values
    Field<Type> internalField(vf.internalField(), cellMap);

    // Create and map the patch field values
    PtrList<fvPatchField<Type> > patchFields(patchMap.size());

    forAll(patchFields, patchI)
    {
        // Set the first one by hand as it corresponds to the
        // exposed internal faces.  Additional interpolation can be put here
        // as necessary.
        if (patchMap[patchI] == -1)
        {
            patchFields.set
            (
                patchI,
                new emptyFvPatchField<Type>
                (
                    sMesh.boundary()[patchI],
                    DimensionedField<Type, volMesh>::null()
                )
            );
        }
        else
        {
            // Construct addressing
            const fvPatch& subPatch = sMesh.boundary()[patchI];
            const fvPatch& basePatch = vf.mesh().boundary()[patchMap[patchI]];
            const label baseStart = basePatch.start();
            const label baseSize = basePatch.size();

            labelList directAddressing(subPatch.size());

            forAll(directAddressing, i)
            {
                label baseFaceI = faceMap[subPatch.start()+i];

                if (baseFaceI >= baseStart && baseFaceI < baseStart+baseSize)
                {
                    directAddressing[i] = baseFaceI-baseStart;
                }
                else
                {
                    // Mapped from internal face. Do what? Map from element
                    // 0 for now.
                    directAddressing[i] = 0;
                }
            }

            patchFields.set
            (
                patchI,
                fvPatchField<Type>::New
                (
                    vf.boundaryField()[patchMap[patchI]],
                    sMesh.boundary()[patchI],
                    DimensionedField<Type, volMesh>::null(),
                    patchFieldSubset(directAddressing)
                )
            );

            // What to do with exposed internal faces if put into this patch?
        }
    }


    // Create the complete field from the pieces
    tmp<GeometricField<Type, fvPatchField, volMesh> > tresF
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "subset"+vf.name(),
                sMesh.time().timeName(),
                sMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            sMesh,
            vf.dimensions(),
            internalField,
            patchFields
        )
    );

    return tresF;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> > fvMeshSubset::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    return interpolate
    (
        vf,
        subMesh(),
        patchMap(),
        cellMap(),
        faceMap()
    );
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > fvMeshSubset::interpolate
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& vf,
    const fvMesh& sMesh,
    const labelList& patchMap,
    const labelList& faceMap
)
{
    // Create and map the internal-field values
    Field<Type> internalField
    (
        vf.internalField(),
        SubList<label>
        (
            faceMap,
            sMesh.nInternalFaces()
        )
    );

    // Create and map the patch field values
    PtrList<fvsPatchField<Type> > patchFields(patchMap.size());

    forAll(patchFields, patchI)
    {
        // Set the first one by hand as it corresponds to the
        // exposed internal faces.  Additional interpolation can be put here
        // as necessary.
        if (patchMap[patchI] == -1)
        {
            patchFields.set
            (
                patchI,
                new emptyFvsPatchField<Type>
                (
                    sMesh.boundary()[patchI],
                    DimensionedField<Type, surfaceMesh>::null()
                )
            );
        }
        else
        {
            // Construct addressing
            const fvPatch& subPatch = sMesh.boundary()[patchI];
            const fvPatch& basePatch = vf.mesh().boundary()[patchMap[patchI]];
            const label baseStart = basePatch.start();
            const label baseSize = basePatch.size();

            labelList directAddressing(subPatch.size());

            forAll(directAddressing, i)
            {
                label baseFaceI = faceMap[subPatch.start()+i];

                if (baseFaceI >= baseStart && baseFaceI < baseStart+baseSize)
                {
                    directAddressing[i] = baseFaceI-baseStart;
                }
                else
                {
                    // Mapped from internal face. Do what? Map from element
                    // 0 for now.
                    directAddressing[i] = 0;
                }
            }

            patchFields.set
            (
                patchI,
                fvsPatchField<Type>::New
                (
                    vf.boundaryField()[patchMap[patchI]],
                    sMesh.boundary()[patchI],
                    DimensionedField<Type, surfaceMesh>::null(),
                    patchFieldSubset(directAddressing)
                )
            );
        }
    }


    // Map exposed internal faces. Note: Only necessary if exposed faces added
    // into existing patch but since we don't know that at this point...
    forAll(patchFields, patchI)
    {
        fvsPatchField<Type>& pfld = patchFields[patchI];

        label meshFaceI = pfld.patch().start();

        forAll(pfld, i)
        {
            label oldFaceI = faceMap[meshFaceI++];

            if (oldFaceI < vf.internalField().size())
            {
                pfld[i] = vf.internalField()[oldFaceI];
            }
        }
    }

    // Create the complete field from the pieces
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tresF
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "subset"+vf.name(),
                sMesh.time().timeName(),
                sMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            sMesh,
            vf.dimensions(),
            internalField,
            patchFields
        )
    );

    return tresF;
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > fvMeshSubset::interpolate
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sf
) const
{
    return interpolate
    (
        sf,
        subMesh(),
        patchMap(),
        faceMap()
    );
}


template<class Type>
tmp<GeometricField<Type, pointPatchField, pointMesh> >
fvMeshSubset::interpolate
(
    const GeometricField<Type, pointPatchField, pointMesh>& vf,
    const pointMesh& sMesh,
    const labelList& patchMap,
    const labelList& pointMap
)
{
    // Create and map the internal-field values
    Field<Type> internalField(vf.internalField(), pointMap);

    // Create and map the patch field values
    PtrList<pointPatchField<Type> > patchFields(patchMap.size());

    forAll(patchFields, patchI)
    {
        // Set the first one by hand as it corresponds to the
        // exposed internal faces.  Additional interpolation can be put here
        // as necessary.
        if (patchMap[patchI] == -1)
        {
            patchFields.set
            (
                patchI,
                new emptyPointPatchField<Type>
                (
                    sMesh.boundary()[patchI],
                    DimensionedField<Type, pointMesh>::null()
                )
            );
        }
        else
        {
            // Construct addressing
            const pointPatch& basePatch =
                vf.mesh().boundary()[patchMap[patchI]];

            const labelList& meshPoints = basePatch.meshPoints();

            // Make addressing from mesh to patch point
            Map<label> meshPointMap(2*meshPoints.size());
            forAll(meshPoints, localI)
            {
                meshPointMap.insert(meshPoints[localI], localI);
            }

            // Find which subpatch points originate from which patch point
            const pointPatch& subPatch = sMesh.boundary()[patchI];
            const labelList& subMeshPoints = subPatch.meshPoints();

            // If mapped from outside patch use point 0 for lack of better.
            labelList directAddressing(subPatch.size(), 0);

            forAll(subMeshPoints, localI)
            {
                // Get mesh point on original mesh.
                label meshPointI = pointMap[subMeshPoints[localI]];

                Map<label>::const_iterator iter = meshPointMap.find(meshPointI);

                if (iter != meshPointMap.end())
                {
                    directAddressing[localI] = iter();
                }
            }

            patchFields.set
            (
                patchI,
                pointPatchField<Type>::New
                (
                    vf.boundaryField()[patchMap[patchI]],
                    subPatch,
                    DimensionedField<Type, pointMesh>::null(),
                    pointPatchFieldSubset(directAddressing)
                )
            );
        }
    }

    // Create the complete field from the pieces
    tmp<GeometricField<Type, pointPatchField, pointMesh> > tresF
    (
        new GeometricField<Type, pointPatchField, pointMesh>
        (
            IOobject
            (
                "subset"+vf.name(),
                vf.time().timeName(),
                sMesh.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            sMesh,
            vf.dimensions(),
            internalField,
            patchFields
        )
    );

    return tresF;
}


template<class Type>
tmp<GeometricField<Type, pointPatchField, pointMesh> > fvMeshSubset::interpolate
(
    const GeometricField<Type, pointPatchField, pointMesh>& sf
) const
{
    return interpolate
    (
        sf,
        pointMesh::New(subMesh()),     // subsetted point mesh
        patchMap(),
        pointMap()
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
