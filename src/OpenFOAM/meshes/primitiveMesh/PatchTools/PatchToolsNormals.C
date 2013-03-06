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

#include "PatchTools.H"
#include "polyMesh.H"
#include "indirectPrimitivePatch.H"
#include "globalMeshData.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    //- Transformation
    class listTransform
    {
    public:

        void operator()
        (
            const vectorTensorTransform& vt,
            const bool forward,
            List<List<point> >& fld
        ) const
        {
            const tensor T
            (
                forward
              ? vt.R()
              : vt.R().T()
            );

            forAll(fld, i)
            {
                List<point>& elems = fld[i];
                forAll(elems, elemI)
                {
                    elems[elemI] = transform(T, elems[elemI]);
                }
            }
        }
    };
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>

Foam::tmp<Foam::pointField>
Foam::PatchTools::pointNormals
(
    const polyMesh& mesh,
    const PrimitivePatch<Face, FaceList, PointField, PointType>& p,
    const labelList& meshFaces
)
{
    // Assume patch is smaller than the globalData().coupledPatch() (?) so
    // loop over patch meshPoints.

    const globalMeshData& globalData = mesh.globalData();
    const indirectPrimitivePatch& coupledPatch = globalData.coupledPatch();
    const Map<label>& coupledPatchMP = coupledPatch.meshPointMap();
    const mapDistribute& map = globalData.globalPointSlavesMap();
    const globalIndexAndTransform& transforms =
        globalData.globalTransforms();


    // 1. Start off with local normals (note:without calculating pointNormals
    //    to avoid them being stored)

    tmp<pointField> textrudeN(new pointField(p.nPoints(), vector::zero));
    pointField& extrudeN = textrudeN();
    {
        const faceList& localFaces = p.localFaces();
        const vectorField& faceNormals = p.faceNormals();

        forAll(localFaces, faceI)
        {
            const face& f = localFaces[faceI];
            const vector& n = faceNormals[faceI];
            forAll(f, fp)
            {
                extrudeN[f[fp]] += n;
            }
        }
        extrudeN /= mag(extrudeN)+VSMALL;
    }


    // Collect local pointFaces
    List<List<point> > pointFaceNormals(map.constructSize());
    forAll(p.meshPoints(), patchPointI)
    {
        label meshPointI = p.meshPoints()[patchPointI];
        Map<label>::const_iterator fnd = coupledPatchMP.find(meshPointI);
        if (fnd != coupledPatchMP.end())
        {
            label coupledPointI = fnd();

            List<point>& pNormals = pointFaceNormals[coupledPointI];
            const labelList& pFaces = p.pointFaces()[patchPointI];
            pNormals.setSize(pFaces.size());
            forAll(pFaces, i)
            {
                pNormals[i] = p.faceNormals()[pFaces[i]];
            }
        }
    }


    // Pull remote data into local slots
    map.distribute
    (
        transforms,
        pointFaceNormals,
        listTransform()
    );


    // Combine normals
    const labelListList& slaves = globalData.globalPointSlaves();
    const labelListList& transformedSlaves =
        globalData.globalPointTransformedSlaves();


    pointField coupledPointNormals(map.constructSize(), vector::zero);

    forAll(p.meshPoints(), patchPointI)
    {
        label meshPointI = p.meshPoints()[patchPointI];
        Map<label>::const_iterator fnd = coupledPatchMP.find(meshPointI);
        if (fnd != coupledPatchMP.end())
        {
            label coupledPointI = fnd();
            const labelList& slaveSlots =
                slaves[coupledPointI];
            const labelList& transformedSlaveSlots =
                transformedSlaves[coupledPointI];

            label nFaces = slaveSlots.size()+transformedSlaveSlots.size();
            if (nFaces > 0)
            {
                // Combine
                point& n = coupledPointNormals[coupledPointI];

                n += sum(pointFaceNormals[coupledPointI]);

                forAll(slaveSlots, i)
                {
                    n += sum(pointFaceNormals[slaveSlots[i]]);
                }
                forAll(transformedSlaveSlots, i)
                {
                    n += sum(pointFaceNormals[transformedSlaveSlots[i]]);
                }
                n /= mag(n)+VSMALL;

                // Put back into slave slots
                forAll(slaveSlots, i)
                {
                    coupledPointNormals[slaveSlots[i]] = n;
                }
                forAll(transformedSlaveSlots, i)
                {
                    coupledPointNormals[transformedSlaveSlots[i]] = n;
                }
            }
        }
    }


    // Send back
    map.reverseDistribute
    (
        transforms,
        coupledPointNormals.size(),
        coupledPointNormals,
        mapDistribute::transform()
    );


    // Override patch normals
    forAll(p.meshPoints(), patchPointI)
    {
        label meshPointI = p.meshPoints()[patchPointI];
        Map<label>::const_iterator fnd = coupledPatchMP.find(meshPointI);
        if (fnd != coupledPatchMP.end())
        {
            label coupledPointI = fnd();
            extrudeN[patchPointI] = coupledPointNormals[coupledPointI];
        }
    }

    return textrudeN;
}


// ************************************************************************* //
