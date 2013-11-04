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

\*---------------------------------------------------------------------------*/

#include "polyMesh.H"
#include "globalIndex.H"
#include "transform.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline findCellData::findCellData()
:
    data_()
{}


inline findCellData::findCellData(const List<faceDataType>& data)
:
    data_(data)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TrackingData>
inline bool findCellData::valid(TrackingData& td) const
{
    return data_.size();
}


template<class TrackingData>
inline bool findCellData::sameGeometry
(
    const polyMesh&,
    const findCellData& w2,
    const scalar tol,
    TrackingData& td
)
 const
{
    // Note: maybe add checking?
    return true;
}


template<class TrackingData>
inline void findCellData::leaveDomain
(
    const polyMesh&,
    const polyPatch&,
    const label,
    const point& faceCentre,
    TrackingData& td
)
{
    forAll(data_, i)
    {
        data_[i].second() -= faceCentre;
    }
}


template<class TrackingData>
inline void findCellData::transform
(
    const polyMesh&,
    const tensor& rotTensor,
    TrackingData& td
)
{
    forAll(data_, i)
    {
        const point& pt = data_[i].second();
        data_[i].second() = Foam::transform(rotTensor, pt);
    }
}


template<class TrackingData>
inline void findCellData::enterDomain
(
    const polyMesh&,
    const polyPatch&,
    const label,
    const point& faceCentre,
    TrackingData& td
)
{
    // back to absolute form
    forAll(data_, i)
    {
        data_[i].second() += faceCentre;
    }
}


template<class TrackingData>
inline bool findCellData::updateCell
(
    const polyMesh& mesh,
    const label thisCellI,
    const label neighbourFaceI,
    const findCellData& neighbourInfo,
    const scalar tol,
    TrackingData& td
)
{
    const point& thisCc = mesh.cellCentres()[thisCellI];
    const List<faceDataType>& nbrData = neighbourInfo.data();

    // Is any data changed in this cell?
    bool changed = false;

    forAll(nbrData, i)
    {
        const point& samplePt = nbrData[i].second();

        if (mesh.pointInCell(samplePt, thisCellI))
        {
            td.found_[thisCellI].append(nbrData[i].first());
        }
        else
        {
            // See if needs transferring:
            // - not yet present
            // - not too far
            if (findIndex(data_, nbrData[i]) == -1)
            {
                if (magSqr(samplePt-thisCc) < td.maxDistSqr_)
                {
                    data_.append(nbrData[i]);
                    changed = true;
                }
            }
        }
    }

    return changed;
}


// Update this with w2 if w2 nearer to pt.
template<class TrackingData>
inline bool findCellData::updateFace
(
    const polyMesh& mesh,
    const label thisFaceI,
    const label neighbourCellI,
    const findCellData& neighbourInfo,
    const scalar tol,
    TrackingData& td
)
{
    const point& thisFc = mesh.faceCentres()[thisFaceI];
    const List<faceDataType>& nbrData = neighbourInfo.data();

    // Is any data changed in this cell?
    bool changed = false;

    forAll(nbrData, i)
    {
        // See if needs transferring:
        // - not yet present
        // - not too far
        if (findIndex(data_, nbrData[i]) == -1)
        {
            const point& samplePt = nbrData[i].second();

            if (magSqr(samplePt-thisFc) < td.maxDistSqr_)
            {
                data_.append(nbrData[i]);
                changed = true;
            }
        }
    }

    return changed;
}


// Update this with w2 if w2 nearer to pt.
template<class TrackingData>
inline bool findCellData::updateFace
(
    const polyMesh& mesh,
    const label thisFaceI,
    const findCellData& neighbourInfo,
    const scalar tol,
    TrackingData& td
)
{
    const List<faceDataType>& nbrData = neighbourInfo.data();

    bool changed = false;

    forAll(nbrData, i)
    {
        if (findIndex(data_, nbrData[i]) == -1)
        {
            data_.append(nbrData[i]);
            changed = true;
        }
    }

    return changed;
}


template<class TrackingData>
inline bool findCellData::equal
(
    const findCellData& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool findCellData::operator==(const findCellData& rhs) const
{
    const List<faceDataType>& rhsData = rhs.data();

    if (data_.size() != rhsData.size())
    {
        return false;
    }

    forAll(data_, i)
    {
        if (data_[i].first() != rhsData[i].first())
        {
            return false;
        }
    }
    return true;
}


inline bool findCellData::operator!=(const findCellData& rhs) const
{
    return !(*this == rhs);
}


// ************************************************************************* //

} // End namespace Foam

// ************************************************************************* //
