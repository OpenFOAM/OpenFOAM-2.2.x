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
    Foam::findCellData

Description

SourceFiles
    findCellDataI.H
    findCellData.C

\*---------------------------------------------------------------------------*/

#ifndef findCellData_H
#define findCellData_H

#include "Tuple2.H"
#include "point.H"
#include "tensor.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;
class polyPatch;
class globalIndex;

/*---------------------------------------------------------------------------*\
                         Class findCellData Declaration
\*---------------------------------------------------------------------------*/

typedef Tuple2<label, point> faceDataType;

class findCellData
{
public:

    //- Class used to pass additional data in
    class trackData
    {
    public:

        const globalIndex& globalCells_;
        const scalar maxDistSqr_;
        labelListList& found_;

        trackData
        (
            const globalIndex& globalCells,
            const scalar maxDistSqr,
            labelListList& found
        )
        :
            globalCells_(globalCells),
            maxDistSqr_(maxDistSqr),
            found_(found)
        {}

    };

private:

    // Private data

        List<faceDataType> data_;

public:

    // Public data


    // Constructors

        //- Construct null
        inline findCellData();

        //- Construct from components
        inline findCellData(const List<faceDataType>& data);


    // Member Functions

        const List<faceDataType>& data() const
        {
            return data_;
        }


        // Needed by FaceCellWave

            //- Check whether origin has been changed at all or
            //  still contains original (invalid) value.
            template<class TrackingData>
            inline bool valid(TrackingData& td) const;

            //- Check for identical geometrical data. Used for cyclics checking.
            template<class TrackingData>
            inline bool sameGeometry
            (
                const polyMesh&,
                const findCellData&,
                const scalar,
                TrackingData& td
            ) const;

            //- Convert any absolute coordinates into relative to (patch)face
            //  centre
            template<class TrackingData>
            inline void leaveDomain
            (
                const polyMesh&,
                const polyPatch&,
                const label patchFaceI,
                const point& faceCentre,
                TrackingData& td
            );

            //- Reverse of leaveDomain
            template<class TrackingData>
            inline void enterDomain
            (
                const polyMesh&,
                const polyPatch&,
                const label patchFaceI,
                const point& faceCentre,
                TrackingData& td
            );

            //- Apply rotation matrix to any coordinates
            template<class TrackingData>
            inline void transform
            (
                const polyMesh&,
                const tensor&,
                TrackingData& td
            );

            //- Influence of neighbouring face.
            //  Calls update(...) with cellCentre of cellI
            template<class TrackingData>
            inline bool updateCell
            (
                const polyMesh& mesh,
                const label thisCellI,
                const label neighbourFaceI,
                const findCellData& neighbourWallInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of neighbouring cell.
            //  Calls update(...) with faceCentre of faceI
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh& mesh,
                const label thisFaceI,
                const label neighbourCellI,
                const findCellData& neighbourWallInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same face.
            //  Merge new and old info.
            //  Calls update(...) with faceCentre of faceI
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh& mesh,
                const label thisFaceI,
                const findCellData& neighbourWallInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Same (like operator==)
            template<class TrackingData>
            inline bool equal(const findCellData&, TrackingData& td) const;


    // Member Operators

        // Needed for List IO
        inline bool operator==(const findCellData&) const;
        inline bool operator!=(const findCellData&) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const findCellData&);
        friend Istream& operator>>(Istream&, findCellData&);
};


//- Data associated with findCellData not contiguous
template<>
inline bool contiguous<findCellData>()
{
    return false;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "findCellDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
