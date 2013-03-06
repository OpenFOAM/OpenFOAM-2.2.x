/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

#include "meshToMeshNew.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline const Foam::labelListList&
Foam::meshToMeshNew::srcToTgtCellAddr() const
{
    return srcToTgtCellAddr_;
}


inline const Foam::labelListList&
Foam::meshToMeshNew::tgtToSrcCellAddr() const
{
    return tgtToSrcCellAddr_;
}


inline const Foam::scalarListList&
Foam::meshToMeshNew::srcToTgtCellWght() const
{
    return srcToTgtCellWght_;
}


inline const Foam::scalarListList&
Foam::meshToMeshNew::tgtToSrcCellWght() const
{
    return tgtToSrcCellWght_;
}


inline Foam::scalar Foam::meshToMeshNew::V() const
{
    return V_;
}


// ************************************************************************* //
