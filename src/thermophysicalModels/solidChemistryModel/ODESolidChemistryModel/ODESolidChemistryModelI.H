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

#include "volFields.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CompType, class SolidThermo, class GasThermo>
inline Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >&
Foam::ODESolidChemistryModel<CompType, SolidThermo, GasThermo>::RRs()
{
    return RRs_;
}


template<class CompType, class SolidThermo, class GasThermo>
inline Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >&
Foam::ODESolidChemistryModel<CompType, SolidThermo, GasThermo>::RRg()
{
    return RRg_;
}


template<class CompType, class SolidThermo, class GasThermo>
inline const Foam::PtrList<Foam::Reaction<SolidThermo> >&
Foam::ODESolidChemistryModel<CompType, SolidThermo,GasThermo>::reactions() const
{
    return reactions_;
}


template<class CompType, class SolidThermo, class GasThermo>
inline const Foam::PtrList<GasThermo>&
Foam::ODESolidChemistryModel<CompType, SolidThermo, GasThermo>::
gasThermo() const
{
    return gasThermo_;
}


template<class CompType, class SolidThermo, class GasThermo>
inline const Foam::speciesTable&
Foam::ODESolidChemistryModel<CompType, SolidThermo, GasThermo>::gasTable() const
{
    return pyrolisisGases_;
}


template<class CompType, class SolidThermo, class GasThermo>
inline Foam::label
Foam::ODESolidChemistryModel<CompType, SolidThermo, GasThermo>::nSpecie() const
{
    return nSpecie_;
}


template<class CompType, class SolidThermo, class GasThermo>
inline Foam::label
Foam::ODESolidChemistryModel<CompType, SolidThermo, GasThermo>::
nReaction() const
{
    return nReaction_;
}


template<class CompType, class SolidThermo, class GasThermo>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ODESolidChemistryModel<CompType, SolidThermo, GasThermo>::RRs
(
    const label i
) const
{
    return RRs_[i];
}


template<class CompType, class SolidThermo, class GasThermo>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ODESolidChemistryModel<CompType, SolidThermo, GasThermo>::RRg
(
    const label i
) const
{
    return RRg_[i];
}


template<class CompType, class SolidThermo, class GasThermo>
inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::ODESolidChemistryModel<CompType, SolidThermo, GasThermo>::RRg() const
{
    tmp<DimensionedField<scalar, volMesh> > tRRg
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "RRg",
                this->time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
        )
    );

    if (this->chemistry_)
    {
        DimensionedField<scalar, volMesh>& RRg = tRRg();
        for (label i=0; i < nGases_; i++)
        {
            RRg += RRg_[i];
        }
    }
    return tRRg;
}


template<class CompType, class SolidThermo, class GasThermo>
inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::ODESolidChemistryModel<CompType, SolidThermo, GasThermo>::RRs() const
{
    tmp<DimensionedField<scalar, volMesh> > tRRs
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "RRs",
                this->time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
        )
    );

    if (this->chemistry_)
    {
        DimensionedField<scalar, volMesh>& RRs = tRRs();
        for (label i=0; i < nSolids_; i++)
        {
            RRs += RRs_[i];
        }
    }
    return tRRs;
}


template<class CompType, class SolidThermo, class GasThermo>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ODESolidChemistryModel<CompType, SolidThermo, GasThermo>::RR
(
    const label i
) const
{
    notImplemented("ODESolidChemistryModel::RR(const label)");
    return (DimensionedField<scalar, volMesh>::null());
}


// ************************************************************************* //
