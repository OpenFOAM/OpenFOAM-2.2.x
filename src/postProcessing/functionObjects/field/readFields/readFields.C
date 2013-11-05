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

#include "readFields.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(readFields, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::readFields::readFields
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    active_(true),
    fieldSet_()
{
    // Check if the available mesh is an fvMesh otherise deactivate
    if (isA<fvMesh>(obr_))
    {
        read(dict);
    }
    else
    {
        active_ = false;
        WarningIn
        (
            "readFields::readFields"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating " << name_
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::readFields::~readFields()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::readFields::read(const dictionary& dict)
{
    if (active_)
    {
        dict.lookup("fields") >> fieldSet_;
    }
}


void Foam::readFields::execute()
{
    //Info<< type() << " " << name_ << ":" << nl;

    // Clear out any previously loaded fields
    vsf_.clear();
    vvf_.clear();
    vSpheretf_.clear();
    vSymmtf_.clear();
    vtf_.clear();

    ssf_.clear();
    svf_.clear();
    sSpheretf_.clear();
    sSymmtf_.clear();
    stf_.clear();

    forAll(fieldSet_, fieldI)
    {
        // If necessary load field
        loadField<scalar>(fieldSet_[fieldI], vsf_, ssf_);
        loadField<vector>(fieldSet_[fieldI], vvf_, svf_);
        loadField<sphericalTensor>(fieldSet_[fieldI], vSpheretf_, sSpheretf_);
        loadField<symmTensor>(fieldSet_[fieldI], vSymmtf_, sSymmtf_);
        loadField<tensor>(fieldSet_[fieldI], vtf_, stf_);
    }
}


void Foam::readFields::end()
{
    execute();
}


void Foam::readFields::timeSet()
{
    // Do nothing
}


void Foam::readFields::write()
{
    // Do nothing
}


// ************************************************************************* //
