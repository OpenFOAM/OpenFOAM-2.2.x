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

#include "nearWallFields.H"
#include "wordReList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(nearWallFields, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nearWallFields::nearWallFields
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
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "nearWallFields::nearWallFields"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating."
            << endl;
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nearWallFields::~nearWallFields()
{
    if (debug)
    {
        Info<< "nearWallFields::~nearWallFields()" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::nearWallFields::read(const dictionary& dict)
{
    if (debug)
    {
        Info<< "nearWallFields::read(const dictionary&)" << endl;
    }

    if (active_)
    {
        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        dict.lookup("fields") >> fieldSet_;
        patchSet_ =
            mesh.boundaryMesh().patchSet(wordReList(dict.lookup("patches")));
        distance_ = readScalar(dict.lookup("distance"));


        // Clear out any previously loaded fields
        vsf_.clear();
        vvf_.clear();
        vSpheretf_.clear();
        vSymmtf_.clear();
        vtf_.clear();
        fieldMap_.clear();
        reverseFieldMap_.clear();


        // Generate fields with mappedField boundary condition

        // Convert field to map
        fieldMap_.resize(2*fieldSet_.size());
        reverseFieldMap_.resize(2*fieldSet_.size());
        forAll(fieldSet_, setI)
        {
            const word& fldName = fieldSet_[setI].first();
            const word& sampleFldName = fieldSet_[setI].second();

            fieldMap_.insert(fldName, sampleFldName);
            reverseFieldMap_.insert(sampleFldName, fldName);
        }

        Info<< "Creating " << fieldMap_.size() << " fields" << endl;
        createFields(vsf_);
        createFields(vvf_);
        createFields(vSpheretf_);
        createFields(vSymmtf_);
        createFields(vtf_);
    }
}


void Foam::nearWallFields::execute()
{
    if (debug)
    {
        Info<< "nearWallFields:execute()" << endl;
    }

    //if (active_)
    //{
    //    sampleFields(vsf_);
    //    sampleFields(vvf_);
    //    sampleFields(vSpheretf_);
    //    sampleFields(vSymmtf_);
    //    sampleFields(vtf_);
    //}
}


void Foam::nearWallFields::end()
{
    if (debug)
    {
        Info<< "nearWallFields:end()" << endl;
    }
    // Update fields
    execute();
}


void Foam::nearWallFields::write()
{
    if (debug)
    {
        Info<< "nearWallFields:write()" << endl;
    }

    // Do nothing
    if (active_)
    {
        Info<< "Writing sampled fields to " << obr_.time().timeName()
            << endl;

        sampleFields(vsf_);
        sampleFields(vvf_);
        sampleFields(vSpheretf_);
        sampleFields(vSymmtf_);
        sampleFields(vtf_);

        // Write fields
        forAll(vsf_, i)
        {
            vsf_[i].write();
        }
        forAll(vvf_, i)
        {
            vvf_[i].write();
        }
        forAll(vSpheretf_, i)
        {
            vSpheretf_[i].write();
        }
        forAll(vSymmtf_, i)
        {
            vSymmtf_[i].write();
        }
        forAll(vtf_, i)
        {
            vtf_[i].write();
        }
    }
}


// ************************************************************************* //
