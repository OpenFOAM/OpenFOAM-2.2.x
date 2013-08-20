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

#include "phaseModel.H"
#include "diameterModel.H"
#include "fixedValueFvPatchFields.H"
#include "slipFvPatchFields.H"
#include "partialSlipFvPatchFields.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseModel::phaseModel
(
    const fvMesh& mesh,
    const dictionary& transportProperties,
    const word& phaseName
)
:
    volScalarField
    (
        IOobject
        (
            "alpha" + phaseName,
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("alpha", dimless, 0)
    ),
    name_(phaseName),
    phaseDict_
    (
        transportProperties.subDict
        (
            phaseName == "1" || phaseName == "2"
          ? "phase" + phaseName
          : word(phaseName)
        )
    ),
    thermo_(rhoThermo::New(mesh, phaseName)),
    U_
    (
        IOobject
        (
            "U" + phaseName,
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
{
    thermo_->validate("phaseModel " + phaseName, "h", "e");

    const word phiName = "phi" + phaseName;

    IOobject phiHeader
    (
        phiName,
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ
    );

    if (phiHeader.headerOk())
    {
        Info<< "Reading face flux field " << phiName << endl;

        phiPtr_.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    phiName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
    }
    else
    {
        Info<< "Calculating face flux field " << phiName << endl;

        wordList phiTypes
        (
            U_.boundaryField().size(),
            calculatedFvPatchScalarField::typeName
        );

        forAll(U_.boundaryField(), i)
        {
            if
            (
                isA<fixedValueFvPatchVectorField>(U_.boundaryField()[i])
             || isA<slipFvPatchVectorField>(U_.boundaryField()[i])
             || isA<partialSlipFvPatchVectorField>(U_.boundaryField()[i])
            )
            {
                phiTypes[i] = fixedValueFvPatchScalarField::typeName;
            }
        }

        phiPtr_.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    phiName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fvc::interpolate(U_) & mesh.Sf(),
                phiTypes
            )
        );
    }

    dPtr_ = diameterModel::New
    (
        phaseDict_,
        *this
    );
}


Foam::autoPtr<Foam::phaseModel> Foam::phaseModel::New
(
    const fvMesh& mesh,
    const dictionary& transportProperties,
    const word& phaseName
)
{
    return autoPtr<phaseModel>
    (
        new phaseModel(mesh, transportProperties, phaseName)
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseModel::~phaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::phaseModel::d() const
{
    return dPtr_().d();
}


// ************************************************************************* //
