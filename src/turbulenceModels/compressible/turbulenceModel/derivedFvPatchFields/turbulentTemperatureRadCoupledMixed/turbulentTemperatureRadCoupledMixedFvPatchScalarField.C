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

#include "turbulentTemperatureRadCoupledMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentTemperatureRadCoupledMixedFvPatchScalarField::
turbulentTemperatureRadCoupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), "undefined", "undefined-K"),
    TnbrName_("undefined-Tnbr"),
    QrNbrName_("undefined-QrNbr"),
    QrName_("undefined-Qr")
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


turbulentTemperatureRadCoupledMixedFvPatchScalarField::
turbulentTemperatureRadCoupledMixedFvPatchScalarField
(
    const turbulentTemperatureRadCoupledMixedFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(psf, p, iF, mapper),
    temperatureCoupledBase(patch(), psf.KMethod(), psf.kappaName()),
    TnbrName_(psf.TnbrName_),
    QrNbrName_(psf.QrNbrName_),
    QrName_(psf.QrName_)
{}


turbulentTemperatureRadCoupledMixedFvPatchScalarField::
turbulentTemperatureRadCoupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), dict),
    TnbrName_(dict.lookup("Tnbr")),
    QrNbrName_(dict.lookup("QrNbr")),
    QrName_(dict.lookup("Qr"))
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorIn
        (
            "turbulentTemperatureRadCoupledMixedFvPatchScalarField::"
            "turbulentTemperatureRadCoupledMixedFvPatchScalarField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<scalar, volMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << dimensionedInternalField().name()
            << " in file " << dimensionedInternalField().objectPath()
            << exit(FatalError);
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 1.0;
    }
}


turbulentTemperatureRadCoupledMixedFvPatchScalarField::
turbulentTemperatureRadCoupledMixedFvPatchScalarField
(
    const turbulentTemperatureRadCoupledMixedFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(psf, iF),
    temperatureCoupledBase(patch(), psf.KMethod(), psf.kappaName()),
    TnbrName_(psf.TnbrName_),
    QrNbrName_(psf.QrNbrName_),
    QrName_(psf.QrName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbulentTemperatureRadCoupledMixedFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());
    const polyMesh& nbrMesh = mpp.sampleMesh();
    const label samplePatchI = mpp.samplePolyPatch().index();
    const fvPatch& nbrPatch =
        refCast<const fvMesh>(nbrMesh).boundary()[samplePatchI];


    scalarField Tc(patchInternalField());
    scalarField& Tp = *this;

    const turbulentTemperatureRadCoupledMixedFvPatchScalarField&
        nbrField = refCast
            <const turbulentTemperatureRadCoupledMixedFvPatchScalarField>
            (
                nbrPatch.lookupPatchField<volScalarField, scalar>(TnbrName_)
            );

    // Swap to obtain full local values of neighbour internal field
    scalarField TcNbr(nbrField.patchInternalField());
    mpp.distribute(TcNbr);


    // Swap to obtain full local values of neighbour K*delta
    scalarField KDeltaNbr(nbrField.kappa(nbrField)*nbrPatch.deltaCoeffs());
    mpp.distribute(KDeltaNbr);


    scalarField KDelta(kappa(*this)*patch().deltaCoeffs());

    scalarField Qr(Tp.size(), 0.0);
    if (QrName_ != "none")
    {
        Qr = patch().lookupPatchField<volScalarField, scalar>(QrName_);
    }

    scalarField QrNbr(Tp.size(), 0.0);
    if (QrNbrName_ != "none")
    {
        QrNbr = nbrPatch.lookupPatchField<volScalarField, scalar>(QrNbrName_);
        mpp.distribute(QrNbr);
    }

    scalarField alpha(KDeltaNbr - (Qr + QrNbr)/Tp);

    valueFraction() = alpha/(alpha + KDelta);

    refValue() = (KDeltaNbr*TcNbr)/alpha;

    // Restore tag
    UPstream::msgType() = oldTag;

    mixedFvPatchScalarField::updateCoeffs();
}


void turbulentTemperatureRadCoupledMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("Tnbr")<< TnbrName_ << token::END_STATEMENT << nl;
    os.writeKeyword("QrNbr")<< QrNbrName_ << token::END_STATEMENT << nl;
    os.writeKeyword("Qr")<< QrName_ << token::END_STATEMENT << nl;
    temperatureCoupledBase::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    turbulentTemperatureRadCoupledMixedFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //
