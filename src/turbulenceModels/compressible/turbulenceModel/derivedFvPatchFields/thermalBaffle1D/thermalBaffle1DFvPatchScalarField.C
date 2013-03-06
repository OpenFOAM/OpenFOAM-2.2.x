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
#include "surfaceFields.H"
#include "mappedPatchBase.H"
#include "turbulenceModel.H"
#include "mapDistribute.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class solidType>
thermalBaffle1DFvPatchScalarField<solidType>::
thermalBaffle1DFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    TName_("T"),
    baffleActivated_(true),
    thickness_(p.size()),
    Qs_(p.size()),
    solidPtr_(NULL)
{}


template<class solidType>
thermalBaffle1DFvPatchScalarField<solidType>::
thermalBaffle1DFvPatchScalarField
(
    const thermalBaffle1DFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    TName_(ptf.TName_),
    baffleActivated_(ptf.baffleActivated_),
    thickness_(ptf.thickness_),
    Qs_(ptf.Qs_),
    solidPtr_(ptf.solidPtr_)
{}


template<class solidType>
thermalBaffle1DFvPatchScalarField<solidType>::
thermalBaffle1DFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    TName_("T"),
    baffleActivated_(readBool(dict.lookup("baffleActivated"))),
    thickness_(scalarField("thickness", dict, p.size())),
    Qs_(scalarField("Qs", dict, p.size())),
    solidPtr_(new solidType(dict))
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorIn
        (
            "thermalBaffle1DFvPatchScalarField::"
            "thermalBaffle1DFvPatchScalarField"
            "("
                "const fvPatch&,\n"
                "const DimensionedField<scalar, volMesh>&, "
                "const dictionary&"
            ")"
        )   << "\n    patch type '" << patch().type()
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << patch().name()
            << " of field " << dimensionedInternalField().name()
            << " in file " << dimensionedInternalField().objectPath()
            << exit(FatalError);
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue") && baffleActivated_)
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume zeroGradient.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 0.0;
    }

}


template<class solidType>
thermalBaffle1DFvPatchScalarField<solidType>::
thermalBaffle1DFvPatchScalarField
(
    const thermalBaffle1DFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    TName_(ptf.TName_),
    baffleActivated_(ptf.baffleActivated_),
    thickness_(ptf.thickness_),
    Qs_(ptf.Qs_),
    solidPtr_(ptf.solidPtr_)
{}


template<class solidType>
thermalBaffle1DFvPatchScalarField<solidType>::
thermalBaffle1DFvPatchScalarField
(
    const thermalBaffle1DFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    TName_(ptf.TName_),
    baffleActivated_(ptf.baffleActivated_),
    thickness_(ptf.thickness_),
    Qs_(ptf.Qs_),
    solidPtr_(ptf.solidPtr_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class solidType>
void thermalBaffle1DFvPatchScalarField<solidType>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
}

template<class solidType>
void thermalBaffle1DFvPatchScalarField<solidType>::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);
}


template<class solidType>
void thermalBaffle1DFvPatchScalarField<solidType>::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());

    const label patchI = patch().index();

    const label nbrPatchI = mpp.samplePolyPatch().index();

    if (baffleActivated_)
    {
        const fvPatch& nbrPatch = patch().boundaryMesh()[nbrPatchI];

        const compressible::turbulenceModel& model =
            db().template lookupObject<compressible::turbulenceModel>
            (
                "turbulenceModel"
            );


        // local properties
        const scalarField kappaw = model.kappaEff()().boundaryField()[patchI];

        const fvPatchScalarField& Tp =
            patch().template lookupPatchField<volScalarField, scalar>(TName_);

        const scalarField qDot(kappaw*Tp.snGrad());

        tmp<scalarField> Ti = patchInternalField();

        scalarField myh(patch().deltaCoeffs()*kappaw);

        // nbr properties

        scalarField nbrKappaw =
            model.kappaEff()().boundaryField()[nbrPatchI];
        mpp.map().distribute(nbrKappaw);

        const fvPatchScalarField& nbrTw =
            model.thermo().T().boundaryField()[nbrPatchI];

        scalarField nbrQDot
        (
            model.kappaEff()().boundaryField()[nbrPatchI]*nbrTw.snGrad()
        );
        mpp.map().distribute(nbrQDot);

        const thermalBaffle1DFvPatchScalarField& nbrField =
        refCast<const thermalBaffle1DFvPatchScalarField>
        (
            nbrPatch.template lookupPatchField<volScalarField, scalar>(TName_)
        );

        scalarField nbrTi(nbrField.patchInternalField());
        mpp.map().distribute(nbrTi);

        scalarField nbrTp =
           nbrPatch.template lookupPatchField<volScalarField, scalar>(TName_);
        mpp.map().distribute(nbrTp);

        scalarField nbrh
        (
            nbrPatch.deltaCoeffs()
           *model.kappaEff()().boundaryField()[nbrPatchI]
        );
        mpp.map().distribute(nbrh);


        // heat source
        const scalarField Q(Qs_/thickness_);

        tmp<scalarField> tKDeltaw(new scalarField(patch().size()));
        scalarField KDeltaw = tKDeltaw();

        // Create fields for solid properties (p paramater not used)
        forAll(KDeltaw, i)
        {
            KDeltaw[i] =
                solidPtr_().kappa(0.0, (Tp[i] + nbrTw[i])/2.0)/thickness_[i];
        }

        const scalarField q
        (
            (Ti() - nbrTi)/(1.0/KDeltaw + 1.0/nbrh + 1.0/myh)
        );

        forAll(qDot, i)
        {
            if (Qs_[i] == 0)
            {
                this->refValue()[i] = Ti()[i] - q[i]/myh[i];
                this->refGrad()[i] = 0.0;
                this->valueFraction()[i] = 1.0;
            }
            else
            {
                if (q[i] > 0)
                {
                    this->refValue()[i] =
                        nbrTp[i]
                      - Q[i]*thickness_[i]/(2*KDeltaw[i]);

                    this->refGrad()[i] = 0.0;
                    this->valueFraction()[i] =
                        1.0
                        /
                        (
                            1.0
                          + patch().deltaCoeffs()[i]*kappaw[i]/KDeltaw[i]
                        );
                }
                else if (q[i] < 0)
                {
                    this->refValue()[i] = 0.0;
                    this->refGrad()[i] =
                          (-nbrQDot[i] + Q[i]*thickness_[i])/kappaw[i];
                    this->valueFraction()[i] = 0.0;
                }
                else
                {
                    scalar Qt = Q[i]*thickness_[i];
                    this->refValue()[i] = 0.0;
                    this->refGrad()[i] = Qt/2/kappaw[i];
                    this->valueFraction()[i] = 0.0;
                }
            }
        }

        if (debug)
        {
            scalar Q = gSum(patch().magSf()*qDot);
            Info<< patch().boundaryMesh().mesh().name() << ':'
                << patch().name() << ':'
                << this->dimensionedInternalField().name() << " <- "
                << nbrPatch.name() << ':'
                << this->dimensionedInternalField().name() << " :"
                << " heat[W]:" << Q
                << " walltemperature "
                << " min:" << gMin(*this)
                << " max:" << gMax(*this)
                << " avg:" << gAverage(*this)
                << endl;
        }
    }

    // Restore tag
    UPstream::msgType() = oldTag;

    mixedFvPatchScalarField::updateCoeffs();
}

template <class solidType>
void thermalBaffle1DFvPatchScalarField<solidType>::
write(Ostream& os) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("TName")
        << TName_ << token::END_STATEMENT << nl;
    thickness_.writeEntry("thickness", os);
    os.writeKeyword("baffleActivated")
        << baffleActivated_ << token::END_STATEMENT << nl;
    Qs_.writeEntry("Qs", os);
    solidPtr_->write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
