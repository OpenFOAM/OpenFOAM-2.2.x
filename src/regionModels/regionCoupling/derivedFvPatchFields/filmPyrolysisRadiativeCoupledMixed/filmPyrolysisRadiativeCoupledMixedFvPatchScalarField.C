/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenCFD Ltd.
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

#include "filmPyrolysisRadiativeCoupledMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const filmPyrolysisRadiativeCoupledMixedFvPatchScalarField::filmModelType&
filmPyrolysisRadiativeCoupledMixedFvPatchScalarField::
filmModel() const
{
    const regionModels::regionModel& model =
        db().time().lookupObject<regionModels::regionModel>
        (
            "surfaceFilmProperties"
        );

    return dynamic_cast<const filmModelType&>(model);
}


const filmPyrolysisRadiativeCoupledMixedFvPatchScalarField::
pyrolysisModelType&
filmPyrolysisRadiativeCoupledMixedFvPatchScalarField::
pyrModel() const
{
    const regionModels::regionModel& model =
        db().time().lookupObject<regionModels::regionModel>
        (
            "pyrolysisProperties"
        );

    return dynamic_cast<const pyrolysisModelType&>(model);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

filmPyrolysisRadiativeCoupledMixedFvPatchScalarField::
filmPyrolysisRadiativeCoupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), "undefined", "undefined-K"),
    TnbrName_("undefined-Tnbr"),
    QrNbrName_("undefined-QrNbr"),
    QrName_("undefined-Qr"),
    convectiveScaling_(1.0),
    filmDeltaDry_(0.0),
    filmDeltaWet_(0.0)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


filmPyrolysisRadiativeCoupledMixedFvPatchScalarField::
filmPyrolysisRadiativeCoupledMixedFvPatchScalarField
(
    const filmPyrolysisRadiativeCoupledMixedFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(psf, p, iF, mapper),
    temperatureCoupledBase(patch(), psf.KMethod(), psf.kappaName()),
    TnbrName_(psf.TnbrName_),
    QrNbrName_(psf.QrNbrName_),
    QrName_(psf.QrName_),
    convectiveScaling_(psf.convectiveScaling_),
    filmDeltaDry_(psf.filmDeltaDry_),
    filmDeltaWet_(psf.filmDeltaWet_)
{}


filmPyrolysisRadiativeCoupledMixedFvPatchScalarField::
filmPyrolysisRadiativeCoupledMixedFvPatchScalarField
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
    QrName_(dict.lookup("Qr")),
    convectiveScaling_(dict.lookupOrDefault<scalar>("convectiveScaling", 1.0)),
    filmDeltaDry_(readScalar(dict.lookup("filmDeltaDry"))),
    filmDeltaWet_(readScalar(dict.lookup("filmDeltaWet")))
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorIn
        (
            "filmPyrolysisRadiativeCoupledMixedFvPatchScalarField::"
            "filmPyrolysisRadiativeCoupledMixedFvPatchScalarField\n"
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


filmPyrolysisRadiativeCoupledMixedFvPatchScalarField::
filmPyrolysisRadiativeCoupledMixedFvPatchScalarField
(
    const filmPyrolysisRadiativeCoupledMixedFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(psf, iF),
    temperatureCoupledBase(patch(), psf.KMethod(), psf.kappaName()),
    TnbrName_(psf.TnbrName_),
    QrNbrName_(psf.QrNbrName_),
    QrName_(psf.QrName_),
    convectiveScaling_(psf.convectiveScaling_),
    filmDeltaDry_(psf.filmDeltaDry_),
    filmDeltaWet_(psf.filmDeltaWet_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void filmPyrolysisRadiativeCoupledMixedFvPatchScalarField::
updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());

    const label patchI = patch().index();
    const label nbrPatchI = mpp.samplePolyPatch().index();

    const polyMesh& mesh = patch().boundaryMesh().mesh();
    const polyMesh& nbrMesh = mpp.sampleMesh();
    const fvPatch& nbrPatch =
        refCast<const fvMesh>(nbrMesh).boundary()[nbrPatchI];

    scalarField intFld(patchInternalField());

    const filmPyrolysisRadiativeCoupledMixedFvPatchScalarField&
        nbrField =
        refCast
        <
            const filmPyrolysisRadiativeCoupledMixedFvPatchScalarField
        >
        (
            nbrPatch.lookupPatchField<volScalarField, scalar>(TnbrName_)
        );

    // Swap to obtain full local values of neighbour internal field
    scalarField nbrIntFld(nbrField.patchInternalField());
    mpp.distribute(nbrIntFld);

    scalarField& Tp = *this;

    const scalarField K(this->kappa(*this));
    const scalarField nbrK(nbrField.kappa(*this));

    // Swap to obtain full local values of neighbour K*delta
    scalarField KDeltaNbr(nbrK*nbrPatch.deltaCoeffs());
    mpp.distribute(KDeltaNbr);

    scalarField myKDelta(K*patch().deltaCoeffs());

    scalarList Tfilm(patch().size(), 0.0);
    scalarList htcwfilm(patch().size(), 0.0);
    scalarList filmDelta(patch().size(), 0.0);

    const pyrolysisModelType& pyrolysis = pyrModel();
    const filmModelType& film = filmModel();

    label myPatchINrbPatchI = -1;

    // Obtain Rad heat (Qr)
    scalarField Qr(patch().size(), 0.0);
    if (QrName_ != "none") //region0
    {
        Qr = patch().lookupPatchField<volScalarField, scalar>(QrName_);
        myPatchINrbPatchI = nbrPatch.index();
    }

    if (QrNbrName_ != "none") //pyrolysis
    {
        Qr = nbrPatch.lookupPatchField<volScalarField, scalar>(QrNbrName_);
        mpp.distribute(Qr);
        myPatchINrbPatchI = patchI;
    }

    const label filmPatchI =
        pyrolysis.nbrCoupledPatchID(film, myPatchINrbPatchI);

    const scalarField htcw(film.htcw().h()().boundaryField()[filmPatchI]);

    // Obtain htcw
    htcwfilm =
        const_cast<pyrolysisModelType&>(pyrolysis).mapRegionPatchField
        (
            film,
            myPatchINrbPatchI,
            filmPatchI,
            htcw,
            true
        );


    // Obtain Tfilm at the boundary through Ts.
    // NOTE: Tf is not good as at the boundary it will retrieve Tp
    Tfilm = film.Ts().boundaryField()[filmPatchI];
    film.toPrimary(filmPatchI, Tfilm);

    // Obtain delta
    filmDelta =
        const_cast<pyrolysisModelType&>(pyrolysis).mapRegionPatchField<scalar>
        (
            film,
            "deltaf",
            myPatchINrbPatchI,
            true
        );

     // Estimate wetness of the film (1: wet , 0: dry)
     scalarField ratio
     (
        min
        (
            max
            (
                (filmDelta - filmDeltaDry_)/(filmDeltaWet_ - filmDeltaDry_),
                scalar(0.0)
            ),
            scalar(1.0)
        )
     );

    scalarField qConv(ratio*htcwfilm*(Tfilm - Tp)*convectiveScaling_);

    scalarField qRad((1.0 - ratio)*Qr);

    scalarField alpha(KDeltaNbr - (qRad + qConv)/Tp);

    valueFraction() = alpha/(alpha + (1.0 - ratio)*myKDelta);

    refValue() = ratio*Tfilm + (1.0 - ratio)*(KDeltaNbr*nbrIntFld)/alpha;

    mixedFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        scalar Qc = gSum(qConv*patch().magSf());
        scalar Qr = gSum(qRad*patch().magSf());
        scalar Qt = gSum((qConv + qRad)*patch().magSf());

        Info<< mesh.name() << ':'
            << patch().name() << ':'
            << this->dimensionedInternalField().name() << " <- "
            << nbrMesh.name() << ':'
            << nbrPatch.name() << ':'
            << this->dimensionedInternalField().name() << " :" << nl
            << "     convective heat[W] : " << Qc << nl
            << "     radiative heat [W] : " << Qr << nl
            << "     total heat     [W] : " << Qt << nl
            << "     walltemperature "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }

}


void filmPyrolysisRadiativeCoupledMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("Tnbr")<< TnbrName_ << token::END_STATEMENT << nl;
    os.writeKeyword("QrNbr")<< QrNbrName_ << token::END_STATEMENT << nl;
    os.writeKeyword("Qr")<< QrName_ << token::END_STATEMENT << nl;
    os.writeKeyword("convectiveScaling") << convectiveScaling_
    << token::END_STATEMENT << nl;
    os.writeKeyword("filmDeltaDry") << filmDeltaDry_ <<
    token::END_STATEMENT << nl;
    os.writeKeyword("filmDeltaWet") << filmDeltaWet_ <<
    token::END_STATEMENT << endl;
    temperatureCoupledBase::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    filmPyrolysisRadiativeCoupledMixedFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //
