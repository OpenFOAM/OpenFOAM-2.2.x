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

#include "reactingOneDim.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "surfaceInterpolate.H"
#include "fvm.H"
#include "fvcDiv.H"
#include "fvcVolumeIntegrate.H"
#include "fvMatrices.H"
#include "absorptionEmissionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace pyrolysisModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(reactingOneDim, 0);

addToRunTimeSelectionTable(pyrolysisModel, reactingOneDim, mesh);
addToRunTimeSelectionTable(pyrolysisModel, reactingOneDim, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void reactingOneDim::readReactingOneDimControls()
{
    const dictionary& solution = this->solution().subDict("SIMPLE");
    solution.lookup("nNonOrthCorr") >> nNonOrthCorr_;
    time().controlDict().lookup("maxDi") >> maxDiff_;

    coeffs().lookup("radFluxName") >> primaryRadFluxName_;
    coeffs().lookup("minimumDelta") >> minimumDelta_;
}


bool reactingOneDim::read()
{
    if (pyrolysisModel::read())
    {
        readReactingOneDimControls();
        return true;
    }
    else
    {
        return false;
    }
}


bool reactingOneDim::read(const dictionary& dict)
{
    if (pyrolysisModel::read(dict))
    {
        readReactingOneDimControls();
        return true;
    }
    else
    {
        return false;
    }
}


void reactingOneDim::updatePhiGas()
{
    phiHsGas_ ==  dimensionedScalar("zero", phiHsGas_.dimensions(), 0.0);
    phiGas_ == dimensionedScalar("zero", phiGas_.dimensions(), 0.0);

    const speciesTable& gasTable = solidChemistry_->gasTable();

    forAll(gasTable, gasI)
    {
        tmp<volScalarField> tHsiGas =
            solidChemistry_->gasHs(solidThermo_.p(), solidThermo_.T(), gasI);

        const volScalarField& HsiGas = tHsiGas();

        const DimensionedField<scalar, volMesh>& RRiGas =
            solidChemistry_->RRg(gasI);

        label totalFaceId = 0;
        forAll(intCoupledPatchIDs_, i)
        {
            const label patchI = intCoupledPatchIDs_[i];

            const scalarField& phiGasp = phiHsGas_.boundaryField()[patchI];

            forAll(phiGasp, faceI)
            {
                const labelList& cells = boundaryFaceCells_[totalFaceId];
                scalar massInt = 0.0;
                forAllReverse(cells, k)
                {
                    const label cellI = cells[k];
                    massInt += RRiGas[cellI]*regionMesh().V()[cellI];
                    phiHsGas_[cellI] += massInt*HsiGas[cellI];
                }

                phiGas_.boundaryField()[patchI][faceI] += massInt;

                if (debug)
                {
                    Info<< " Gas : " << gasTable[gasI]
                        << " on patch : " << patchI
                        << " mass produced at face(local) : "
                        <<  faceI
                        << " is : " << massInt
                        << " [kg/s] " << endl;
                }
                totalFaceId ++;
            }
        }
        tHsiGas().clear();
    }
}


void reactingOneDim::updateFields()
{
    updatePhiGas();
}


void reactingOneDim::updateMesh(const scalarField& mass0)
{
    if (!moveMesh_)
    {
        return;
    }

    const scalarField newV(mass0/rho_);

    Info<< "Initial/final volumes = " << gSum(regionMesh().V()) << ", "
        << gSum(newV) << " [m3]" << endl;

    // move the mesh
    const labelList moveMap = moveMesh(regionMesh().V() - newV, minimumDelta_);

    // flag any cells that have not moved as non-reacting
    forAll(moveMap, i)
    {
        if (moveMap[i] == 0)
        {
            solidChemistry_->setCellReacting(i, false);
        }
    }
}


void reactingOneDim::solveContinuity()
{
    if (debug)
    {
        Info<< "reactingOneDim::solveContinuity()" << endl;
    }

    if (moveMesh_)
    {
        const scalarField mass0 = rho_*regionMesh().V();

        fvScalarMatrix rhoEqn
        (
            fvm::ddt(rho_)
         ==
          - solidChemistry_->RRg()
        );

        rhoEqn.solve();

        updateMesh(mass0);

    }
}


void reactingOneDim::solveSpeciesMass()
{
    if (debug)
    {
        Info<< "reactingOneDim::solveSpeciesMass()" << endl;
    }

    volScalarField Yt(0.0*Ys_[0]);

    for (label i=0; i<Ys_.size()-1; i++)
    {
        volScalarField& Yi = Ys_[i];

        fvScalarMatrix YiEqn
        (
            fvm::ddt(rho_, Yi)
         ==
            solidChemistry_->RRs(i)
        );

        if (regionMesh().moving())
        {
            surfaceScalarField phiYiRhoMesh
            (
                fvc::interpolate(Yi*rho_)*regionMesh().phi()
            );

            YiEqn += fvc::div(phiYiRhoMesh);

        }

        YiEqn.solve(regionMesh().solver("Yi"));
        Yi.max(0.0);
        Yt += Yi;
    }

    Ys_[Ys_.size() - 1] = 1.0 - Yt;
}


void reactingOneDim::solveEnergy()
{
    if (debug)
    {
        Info<< "reactingOneDim::solveEnergy()" << endl;
    }

    tmp<volScalarField> alpha(solidThermo_.alpha());

    fvScalarMatrix hEqn
    (
        fvm::ddt(rho_, h_)
      - fvm::laplacian(alpha, h_)
     ==
        chemistrySh_
    );

    if (regionMesh().moving())
    {
        surfaceScalarField phihMesh
        (
            fvc::interpolate(rho_*h_)*regionMesh().phi()
        );

        hEqn += fvc::div(phihMesh);
    }

    hEqn.relax();
    hEqn.solve();

    Info<< "pyrolysis min/max(T) = " << min(solidThermo_.T()) << ", "
        << max(solidThermo_.T()) << endl;
}


void reactingOneDim::calculateMassTransfer()
{
    totalGasMassFlux_ = 0;
    forAll(intCoupledPatchIDs_, i)
    {
        const label patchI = intCoupledPatchIDs_[i];
        totalGasMassFlux_ += gSum(phiGas_.boundaryField()[patchI]);
    }

    if (infoOutput_)
    {
        totalHeatRR_ = fvc::domainIntegrate(chemistrySh_);

        addedGasMass_ +=
            fvc::domainIntegrate(solidChemistry_->RRg())*time_.deltaT();
        lostSolidMass_ +=
            fvc::domainIntegrate(solidChemistry_->RRs())*time_.deltaT();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

reactingOneDim::reactingOneDim(const word& modelType, const fvMesh& mesh)
:
    pyrolysisModel(modelType, mesh),
    solidChemistry_(solidChemistryModel::New(regionMesh())),
    solidThermo_(solidChemistry_->solidThermo()),
    radiation_(radiation::radiationModel::New(solidThermo_.T())),
    rho_
    (
        IOobject
        (
            "rho",
            regionMesh().time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        solidThermo_.rho()
    ),
    Ys_(solidThermo_.composition().Y()),
    h_(solidThermo_.he()),
    nNonOrthCorr_(-1),
    maxDiff_(10),
    minimumDelta_(1e-4),

    phiGas_
    (
        IOobject
        (
            "phiGas",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/dimTime, 0.0)
    ),

    phiHsGas_
    (
        IOobject
        (
            "phiHsGas",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime, 0.0)
    ),

    chemistrySh_
    (
        IOobject
        (
            "chemistrySh",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
    ),

    lostSolidMass_(dimensionedScalar("zero", dimMass, 0.0)),
    addedGasMass_(dimensionedScalar("zero", dimMass, 0.0)),
    totalGasMassFlux_(0.0),
    totalHeatRR_(dimensionedScalar("zero", dimEnergy/dimTime, 0.0))
{
    if (active_)
    {
        read();
    }
}


reactingOneDim::reactingOneDim
(
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    pyrolysisModel(modelType, mesh, dict),
    solidChemistry_(solidChemistryModel::New(regionMesh())),
    solidThermo_(solidChemistry_->solidThermo()),
    radiation_(radiation::radiationModel::New(solidThermo_.T())),
    rho_
    (
        IOobject
        (
            "rho",
            regionMesh().time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        solidThermo_.rho()
    ),
    Ys_(solidThermo_.composition().Y()),
    h_(solidThermo_.he()),
    nNonOrthCorr_(-1),
    maxDiff_(10),
    minimumDelta_(1e-4),

    phiGas_
    (
        IOobject
        (
            "phiGas",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/dimTime, 0.0)
    ),

    phiHsGas_
    (
        IOobject
        (
            "phiHsGas",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime, 0.0)
    ),

    chemistrySh_
    (
        IOobject
        (
            "chemistrySh",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
    ),

    lostSolidMass_(dimensionedScalar("zero", dimMass, 0.0)),
    addedGasMass_(dimensionedScalar("zero", dimMass, 0.0)),
    totalGasMassFlux_(0.0),
    totalHeatRR_(dimensionedScalar("zero", dimEnergy/dimTime, 0.0))
{
    if (active_)
    {
        read(dict);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

reactingOneDim::~reactingOneDim()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar reactingOneDim::addMassSources(const label patchI, const label faceI)
{
    label index = 0;
    forAll(primaryPatchIDs_, i)
    {
        if (primaryPatchIDs_[i] == patchI)
        {
            index = i;
            break;
        }
    }

    const label localPatchId =  intCoupledPatchIDs_[index];

    const scalar massAdded = phiGas_.boundaryField()[localPatchId][faceI];

    if (debug)
    {
        Info<< "\nPyrolysis region: " << type() << "added mass : "
            << massAdded << endl;
    }

    return massAdded;
}


scalar reactingOneDim::solidRegionDiffNo() const
{
    scalar DiNum = -GREAT;

    if (regionMesh().nInternalFaces() > 0)
    {
        surfaceScalarField KrhoCpbyDelta
        (
            regionMesh().surfaceInterpolation::deltaCoeffs()
          * fvc::interpolate(kappa())
          / fvc::interpolate(Cp()*rho_)
        );

        DiNum = max(KrhoCpbyDelta.internalField())*time().deltaTValue();
    }

    return DiNum;
}


scalar reactingOneDim::maxDiff() const
{
    return maxDiff_;
}


const volScalarField& reactingOneDim::rho() const
{
    return rho_;
}


const volScalarField& reactingOneDim::T() const
{
    return solidThermo_.T();
}


const tmp<volScalarField> reactingOneDim::Cp() const
{
    return solidThermo_.Cp();
}


tmp<volScalarField> reactingOneDim::kappaRad() const
{
    return radiation_->absorptionEmission().a();
}


tmp<volScalarField> reactingOneDim::kappa() const
{
    return solidThermo_.kappa();
}


const surfaceScalarField& reactingOneDim::phiGas() const
{
    return phiGas_;
}


void reactingOneDim::preEvolveRegion()
{
    pyrolysisModel::preEvolveRegion();

    // Initialise all cells as able to react
    forAll(h_, cellI)
    {
        solidChemistry_->setCellReacting(cellI, true);
    }
}


void reactingOneDim::evolveRegion()
{
    Info<< "\nEvolving pyrolysis in region: " << regionMesh().name() << endl;

    solidChemistry_->solve
    (
        time().value() - time().deltaTValue(),
        time().deltaTValue()
    );

    calculateMassTransfer();

    solveContinuity();

    chemistrySh_ = solidChemistry_->Sh()();

    updateFields();

    solveSpeciesMass();

    for (int nonOrth=0; nonOrth<=nNonOrthCorr_; nonOrth++)
    {
        solveEnergy();
    }

    solidThermo_.correct();

    rho_ = solidThermo_.rho();
}


void reactingOneDim::info() const
{
    Info<< "\nPyrolysis in region: " << regionMesh().name() << endl;

    Info<< indent << "Total gas mass produced  [kg] = "
        << addedGasMass_.value() << nl
        << indent << "Total solid mass lost    [kg] = "
        << lostSolidMass_.value() << nl
        << indent << "Total pyrolysis gases  [kg/s] = "
        << totalGasMassFlux_ << nl
        << indent << "Total heat release rate [J/s] = "
        << totalHeatRR_.value() << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace regionModels
} // End namespace pyrolysisModels

// ************************************************************************* //
