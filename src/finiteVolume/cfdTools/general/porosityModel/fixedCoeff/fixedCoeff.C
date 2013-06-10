/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "addToRunTimeSelectionTable.H"
#include "fixedCoeff.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace porosityModels
    {
        defineTypeNameAndDebug(fixedCoeff, 0);
        addToRunTimeSelectionTable(porosityModel, fixedCoeff, mesh);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::porosityModels::fixedCoeff::apply
(
    scalarField& Udiag,
    vectorField& Usource,
    const scalarField& V,
    const vectorField& U,
    const scalar rho
) const
{
    const tensor& alpha = alpha_.value();
    const tensor& beta = beta_.value();

    forAll(cellZoneIds_, zoneI)
    {
        const labelList& cells = mesh_.cellZones()[cellZoneIds_[zoneI]];

        forAll(cells, i)
        {
            const label cellI = cells[i];

            const tensor Cd = rho*(alpha + beta*mag(U[cellI]));

            const scalar isoCd = tr(Cd);

            Udiag[cellI] += V[cellI]*isoCd;
            Usource[cellI] -= V[cellI]*((Cd - I*isoCd) & U[cellI]);
        }
    }
}


void Foam::porosityModels::fixedCoeff::apply
(
    tensorField& AU,
    const vectorField& U,
    const scalar rho
) const
{
    const tensor& alpha = alpha_.value();
    const tensor& beta = beta_.value();

    forAll(cellZoneIds_, zoneI)
    {
        const labelList& cells = mesh_.cellZones()[cellZoneIds_[zoneI]];

        forAll(cells, i)
        {
            const label cellI = cells[i];
            AU[cellI] += rho*(alpha + beta*mag(U[cellI]));
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porosityModels::fixedCoeff::fixedCoeff
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& cellZoneName
)
:
    porosityModel(name, modelType, mesh, dict, cellZoneName),
    coordSys_(coeffs_, mesh),
    alpha_("alpha", dimless/dimTime, tensor::zero),
    beta_("beta", dimless/dimLength, tensor::zero)
{
    // local-to-global transformation tensor
    const tensor& E = coordSys_.R();

    dimensionedVector alpha(coeffs_.lookup("alpha"));
    if (alpha_.dimensions() != alpha.dimensions())
    {
        FatalIOErrorIn
        (
            "Foam::porosityModels::fixedCoeff::fixedCoeff"
            "("
                "const word&, "
                "const word&, "
                "const fvMesh&, "
                "const dictionary&"
            ")",
            coeffs_
        )   << "incorrect dimensions for alpha: " << alpha.dimensions()
            << " should be " << alpha_.dimensions()
            << exit(FatalIOError);
    }

    adjustNegativeResistance(alpha);

    alpha_.value().xx() = alpha.value().x();
    alpha_.value().yy() = alpha.value().y();
    alpha_.value().zz() = alpha.value().z();
    alpha_.value() = (E & alpha_ & E.T()).value();

    dimensionedVector beta(coeffs_.lookup("beta"));
    if (beta_.dimensions() != beta.dimensions())
    {
        FatalIOErrorIn
        (
            "Foam::porosityModels::fixedCoeff::fixedCoeff"
            "("
                "const word&, "
                "const word&, "
                "const fvMesh&, "
                "const dictionary&"
            ")",
            coeffs_
        )   << "incorrect dimensions for beta: " << beta.dimensions()
            << " should be " << beta_.dimensions()
            << exit(FatalIOError);
    }

    adjustNegativeResistance(beta);

    beta_.value().xx() = beta.value().x();
    beta_.value().yy() = beta.value().y();
    beta_.value().zz() = beta.value().z();
    beta_.value() = (E & beta_ & E.T()).value();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::porosityModels::fixedCoeff::~fixedCoeff()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porosityModels::fixedCoeff::calcForce
(
    const volVectorField& U,
    const volScalarField& rho,
    const volScalarField& mu,
    vectorField& force
) const
{
    scalarField Udiag(U.size(), 0.0);
    vectorField Usource(U.size(), vector::zero);
    const scalarField& V = mesh_.V();
    scalar rhoRef = readScalar(coeffs_.lookup("rhoRef"));

    apply(Udiag, Usource, V, U, rhoRef);

    force = Udiag*U - Usource;
}


void Foam::porosityModels::fixedCoeff::correct
(
    fvVectorMatrix& UEqn
) const
{
    const vectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();
 
    scalar rho = 1.0;
    if (UEqn.dimensions() == dimForce)
    {
        coeffs_.lookup("rhoRef") >> rho;
    }

    apply(Udiag, Usource, V, U, rho);
}


void Foam::porosityModels::fixedCoeff::correct
(
    fvVectorMatrix& UEqn,
    const volScalarField&,
    const volScalarField&
) const
{
    const vectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();

    scalar rho = 1.0;
    if (UEqn.dimensions() == dimForce)
    {
        coeffs_.lookup("rhoRef") >> rho;
    }
 
    apply(Udiag, Usource, V, U, rho);
}


void Foam::porosityModels::fixedCoeff::correct
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU
) const
{
    const vectorField& U = UEqn.psi();

    scalar rho = 1.0;
    if (UEqn.dimensions() == dimForce)
    {
        coeffs_.lookup("rhoRef") >> rho;
    }

    apply(AU, U, rho);
}


bool Foam::porosityModels::fixedCoeff::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);

    return true;
}


// ************************************************************************* //
