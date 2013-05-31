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
#include "DarcyForchheimer.H"
#include "geometricOneField.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace porosityModels
    {
        defineTypeNameAndDebug(DarcyForchheimer, 0);
        addToRunTimeSelectionTable(porosityModel, DarcyForchheimer, mesh);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porosityModels::DarcyForchheimer::DarcyForchheimer
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
    D_("D", dimless/sqr(dimLength), tensor::zero),
    F_("F", dimless/dimLength, tensor::zero),
    rhoName_(coeffs_.lookupOrDefault<word>("rho", "rho")),
    muName_(coeffs_.lookupOrDefault<word>("mu", "thermo:mu")),
    nuName_(coeffs_.lookupOrDefault<word>("nu", "nu"))
{
    // local-to-global transformation tensor
    const tensor& E = coordSys_.R();

    dimensionedVector d(coeffs_.lookup("d"));
    if (D_.dimensions() != d.dimensions())
    {
        FatalIOErrorIn
        (
            "Foam::porosityModels::DarcyForchheimer::DarcyForchheimer"
            "("
                "const word&, "
                "const word&, "
                "const fvMesh&, "
                "const dictionary&"
            ")",
            coeffs_
        )   << "incorrect dimensions for d: " << d.dimensions()
            << " should be " << D_.dimensions()
            << exit(FatalIOError);
    }

    adjustNegativeResistance(d);

    D_.value().xx() = d.value().x();
    D_.value().yy() = d.value().y();
    D_.value().zz() = d.value().z();
    D_.value() = (E & D_ & E.T()).value();

    dimensionedVector f(coeffs_.lookup("f"));
    if (F_.dimensions() != f.dimensions())
    {
        FatalIOErrorIn
        (
            "Foam::porosityModels::DarcyForchheimer::DarcyForchheimer"
            "("
                "const word&, "
                "const word&, "
                "const fvMesh&, "
                "const dictionary&"
            ")",
            coeffs_
        )   << "incorrect dimensions for f: " << f.dimensions()
            << " should be " << F_.dimensions()
            << exit(FatalIOError);
    }

    adjustNegativeResistance(f);

    // leading 0.5 is from 1/2*rho
    F_.value().xx() = 0.5*f.value().x();
    F_.value().yy() = 0.5*f.value().y();
    F_.value().zz() = 0.5*f.value().z();
    F_.value() = (E & F_ & E.T()).value();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::porosityModels::DarcyForchheimer::~DarcyForchheimer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porosityModels::DarcyForchheimer::calcForce
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

    apply(Udiag, Usource, V, rho, mu, U);

    force = Udiag*U - Usource;
}


void Foam::porosityModels::DarcyForchheimer::correct
(
    fvVectorMatrix& UEqn
) const
{
    const vectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();
 
    if (UEqn.dimensions() == dimForce)
    {
        const volScalarField& rho =
            mesh_.lookupObject<volScalarField>(rhoName_);
        const volScalarField& mu =
            mesh_.lookupObject<volScalarField>(muName_);

        apply(Udiag, Usource, V, rho, mu, U);
    }
    else
    {
        const volScalarField& nu =
            mesh_.lookupObject<volScalarField>(nuName_);

        apply(Udiag, Usource, V, geometricOneField(), nu, U);
    }
}


void Foam::porosityModels::DarcyForchheimer::correct
(
    fvVectorMatrix& UEqn,
    const volScalarField& rho,
    const volScalarField& mu
) const
{
    const vectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();
 
    apply(Udiag, Usource, V, rho, mu, U);
}


void Foam::porosityModels::DarcyForchheimer::correct
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU
) const
{
    const vectorField& U = UEqn.psi();

    if (UEqn.dimensions() == dimForce)
    {
        const volScalarField& rho =
            mesh_.lookupObject<volScalarField>(rhoName_);
        const volScalarField& mu =
            mesh_.lookupObject<volScalarField>(muName_);

        apply(AU, rho, mu, U);
    }
    else
    {
        const volScalarField& nu =
            mesh_.lookupObject<volScalarField>(nuName_);

        apply(AU, geometricOneField(), nu, U);
    }
}


bool Foam::porosityModels::DarcyForchheimer::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);

    return true;
}


// ************************************************************************* //
