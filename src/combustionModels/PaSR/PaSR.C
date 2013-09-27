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

#include "PaSR.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::combustionModels::PaSR<Type>::PaSR
(
    const word& modelType,
    const fvMesh& mesh
)
:
    Type(modelType, mesh),
    Cmix_(readScalar(this->coeffs().lookup("Cmix"))),
    turbulentReaction_(this->coeffs().lookup("turbulentReaction")),
    kappa_
    (
        IOobject
        (
            "PaSR::kappa",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("kappa", dimless, 0.0)
    ),
    useReactionRate_(this->coeffs().lookupOrDefault("useReactionRate", false))
{
    if (useReactionRate_)
    {
        Info<< "    using reaction rate" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::combustionModels::PaSR<Type>::~PaSR()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::volScalarField> Foam::combustionModels::PaSR<Type>::tc() const
{
    return this->chemistryPtr_->tc();
}


template<class Type>
void Foam::combustionModels::PaSR<Type>::correct()
{
    if (this->active())
    {
        const scalar t = this->mesh().time().value();
        const scalar dt = this->mesh().time().deltaTValue();

        if (!useReactionRate_)
        {
            this->chemistryPtr_->solve(t - dt, dt);
        }
        else
        {
            this->chemistryPtr_->calculate();
        }

        if (turbulentReaction_)
        {
            tmp<volScalarField> trho(this->rho());
            const volScalarField& rho = trho();
            tmp<volScalarField> tepsilon(this->turbulence().epsilon());
            const volScalarField& epsilon = tepsilon();
            tmp<volScalarField> tmuEff(this->turbulence().muEff());
            const volScalarField& muEff = tmuEff();

            tmp<volScalarField> ttc(tc());
            const volScalarField& tc = ttc();

            forAll(epsilon, i)
            {
                scalar tk =
                    Cmix_*Foam::sqrt(muEff[i]/rho[i]/(epsilon[i] + SMALL));

                if (tk > SMALL)
                {
                    kappa_[i] = tc[i]/(tc[i] + tk);
                }
                else
                {
                    kappa_[i] = 1.0;
                }
            }
        }
        else
        {
            kappa_ = 1.0;
        }
    }
}


template<class Type>
Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::PaSR<Type>::R(volScalarField& Y) const
{
    tmp<fvScalarMatrix> tSu(new fvScalarMatrix(Y, dimMass/dimTime));

    fvScalarMatrix& Su = tSu();

    if (this->active())
    {
        const label specieI = this->thermo().composition().species()[Y.name()];

        Su += kappa_*this->chemistryPtr_->RR(specieI);
    }

    return tSu;
}


template<class Type>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::PaSR<Type>::dQ() const
{
    tmp<volScalarField> tdQ
    (
        new volScalarField
        (
            IOobject
            (
                "dQ",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("dQ", dimEnergy/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (this->active())
    {
        volScalarField& dQ = tdQ();
        dQ = kappa_*this->chemistryPtr_->dQ();
    }

    return tdQ;
}


template<class Type>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::PaSR<Type>::Sh() const
{
    tmp<volScalarField> tSh
    (
        new volScalarField
        (
            IOobject
            (
                "Sh",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (this->active())
    {
        scalarField& Sh = tSh();
        Sh = kappa_*this->chemistryPtr_->Sh();
    }

    return tSh;
}


template<class Type>
bool Foam::combustionModels::PaSR<Type>::read()
{
    if (Type::read())
    {
        this->coeffs().lookup("Cmix") >> Cmix_;
        this->coeffs().lookup("turbulentReaction") >> turbulentReaction_;
        this->coeffs().lookup("useReactionRate") >> useReactionRate_;
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
