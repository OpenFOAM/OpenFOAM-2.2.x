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

Class
    Foam::phaseModel

SourceFiles
    phaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef phaseModel_H
#define phaseModel_H

#include "dictionary.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "rhoThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class diameterModel;

/*---------------------------------------------------------------------------*\
                           Class phaseModel Declaration
\*---------------------------------------------------------------------------*/

class phaseModel
:
    public volScalarField
{
    // Private data

        //- Name of phase
        word name_;

        dictionary phaseDict_;

        //- Thermophysical properties
        autoPtr<rhoThermo> thermo_;

        //- Velocity
        volVectorField U_;

        //- Fluxes
        autoPtr<surfaceScalarField> phiPtr_;

        //- Diameter model
        autoPtr<diameterModel> dPtr_;


public:

    // Constructors

        phaseModel
        (
            const fvMesh& mesh,
            const dictionary& transportProperties,
            const word& phaseName
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<phaseModel> New
        (
            const fvMesh& mesh,
            const dictionary& transportProperties,
            const word& phaseName
        );


    //- Destructor
    virtual ~phaseModel();


    // Member Functions

        const word& name() const
        {
            return name_;
        }

        tmp<volScalarField> d() const;

        tmp<volScalarField> nu() const
        {
            return thermo_->mu()/thermo_->rho();
        }

        tmp<volScalarField> kappa() const
        {
            return thermo_->kappa();
        }

        tmp<volScalarField> Cp() const
        {
            return thermo_->Cp();
        }

        const volScalarField& rho() const
        {
            return thermo_->rho();
        }

        const rhoThermo& thermo() const
        {
            return thermo_();
        }

        rhoThermo& thermo()
        {
            return thermo_();
        }

        const volVectorField& U() const
        {
            return U_;
        }

        volVectorField& U()
        {
            return U_;
        }

        const surfaceScalarField& phi() const
        {
            return phiPtr_();
        }

        surfaceScalarField& phi()
        {
            return phiPtr_();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
