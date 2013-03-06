/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "singlePhaseTransportModel.H"
#include "viscosityModel.H"
#include "volFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::singlePhaseTransportModel::singlePhaseTransportModel
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    transportModel(U, phi),
    viscosityModelPtr_(viscosityModel::New("nu", *this, U, phi))
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::singlePhaseTransportModel::~singlePhaseTransportModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::singlePhaseTransportModel::nu() const
{
    return viscosityModelPtr_->nu();
}


void Foam::singlePhaseTransportModel::correct()
{
    viscosityModelPtr_->correct();
}


bool Foam::singlePhaseTransportModel::read()
{
    if (regIOobject::read())
    {
        return viscosityModelPtr_->read(*this);
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
