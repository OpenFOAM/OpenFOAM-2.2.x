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

#include "MRFZoneList.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MRFZoneList::MRFZoneList
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    PtrList<MRFZone>(),
    mesh_(mesh)
{
    reset(dict);

    active(true);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::MRFZoneList::~MRFZoneList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::MRFZoneList::active(const bool warn) const
{
    bool a = false;
    forAll(*this, i)
    {
        a = a || this->operator[](i).active();
    }

    if (warn && this->size() && !a)
    {
        Info<< "    No MRF zones active" << endl;
    }

    return a;
}


void Foam::MRFZoneList::reset(const dictionary& dict)
{
    label count = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        if (iter().isDict())
        {
            count++;
        }
    }

    this->setSize(count);
    label i = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        if (iter().isDict())
        {
            const word& name = iter().keyword();
            const dictionary& modelDict = iter().dict();

            Info<< "    creating MRF zone: " << name << endl;

            this->set
            (
                i++,
                new MRFZone(name, mesh_, modelDict)
            );
        }
    }
}


bool Foam::MRFZoneList::read(const dictionary& dict)
{
    bool allOk = true;
    forAll(*this, i)
    {
        MRFZone& pm = this->operator[](i);
        bool ok = pm.read(dict.subDict(pm.name()));
        allOk = (allOk && ok);
    }
    return allOk;
}


bool Foam::MRFZoneList::writeData(Ostream& os) const
{
    forAll(*this, i)
    {
        os  << nl;
        this->operator[](i).writeData(os);
    }

    return os.good();
}


void Foam::MRFZoneList::addCoriolis
(
    const volVectorField& U,
    volVectorField& ddtU
) const
{
    forAll(*this, i)
    {
        operator[](i).addCoriolis(U, ddtU);
    }
}


void Foam::MRFZoneList::addCoriolis(fvVectorMatrix& UEqn) const
{
    forAll(*this, i)
    {
        operator[](i).addCoriolis(UEqn);
    }
}


void Foam::MRFZoneList::addCoriolis
(
    const volScalarField& rho,
    fvVectorMatrix& UEqn
) const
{
    forAll(*this, i)
    {
        operator[](i).addCoriolis(rho, UEqn);
    }
}


void Foam::MRFZoneList::relativeVelocity(volVectorField& U) const
{
    forAll(*this, i)
    {
        operator[](i).relativeVelocity(U);
    }
}


void Foam::MRFZoneList::absoluteVelocity(volVectorField& U) const
{
    forAll(*this, i)
    {
        operator[](i).absoluteVelocity(U);
    }
}


void Foam::MRFZoneList::relativeFlux(surfaceScalarField& phi) const
{
    forAll(*this, i)
    {
        operator[](i).relativeFlux(phi);
    }
}


void Foam::MRFZoneList::relativeFlux
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    forAll(*this, i)
    {
        operator[](i).relativeFlux(rho, phi);
    }
}


void Foam::MRFZoneList::absoluteFlux(surfaceScalarField& phi) const
{
    forAll(*this, i)
    {
        operator[](i).absoluteFlux(phi);
    }
}


void Foam::MRFZoneList::absoluteFlux
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    forAll(*this, i)
    {
        operator[](i).absoluteFlux(rho, phi);
    }
}


void Foam::MRFZoneList::correctBoundaryVelocity(volVectorField& U) const
{
    forAll(*this, i)
    {
        operator[](i).correctBoundaryVelocity(U);
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const MRFZoneList& models
)
{
    models.writeData(os);
    return os;
}


// ************************************************************************* //
