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

#include "toroidalCS.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(toroidalCS, 0);
    addToRunTimeSelectionTable(coordinateSystem, toroidalCS, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::toroidalCS::toroidalCS
(
    const word& name,
    const point& origin,
    const coordinateRotation& cr,
    const scalar radius
)
:
    coordinateSystem(name, origin, cr),
    radius_(radius)
{}


Foam::toroidalCS::toroidalCS
(
    const word& name,
    const dictionary& dict
)
:
    coordinateSystem(name, dict),
    radius_(readScalar(dict.lookup("radius")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::toroidalCS::localToGlobal
(
    const vector& local,
    bool translate
) const
{
    // Notation: r = local.x()
    scalar theta = degToRad(local.y());
    scalar phi = degToRad(local.z());

    scalar rprime = radius_ + local.x()*sin(phi);

    if ((local.x()*sin(phi)) > (radius_))
    {
        FatalErrorIn("toroidalCS::toGlobal(vector) const")
            << "Badly defined toroidal coordinates"
            << abort(FatalError);
    }

    return coordinateSystem::localToGlobal
    (
        vector(rprime*cos(theta), rprime*sin(theta), local.x()*cos(phi)),
        translate
    );
}


Foam::tmp<Foam::vectorField> Foam::toroidalCS::localToGlobal
(
    const vectorField& local,
    bool translate
) const
{
    const scalarField r
    (
        local.component(vector::X)
    );

    const scalarField theta
    (
        local.component(vector::Y)*constant::mathematical::pi/180.0
    );

    const scalarField phi
    (
        local.component(vector::Z)*constant::mathematical::pi/180.0
    );

    const scalarField rprime
    (
        radius_ + r*sin(phi)
    );

    vectorField lc(local.size());
    lc.replace(vector::X, rprime*cos(theta));
    lc.replace(vector::Y, rprime*sin(theta));
    lc.replace(vector::Z, r*cos(phi));

    return coordinateSystem::localToGlobal(lc, translate);
}


Foam::vector Foam::toroidalCS::globalToLocal
(
    const vector& global,
    bool translate
) const
{
    notImplemented
    (
        "toroidalCS::globalToLocal(const vector&, bool) const"
    );

    return vector::zero;
}


Foam::tmp<Foam::vectorField> Foam::toroidalCS::globalToLocal
(
    const vectorField& global,
    bool translate
) const
{
    notImplemented
    (
        "toroidalCS::globalToLocal(const vectorField&, bool) const"
    );

    return tmp<vectorField>(vectorField::null());
}


void Foam::toroidalCS::write(Ostream& os) const
{
    coordinateSystem::write(os);
    os << "radius: " << radius() << endl;
}


void Foam::toroidalCS::writeDict(Ostream& os, bool subDict) const
{
    if (subDict)
    {
        os  << indent << name() << nl
            << indent << token::BEGIN_BLOCK << incrIndent << nl;
    }

    coordinateSystem::writeDict(os, false);
    os.writeKeyword("radius") << radius() << token::END_STATEMENT << nl;

    if (subDict)
    {
        os << decrIndent << indent << token::END_BLOCK << endl;
    }
}


// ************************************************************************* //
