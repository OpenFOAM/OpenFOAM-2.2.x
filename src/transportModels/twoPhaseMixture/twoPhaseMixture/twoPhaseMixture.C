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

#include "twoPhaseMixture.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseMixture::twoPhaseMixture
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& alpha1Name,
    const word& alpha2Name
)
:
    phase1Name_
    (
        dict.found("phases")
      ? wordList(dict.lookup("phases"))[0]
      : "1"
    ),
    phase2Name_
    (
        dict.found("phases")
      ? wordList(dict.lookup("phases"))[1]
      : "2"
    ),

    alpha1_
    (
        IOobject
        (
            dict.found("phases") ? word("alpha" + phase1Name_) : alpha1Name,
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    alpha2_
    (
        IOobject
        (
            dict.found("phases") ? word("alpha" + phase2Name_) : alpha2Name,
            mesh.time().timeName(),
            mesh
        ),
        1.0 - alpha1_
    )
{}


// ************************************************************************* //
