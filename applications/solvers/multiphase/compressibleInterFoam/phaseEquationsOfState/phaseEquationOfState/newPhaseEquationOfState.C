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

#include "phaseEquationOfState.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::phaseEquationOfState> Foam::phaseEquationOfState::New
(
    const dictionary& dict
)
{
    word phaseEquationOfStateType
    (
        dict.subDict("equationOfState").lookup("type")
    );

    Info<< "Selecting phaseEquationOfState "
        << phaseEquationOfStateType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(phaseEquationOfStateType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("phaseEquationOfState::New")
           << "Unknown phaseEquationOfStateType type "
           << phaseEquationOfStateType << endl << endl
           << "Valid phaseEquationOfState types are : " << endl
           << dictionaryConstructorTablePtr_->sortedToc()
           << exit(FatalError);
    }

    return cstrIter()(dict.subDict("equationOfState"));
}


// ************************************************************************* //
