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

Description
    Implementation of parser: test the state of either an istream or an
    ostream. Report an error if there is one.

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "token.H"
#include "int.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Print state of ostream.

void Foam::state(ostream& to, const string& s)
{
    state_value osState = state_value(to.rdstate());

    switch (osState)
    {
        case _good:                // Do not anything 'unusual'.
            break;

        case _eof:
            Info
                << "Output stream: premature end of stream", s << endl;
        break;

        case _fail:
            SeriousErrorIn("state(ostream& to, const string& s)")
                << "Output stream failure (bad format?)", s << endl;
        break;

        case (_fail + _eof) :
         SeriousErrorIn("state(ostream& to, const string& s)")
             << "Output stream failure and end of stream", s << endl;
        break;

        case _bad:
            SeriousErrorIn("state(ostream& to, const string& s)")
                << "Serious output stream failure", s << endl;
        break;

        default:
            SeriousErrorIn("state(ostream& to, const string& s)")
                << "Output stream failure of unknown type", s << endl
                << "Stream state value = ", osState << endl;
        break;
    }

    return;
}


//  Print state of istream.
void Foam::state(istream& from, const string& s)
{
    state_value isState = state_value(from.rdstate());

    switch (isState)
    {
        case _good:                // Do not anything 'unusual'.
            break;

        case _eof:
            Info
                << "Input stream: premature end of stream", s << endl;
            Info<< "If all else well, possibly a quote mark missing" << endl;
        break;

        case _fail:
            SeriousErrorIn("state(istream& from, const string& s)")
                << "Input stream failure (bad format?)", s << endl;
            Info<< "If all else well, possibly a quote mark missing" << endl;
        break;

        case (_fail + _eof) :
            SeriousErrorIn("state(istream& from, const string& s)")
                << "Input stream failure and end of stream", s << endl;
            Info<< "If all else well, possibly a quote mark missing" << endl;
        break;

        case _bad:
            SeriousErrorIn("state(istream& from, const string& s)")
                << "Serious input stream failure", s << endl;
        break;

        default:
            SeriousErrorIn("state(istream& from, const string& s)")
                << "Input stream failure of unknown type", s << endl;
            SeriousErrorIn("state(istream& from, const string& s)")
                << "Stream state value = ", isState << endl;
        break;
    }

    return;
}


// ************************************************************************* //
