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

#include "SquareMatrix.H"
#include "vector.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    SquareMatrix<scalar> hmm(3);

    hmm[0][0] = -3.0;
    hmm[0][1] = 10.0;
    hmm[0][2] = -4.0;
    hmm[1][0] = 2.0;
    hmm[1][1] = 3.0;
    hmm[1][2] = 10.0;
    hmm[2][0] = 2.0;
    hmm[2][1] = 6.0;
    hmm[2][2] = 1.0;

    //Info<< hmm << endl << hmm - 2.0*(-hmm) << endl;
    Info<< max(hmm) << endl;
    Info<< min(hmm) << endl;

    SquareMatrix<scalar> hmm2(3, 1.0);

    hmm = hmm2;

    Info<< hmm << endl;

    SquareMatrix<scalar> hmm3(Sin);

    Info<< hmm3 << endl;

    SquareMatrix<scalar> hmm4;

    hmm4 = hmm2;

    Info<< hmm4 << endl;

    SquareMatrix<scalar> hmm5;

    hmm4 = hmm5;
    Info<< hmm5 << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
