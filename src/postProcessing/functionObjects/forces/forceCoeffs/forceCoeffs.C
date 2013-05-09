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

#include "forceCoeffs.H"
#include "dictionary.H"
#include "Time.H"
#include "Pstream.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(forceCoeffs, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::forceCoeffs::forceCoeffs
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    forces(name, obr, dict, loadFromFiles),
    liftDir_(vector::zero),
    dragDir_(vector::zero),
    pitchAxis_(vector::zero),
    magUInf_(0.0),
    lRef_(0.0),
    Aref_(0.0)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::forceCoeffs::~forceCoeffs()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::forceCoeffs::read(const dictionary& dict)
{
    if (active_)
    {
        forces::read(dict);

        // Directions for lift and drag forces, and pitch moment
        dict.lookup("liftDir") >> liftDir_;
        dict.lookup("dragDir") >> dragDir_;
        dict.lookup("pitchAxis") >> pitchAxis_;

        // Free stream velocity magnitude
        dict.lookup("magUInf") >> magUInf_;

        // Reference length and area scales
        dict.lookup("lRef") >> lRef_;
        dict.lookup("Aref") >> Aref_;
    }
}


void Foam::forceCoeffs::writeFileHeader(const label i)
{
    file()
        << "# Time" << tab << "Cm" << tab << "Cd" << tab << "Cl" << tab
        << "Cl(f)" << "Cl(r)" << endl;
}


void Foam::forceCoeffs::execute()
{
    // Do nothing - only valid on write
}


void Foam::forceCoeffs::end()
{
    // Do nothing - only valid on write
}


void Foam::forceCoeffs::write()
{
    if (active_)
    {
        forces::calcForcesMoment();

        if (Pstream::master())
        {
            functionObjectFile::write();

            scalar pDyn = 0.5*rhoRef_*magUInf_*magUInf_;

            Field<vector> totForce(force_[0] + force_[1] + force_[2]);
            Field<vector> totMoment(moment_[0] + moment_[1] + moment_[2]);

            List<Field<scalar> > coeffs(3);
            coeffs[0].setSize(nBin_);
            coeffs[1].setSize(nBin_);
            coeffs[2].setSize(nBin_);

            // lift, drag and moment
            coeffs[0] = (totForce & liftDir_)/(Aref_*pDyn);
            coeffs[1] = (totForce & dragDir_)/(Aref_*pDyn);
            coeffs[2] = (totMoment & pitchAxis_)/(Aref_*lRef_*pDyn);

            scalar Cl = sum(coeffs[0]);
            scalar Cd = sum(coeffs[1]);
            scalar Cm = sum(coeffs[2]);

            scalar Clf = Cl/2.0 + Cm;
            scalar Clr = Cl/2.0 - Cm;

            file()
                << obr_.time().value() << tab
                << Cm << tab << Cd << tab << Cl << tab << Clf << tab << Clr
                << endl;

            if (log_)
            {
                Info<< type() << " output:" << nl
                    << "    Cm    = " << Cm << nl
                    << "    Cd    = " << Cd << nl
                    << "    Cl    = " << Cl << nl
                    << "    Cl(f) = " << Clf << nl
                    << "    Cl(r) = " << Clr << endl;
            }

            if (nBin_ > 1)
            {
                autoPtr<writer<scalar> >
                    binWriterPtr(writer<scalar>::New(binFormat_));
                wordList fieldNames(IStringStream("(lift drag moment)")());

                coordSet axis
                (
                    "forceCoeffs",
                    "distance",
                    binPoints_,
                    mag(binPoints_)
                );

                fileName forcesDir = baseTimeDir();
                mkDir(forcesDir);

                if (log_)
                {
                    Info<< "    Writing bins to " << forcesDir << endl;
                }

                OFstream osCoeffs(forcesDir/"forceCoeffs_bins");

                if (binCumulative_)
                {
                    for (label i = 1; i < coeffs[0].size(); i++)
                    {
                        coeffs[0][i] += coeffs[0][i-1];
                        coeffs[1][i] += coeffs[1][i-1];
                        coeffs[2][i] += coeffs[2][i-1];
                    }
                }

                binWriterPtr->write(axis, fieldNames, coeffs, osCoeffs);
            }

            if (log_)
            {
                Info<< endl;
            }
        }
    }
}


// ************************************************************************* //
