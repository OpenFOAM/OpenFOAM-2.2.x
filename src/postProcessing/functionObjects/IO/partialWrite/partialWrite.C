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

#include "partialWrite.H"
#include "dictionary.H"
#include "Time.H"
#include "IOobjectList.H"
#include "polyMesh.H"
#include "cloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(partialWrite, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::partialWrite::partialWrite
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::partialWrite::~partialWrite()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::partialWrite::read(const dictionary& dict)
{
    dict.lookup("objectNames") >> objectNames_;
    dict.lookup("writeInterval") >> writeInterval_;
    writeInstance_ = 0;

    Info<< type() << " " << name() << ":" << nl
        << "    dumping every outputTime :";
    forAllConstIter(HashSet<word>, objectNames_, iter)
    {
        Info<< ' ' << iter.key();
    }
    Info<< nl
        << "    dumping all other fields every "
        << writeInterval_ << "th outputTime" << nl
        << endl;

    if (writeInterval_ < 1)
    {
        FatalIOErrorIn("partialWrite::read(const dictionary&)", dict)
            << "Illegal value for writeInterval " << writeInterval_
            << ". It should be >= 1."
            << exit(FatalIOError);
    }
}


void Foam::partialWrite::execute()
{
    //Pout<< "execute at time " << obr_.time().timeName()
    //    << " index:" << obr_.time().timeIndex() << endl;
}


void Foam::partialWrite::end()
{
    //Pout<< "end at time " << obr_.time().timeName() << endl;
    // Do nothing - only valid on write
}


void Foam::partialWrite::write()
{
    //Pout<< "write at time " << obr_.time().timeName() << endl;
    if (obr_.time().outputTime())
    {
        // Above check so it can be used both with
        //  outputControl   timeStep;
        //  outputInterval  1;
        // or with
        //  outputControl   outputTime;

        writeInstance_++;

        if (writeInstance_ == writeInterval_)
        {
            // Normal dump
            writeInstance_ = 0;
        }
        else
        {
            // Delete all but marked objects
            fileName dbDir;
            if (isA<polyMesh>(obr_))
            {
                dbDir = dynamic_cast<const polyMesh&>(obr_).dbDir();
            }

            IOobjectList objects(obr_, obr_.time().timeName());

            if (debug)
            {
                Pout<< "For region:" << obr_.name() << endl;
            }

            forAllConstIter(IOobjectList, objects, iter)
            {
                if (!objectNames_.found(iter()->name()))
                {
                    const fileName f =
                        obr_.time().timePath()
                       /dbDir
                       /iter()->name();
                    if (debug)
                    {
                        Pout<< "   rm " << f << endl;
                    }
                    rm(f);
                }
            }

            // Do the lagrangian files as well.
            fileNameList cloudDirs
            (
                readDir
                (
                    obr_.time().timePath()/dbDir/cloud::prefix,
                    fileName::DIRECTORY
                )
            );
            forAll(cloudDirs, i)
            {
                if (debug)
                {
                    Pout<< "For cloud:" << cloudDirs[i] << endl;
                }

                IOobjectList sprayObjs
                (
                    obr_,
                    obr_.time().timeName(),
                    cloud::prefix/cloudDirs[i]
                );
                forAllConstIter(IOobjectList, sprayObjs, iter)
                {
                    if (!objectNames_.found(iter()->name()))
                    {
                        const fileName f =
                            obr_.time().timePath()
                           /dbDir
                           /cloud::prefix
                           /cloudDirs[i]
                           /iter()->name();
                        if (debug)
                        {
                            Pout<< "   rm " << f << endl;
                        }
                        rm(f);
                    }
                }
            }
        }
    }
}


// ************************************************************************* //
