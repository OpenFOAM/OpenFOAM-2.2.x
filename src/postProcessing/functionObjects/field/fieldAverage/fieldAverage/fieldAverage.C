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

#include "fieldAverage.H"
#include "volFields.H"
#include "Time.H"

#include "fieldAverageItem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(fieldAverage, 0);

const word fieldAverage::EXT_MEAN = "Mean";
const word fieldAverage::EXT_PRIME2MEAN = "Prime2Mean";
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fieldAverage::resetFields(wordList& names)
{
    forAll(names, fieldI)
    {
        if (names[fieldI].size())
        {
            obr_.checkOut(*obr_[names[fieldI]]);
        }
    }

    names.clear();
    names.setSize(faItems_.size());
}


void Foam::fieldAverage::initialize()
{
    resetFields(meanScalarFields_);
    resetFields(meanVectorFields_);
    resetFields(meanSphericalTensorFields_);
    resetFields(meanSymmTensorFields_);
    resetFields(meanTensorFields_);

    resetFields(prime2MeanScalarFields_);
    resetFields(prime2MeanSymmTensorFields_);


    // Add mean fields to the field lists
    forAll(faItems_, fieldI)
    {
        const word& fieldName = faItems_[fieldI].fieldName();
        if (obr_.foundObject<volScalarField>(fieldName))
        {
            addMeanField<scalar>(fieldI, meanScalarFields_);
        }
        else if (obr_.foundObject<volVectorField>(fieldName))
        {
            addMeanField<vector>(fieldI, meanVectorFields_);
        }
        else if (obr_.foundObject<volSphericalTensorField>(fieldName))
        {
            addMeanField<sphericalTensor>(fieldI, meanSphericalTensorFields_);
        }
        else if (obr_.foundObject<volSymmTensorField>(fieldName))
        {
            addMeanField<symmTensor>(fieldI, meanSymmTensorFields_);
        }
        else if (obr_.foundObject<volTensorField>(fieldName))
        {
            addMeanField<tensor>(fieldI, meanTensorFields_);
        }
        else
        {
            FatalErrorIn("Foam::fieldAverage::initialize()")
                << "Requested field " << faItems_[fieldI].fieldName()
                << " does not exist in the database" << nl
                << exit(FatalError);
        }
    }

    // Add prime-squared mean fields to the field lists
    forAll(faItems_, fieldI)
    {
        if (faItems_[fieldI].prime2Mean())
        {
            const word& fieldName = faItems_[fieldI].fieldName();
            if (!faItems_[fieldI].mean())
            {
                FatalErrorIn("Foam::fieldAverage::initialize()")
                    << "To calculate the prime-squared average, the "
                    << "mean average must also be selected for field "
                    << fieldName << nl << exit(FatalError);
            }

            if (obr_.foundObject<volScalarField>(fieldName))
            {
                addPrime2MeanField<scalar, scalar>
                (
                    fieldI,
                    meanScalarFields_,
                    prime2MeanScalarFields_
                );
            }
            else if (obr_.foundObject<volVectorField>(fieldName))
            {
                addPrime2MeanField<vector, symmTensor>
                (
                    fieldI,
                    meanVectorFields_,
                    prime2MeanSymmTensorFields_
                );
            }
            else
            {
                FatalErrorIn("Foam::fieldAverage::initialize()")
                    << "prime2Mean average can only be applied to "
                    << "volScalarFields and volVectorFields"
                    << nl << "    Field: " << fieldName << nl
                    << exit(FatalError);
            }
        }
    }

    // ensure first averaging works unconditionally
    prevTimeIndex_ = -1;

    initialised_ = true;
}


void Foam::fieldAverage::calcAverages()
{
    if (!initialised_)
    {
        initialize();
    }

    const label currentTimeIndex =
        static_cast<const fvMesh&>(obr_).time().timeIndex();

    if (prevTimeIndex_ == currentTimeIndex)
    {
        return;
    }
    else
    {
        prevTimeIndex_ = currentTimeIndex;
    }

    Info<< "Calculating averages" << nl << endl;

    addMeanSqrToPrime2Mean<scalar, scalar>
    (
        meanScalarFields_,
        prime2MeanScalarFields_
    );
    addMeanSqrToPrime2Mean<vector, symmTensor>
    (
        meanVectorFields_,
        prime2MeanSymmTensorFields_
    );

    calculateMeanFields<scalar>(meanScalarFields_);
    calculateMeanFields<vector>(meanVectorFields_);
    calculateMeanFields<sphericalTensor>(meanSphericalTensorFields_);
    calculateMeanFields<symmTensor>(meanSymmTensorFields_);
    calculateMeanFields<tensor>(meanTensorFields_);

    calculatePrime2MeanFields<scalar, scalar>
    (
        meanScalarFields_,
        prime2MeanScalarFields_
    );
    calculatePrime2MeanFields<vector, symmTensor>
    (
        meanVectorFields_,
        prime2MeanSymmTensorFields_
    );

    forAll(faItems_, fieldI)
    {
        totalIter_[fieldI]++;
        totalTime_[fieldI] += obr_.time().deltaTValue();
    }
}


void Foam::fieldAverage::writeAverages() const
{
    writeFieldList<scalar>(meanScalarFields_);
    writeFieldList<vector>(meanVectorFields_);
    writeFieldList<sphericalTensor>(meanSphericalTensorFields_);
    writeFieldList<symmTensor>(meanSymmTensorFields_);
    writeFieldList<tensor>(meanTensorFields_);

    writeFieldList<scalar>(prime2MeanScalarFields_);
    writeFieldList<symmTensor>(prime2MeanSymmTensorFields_);
}


void Foam::fieldAverage::writeAveragingProperties() const
{
    IOdictionary propsDict
    (
        IOobject
        (
            "fieldAveragingProperties",
            obr_.time().timeName(),
            "uniform",
            obr_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    forAll(faItems_, fieldI)
    {
        const word& fieldName = faItems_[fieldI].fieldName();
        propsDict.add(fieldName, dictionary());
        propsDict.subDict(fieldName).add("totalIter", totalIter_[fieldI]);
        propsDict.subDict(fieldName).add("totalTime", totalTime_[fieldI]);
    }

    propsDict.regIOobject::write();
}


void Foam::fieldAverage::readAveragingProperties()
{
    totalIter_.clear();
    totalIter_.setSize(faItems_.size(), 1);

    totalTime_.clear();
    totalTime_.setSize(faItems_.size(), obr_.time().deltaTValue());

    if (resetOnRestart_)
    {
        Info<< "fieldAverage: starting averaging at time "
            << obr_.time().timeName() << nl << endl;
    }
    else
    {
        IOobject propsDictHeader
        (
            "fieldAveragingProperties",
            obr_.time().timeName(obr_.time().startTime().value()),
            "uniform",
            obr_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        );

        if (!propsDictHeader.headerOk())
        {
            Info<< "fieldAverage: starting averaging at time "
                << obr_.time().timeName() << nl << endl;
            return;
        }

        IOdictionary propsDict(propsDictHeader);

        Info<< "fieldAverage: restarting averaging for fields:" << endl;
        forAll(faItems_, fieldI)
        {
            const word& fieldName = faItems_[fieldI].fieldName();
            if (propsDict.found(fieldName))
            {
                dictionary fieldDict(propsDict.subDict(fieldName));

                totalIter_[fieldI] = readLabel(fieldDict.lookup("totalIter"));
                totalTime_[fieldI] = readScalar(fieldDict.lookup("totalTime"));
                Info<< "    " << fieldName
                    << " iters = " << totalIter_[fieldI]
                    << " time = " << totalTime_[fieldI] << endl;
            }
        }

        Info<< endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldAverage::fieldAverage
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    active_(true),
    prevTimeIndex_(-1),
    resetOnRestart_(false),
    resetOnOutput_(false),
    initialised_(false),
    faItems_(),
    meanScalarFields_(),
    meanVectorFields_(),
    meanSphericalTensorFields_(),
    meanSymmTensorFields_(),
    meanTensorFields_(),
    prime2MeanScalarFields_(),
    prime2MeanSymmTensorFields_(),
    totalIter_(),
    totalTime_()
{
    // Only active if a fvMesh is available
    if (isA<fvMesh>(obr_))
    {
        read(dict);
    }
    else
    {
        active_ = false;
        WarningIn
        (
            "fieldAverage::fieldAverage\n"
            "(\n"
                "const word&,\n"
                "const objectRegistry&,\n"
                "const dictionary&,\n"
                "const bool\n"
            ")"
        )   << "No fvMesh available, deactivating."
            << nl << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fieldAverage::~fieldAverage()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fieldAverage::read(const dictionary& dict)
{
    if (active_)
    {
        initialised_ = false;

        dict.readIfPresent("resetOnRestart", resetOnRestart_);
        dict.readIfPresent("resetOnOutput", resetOnOutput_);
        dict.lookup("fields") >> faItems_;

        readAveragingProperties();
    }
}


void Foam::fieldAverage::execute()
{
    if (active_)
    {
        calcAverages();
    }
}


void Foam::fieldAverage::end()
{}


void Foam::fieldAverage::timeSet()
{
    // Do nothing
}


void Foam::fieldAverage::write()
{
    if (active_)
    {
        calcAverages();
        writeAverages();
        writeAveragingProperties();

        if (resetOnOutput_)
        {
            Info<< "fieldAverage: restarting averaging at time "
                << obr_.time().timeName() << nl << endl;

            initialize();

            // ensure first averaging works unconditionally
            prevTimeIndex_ = -1;
        }
    }
}


void Foam::fieldAverage::updateMesh(const mapPolyMesh&)
{
    // Do nothing
}


void Foam::fieldAverage::movePoints(const polyMesh&)
{
    // Do nothing
}


// ************************************************************************* //
