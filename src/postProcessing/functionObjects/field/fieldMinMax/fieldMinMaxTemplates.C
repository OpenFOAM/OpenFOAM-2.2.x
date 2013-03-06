/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "fieldMinMax.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::fieldMinMax::calcMinMaxFields
(
    const word& fieldName,
    const modeType& mode
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (obr_.foundObject<fieldType>(fieldName))
    {
        const label procI = Pstream::myProcNo();

        const fieldType& field = obr_.lookupObject<fieldType>(fieldName);
        const fvMesh& mesh = field.mesh();

        const volVectorField::GeometricBoundaryField& CfBoundary =
            mesh.C().boundaryField();

        switch (mode)
        {
            case mdMag:
            {
                const volScalarField magField(mag(field));
                const volScalarField::GeometricBoundaryField& magFieldBoundary =
                    magField.boundaryField();

                scalarList minVs(Pstream::nProcs());
                List<vector> minCs(Pstream::nProcs());
                label minProcI = findMin(magField);
                minVs[procI] = magField[minProcI];
                minCs[procI] = field.mesh().C()[minProcI];


                labelList maxIs(Pstream::nProcs());
                scalarList maxVs(Pstream::nProcs());
                List<vector> maxCs(Pstream::nProcs());
                label maxProcI = findMax(magField);
                maxVs[procI] = magField[maxProcI];
                maxCs[procI] = field.mesh().C()[maxProcI];

                forAll(magFieldBoundary, patchI)
                {
                    const scalarField& mfp = magFieldBoundary[patchI];
                    if (mfp.size())
                    {
                        const vectorField& Cfp = CfBoundary[patchI];

                        label minPI = findMin(mfp);
                        if (mfp[minPI] < minVs[procI])
                        {
                            minVs[procI] = mfp[minPI];
                            minCs[procI] = Cfp[minPI];
                        }

                        label maxPI = findMax(mfp);
                        if (mfp[maxPI] > maxVs[procI])
                        {
                            maxVs[procI] = mfp[maxPI];
                            maxCs[procI] = Cfp[maxPI];
                        }
                    }
                }

                Pstream::gatherList(minVs);
                Pstream::gatherList(minCs);

                Pstream::gatherList(maxVs);
                Pstream::gatherList(maxCs);

                if (Pstream::master())
                {
                    label minI = findMin(minVs);
                    scalar minValue = minVs[minI];
                    const vector& minC = minCs[minI];

                    label maxI = findMax(maxVs);
                    scalar maxValue = maxVs[maxI];
                    const vector& maxC = maxCs[maxI];

                    file()
                        << obr_.time().value() << token::TAB
                        << fieldName << token::TAB
                        << minValue << token::TAB << minC;

                    if (Pstream::parRun())
                    {
                        file() << token::TAB << minI;
                    }

                    file() << token::TAB << maxValue << token::TAB << maxC;

                    if (Pstream::parRun())
                    {
                        file() << token::TAB << maxI;
                    }

                    file() << endl;

                    if (log_)
                    {
                        Info<< "    min(mag(" << fieldName << ")) = "
                            << minValue << " at position " << minC;

                        if (Pstream::parRun())
                        {
                            Info<< " on processor " << minI;
                        }

                        Info<< nl << "    max(mag(" << fieldName << ")) = "
                            << maxValue << " at position " << maxC;

                        if (Pstream::parRun())
                        {
                            Info<< " on processor " << maxI;
                        }

                        Info<< endl;
                    }
                }
                break;
            }
            case mdCmpt:
            {
                const typename fieldType::GeometricBoundaryField&
                    fieldBoundary = field.boundaryField();

                List<Type> minVs(Pstream::nProcs());
                List<vector> minCs(Pstream::nProcs());
                label minProcI = findMin(field);
                minVs[procI] = field[minProcI];
                minCs[procI] = field.mesh().C()[minProcI];

                Pstream::gatherList(minVs);
                Pstream::gatherList(minCs);

                List<Type> maxVs(Pstream::nProcs());
                List<vector> maxCs(Pstream::nProcs());
                label maxProcI = findMax(field);
                maxVs[procI] = field[maxProcI];
                maxCs[procI] = field.mesh().C()[maxProcI];

                forAll(fieldBoundary, patchI)
                {
                    const Field<Type>& fp = fieldBoundary[patchI];
                    if (fp.size())
                    {
                        const vectorField& Cfp = CfBoundary[patchI];

                        label minPI = findMin(fp);
                        if (fp[minPI] < minVs[procI])
                        {
                            minVs[procI] = fp[minPI];
                            minCs[procI] = Cfp[minPI];
                        }

                        label maxPI = findMax(fp);
                        if (fp[maxPI] > maxVs[procI])
                        {
                            maxVs[procI] = fp[maxPI];
                            maxCs[procI] = Cfp[maxPI];
                        }
                    }
                }

                Pstream::gatherList(minVs);
                Pstream::gatherList(minCs);

                Pstream::gatherList(maxVs);
                Pstream::gatherList(maxCs);

                if (Pstream::master())
                {
                    label minI = findMin(minVs);
                    Type minValue = minVs[minI];
                    const vector& minC = minCs[minI];

                    label maxI = findMax(maxVs);
                    Type maxValue = maxVs[maxI];
                    const vector& maxC = maxCs[maxI];

                    file()
                        << obr_.time().value() << token::TAB
                        << fieldName << token::TAB
                        << minValue << token::TAB << minC;

                    if (Pstream::parRun())
                    {
                        file() << token::TAB << minI;
                    }

                    file() << token::TAB << maxValue << token::TAB << maxC;

                    if (Pstream::parRun())
                    {
                        file() << token::TAB << maxI;
                    }

                    file() << endl;

                    if (log_)
                    {
                        Info<< "    min(" << fieldName << ") = "
                            << minValue << " at position " << minC;

                        if (Pstream::parRun())
                        {
                            Info<< " on processor " << minI;
                        }

                        Info<< nl << "    max(" << fieldName << ") = "
                            << maxValue << " at position " << maxC;

                        if (Pstream::parRun())
                        {
                            Info<< " on processor " << maxI;
                        }

                        Info<< endl;
                    }
                }
                break;
            }
            default:
            {
                FatalErrorIn
                (
                    "Foam::fieldMinMax::calcMinMaxFields"
                    "("
                        "const word&, "
                        "const modeType&"
                    ")"
                )
                    << "Unknown min/max mode: " << modeTypeNames_[mode_]
                    << exit(FatalError);
            }
        }
    }
}


// ************************************************************************* //
