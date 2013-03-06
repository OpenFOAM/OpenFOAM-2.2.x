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

#include "codedBasicSourceTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
//{{{ begin codeInclude
${codeInclude}
//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode
${localCode}
//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = ${SHA1sum}
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void ${typeName}_${SHA1sum}(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

//makeRemovablePatchTypeField
//(
//    fvPatch${FieldType},
//    ${typeName}CodedBasic${SourceType}
//);
defineTypeNameAndDebug(${typeName}CodedBasic${SourceType}, 0);
addRemovableToRunTimeSelectionTable
(
    basicSource,
    ${typeName}CodedBasic${SourceType},
    dictionary
);


const char* const ${typeName}CodedBasic${SourceType}::SHA1sum =
    "${SHA1sum}";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

${typeName}CodedBasic${SourceType}::
${typeName}CodedBasic${SourceType}
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    basicSource(name, modelType, dict, mesh)
{
    if (${verbose:-false})
    {
        Info<<"construct ${typeName} sha1: ${SHA1sum}"
            " from components\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

${typeName}CodedBasic${SourceType}::
~${typeName}CodedBasic${SourceType}()
{
    if (${verbose:-false})
    {
        Info<<"destroy ${typeName} sha1: ${SHA1sum}\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ${typeName}CodedBasic${SourceType}::correct
(
    GeometricField<${TemplateType}, fvPatchField, volMesh>& fld
)
{
    if (${verbose:-false})
    {
        Info<<"${typeName}CodedBasic${SourceType}::correct()\n";
    }

//{{{ begin code
    ${codeCorrect}
//}}} end code
}


void ${typeName}CodedBasic${SourceType}::addSup
(
    fvMatrix<${TemplateType}>& eqn,
    const label fieldI
)
{
    if (${verbose:-false})
    {
        Info<<"${typeName}CodedBasic${SourceType}::addSup()\n";
    }

//{{{ begin code
    ${codeAddSup}
//}}} end code
}


void ${typeName}CodedBasic${SourceType}::setValue
(
    fvMatrix<${TemplateType}>& eqn,
    const label fieldI
)
{
    if (${verbose:-false})
    {
        Info<<"${typeName}CodedBasic${SourceType}::setValue()\n";
    }

//{{{ begin code
    ${codeSetValue}
//}}} end code
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
