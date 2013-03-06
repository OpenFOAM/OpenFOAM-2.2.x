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

Description
    Template for use with dynamic code generation of a basicSource

    - without state

SourceFiles
    codedBasicSourceTemplate.C

\*---------------------------------------------------------------------------*/

#ifndef codedBasicSourceTemplate_H
#define codedBasicSourceTemplate_H

#include "basicSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         A templated CodedBasicSource
\*---------------------------------------------------------------------------*/

class ${typeName}CodedBasic${SourceType}
:
    public basicSource
{
public:

    //- Information about the SHA1 of the code itself
    static const char* const SHA1sum;

    //- Runtime type information
    TypeName("${typeName}");


    // Constructors

        //- Construct from patch and internal field
        ${typeName}CodedBasic${SourceType}
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

    //- Destructor
    virtual ~${typeName}CodedBasic${SourceType}();


    // Member functions

            //- Correct field
            virtual void correct
            (
                GeometricField<${TemplateType}, fvPatchField, volMesh>&
            );

            //- Explicit and implicit matrix contributions
            virtual void addSup
            (
                fvMatrix<${TemplateType}>& eqn,
                const label fieldI
            );

            //- Set value
            virtual void setValue
            (
                fvMatrix<${TemplateType}>& eqn,
                const label fieldI
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
