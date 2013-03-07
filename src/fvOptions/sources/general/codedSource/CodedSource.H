/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

Class
    Foam::fv::codedSource

Description
    Constructs on-the-fly source

    \heading Source usage

    Example usage:
    \verbatim
    vectorCodedSourceCoeffs
    {
        fieldNames  (U);
        redirectType    ramp;

        codeCorrect
        #{
            Pout<< "**codeCorrect**" << endl;
        #};

        codeAddSup
        #{
            Pout<< "**codeAddSup**" << endl;
        #};

        codeSetValue
        #{
            Pout<< "**codeSetValue**" << endl;
        #};

        // Dummy entry. Make dependent on above to trigger recompilation
        code
        #{
            $codeCorrect
            $codeAddSup
            $codeSetValue
        #};
    }

    // Dummy entry
    rampCoeffs
    {}
    \endverbatim


SourceFiles
    codedSource.C

\*---------------------------------------------------------------------------*/

#ifndef CodedSource_H
#define CodedSource_H

#include "fvOption.H"
#include "codedBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                         Class codedSource Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class CodedSource
:
    public option,
    public codedBase
{

protected:

    // Protected data

        word redirectType_;

        string codeCorrect_;
        string codeAddSup_;
        string codeSetValue_;

        //- Underlying functionObject
        mutable autoPtr<option> redirectFvOptionPtr_;


    // Protected Member Functions

        //- get the loaded dynamic libraries
        virtual dlLibraryTable& libs() const;

        //- adapt the context for the current object
        virtual void prepare(dynamicCode&, const dynamicCodeContext&) const;

        // Return a description (type + name) for the output
        virtual string description() const;

        // Clear any redirected objects
        virtual void clearRedirect() const;

        // Get the dictionary to initialize the codeContext
        virtual const dictionary& codeDict() const;


public:

    //- Runtime type information
    TypeName("coded");


    // Constructors

        //- Construct from components
        CodedSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Member Functions

        //- Dynamically compiled fvOption
        option& redirectFvOption() const;

        // Evaluation

            //- Correct field
            virtual void correct
            (
                GeometricField<Type, fvPatchField, volMesh>&
            );

            //- Explicit and implicit matrix contributions
            virtual void addSup
            (
                fvMatrix<Type>& eqn,
                const label fieldI
            );

            //- Set value
            virtual void setValue
            (
                fvMatrix<Type>& eqn,
                const label fieldI
            );


        // I-O

            //- Write the source properties
            virtual void writeData(Ostream&) const;

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "CodedSource.C"
#   include "CodedSourceIO.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
