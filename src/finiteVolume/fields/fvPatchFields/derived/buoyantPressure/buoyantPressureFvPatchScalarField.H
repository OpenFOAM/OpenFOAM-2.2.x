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

Class
    Foam::buoyantPressureFvPatchScalarField

Group
    grpGenericBoundaryConditions grpWallBoundaryConditions

Description
    This boundary condition sets the pressure gradient appropriately for
    buoyant flow.

    If the variable name is one of:
    - \c pd
    - \c p_rgh
    - \c ph_rgh
    we assume that the pressure variable is \f$p - \rho(g \cdot h)\f$ and the
    gradient set using:
        \f[
            \nabla(p) = -\nabla_\perp(\rho)(g \cdot h)
        \f]

    where

    \vartable
        \rho   | density [kg/m3]
        g      | acceleration due to gravity [m/s2]
        h      | patch face centres [m]
    \endvartable

    Otherwise we assume that it is the static pressure, and the gradient
    calculated using:
        \f[
            \nabla(p) = \rho(g \cdot n)
        \f]

    where
    \vartable
        n      | patch face normal vectors
    \endvartable

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        rho          | density field name      | no          | rho
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            buoyantPressure;
        rho             rho;
        value           uniform 0;
    }
    \endverbatim

SourceFiles
    buoyantPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef buoyantPressureFvPatchScalarFields_H
#define buoyantPressureFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class buoyantPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class buoyantPressureFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    // Private data

        //- Name of the density field used to calculate the buoyancy force
        word rhoName_;


public:

    //- Runtime type information
    TypeName("buoyantPressure");


    // Constructors

        //- Construct from patch and internal field
        buoyantPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        buoyantPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  buoyantPressureFvPatchScalarField onto a new patch
        buoyantPressureFvPatchScalarField
        (
            const buoyantPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        buoyantPressureFvPatchScalarField
        (
            const buoyantPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new buoyantPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        buoyantPressureFvPatchScalarField
        (
            const buoyantPressureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new buoyantPressureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
