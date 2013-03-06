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
    Foam::multiphaseFixedFluxPressureFvPatchScalarField

Group
    grpWallBoundaryConditions grpGenericBoundaryConditions

Description
    This boundary condition adjusts the pressure gradient such that the flux
    on the boundary is that specified by the velocity boundary condition.

    The predicted flux to be compensated by the pressure gradient is evaluated
    as \f$(\phi - \phi_{H/A})\f$, both of which are looked-up from the database,
    as is the pressure diffusivity Dp used to calculate the gradient using:

        \f[
            \nabla(p) = \frac{\phi_{H/A} - \phi}{|Sf| Dp}
        \f]

    where

    \vartable
        \phi    | flux
        Dp      | pressure diffusivity
        Sf      | patch face areas [m2]
    \endvartable

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        phiHbyA      | name of predicted flux field | no | phiHbyA
        phi          | name of flux field      | no |   phi
        rho          | name of density field   | no |   rho
        Dp           | name of pressure diffusivity field | no |   Dp
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            multiphaseFixedFluxPressure;
        phiHbyA         phiHbyA;
        phi             phi;
        rho             rho;
        Dp              Dp;
    }
    \endverbatim

SourceFiles
    multiphaseFixedFluxPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef multiphaseFixedFluxPressureFvPatchScalarFields_H
#define multiphaseFixedFluxPressureFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class multiphaseFixedFluxPressureFvPatch Declaration
\*---------------------------------------------------------------------------*/

class multiphaseFixedFluxPressureFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    // Private data

        //- Name of the predicted flux transporting the field
        word phiHbyAName_;

        //- Name of the flux transporting the field
        word phiName_;

        //- Name of the density field used to normalise the mass flux
        //  if neccessary
        word rhoName_;

        //- Name of the pressure diffusivity field
        word DpName_;


public:

    //- Runtime type information
    TypeName("multiphaseFixedFluxPressure");


    // Constructors

        //- Construct from patch and internal field
        multiphaseFixedFluxPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        multiphaseFixedFluxPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  multiphaseFixedFluxPressureFvPatchScalarField onto a new patch
        multiphaseFixedFluxPressureFvPatchScalarField
        (
            const multiphaseFixedFluxPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        multiphaseFixedFluxPressureFvPatchScalarField
        (
            const multiphaseFixedFluxPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new multiphaseFixedFluxPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        multiphaseFixedFluxPressureFvPatchScalarField
        (
            const multiphaseFixedFluxPressureFvPatchScalarField&,
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
                new multiphaseFixedFluxPressureFvPatchScalarField(*this, iF)
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
