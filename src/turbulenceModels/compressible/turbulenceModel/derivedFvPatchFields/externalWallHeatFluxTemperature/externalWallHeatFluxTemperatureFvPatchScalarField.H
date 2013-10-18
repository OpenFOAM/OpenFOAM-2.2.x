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

Class
    externalWallHeatFluxTemperatureFvPatchScalarField

Group
    grpThermoBoundaryConditions grpWallBoundaryConditions

Description
    This boundary condition supplies a heat flux condition for temperature
    on an external wall. Optional thin thermal layer resistances can be
    specified through thicknessLayers and kappaLayers entries for the
    fixed heat transfer coefficient mode.

    The condition can operate in two modes:
        \li fixed heat transfer coefficient: supply h and Ta
        \li fixed heat flux: supply q

    where
        \var h  = heat transfer coefficient [W/m2/K]
        \var Ta = ambient temperature [K]
        \var q  = heat flux [W/m2]

    The thermal conductivity, \c kappa, can either be retrieved from the
    mesh database using the \c lookup option, or from a \c solidThermo
    thermophysical package.

    Example of the boundary condition specification:
    \verbatim
        myPatch
        {
            type            externalWallHeatFluxTemperature;
            kappa           fluidThermo;       // fluidThermo, solidThermo or
                                               // lookup
            q               uniform 1000;      // heat flux / [W/m2]
            Ta              uniform 300.0;     // ambient temperature /[K]
            h               uniform 10.0;      // heat transfer coeff /[W/K/m2]
            thicknessLayers (0.1 0.2 0.3 0.4); // thickness of layer [m]
            kappaLayers     (1 2 3 4)          // thermal conductivity of
                                               // layer [W/m/K]
            value           uniform 300.0;     // initial temperature / [K]
            kappaName       none;
        }
    \endverbatim

Note
    Only supply \c h and \c Ta, or \c q in the dictionary (see above)

SourceFiles
    externalWallHeatFluxTemperatureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef solidWallHeatFluxTemperatureFvPatchScalarField_H
#define solidWallHeatFluxTemperatureFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "temperatureCoupledBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class externalWallHeatFluxTemperatureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class externalWallHeatFluxTemperatureFvPatchScalarField
:
    public mixedFvPatchScalarField,
    public temperatureCoupledBase
{
public:

    // Public data


        //- Operation mode enumeration
        enum operationMode
        {
            fixedHeatFlux,
            fixedHeatTransferCoeff,
            unknown
        };
        static const NamedEnum<operationMode, 3> operationModeNames;


private:

    // Private data

        //- Operation mode
        operationMode oldMode_;

        //- Heat flux / [W/m2]
        scalarField q_;

        //- Heat transfer coefficient / [W/m2K]
        scalarField h_;

        //- Ambient temperature / [K]
        scalarField Ta_;

        //- Thickness of layers
        scalarList thicknessLayers_;

        //- Conductivity of layers
        scalarList kappaLayers_;


public:

    //- Runtime type information
    TypeName("externalWallHeatFluxTemperature");


    // Constructors

        //- Construct from patch and internal field
        externalWallHeatFluxTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        externalWallHeatFluxTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // externalWallHeatFluxTemperatureFvPatchScalarField
        // onto a new patch
        externalWallHeatFluxTemperatureFvPatchScalarField
        (
            const externalWallHeatFluxTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        externalWallHeatFluxTemperatureFvPatchScalarField
        (
            const externalWallHeatFluxTemperatureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new externalWallHeatFluxTemperatureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        externalWallHeatFluxTemperatureFvPatchScalarField
        (
            const externalWallHeatFluxTemperatureFvPatchScalarField&,
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
                new externalWallHeatFluxTemperatureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
