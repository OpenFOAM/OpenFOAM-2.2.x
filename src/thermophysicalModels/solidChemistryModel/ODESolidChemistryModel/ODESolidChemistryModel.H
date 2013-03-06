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
    Foam::ODESolidChemistryModel

Description
    Extends base chemistry model by adding a thermo package, and ODE functions.
    Introduces chemistry equation system and evaluation of chemical source
    terms.

SourceFiles
    ODESolidChemistryModelI.H
    ODESolidChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef ODESolidChemistryModel_H
#define ODESolidChemistryModel_H

#include "Reaction.H"
#include "ODE.H"
#include "volFieldsFwd.H"
#include "DimensionedField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                   Class ODESolidChemistryModel Declaration
\*---------------------------------------------------------------------------*/

template<class CompType, class SolidThermo, class GasThermo>
class ODESolidChemistryModel
:
    public CompType,
    public ODE
{
    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const ODESolidChemistryModel&);


protected:

        //- Reference to solid mass fractions
        PtrList<volScalarField>& Ys_;

        //- Reactions
        const PtrList<Reaction<SolidThermo> >& reactions_;

        //- List of gas species present in reaction system
        speciesTable pyrolisisGases_;

        //- Thermodynamic data of solids
        const PtrList<SolidThermo>& solidThermo_;

        //- Thermodynamic data of gases
        PtrList<GasThermo> gasThermo_;

         //- Number of gas species
        label nGases_;

        //- Number of components being solved by ODE
        label nSpecie_;

        //- Number of solid components
        label nSolids_;

        //- Number of solid reactions
        label nReaction_;

        //- List of reaction rate per solid [kg/m3/s]
        PtrList<DimensionedField<scalar, volMesh> > RRs_;

        //- List of reaction rate per gas [kg/m3/s]
        PtrList<DimensionedField<scalar, volMesh> > RRg_;


    // Protected Member Functions

        //- Write access to source terms for solids
        inline PtrList<DimensionedField<scalar, volMesh> >& RRs();

        //- Write access to source terms for gases
        inline PtrList<DimensionedField<scalar, volMesh> >& RRg();


private:

        //- List of accumulative solid concentrations
        mutable PtrList<volScalarField> Ys0_;

        //- Cell counter
        label cellCounter_;

        //- List of active reacting cells
        List<bool> reactingCells_;


        // Private members

            //- Set reacting status of cell, cellI
            void setCellReacting(const label cellI, const bool active);


public:

    //- Runtime type information
    TypeName("ODESolidChemistryModel");


    // Constructors

        //- Construct from mesh
        ODESolidChemistryModel(const fvMesh& mesh);


    //- Destructor
    virtual ~ODESolidChemistryModel();


    // Member Functions

        //- The reactions
        inline const PtrList<Reaction<SolidThermo> >& reactions() const;

        //- Thermodynamic data of gases
        inline const PtrList<GasThermo>& gasThermo() const;

        //- Gases table
        inline const speciesTable& gasTable() const;

        //- The number of solids
        inline label nSpecie() const;

        //- The number of solids
        inline label nGases() const;

        //- The number of reactions
        inline label nReaction() const;


        //- dc/dt = omega, rate of change in concentration, for each species
        virtual scalarField omega
        (
            const scalarField& c,
            const scalar T,
            const scalar p,
            const bool updateC0 = false
        ) const;

        //- Return the reaction rate for reaction r and the reference
        //  species and charateristic times
        virtual scalar omega
        (
            const Reaction<SolidThermo>& r,
            const scalarField& c,
            const scalar T,
            const scalar p,
            scalar& pf,
            scalar& cf,
            label& lRef,
            scalar& pr,
            scalar& cr,
            label& rRef
        ) const;

        //- Calculates the reaction rates
        virtual void calculate();


        // Chemistry model functions

            //- Return const access to the chemical source terms for solids
            inline const DimensionedField<scalar, volMesh>& RRs
            (
                const label i
            ) const;

            //- Return const access to the chemical source terms for gases
            inline const DimensionedField<scalar, volMesh>& RRg
            (
                const label i
            ) const;

            //- Return total gas source term
            inline tmp<DimensionedField<scalar, volMesh> > RRg() const;

            //- Return total solid source term
            inline tmp<DimensionedField<scalar, volMesh> > RRs() const;

            //- Return const access to the total source terms
            inline const DimensionedField<scalar, volMesh>& RR
            (
                const label i
            ) const;

            //- Return sensible enthalpy for gas i [J/Kg]
            virtual tmp<volScalarField> gasHs
            (
                const volScalarField& p,
                const volScalarField& T,
                const label i
            ) const;

            //- Solve the reaction system for the given start time and time
            //  step and return the characteristic time
            virtual scalar solve(const scalar t0, const scalar deltaT);

            //- Return the chemical time scale
            virtual tmp<volScalarField> tc() const;

            //- Return source for enthalpy equation [kg/m/s3]
            virtual tmp<volScalarField> Sh() const;

            //- Return the heat release, i.e. enthalpy/sec [m2/s3]
            virtual tmp<volScalarField> dQ() const;


        // ODE functions (overriding abstract functions in ODE.H)

            //- Number of ODE's to solve
            virtual label nEqns() const;

            virtual void derivatives
            (
                const scalar t,
                const scalarField& c,
                scalarField& dcdt
            ) const;

            virtual void jacobian
            (
                const scalar t,
                const scalarField& c,
                scalarField& dcdt,
                scalarSquareMatrix& dfdc
            ) const;

            virtual scalar solve
            (
                scalarField &c,
                const scalar T,
                const scalar p,
                const scalar t0,
                const scalar dt
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "ODESolidChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ODESolidChemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
