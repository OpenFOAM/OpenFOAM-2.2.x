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

Class
    Foam::fv::boundedBackwardDdtScheme

Description
    Second-order bounded-backward-differencing ddt using the current and
    two previous time-step values.

SourceFiles
    boundedBackwardDdtScheme.C

\*---------------------------------------------------------------------------*/

#ifndef boundedBackwardDdtScheme_H
#define boundedBackwardDdtScheme_H

#include "ddtScheme.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class boundedBackwardDdtScheme Declaration
\*---------------------------------------------------------------------------*/

class boundedBackwardDdtScheme
:
    public fv::ddtScheme<scalar>
{
    // Private Member Functions

        //- Return the current time-step
        scalar deltaT_() const;

        //- Return the previous time-step
        scalar deltaT0_() const;

        //- Return the previous time-step or GREAT if the old timestep field
        //  wasn't available in which case Euler ddt is used
        template<class GeoField>
        scalar deltaT0_(const GeoField& vf) const
        {
            if (vf.oldTime().timeIndex() == vf.oldTime().oldTime().timeIndex())
            {
                return GREAT;
            }
            else
            {
                return deltaT0_();
            }
        }


        //- Disallow default bitwise copy construct
        boundedBackwardDdtScheme(const boundedBackwardDdtScheme&);

        //- Disallow default bitwise assignment
        void operator=(const boundedBackwardDdtScheme&);


public:

    //- Runtime type information
    TypeName("boundedBackward");


    // Constructors

        //- Construct from mesh
        boundedBackwardDdtScheme(const fvMesh& mesh)
        :
            ddtScheme<scalar>(mesh)
        {}

        //- Construct from mesh and Istream
        boundedBackwardDdtScheme(const fvMesh& mesh, Istream& is)
        :
            ddtScheme<scalar>(mesh, is)
        {}


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return fv::ddtScheme<scalar>::mesh();
        }

        tmp<volScalarField> fvcDdt
        (
            const dimensionedScalar&
        );

        tmp<volScalarField> fvcDdt
        (
            const volScalarField&
        );

        tmp<volScalarField> fvcDdt
        (
            const dimensionedScalar&,
            const volScalarField&
        );

        tmp<volScalarField> fvcDdt
        (
            const volScalarField&,
            const volScalarField&
        );

        tmp<fvScalarMatrix> fvmDdt
        (
            const volScalarField&
        );

        tmp<fvScalarMatrix> fvmDdt
        (
            const dimensionedScalar&,
            const volScalarField&
        );

        tmp<fvScalarMatrix> fvmDdt
        (
            const volScalarField&,
            const volScalarField&
        );

        tmp<surfaceScalarField> fvcDdtPhiCorr
        (
            const volScalarField& rA,
            const volScalarField& U,
            const surfaceScalarField& phi
        );

        tmp<surfaceScalarField> fvcDdtPhiCorr
        (
            const volScalarField& rA,
            const volScalarField& rho,
            const volScalarField& U,
            const surfaceScalarField& phi
        );

        tmp<surfaceScalarField> meshPhi
        (
            const volScalarField&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
