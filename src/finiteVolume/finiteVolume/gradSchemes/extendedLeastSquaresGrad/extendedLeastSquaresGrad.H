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
    Foam::fv::extendedLeastSquaresGrad

Description
    Second-order gradient scheme using least-squares.

SourceFiles
    extendedLeastSquaresGrad.C

\*---------------------------------------------------------------------------*/

#ifndef extendedLeastSquaresGrad_H
#define extendedLeastSquaresGrad_H

#include "gradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class extendedLeastSquaresGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class extendedLeastSquaresGrad
:
    public fv::gradScheme<Type>
{
    // Private Data

        //- Minimum determinant criterion to choose extra cells
        scalar minDet_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        extendedLeastSquaresGrad(const extendedLeastSquaresGrad&);

        //- Disallow default bitwise assignment
        void operator=(const extendedLeastSquaresGrad&);


public:

    //- Runtime type information
    TypeName("extendedLeastSquares");


    // Constructors

        //- Construct from Istream
        extendedLeastSquaresGrad(const fvMesh& mesh, Istream& schemeData)
        :
            gradScheme<Type>(mesh),
            minDet_(readScalar(schemeData))
        {
            if (minDet_ < 0)    //-for facearea weighted: || minDet_ > 8)
            {
                FatalIOErrorIn
                (
                    "extendedLeastSquaresGrad"
                    "(const fvMesh&, Istream& schemeData)",
                    schemeData
                )   << "Minimum determinant = " << minDet_
                    << " should be >= 0"    // and <= 8"
                    << exit(FatalIOError);
            }
        }


    // Member Functions

        //- Return the gradient of the given field to the gradScheme::grad
        //  for optional caching
        virtual tmp
        <
            GeometricField
            <typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
        > calcGrad
        (
            const GeometricField<Type, fvPatchField, volMesh>& vsf,
            const word& name
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "extendedLeastSquaresGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
