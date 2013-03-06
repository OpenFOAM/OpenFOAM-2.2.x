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

Class
    Foam::CoBlended

Description
    Two-scheme Courant number based blending differencing scheme.

    Similar to localBlended but uses a blending factor computed from the
    face-based Courant number and the alpha factor supplied:

    weight = 1 - Co/alpha

    The weight applies to the first scheme and 1-factor to the second scheme.

SourceFiles
    CoBlended.C

\*---------------------------------------------------------------------------*/

#ifndef CoBlended_H
#define CoBlended_H

#include "surfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class CoBlended Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class CoBlended
:
    public surfaceInterpolationScheme<Type>
{
    // Private data

        const scalar alpha_;

    // Private Member Functions

        //- Scheme 1
        tmp<surfaceInterpolationScheme<Type> > tScheme1_;

        //- Scheme 2
        tmp<surfaceInterpolationScheme<Type> > tScheme2_;

        //- The face-flux used to compute the face Courant number
        const surfaceScalarField& faceFlux_;


        //- Disallow default bitwise copy construct
        CoBlended(const CoBlended&);

        //- Disallow default bitwise assignment
        void operator=(const CoBlended&);


public:

    //- Runtime type information
    TypeName("CoBlended");


    // Constructors

        //- Construct from mesh and Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        CoBlended
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            alpha_(readScalar(is)),
            tScheme1_
            (
                surfaceInterpolationScheme<Type>::New(mesh, is)
            ),
            tScheme2_
            (
                surfaceInterpolationScheme<Type>::New(mesh, is)
            ),
            faceFlux_
            (
                mesh.lookupObject<surfaceScalarField>
                (
                    word(is)
                )
            )
        {
            if (alpha_ <= 0 )
            {
                FatalIOErrorIn("CoBlended(const fvMesh&, Istream&)", is)
                    << "coefficient = " << alpha_
                    << " should be > 0"
                    << exit(FatalIOError);
            }
        }


        //- Construct from mesh, faceFlux and Istream
        CoBlended
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            alpha_(readScalar(is)),
            tScheme1_
            (
                surfaceInterpolationScheme<Type>::New(mesh, faceFlux, is)
            ),
            tScheme2_
            (
                surfaceInterpolationScheme<Type>::New(mesh, faceFlux, is)
            ),
             faceFlux_(faceFlux)
        {
            if (alpha_ <= 0)
            {
                FatalIOErrorIn("CoBlended(const fvMesh&, Istream&)", is)
                    << "coefficient = " << alpha_
                    << " should be > 0"
                    << exit(FatalIOError);
            }
        }


    // Member Functions

        //- Return the face-based Courant number blending factor
        tmp<surfaceScalarField> blendingFactor() const
        {
            const fvMesh& mesh = faceFlux_.mesh();

            return
            (
                max
                (
                    scalar(1)
                  - mesh.time().deltaT()*mesh.deltaCoeffs()*mag(faceFlux_)
                   /(mesh.magSf()*alpha_),
                    scalar(0)
                )
            );
        }


        //- Return the interpolation weighting factors
        tmp<surfaceScalarField>
        weights
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const
        {
            surfaceScalarField bf(blendingFactor());

            return
                bf*tScheme1_().weights(vf)
              + (scalar(1.0) - bf)*tScheme2_().weights(vf);
        }


        //- Return the face-interpolate of the given cell field
        //  with explicit correction
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
        interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const
        {
            surfaceScalarField bf(blendingFactor());

            return
                bf*tScheme1_().interpolate(vf)
              + (scalar(1.0) - bf)*tScheme2_().interpolate(vf);
        }


        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return tScheme1_().corrected() || tScheme2_().corrected();
        }


        //- Return the explicit correction to the face-interpolate
        //  for the given field
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
        correction
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const
        {
            surfaceScalarField bf(blendingFactor());

            if (tScheme1_().corrected())
            {
                if (tScheme2_().corrected())
                {
                    return
                    (
                        bf
                      * tScheme1_().correction(vf)
                      + (scalar(1.0) - bf)
                      * tScheme2_().correction(vf)
                    );
                }
                else
                {
                    return
                    (
                        bf
                      * tScheme1_().correction(vf)
                    );
                }
            }
            else if (tScheme2_().corrected())
            {
                return
                (
                    (scalar(1.0) - bf)
                  * tScheme2_().correction(vf)
                );
            }
            else
            {
                return tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
                (
                    NULL
                );
            }
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
