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

\*---------------------------------------------------------------------------*/

#include "coupledFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
coupledFvPatchField<Type>::coupledFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    LduInterfaceField<Type>(refCast<const lduInterface>(p)),
    fvPatchField<Type>(p, iF)
{}


template<class Type>
coupledFvPatchField<Type>::coupledFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Field<Type>& f
)
:
    LduInterfaceField<Type>(refCast<const lduInterface>(p)),
    fvPatchField<Type>(p, iF, f)
{}


template<class Type>
coupledFvPatchField<Type>::coupledFvPatchField
(
    const coupledFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    LduInterfaceField<Type>(refCast<const lduInterface>(p)),
    fvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
coupledFvPatchField<Type>::coupledFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    LduInterfaceField<Type>(refCast<const lduInterface>(p)),
    fvPatchField<Type>(p, iF, dict)
{}


template<class Type>
coupledFvPatchField<Type>::coupledFvPatchField
(
    const coupledFvPatchField<Type>& ptf
)
:
    LduInterfaceField<Type>(refCast<const lduInterface>(ptf.patch())),
    fvPatchField<Type>(ptf)
{}


template<class Type>
coupledFvPatchField<Type>::coupledFvPatchField
(
    const coupledFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    LduInterfaceField<Type>(refCast<const lduInterface>(ptf.patch())),
    fvPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<Field<Type> > coupledFvPatchField<Type>::snGrad() const
{
    return
        (this->patchNeighbourField() - this->patchInternalField())
       *this->patch().deltaCoeffs();
}


template<class Type>
void coupledFvPatchField<Type>::initEvaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }
}


template<class Type>
void coupledFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    Field<Type>::operator=
    (
        this->patch().weights()*this->patchInternalField()
      + (1.0 - this->patch().weights())*this->patchNeighbourField()
    );

    fvPatchField<Type>::evaluate();
}


template<class Type>
tmp<Field<Type> > coupledFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>& w
) const
{
    return Type(pTraits<Type>::one)*w;
}

template<class Type>
tmp<Field<Type> > coupledFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>& w
) const
{
    return Type(pTraits<Type>::one)*(1.0 - w);
}

template<class Type>
tmp<Field<Type> > coupledFvPatchField<Type>::gradientInternalCoeffs() const
{
    return -Type(pTraits<Type>::one)*this->patch().deltaCoeffs();
}


template<class Type>
tmp<Field<Type> > coupledFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    return -this->gradientInternalCoeffs();
}


template<class Type>
void coupledFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
