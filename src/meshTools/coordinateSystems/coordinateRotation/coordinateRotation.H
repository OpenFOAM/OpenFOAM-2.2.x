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
    Foam::coordinateRotation

Description
    A coordinate rotation specified per local axes and the base class for
    other rotation specifications

    The rotation is defined by a combination of local vectors (e1/e2), (e2/e3)
    or (e3/e1). Any nonorthogonality will be absorbed into the second vector.

    For convenience, the dictionary constructor forms allow a few shortcuts:
    - if the \c type is not otherwise specified, the type \c axes
      is implicit
    - if an axes specification (eg, e3/e1) is used, the coordinateRotation
      sub-dictionary can be dropped.

    Specifying the rotation by an EulerCoordinateRotation
    (type "EulerRotation") or by a STARCDCoordinateRotation
    (type "STARCDRotation") requires the coordinateRotation sub-dictionary.

    \verbatim
        coordinateRotation
        {
            type        STARCDRotation
            rotation    (0 0 90);
        }
    \endverbatim

    - the rotation angles are in degrees, unless otherwise explictly specified:

    \verbatim
        coordinateRotation
        {
            type        STARCDRotation
            degrees     false;
            rotation    (0 0 3.141592654);
        }
    \endverbatim

Deprecated
    Specifying the local vectors as an \c axis (corresponding to e3) and a
    \c direction (corresponding to e1), is allowed for backwards
    compatibility, but this terminology is generally a bit confusing.
    (deprecated Apr 2008)

\*---------------------------------------------------------------------------*/

#ifndef coordinateRotation_H
#define coordinateRotation_H

#include "vector.H"
#include "tensor.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class coordinateRotation Declaration
\*---------------------------------------------------------------------------*/

class coordinateRotation
:
    public tensor
{
    // Private data

        //- the combination of local axes to be used
        enum axisOrder
        {
            e1e2,
            e2e3,
            e3e1
        };

    // Private Member Functions

        //- Calculate transformation tensor
        void calcTransform
        (
            const vector& axis1,
            const vector& axis2,
            const axisOrder& order = e3e1
        );

public:

    //- Runtime type information
    TypeName("coordinateRotation");

    // Constructors

        //- Construct null
        coordinateRotation();

        //- Construct from 2 axes
        coordinateRotation
        (
            const vector& axis,
            const vector& dir
        );

        //- Construct from dictionary
        coordinateRotation(const dictionary&);

        //- Return clone
        autoPtr<coordinateRotation> clone() const
        {
            return autoPtr<coordinateRotation>(new coordinateRotation(*this));
        }

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            coordinateRotation,
            dictionary,
            (
                const dictionary& dict
            ),
            (dict)
        );


    // Selectors

        //- Select constructed from Istream
        static autoPtr<coordinateRotation> New
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~coordinateRotation()
    {}


    // Member Functions

        //- Reset rotation to an identity rotation
        virtual void clear();

        //- Return local-to-global transformation tensor
        const tensor& R() const
        {
            return (*this);
        }

        //- Return local Cartesian x-axis
        const vector e1() const
        {
            return tensor::T().x();
        }

        //- Return local Cartesian y-axis
        const vector e2() const
        {
            return tensor::T().y();
        }

        //- Return local Cartesian z-axis
        const vector e3() const
        {
            return tensor::T().z();
        }


    // Member Operators

        //- assign from dictionary
        void operator=(const dictionary&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
