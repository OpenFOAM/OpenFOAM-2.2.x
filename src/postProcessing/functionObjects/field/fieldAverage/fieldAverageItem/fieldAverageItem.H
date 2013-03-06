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
    Foam::fieldAverageItem

Description
    Helper class to describe what form of averaging to apply.  A set will be
    applied to each base field in Foam::fieldAverage, of the form:

    \verbatim
    {
        mean            on;
        prime2Mean      on;
        base            time; // iteration
        window          200;  // optional averaging window
        windowName      w1;   // optional window name (default = "")
    }
    \endverbatim

    The averaging window corresponds to the averaging interval (iters or time)
    If not specified, the averaging is over 'all iters/time'

SourceFiles
    fieldAverageItem.C
    fieldAverageItemIO.C

\*---------------------------------------------------------------------------*/

#ifndef fieldAverageItem_H
#define fieldAverageItem_H

#include "NamedEnum.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class fieldAverageItem;
Istream& operator>>(Istream&, fieldAverageItem&);
Ostream& operator<<(Ostream&, const fieldAverageItem&);


/*---------------------------------------------------------------------------*\
                      Class fieldAverageItem Declaration
\*---------------------------------------------------------------------------*/

class fieldAverageItem
{
public:

    // Public data

        //- Enumeration defining the averaging base type
        enum baseType
        {
            ITER,
            TIME
        };


private:

    // Private data

        //- Field name
        word fieldName_;

        //- Compute mean flag
        Switch mean_;

        //- Compute prime-squared mean flag
        Switch prime2Mean_;

        //- Averaging base type names
        static const NamedEnum<baseType, 2> baseTypeNames_;

        //- Averaging base type
        baseType base_;

        //- Averaging window - defaults to -1 for 'all iters/time'
        scalar window_;

        //- Averaging window name - defaults to 'window'
        word windowName_;


public:

    // Constructors

        //- Construct null
        fieldAverageItem();

        //- Construct from Istream
        fieldAverageItem(Istream&);

        //- Construct as copy
        fieldAverageItem(const fieldAverageItem&);


    //- Destructor
    ~fieldAverageItem();


    // Member Functions

        // Access

            //- Return const access to the field name
            const word& fieldName() const
            {
                return fieldName_;
            }

            //- Return const access to the mean flag
            const Switch& mean() const
            {
                return mean_;
            }

            //- Return const access to the prime-squared mean flag
            const Switch& prime2Mean() const
            {
                return prime2Mean_;
            }

            //- Return averaging base type name
            const word base() const
            {
                return baseTypeNames_[base_];
            }

            //- Return true if base is ITER
            Switch iterBase() const
            {
                return base_ == ITER;
            }

            //- Return true if base is time
            Switch timeBase() const
            {
                return base_ == TIME;
            }

            scalar window() const
            {
                return window_;
            }

            const word& windowName() const
            {
                return windowName_;
            }


    // Member Operators

        void operator=(const fieldAverageItem&);


    // Friend Operators

        friend bool operator==
        (
            const fieldAverageItem& a,
            const fieldAverageItem& b
        )
        {
            return
                a.fieldName_ == b.fieldName_
             && a.mean_ == b.mean_
             && a.prime2Mean_ == b.prime2Mean_
             && a.base_ == b.base_
             && a.window_ == b.window_
             && a.windowName_ == b.windowName_;
        }

        friend bool operator!=
        (
            const fieldAverageItem& a,
            const fieldAverageItem& b
        )
        {
            return !(a == b);
        }


    // IOstream Operators

        friend Istream& operator>>(Istream&, fieldAverageItem&);
        friend Ostream& operator<<(Ostream&, const fieldAverageItem&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
