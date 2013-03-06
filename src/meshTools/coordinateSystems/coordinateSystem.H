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
    Foam::coordinateSystem

Description
    A cartesian coordinate system and the base class for other coordinate
    system specifications.

    All systems are defined by an origin point and a coordinateRotation.
    For convenience, the dictionary constructor forms allow a few shortcuts:
    - the default origin corresponds to <em>(0 0 0)</em>
    - if the \c type is not otherwise specified, a Cartesian coordinateSystem
      is implicit

    \verbatim
        flipped
        {
            origin  (0 0 0);
            coordinateRotation
            {
                type        STARCDRotation;
                rotation    (0 0 90);
            }
        }
    \endverbatim

    - if an axes specification (eg, e3/e1) is used, the coordinateRotation
      sub-dictionary can be dropped.

    \verbatim
        flipped     // the same, specified as axes
        {
            origin  (0 0 0);
            coordinateRotation
            {
                type    axes;
                e3      (1 0 0);
                e1      (0 0 -1);
            }
        }
        flipped     // the same, using all the shortcuts
        {
            e3      (1 0 0);
            e1      (0 0 -1);
        }
    \endverbatim

    - if a sub-dictionary coordinateSystem is found within the dictionary, it
      will be used. This provides a convenient means of embedding
      coordinateSystem information in another dictionary.
      This is used, for example, in the porousZones:

    \verbatim
        1
        (
        cat1
        {
            coordinateSystem
            {
                origin  (0 0 0);
                coordinateRotation
                {
                    type        STARCDRotation;
                    rotation    (0 0 90);
                }
            }
            porosity        0.781;
            Darcy
            {
                d   d [0 -2 0 0 0]  (-1000 -1000 0.50753e+08);
                f   f [0 -1 0 0 0]  (-1000 -1000 12.83);
            }
        }
        )
    \endverbatim

    - additionally, if the coordinateSystem points to a plain entry,
      it can be used to reference one of the global coordinateSystems

    \verbatim
        1
        (
        cat1
        {
            coordinateSystem  system_10;
            porosity        0.781;
            Darcy
            {
                d   d [0 -2 0 0 0]  (-1000 -1000 0.50753e+08);
                f   f [0 -1 0 0 0]  (-1000 -1000 12.83);
            }
        }
        )
    \endverbatim
    For this to work correctly, the coordinateSystem constructor must be
    supplied with both a dictionary and an objectRegistry.

See Also
    coordinateSystems and coordinateSystems::New

SourceFiles
    coordinateSystem.C
    newCoordinateSystem.C
\*---------------------------------------------------------------------------*/

#ifndef coordinateSystem_H
#define coordinateSystem_H

#include "vector.H"
#include "point.H"
#include "tensor.H"
#include "vectorField.H"
#include "pointField.H"
#include "tmp.H"
#include "coordinateRotation.H"
#include "objectRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class coordinateSystem Declaration
\*---------------------------------------------------------------------------*/

class coordinateSystem
{
    // Private data

        //- Name of coordinate system
        mutable word name_;

        //- Optional note
        mutable string note_;

        //- Origin
        mutable point origin_;

        //- Local-to-Global transformation tensor
        coordinateRotation R_;

        //- Global-to-Local transformation tensor
        tensor Rtr_;


protected:

    // Protected Member Functions

        //- Convert from local coordinate system to the global Cartesian system
        //  with optional translation for the origin
        virtual vector localToGlobal(const vector&, bool translate) const;

        //- Convert from local coordinate system to the global Cartesian system
        //  with optional translation for the origin
        virtual tmp<vectorField> localToGlobal
        (
            const vectorField&,
            bool translate
        ) const;

        //- Convert from global Cartesian system to the local coordinate system
        //  with optional translation for the origin
        virtual vector globalToLocal(const vector&, bool translate) const;

        //- Convert from global Cartesian system to the local coordinate system
        //  with optional translation for the origin
        virtual tmp<vectorField> globalToLocal
        (
            const vectorField&,
            bool translate
        ) const;


public:

    //- Runtime type information
    TypeName("coordinateSystem");


    // Constructors

        //- Construct null. This is equivalent to an identity coordinateSystem
        coordinateSystem();

        //- Construct copy with a different name
        coordinateSystem
        (
            const word& name,
            const coordinateSystem&
        );

        //- Construct from origin and rotation
        coordinateSystem
        (
            const word& name,
            const point& origin,
            const coordinateRotation&
        );

        //- Construct from origin and 2 axes
        coordinateSystem
        (
            const word& name,
            const point& origin,
            const vector& axis,
            const vector& dirn
        );

        //- Construct from dictionary with a given name
        coordinateSystem(const word& name, const dictionary&);

        //- Construct from dictionary with default name
        coordinateSystem(const dictionary&);

        //- Construct from dictionary (default name)
        //  With the ability to reference global coordinateSystems
        coordinateSystem(const dictionary&, const objectRegistry&);

        //- Construct from Istream
        //  The Istream contains a word followed by a dictionary
        coordinateSystem(Istream&);

        //- Return clone
        autoPtr<coordinateSystem> clone() const
        {
            return autoPtr<coordinateSystem>(new coordinateSystem(*this));
        }


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            coordinateSystem,
            dictionary,
            (
                const word& name,
                const dictionary& dict
            ),
            (name, dict)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            coordinateSystem,
            origRotation,
            (
                const word& name,
                const point& origin,
                const coordinateRotation& cr
            ),
            (name, origin, cr)
        );


    // Selectors

        //- Select constructed from dictionary
        static autoPtr<coordinateSystem> New
        (
            const word& name,
            const dictionary&
        );

        //- Select constructed from origin and rotation
        static autoPtr<coordinateSystem> New
        (
            const word& coordType,
            const word& name,
            const point& origin,
            const coordinateRotation&
        );

        //- Select constructed from Istream
        static autoPtr<coordinateSystem> New(Istream& is);


    //- Destructor
    virtual ~coordinateSystem();


    // Member Functions

        // Access

            //- Return name
            const word& name() const
            {
                return name_;
            }

            //- Return non-constant access to the optional note
            string& note()
            {
                return note_;
            }

            //- Return the optional note
            const string& note() const
            {
                return note_;
            }

            //- Return origin
            const point& origin() const
            {
                return origin_;
            }

            //- Return coordinate rotation
            const coordinateRotation& rotation() const
            {
                return R_;
            }

            //- Return local-to-global transformation tensor
            const tensor& R() const
            {
                return R_;
            }

            //- Return local Cartesian x-axis
            const vector e1() const
            {
                return Rtr_.x();
            }

            //- Return local Cartesian y-axis
            const vector e2() const
            {
                return Rtr_.y();
            }

            //- Return local Cartesian z-axis
            const vector e3() const
            {
                return Rtr_.z();
            }

            //- Return axis (e3: local Cartesian z-axis)
            // \deprecated method e3 is preferred (deprecated Apr 2008)
            const vector axis() const
            {
                return Rtr_.z();
            }

            //- Return direction (e1: local Cartesian x-axis)
            // \deprecated method e1 is preferred (deprecated Apr 2008)
            const vector direction() const
            {
                return Rtr_.x();
            }

            //- Return as dictionary of entries
            //  \param [in] ignoreType drop type (cartesian, cylindrical, etc)
            //  when generating the dictionary
            virtual dictionary dict(bool ignoreType=false) const;


        // Edit

            //- Rename
            virtual void rename(const word& newName)
            {
                name_ = newName;
            }

            //- Edit access to origin
            point& origin()
            {
                return origin_;
            }

            //- Reset origin and rotation to an identity coordinateSystem
            //  Also resets the note
            virtual void clear();


        // Write

            //- Write
            virtual void write(Ostream&) const;

            //- Write dictionary
            virtual void writeDict(Ostream&, bool subDict=true) const;


        // Transformations

            //- Convert from position in local coordinate system to global
            //  Cartesian position
            point globalPosition(const point& local) const
            {
                return localToGlobal(local, true);
            }

            //- Convert from position in local coordinate system to global
            //  Cartesian position
            tmp<pointField> globalPosition(const pointField& local) const
            {
                return localToGlobal(local, true);
            }

            //- Convert from vector components in local coordinate system to
            //  global Cartesian vector
            vector globalVector(const vector& local) const
            {
                return localToGlobal(local, false);
            }

            //- Convert from vector components in local coordinate system to
            //  global Cartesian vector
            tmp<vectorField> globalVector(const vectorField& local) const
            {
                return localToGlobal(local, false);
            }

            //- Convert from global Cartesian position to position in local
            //  coordinate system
            point localPosition(const point& global) const
            {
                return globalToLocal(global, true);
            }

            //- Convert from global Cartesian position to position in local
            //  coordinate system
            tmp<pointField> localPosition(const pointField& global) const
            {
                return globalToLocal(global, true);
            }

            //- Convert from global Cartesian vector to components in local
            //  coordinate system
            vector localVector(const vector& global) const
            {
                return globalToLocal(global, false);
            }

            //- Convert from global Cartesian vector to components in local
            //  coordinate system
            tmp<vectorField> localVector(const vectorField& global) const
            {
                return globalToLocal(global, false);
            }


    // Member Operators

        //- assign from dictionary
        void operator=(const dictionary&);


    // friend Operators

        friend bool operator!=
        (
            const coordinateSystem&,
            const coordinateSystem&
        );


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const coordinateSystem&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
