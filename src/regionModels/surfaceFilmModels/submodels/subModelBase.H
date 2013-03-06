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
    Foam::subModelBase

Description
    Base class for surface film sub-models

SourceFiles
    subModelBaseI.H
    subModelBase.C

\*---------------------------------------------------------------------------*/

#ifndef subModelBase_H
#define subModelBase_H

#include "surfaceFilmModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                      Class subModelBase Declaration
\*---------------------------------------------------------------------------*/

class subModelBase
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        subModelBase(const subModelBase&);

        //- Disallow default bitwise assignment
        void operator=(const subModelBase&);


protected:

    // Protected data

        //- Reference to the owner surface film model
        const surfaceFilmModel& owner_;

        //- Model coefficients dictionary
        dictionary coeffs_;


public:

    // Constructors

        //- Construct null
        subModelBase(const surfaceFilmModel& owner);

        //- Construct from type name, dictionary and surface film model
        subModelBase
        (
            const word& type,
            const surfaceFilmModel& owner,
            const dictionary& dict
        );


    //- Destructor
    virtual ~subModelBase();


    // Member Functions

        // Access

            //- Return the reference to the owner surface film model
            inline const surfaceFilmModel& owner() const;

            //- Return the model coefficients dictionary
            inline const dictionary& coeffs() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "subModelBaseI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
