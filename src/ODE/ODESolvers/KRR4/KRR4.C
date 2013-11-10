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

#include "KRR4.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(KRR4, 0);

    addToRunTimeSelectionTable(ODESolver, KRR4, ODE);

const scalar
    KRR4::safety = 0.9, KRR4::grow = 1.5, KRR4::pgrow = -0.25,
    KRR4::shrink = 0.5, KRR4::pshrink = (-1.0/3.0), KRR4::errcon = 0.1296;

const scalar
    KRR4::gamma = 1.0/2.0,
    KRR4::a21 = 2.0, KRR4::a31 = 48.0/25.0, KRR4::a32 = 6.0/25.0,
    KRR4::c21 = -8.0, KRR4::c31 = 372.0/25.0, KRR4::c32 = 12.0/5.0,
    KRR4::c41 = -112.0/125.0, KRR4::c42 = -54.0/125.0, KRR4::c43 = -2.0/5.0,
    KRR4::b1 = 19.0/9.0, KRR4::b2 = 1.0/2.0, KRR4::b3 = 25.0/108.0,
    KRR4::b4 = 125.0/108.0,
    KRR4::e1 = 17.0/54.0, KRR4::e2 = 7.0/36.0, KRR4::e3 = 0.0,
    KRR4::e4 = 125.0/108.0,
    KRR4::c1X = 1.0/2.0, KRR4::c2X = -3.0/2.0, KRR4::c3X = 121.0/50.0,
    KRR4::c4X = 29.0/250.0,
    KRR4::a2X = 1.0, KRR4::a3X = 3.0/5.0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::KRR4::KRR4(const ODE& ode)
:
    ODESolver(ode),
    yTemp_(n_, 0.0),
    dydxTemp_(n_, 0.0),
    g1_(n_, 0.0),
    g2_(n_, 0.0),
    g3_(n_, 0.0),
    g4_(n_, 0.0),
    yErr_(n_, 0.0),
    dfdx_(n_, 0.0),
    dfdy_(n_, n_, 0.0),
    a_(n_, n_, 0.0),
    pivotIndices_(n_, 0.0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::KRR4::solve
(
    const ODE& ode,
    scalar& x,
    scalarField& y,
    scalarField& dydx,
    const scalar eps,
    const scalarField& yScale,
    const scalar hTry,
    scalar& hDid,
    scalar& hNext
) const
{
    scalar xTemp = x;
    yTemp_ = y;
    dydxTemp_ = dydx;

    ode.jacobian(xTemp, yTemp_, dfdx_, dfdy_);

    scalar h = hTry;

    for (register label jtry=0; jtry<maxtry; jtry++)
    {
        for (register label i=0; i<n_; i++)
        {
            for (register label j=0; j<n_; j++)
            {
                a_[i][j] = -dfdy_[i][j];
            }

            a_[i][i] += 1.0/(gamma*h);
        }

        LUDecompose(a_, pivotIndices_);

        for (register label i=0; i<n_; i++)
        {
            g1_[i] = dydxTemp_[i] + h*c1X*dfdx_[i];
        }

        LUBacksubstitute(a_, pivotIndices_, g1_);

        for (register label i=0; i<n_; i++)
        {
            y[i] = yTemp_[i] + a21*g1_[i];
        }

        x = xTemp + a2X*h;
        ode.derivatives(x, y, dydx_);

        for (register label i=0; i<n_; i++)
        {
            g2_[i] = dydx_[i] + h*c2X*dfdx_[i] + c21*g1_[i]/h;
        }

        LUBacksubstitute(a_, pivotIndices_, g2_);

        for (register label i=0; i<n_; i++)
        {
            y[i] = yTemp_[i] + a31*g1_[i] + a32*g2_[i];
        }

        x = xTemp + a3X*h;
        ode.derivatives(x, y, dydx_);

        for (register label i=0; i<n_; i++)
        {
            g3_[i] = dydx_[i] + h*c3X*dfdx_[i] + (c31*g1_[i] + c32*g2_[i])/h;
        }

        LUBacksubstitute(a_, pivotIndices_, g3_);

        for (register label i=0; i<n_; i++)
        {
            g4_[i] = dydx_[i] + h*c4X*dfdx_[i]
                + (c41*g1_[i] + c42*g2_[i] + c43*g3_[i])/h;
        }

        LUBacksubstitute(a_, pivotIndices_, g4_);

        for (register label i=0; i<n_; i++)
        {
            y[i] = yTemp_[i] + b1*g1_[i] + b2*g2_[i] + b3*g3_[i] + b4*g4_[i];
            yErr_[i] = e1*g1_[i] + e2*g2_[i] + e3*g3_[i] + e4*g4_[i];
        }

        x = xTemp + h;

        if (x == xTemp)
        {
            FatalErrorIn
            (
                "void Foam::KRR4::solve"
                "("
                    "const ODE&, "
                    "scalar&, "
                    "scalarField&, "
                    "scalarField&, "
                    "const scalar, "
                    "const scalarField&, "
                    "const scalar, "
                    "scalar&, "
                    "scalar&"
                ") const"
            )   << "solver stalled: step size = 0"
                << exit(FatalError);
        }

        scalar maxErr = 0.0;
        for (register label i=0; i<n_; i++)
        {
            maxErr = max(maxErr, mag(yErr_[i]/yScale[i]));
        }
        maxErr /= eps;

        if (maxErr <= 1.0)
        {
            hDid = h;
            hNext = (maxErr > errcon ? safety*h*pow(maxErr, pgrow) : grow*h);
            return;
        }
        else
        {
            hNext = safety*h*pow(maxErr, pshrink);
            h = (h >= 0.0 ? max(hNext, shrink*h) : min(hNext, shrink*h));
        }
    }

    FatalErrorIn
    (
        "void Foam::KRR4::solve"
        "("
            "const ODE&, "
            "scalar&, "
            "scalarField&, "
            "scalarField&, "
            "const scalar, "
            "const scalarField&, "
            "const scalar, "
            "scalar&, "
            "scalar&"
        ") const"
    )   << "Maximum number of solver iterations exceeded"
        << exit(FatalError);
}


// ************************************************************************* //
