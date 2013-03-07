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

\*---------------------------------------------------------------------------*/

#include "volPointInterpolation.H"
#include "fvMesh.H"
#include "volFields.H"
#include "pointFields.H"
#include "demandDrivenData.H"
#include "coupledPointPatchFields.H"
#include "pointConstraint.H"

#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(volPointInterpolation, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void volPointInterpolation::calcBoundaryAddressing()
{
    if (debug)
    {
        Pout<< "volPointInterpolation::calcBoundaryAddressing() : "
            << "constructing boundary addressing"
            << endl;
    }

    boundaryPtr_.reset
    (
        new primitivePatch
        (
            SubList<face>
            (
                mesh().faces(),
                mesh().nFaces()-mesh().nInternalFaces(),
                mesh().nInternalFaces()
            ),
            mesh().points()
        )
    );
    const primitivePatch& boundary = boundaryPtr_();

    boundaryIsPatchFace_.setSize(boundary.size());
    boundaryIsPatchFace_ = false;

    isPatchPoint_.setSize(mesh().nPoints());
    isPatchPoint_ = false;

    const polyBoundaryMesh& pbm = mesh().boundaryMesh();

    // Get precalculated volField only so we can use coupled() tests for
    // cyclicAMI
    const surfaceScalarField& magSf = mesh().magSf();

    forAll(pbm, patchI)
    {
        const polyPatch& pp = pbm[patchI];

        if
        (
            !isA<emptyPolyPatch>(pp)
         && !magSf.boundaryField()[patchI].coupled()
        )
        {
            label bFaceI = pp.start()-mesh().nInternalFaces();

            forAll(pp, i)
            {
                boundaryIsPatchFace_[bFaceI] = true;

                const face& f = boundary[bFaceI++];

                forAll(f, fp)
                {
                    isPatchPoint_[f[fp]] = true;
                }
            }
        }
    }

    // Make sure point status is synchronised so even processor that holds
    // no face of a certain patch still can have boundary points marked.
    if (debug)
    {
        boolList oldData(isPatchPoint_);

        syncUntransformedData(isPatchPoint_, orEqOp<bool>());

        forAll(isPatchPoint_, pointI)
        {
            if (isPatchPoint_[pointI] != oldData[pointI])
            {
                Pout<< "volPointInterpolation::calcBoundaryAddressing():"
                    << " added dangling mesh point:" << pointI
                    << " at:" << mesh().points()[pointI]
                    << endl;
            }
        }

        label nPatchFace = 0;
        forAll(boundaryIsPatchFace_, i)
        {
            if (boundaryIsPatchFace_[i])
            {
                nPatchFace++;
            }
        }
        label nPatchPoint = 0;
        forAll(isPatchPoint_, i)
        {
            if (isPatchPoint_[i])
            {
                nPatchPoint++;
            }
        }
        Pout<< "boundary:" << nl
            << "    faces :" << boundary.size() << nl
            << "    of which on proper patch:" << nPatchFace << nl
            << "    points:" << boundary.nPoints() << nl
            << "    of which on proper patch:" << nPatchPoint << endl;
    }
}


void volPointInterpolation::makeInternalWeights(scalarField& sumWeights)
{
    if (debug)
    {
        Pout<< "volPointInterpolation::makeInternalWeights() : "
            << "constructing weighting factors for internal and non-coupled"
            << " points." << endl;
    }

    const pointField& points = mesh().points();
    const labelListList& pointCells = mesh().pointCells();
    const vectorField& cellCentres = mesh().cellCentres();

    // Allocate storage for weighting factors
    pointWeights_.clear();
    pointWeights_.setSize(points.size());

    // Calculate inverse distances between cell centres and points
    // and store in weighting factor array
    forAll(points, pointi)
    {
        if (!isPatchPoint_[pointi])
        {
            const labelList& pcp = pointCells[pointi];

            scalarList& pw = pointWeights_[pointi];
            pw.setSize(pcp.size());

            forAll(pcp, pointCelli)
            {
                pw[pointCelli] =
                    1.0/mag(points[pointi] - cellCentres[pcp[pointCelli]]);

                sumWeights[pointi] += pw[pointCelli];
            }
        }
    }
}


void volPointInterpolation::makeBoundaryWeights(scalarField& sumWeights)
{
    if (debug)
    {
        Pout<< "volPointInterpolation::makeBoundaryWeights() : "
            << "constructing weighting factors for boundary points." << endl;
    }

    const pointField& points = mesh().points();
    const pointField& faceCentres = mesh().faceCentres();

    const primitivePatch& boundary = boundaryPtr_();

    boundaryPointWeights_.clear();
    boundaryPointWeights_.setSize(boundary.meshPoints().size());

    forAll(boundary.meshPoints(), i)
    {
        label pointI = boundary.meshPoints()[i];

        if (isPatchPoint_[pointI])
        {
            const labelList& pFaces = boundary.pointFaces()[i];

            scalarList& pw = boundaryPointWeights_[i];
            pw.setSize(pFaces.size());

            sumWeights[pointI] = 0.0;

            forAll(pFaces, i)
            {
                if (boundaryIsPatchFace_[pFaces[i]])
                {
                    label faceI = mesh().nInternalFaces() + pFaces[i];

                    pw[i] = 1.0/mag(points[pointI] - faceCentres[faceI]);
                    sumWeights[pointI] += pw[i];
                }
                else
                {
                    pw[i] = 0.0;
                }
            }
        }
    }
}


void volPointInterpolation::makeWeights()
{
    if (debug)
    {
        Pout<< "volPointInterpolation::makeWeights() : "
            << "constructing weighting factors"
            << endl;
    }

    // Update addressing over all boundary faces
    calcBoundaryAddressing();


    // Running sum of weights
    pointScalarField sumWeights
    (
        IOobject
        (
            "volPointSumWeights",
            mesh().polyMesh::instance(),
            mesh()
        ),
        pointMesh::New(mesh()),
        dimensionedScalar("zero", dimless, 0)
    );


    // Create internal weights; add to sumWeights
    makeInternalWeights(sumWeights);


    // Create boundary weights; override sumWeights
    makeBoundaryWeights(sumWeights);


    //forAll(boundary.meshPoints(), i)
    //{
    //    label pointI = boundary.meshPoints()[i];
    //
    //    if (isPatchPoint_[pointI])
    //    {
    //        Pout<< "Calculated Weight at boundary point:" << i
    //            << " at:" << mesh().points()[pointI]
    //            << " sumWeight:" << sumWeights[pointI]
    //            << " from:" << boundaryPointWeights_[i]
    //            << endl;
    //    }
    //}


    // Sum collocated contributions
    syncUntransformedData(sumWeights, plusEqOp<scalar>());

    // And add separated contributions
    addSeparated(sumWeights);

    // Push master data to slaves. It is possible (not sure how often) for
    // a coupled point to have its master on a different patch so
    // to make sure just push master data to slaves. Reuse the syncPointData
    // structure.
    pushUntransformedData(sumWeights);


    // Normalise internal weights
    forAll(pointWeights_, pointI)
    {
        scalarList& pw = pointWeights_[pointI];
        // Note:pw only sized for !isPatchPoint
        forAll(pw, i)
        {
            pw[i] /= sumWeights[pointI];
        }
    }

    // Normalise boundary weights
    const primitivePatch& boundary = boundaryPtr_();

    forAll(boundary.meshPoints(), i)
    {
        label pointI = boundary.meshPoints()[i];

        scalarList& pw = boundaryPointWeights_[i];
        // Note:pw only sized for isPatchPoint
        forAll(pw, i)
        {
            pw[i] /= sumWeights[pointI];
        }
    }


    if (debug)
    {
        Pout<< "volPointInterpolation::makeWeights() : "
            << "finished constructing weighting factors"
            << endl;
    }
}


void volPointInterpolation::makePatchPatchAddressing()
{
    if (debug)
    {
        Pout<< "volPointInterpolation::makePatchPatchAddressing() : "
            << "constructing boundary addressing"
            << endl;
    }

    const fvBoundaryMesh& bm = mesh().boundary();
    const pointBoundaryMesh& pbm = pointMesh::New(mesh()).boundary();


    // first count the total number of patch-patch points

    label nPatchPatchPoints = 0;

    forAll(bm, patchi)
    {
        if (!isA<emptyFvPatch>(bm[patchi]) && !bm[patchi].coupled())
        {
            nPatchPatchPoints += bm[patchi].patch().boundaryPoints().size();

            if (debug)
            {
                Pout<< "On patch:" << bm[patchi].patch()
                    << " nBoundaryPoints:"
                    << bm[patchi].patch().boundaryPoints().size() << endl;
            }
        }
    }

    if (debug)
    {
        Pout<< "Found nPatchPatchPoints:" << nPatchPatchPoints << endl;
    }


    // Go through all patches and mark up the external edge points
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // From meshpoint to index in patchPatchPointConstraints.
    Map<label> patchPatchPointSet(2*nPatchPatchPoints);

    // Constraints (initialised to unconstrained)
    List<pointConstraint> patchPatchPointConstraints(nPatchPatchPoints);

    // From constraint index to mesh point
    labelList patchPatchPoints(nPatchPatchPoints);

    label pppi = 0;

    forAll(bm, patchi)
    {
        if (!isA<emptyFvPatch>(bm[patchi]) && !bm[patchi].coupled())
        {
            const labelList& bp = bm[patchi].patch().boundaryPoints();
            const labelList& meshPoints = bm[patchi].patch().meshPoints();

            forAll(bp, pointi)
            {
                label ppp = meshPoints[bp[pointi]];

                Map<label>::iterator iter = patchPatchPointSet.find(ppp);

                label constraintI = -1;

                if (iter == patchPatchPointSet.end())
                {
                    patchPatchPointSet.insert(ppp, pppi);
                    patchPatchPoints[pppi] = ppp;
                    constraintI = pppi++;
                }
                else
                {
                    constraintI = iter();
                }

                // Apply to patch constraints
                pbm[patchi].applyConstraint
                (
                    bp[pointi],
                    patchPatchPointConstraints[constraintI]
                );
            }
        }
    }

    if (debug)
    {
        Pout<< "Have (local) constrained points:"
            << nPatchPatchPoints << endl;
    }


    // Extend set with constraints across coupled points
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        const globalMeshData& gd = mesh().globalData();
        const labelListList& globalPointSlaves = gd.globalPointSlaves();
        const mapDistribute& globalPointSlavesMap = gd.globalPointSlavesMap();
        const Map<label>& cpPointMap = gd.coupledPatch().meshPointMap();
        const labelList& cpMeshPoints = gd.coupledPatch().meshPoints();

        // Constraints on coupled points
        List<pointConstraint> constraints
        (
            globalPointSlavesMap.constructSize()
        );

        // Copy from patchPatch constraints into coupledConstraints.
        forAll(bm, patchi)
        {
            if (!isA<emptyFvPatch>(bm[patchi]) && !bm[patchi].coupled())
            {
                const labelList& bp = bm[patchi].patch().boundaryPoints();
                const labelList& meshPoints = bm[patchi].patch().meshPoints();

                forAll(bp, pointi)
                {
                    label ppp = meshPoints[bp[pointi]];

                    Map<label>::const_iterator fnd = cpPointMap.find(ppp);
                    if (fnd != cpPointMap.end())
                    {
                        // Can just copy (instead of apply) constraint
                        // will already be consistent across multiple patches.
                        constraints[fnd()] = patchPatchPointConstraints
                        [
                            patchPatchPointSet[ppp]
                        ];
                    }
                }
            }
        }

        // Exchange data
        globalPointSlavesMap.distribute(constraints);

        // Combine master with slave constraints
        forAll(globalPointSlaves, pointI)
        {
            const labelList& slaves = globalPointSlaves[pointI];

            // Combine master constraint with slave constraints
            forAll(slaves, i)
            {
                constraints[pointI].combine(constraints[slaves[i]]);
            }
            // Duplicate master constraint into slave slots
            forAll(slaves, i)
            {
                constraints[slaves[i]] = constraints[pointI];
            }
        }

        // Send back
        globalPointSlavesMap.reverseDistribute
        (
            cpMeshPoints.size(),
            constraints
        );

        // Add back into patchPatch constraints
        forAll(constraints, coupledPointI)
        {
            if (constraints[coupledPointI].first() != 0)
            {
                label meshPointI = cpMeshPoints[coupledPointI];

                Map<label>::iterator iter = patchPatchPointSet.find(meshPointI);

                label constraintI = -1;

                if (iter == patchPatchPointSet.end())
                {
                    //Pout<< "on meshpoint:" << meshPointI
                    //    << " coupled:" << coupledPointI
                    //    << " at:" << mesh().points()[meshPointI]
                    //    << " have new constraint:"
                    //    << constraints[coupledPointI]
                    //    << endl;

                    // Allocate new constraint
                    if (patchPatchPoints.size() <= pppi)
                    {
                        patchPatchPoints.setSize(pppi+100);
                    }
                    patchPatchPointSet.insert(meshPointI, pppi);
                    patchPatchPoints[pppi] = meshPointI;
                    constraintI = pppi++;
                }
                else
                {
                    //Pout<< "on meshpoint:" << meshPointI
                    //    << " coupled:" << coupledPointI
                    //    << " at:" << mesh().points()[meshPointI]
                    //    << " have possibly extended constraint:"
                    //    << constraints[coupledPointI]
                    //    << endl;

                    constraintI = iter();
                }

                // Combine (new or existing) constraint with one
                // on coupled.
                patchPatchPointConstraints[constraintI].combine
                (
                    constraints[coupledPointI]
                );
            }
        }
    }



    nPatchPatchPoints = pppi;
    patchPatchPoints.setSize(nPatchPatchPoints);
    patchPatchPointConstraints.setSize(nPatchPatchPoints);


    if (debug)
    {
        Pout<< "Have (global) constrained points:"
            << nPatchPatchPoints << endl;
    }


    // Copy out all non-trivial constraints
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    patchPatchPointConstraintPoints_.setSize(nPatchPatchPoints);
    patchPatchPointConstraintTensors_.setSize(nPatchPatchPoints);

    label nConstraints = 0;

    forAll(patchPatchPointConstraints, i)
    {
        if (patchPatchPointConstraints[i].first() != 0)
        {
            patchPatchPointConstraintPoints_[nConstraints] =
                patchPatchPoints[i];

            patchPatchPointConstraintTensors_[nConstraints] =
                patchPatchPointConstraints[i].constraintTransformation();

            nConstraints++;
        }
    }

    if (debug)
    {
        Pout<< "Have non-trivial constrained points:"
            << nConstraints << endl;
    }

    patchPatchPointConstraintPoints_.setSize(nConstraints);
    patchPatchPointConstraintTensors_.setSize(nConstraints);


    if (debug)
    {
        Pout<< "volPointInterpolation::makePatchPatchAddressing() : "
            << "finished constructing boundary addressing"
            << endl;
    }
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

volPointInterpolation::volPointInterpolation(const fvMesh& vm)
:
    MeshObject<fvMesh, Foam::UpdateableMeshObject, volPointInterpolation>(vm)
{
    makeWeights();
    makePatchPatchAddressing();
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

volPointInterpolation::~volPointInterpolation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void volPointInterpolation::updateMesh(const mapPolyMesh&)
{
    makeWeights();
    makePatchPatchAddressing();
}


bool volPointInterpolation::movePoints()
{
    makeWeights();

    return true;
}


// Specialisaion of applyCornerConstraints for scalars because
// no constraint need be applied
template<>
void volPointInterpolation::applyCornerConstraints<scalar>
(
    GeometricField<scalar, pointPatchField, pointMesh>& pf
) const
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
