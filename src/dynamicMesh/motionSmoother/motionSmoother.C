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

#include "motionSmoother.H"
#include "twoDPointCorrector.H"
#include "faceSet.H"
#include "pointSet.H"
#include "fixedValuePointPatchFields.H"
#include "pointConstraint.H"
#include "syncTools.H"
#include "meshTools.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(motionSmoother, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::motionSmoother::testSyncPositions
(
    const pointField& fld,
    const scalar maxMag
) const
{
    pointField syncedFld(fld);

    syncTools::syncPointPositions
    (
        mesh_,
        syncedFld,
        minEqOp<point>(),           // combine op
        point(GREAT,GREAT,GREAT)    // null
    );

    forAll(syncedFld, i)
    {
        if (mag(syncedFld[i] - fld[i]) > maxMag)
        {
            FatalErrorIn
            (
                "motionSmoother::testSyncPositions"
                "("
                    "const pointField&, "
                    "const scalar"
                ")"
            )   << "On point " << i << " point:" << fld[i]
                << " synchronised point:" << syncedFld[i]
                << abort(FatalError);
        }
    }
}


//Foam::tmp<Foam::scalarField> Foam::motionSmoother::sumWeights
//(
//    const scalarField& edgeWeight
//) const
//{
//    tmp<scalarField> tsumWeight
//    (
//        new scalarField
//        (
//            mesh_.nPoints(),
//            0.0
//        )
//    );
//    scalarField& sumWeight = tsumWeight();
//
//    const edgeList& edges = mesh_.edges();
//
//    forAll(edges, edgeI)
//    {
//        if (isMasterEdge_.get(edgeI) == 1)
//        {
//            const edge& e = edges[edgeI];
//            const scalar w = edgeWeight[edgeI];
//            sumWeight[e[0]] += w;
//            sumWeight[e[1]] += w;
//        }
//    }
//
//
//    // Add coupled contributions
//    // ~~~~~~~~~~~~~~~~~~~~~~~~~
//    syncTools::syncPointList
//    (
//        mesh,
//        sumWeight,
//        plusEqOp<scalar>(),
//        scalar(0)               // null value
//    );
//
//    return tsumWeight;
//}


// From pointPatchInterpolation
void Foam::motionSmoother::makePatchPatchAddressing()
{
    if (debug)
    {
        Pout<< "motionSmoother::makePatchPatchAddressing() : "
            << "constructing boundary addressing"
            << endl;
    }

    const polyBoundaryMesh& bm = mesh_.boundaryMesh();
    const pointBoundaryMesh& pbm = pMesh_.boundary();

    // first count the total number of patch-patch points

    label nPatchPatchPoints = 0;

    forAll(bm, patchi)
    {
        if (!isA<emptyPolyPatch>(bm[patchi]))
        {
            nPatchPatchPoints += bm[patchi].boundaryPoints().size();
        }
    }


    // Go through all patches and mark up the external edge points
    Map<label> patchPatchPointSet(2*nPatchPatchPoints);

    labelList patchPatchPoints(nPatchPatchPoints);

    List<pointConstraint> patchPatchPointConstraints(nPatchPatchPoints);

    label pppi = 0;

    forAll(bm, patchi)
    {
        if (!isA<emptyPolyPatch>(bm[patchi]))
        {
            const labelList& bp = bm[patchi].boundaryPoints();
            const labelList& meshPoints = bm[patchi].meshPoints();

            forAll(bp, pointi)
            {
                label ppp = meshPoints[bp[pointi]];

                Map<label>::iterator iter = patchPatchPointSet.find(ppp);

                if (iter == patchPatchPointSet.end())
                {
                    patchPatchPointSet.insert(ppp, pppi);
                    patchPatchPoints[pppi] = ppp;
                    pbm[patchi].applyConstraint
                    (
                        bp[pointi],
                        patchPatchPointConstraints[pppi]
                    );
                    pppi++;
                }
                else
                {
                    pbm[patchi].applyConstraint
                    (
                        bp[pointi],
                        patchPatchPointConstraints[iter()]
                    );
                }
            }
        }
    }

    nPatchPatchPoints = pppi;
    patchPatchPoints.setSize(nPatchPatchPoints);
    patchPatchPointConstraints.setSize(nPatchPatchPoints);

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

    patchPatchPointConstraintPoints_.setSize(nConstraints);
    patchPatchPointConstraintTensors_.setSize(nConstraints);


    if (debug)
    {
        OFstream str(mesh_.db().path()/"constraintPoints.obj");

        Pout<< "Dumping " << patchPatchPointConstraintPoints_.size()
            << " constraintPoints to " << str.name() << endl;
        forAll(patchPatchPointConstraintPoints_, i)
        {
            label pointI = patchPatchPointConstraintPoints_[i];

            meshTools::writeOBJ(str, mesh_.points()[pointI]);
        }

        Pout<< "motionSmoother::makePatchPatchAddressing() : "
            << "finished constructing boundary addressing"
            << endl;
    }
}


void Foam::motionSmoother::checkFld(const pointScalarField& fld)
{
    forAll(fld, pointI)
    {
        const scalar val = fld[pointI];

        if ((val > -GREAT) && (val < GREAT))
        {}
        else
        {
            FatalErrorIn("motionSmoother::checkFld(const pointScalarField&)")
                << "Problem : point:" << pointI << " value:" << val
                << abort(FatalError);
        }
    }
}


Foam::labelHashSet Foam::motionSmoother::getPoints
(
    const labelHashSet& faceLabels
) const
{
    labelHashSet usedPoints(mesh_.nPoints()/100);

    forAllConstIter(labelHashSet, faceLabels, iter)
    {
        const face& f = mesh_.faces()[iter.key()];

        forAll(f, fp)
        {
            usedPoints.insert(f[fp]);
        }
    }

    return usedPoints;
}


Foam::tmp<Foam::scalarField> Foam::motionSmoother::calcEdgeWeights
(
    const pointField& points
) const
{
    const edgeList& edges = mesh_.edges();

    tmp<scalarField> twght(new scalarField(edges.size()));
    scalarField& wght = twght();

    forAll(edges, edgeI)
    {
        wght[edgeI] = 1.0/(edges[edgeI].mag(points)+SMALL);
    }
    return twght;
}


// Smooth on selected points (usually patch points)
void Foam::motionSmoother::minSmooth
(
    const scalarField& edgeWeights,
    const PackedBoolList& isAffectedPoint,
    const labelList& meshPoints,
    const pointScalarField& fld,
    pointScalarField& newFld
) const
{
    tmp<pointScalarField> tavgFld = avg
    (
        fld,
        edgeWeights //scalarField(mesh_.nEdges(), 1.0)    // uniform weighting
    );
    const pointScalarField& avgFld = tavgFld();

    forAll(meshPoints, i)
    {
        label pointI = meshPoints[i];
        if (isAffectedPoint.get(pointI) == 1)
        {
            newFld[pointI] = min
            (
                fld[pointI],
                0.5*fld[pointI] + 0.5*avgFld[pointI]
            );
        }
    }

    newFld.correctBoundaryConditions();
    applyCornerConstraints(newFld);
}


// Smooth on all internal points
void Foam::motionSmoother::minSmooth
(
    const scalarField& edgeWeights,
    const PackedBoolList& isAffectedPoint,
    const pointScalarField& fld,
    pointScalarField& newFld
) const
{
    tmp<pointScalarField> tavgFld = avg
    (
        fld,
        edgeWeights //scalarField(mesh_.nEdges(), 1.0)    // uniform weighting
    );
    const pointScalarField& avgFld = tavgFld();

    forAll(fld, pointI)
    {
        if (isAffectedPoint.get(pointI) == 1 && isInternalPoint(pointI))
        {
            newFld[pointI] = min
            (
                fld[pointI],
                0.5*fld[pointI] + 0.5*avgFld[pointI]
            );
        }
    }

    newFld.correctBoundaryConditions();
    applyCornerConstraints(newFld);
}


// Scale on all internal points
void Foam::motionSmoother::scaleField
(
    const labelHashSet& pointLabels,
    const scalar scale,
    pointScalarField& fld
) const
{
    forAllConstIter(labelHashSet, pointLabels, iter)
    {
        if (isInternalPoint(iter.key()))
        {
            fld[iter.key()] *= scale;
        }
    }
    fld.correctBoundaryConditions();
    applyCornerConstraints(fld);
}


// Scale on selected points (usually patch points)
void Foam::motionSmoother::scaleField
(
    const labelList& meshPoints,
    const labelHashSet& pointLabels,
    const scalar scale,
    pointScalarField& fld
) const
{
    forAll(meshPoints, i)
    {
        label pointI = meshPoints[i];

        if (pointLabels.found(pointI))
        {
            fld[pointI] *= scale;
        }
    }
}


// Lower on internal points
void Foam::motionSmoother::subtractField
(
    const labelHashSet& pointLabels,
    const scalar f,
    pointScalarField& fld
) const
{
    forAllConstIter(labelHashSet, pointLabels, iter)
    {
        if (isInternalPoint(iter.key()))
        {
            fld[iter.key()] = max(0.0, fld[iter.key()]-f);
        }
    }
    fld.correctBoundaryConditions();
    applyCornerConstraints(fld);
}


// Scale on selected points (usually patch points)
void Foam::motionSmoother::subtractField
(
    const labelList& meshPoints,
    const labelHashSet& pointLabels,
    const scalar f,
    pointScalarField& fld
) const
{
    forAll(meshPoints, i)
    {
        label pointI = meshPoints[i];

        if (pointLabels.found(pointI))
        {
            fld[pointI] = max(0.0, fld[pointI]-f);
        }
    }
}


bool Foam::motionSmoother::isInternalPoint(const label pointI) const
{
    return isInternalPoint_.get(pointI) == 1;
}


void Foam::motionSmoother::getAffectedFacesAndPoints
(
    const label nPointIter,
    const faceSet& wrongFaces,

    labelList& affectedFaces,
    PackedBoolList& isAffectedPoint
) const
{
    isAffectedPoint.setSize(mesh_.nPoints());
    isAffectedPoint = 0;

    faceSet nbrFaces(mesh_, "checkFaces", wrongFaces);

    // Find possible points influenced by nPointIter iterations of
    // scaling and smoothing by doing pointCellpoint walk.
    // Also update faces-to-be-checked to extend one layer beyond the points
    // that will get updated.

    for (label i = 0; i < nPointIter; i++)
    {
        pointSet nbrPoints(mesh_, "grownPoints", getPoints(nbrFaces.toc()));

        forAllConstIter(pointSet, nbrPoints, iter)
        {
            const labelList& pCells = mesh_.pointCells(iter.key());

            forAll(pCells, pCellI)
            {
                const cell& cFaces = mesh_.cells()[pCells[pCellI]];

                forAll(cFaces, cFaceI)
                {
                    nbrFaces.insert(cFaces[cFaceI]);
                }
            }
        }
        nbrFaces.sync(mesh_);

        if (i == nPointIter - 2)
        {
            forAllConstIter(faceSet, nbrFaces, iter)
            {
                const face& f = mesh_.faces()[iter.key()];
                forAll(f, fp)
                {
                    isAffectedPoint.set(f[fp], 1);
                }
            }
        }
    }

    affectedFaces = nbrFaces.toc();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::motionSmoother::motionSmoother
(
    polyMesh& mesh,
    pointMesh& pMesh,
    indirectPrimitivePatch& pp,
    const labelList& adaptPatchIDs,
    const dictionary& paramDict
)
:
    mesh_(mesh),
    pMesh_(pMesh),
    pp_(pp),
    adaptPatchIDs_(adaptPatchIDs),
    paramDict_(paramDict),
    displacement_
    (
        IOobject
        (
            "displacement",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh_
    ),
    scale_
    (
        IOobject
        (
            "scale",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh_,
        dimensionedScalar("scale", dimless, 1.0)
    ),
    oldPoints_(mesh_.points()),
    isInternalPoint_(mesh_.nPoints(), 1),
    twoDCorrector_(mesh_)
{
    updateMesh();
}


Foam::motionSmoother::motionSmoother
(
    polyMesh& mesh,
    indirectPrimitivePatch& pp,
    const labelList& adaptPatchIDs,
    const pointVectorField& displacement,
    const dictionary& paramDict
)
:
    mesh_(mesh),
    pMesh_(const_cast<pointMesh&>(displacement.mesh())),
    pp_(pp),
    adaptPatchIDs_(adaptPatchIDs),
    paramDict_(paramDict),
    displacement_
    (
        IOobject
        (
            "displacement",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        displacement
    ),
    scale_
    (
        IOobject
        (
            "scale",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh_,
        dimensionedScalar("scale", dimless, 1.0)
    ),
    oldPoints_(mesh_.points()),
    isInternalPoint_(mesh_.nPoints(), 1),
    twoDCorrector_(mesh_)
{
    updateMesh();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::motionSmoother::~motionSmoother()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::polyMesh& Foam::motionSmoother::mesh() const
{
    return mesh_;
}


const Foam::pointMesh& Foam::motionSmoother::pMesh() const
{
    return pMesh_;
}


const Foam::indirectPrimitivePatch& Foam::motionSmoother::patch() const
{
    return pp_;
}


const Foam::labelList& Foam::motionSmoother::adaptPatchIDs() const
{
    return adaptPatchIDs_;
}


const Foam::dictionary& Foam::motionSmoother::paramDict() const
{
    return paramDict_;
}


Foam::pointVectorField& Foam::motionSmoother::displacement()
{
    return displacement_;
}


const Foam::pointVectorField& Foam::motionSmoother::displacement() const
{
    return displacement_;
}


const Foam::pointScalarField& Foam::motionSmoother::scale() const
{
    return scale_;
}


const Foam::pointField& Foam::motionSmoother::oldPoints() const
{
    return oldPoints_;
}


void Foam::motionSmoother::correct()
{
    oldPoints_ = mesh_.points();

    scale_ = 1.0;

    // No need to update twoDmotion corrector since only holds edge labels
    // which will remain the same as before. So unless the mesh was distorted
    // severely outside of motionSmoother there will be no need.
}


void Foam::motionSmoother::setDisplacement(pointField& patchDisp)
{
    // See comment in .H file about shared points.
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            const labelList& meshPoints = pp.meshPoints();

            forAll(meshPoints, i)
            {
                displacement_[meshPoints[i]] = vector::zero;
            }
        }
    }

    const labelList& ppMeshPoints = pp_.meshPoints();

    // Set internal point data from displacement on combined patch points.
    forAll(ppMeshPoints, patchPointI)
    {
        displacement_[ppMeshPoints[patchPointI]] = patchDisp[patchPointI];
    }

    // Adapt the fixedValue bc's (i.e. copy internal point data to
    // boundaryField for all affected patches)
    forAll(adaptPatchIDs_, i)
    {
        label patchI = adaptPatchIDs_[i];

        displacement_.boundaryField()[patchI] ==
            displacement_.boundaryField()[patchI].patchInternalField();
    }

    // Make consistent with non-adapted bc's by evaluating those now and
    // resetting the displacement from the values.
    // Note that we're just doing a correctBoundaryConditions with
    // fixedValue bc's first.
    labelHashSet adaptPatchSet(adaptPatchIDs_);

    const lduSchedule& patchSchedule = mesh_.globalData().patchSchedule();

    forAll(patchSchedule, patchEvalI)
    {
        label patchI = patchSchedule[patchEvalI].patch;

        if (!adaptPatchSet.found(patchI))
        {
            if (patchSchedule[patchEvalI].init)
            {
                displacement_.boundaryField()[patchI]
                    .initEvaluate(Pstream::scheduled);
            }
            else
            {
                displacement_.boundaryField()[patchI]
                    .evaluate(Pstream::scheduled);
            }
        }
    }

    // Multi-patch constraints
    applyCornerConstraints(displacement_);

    // Correct for problems introduced by corner constraints
    syncTools::syncPointList
    (
        mesh_,
        displacement_,
        maxMagEqOp(),   // combine op
        vector::zero    // null value
    );

    // Adapt the fixedValue bc's (i.e. copy internal point data to
    // boundaryField for all affected patches) to take the changes caused
    // by multi-corner constraints into account.
    forAll(adaptPatchIDs_, i)
    {
        label patchI = adaptPatchIDs_[i];

        displacement_.boundaryField()[patchI] ==
            displacement_.boundaryField()[patchI].patchInternalField();
    }

    if (debug)
    {
        OFstream str(mesh_.db().path()/"changedPoints.obj");
        label nVerts = 0;
        forAll(ppMeshPoints, patchPointI)
        {
            const vector& newDisp = displacement_[ppMeshPoints[patchPointI]];

            if (mag(newDisp-patchDisp[patchPointI]) > SMALL)
            {
                const point& pt = mesh_.points()[ppMeshPoints[patchPointI]];

                meshTools::writeOBJ(str, pt);
                nVerts++;
                //Pout<< "Point:" << pt
                //    << " oldDisp:" << patchDisp[patchPointI]
                //    << " newDisp:" << newDisp << endl;
            }
        }
        Pout<< "Written " << nVerts << " points that are changed to file "
            << str.name() << endl;
    }

    // Now reset input displacement
    forAll(ppMeshPoints, patchPointI)
    {
        patchDisp[patchPointI] = displacement_[ppMeshPoints[patchPointI]];
    }
}


// correctBoundaryConditions with fixedValue bc's first.
void Foam::motionSmoother::correctBoundaryConditions
(
    pointVectorField& displacement
) const
{
    labelHashSet adaptPatchSet(adaptPatchIDs_);

    const lduSchedule& patchSchedule = mesh_.globalData().patchSchedule();

    // 1. evaluate on adaptPatches
    forAll(patchSchedule, patchEvalI)
    {
        label patchI = patchSchedule[patchEvalI].patch;

        if (adaptPatchSet.found(patchI))
        {
            if (patchSchedule[patchEvalI].init)
            {
                displacement.boundaryField()[patchI]
                    .initEvaluate(Pstream::blocking);
            }
            else
            {
                displacement.boundaryField()[patchI]
                    .evaluate(Pstream::blocking);
            }
        }
    }


    // 2. evaluate on non-AdaptPatches
    forAll(patchSchedule, patchEvalI)
    {
        label patchI = patchSchedule[patchEvalI].patch;

        if (!adaptPatchSet.found(patchI))
        {
            if (patchSchedule[patchEvalI].init)
            {
                displacement.boundaryField()[patchI]
                    .initEvaluate(Pstream::blocking);
            }
            else
            {
                displacement.boundaryField()[patchI]
                    .evaluate(Pstream::blocking);
            }
        }
    }

    // Multi-patch constraints
    applyCornerConstraints(displacement);

    // Correct for problems introduced by corner constraints
    syncTools::syncPointList
    (
        mesh_,
        displacement,
        maxMagEqOp(),           // combine op
        vector::zero            // null value
    );
}


Foam::tmp<Foam::scalarField> Foam::motionSmoother::movePoints
(
    pointField& newPoints
)
{
    // Correct for 2-D motion
    if (twoDCorrector_.required())
    {
        Info<< "Correcting 2-D mesh motion";

        if (mesh_.globalData().parallel())
        {
            WarningIn("motionSmoother::movePoints(pointField&)")
                << "2D mesh-motion probably not correct in parallel" << endl;
        }

        // We do not want to move 3D planes so project all points onto those
        const pointField& oldPoints = mesh_.points();
        const edgeList& edges = mesh_.edges();

        const labelList& neIndices = twoDCorrector().normalEdgeIndices();
        const vector& pn = twoDCorrector().planeNormal();

        forAll(neIndices, i)
        {
            const edge& e = edges[neIndices[i]];

            point& pStart = newPoints[e.start()];

            pStart += pn*(pn & (oldPoints[e.start()] - pStart));

            point& pEnd = newPoints[e.end()];

            pEnd += pn*(pn & (oldPoints[e.end()] - pEnd));
        }

        // Correct tangentially
        twoDCorrector_.correctPoints(newPoints);
        Info<< " ...done" << endl;
    }

    if (debug)
    {
        Pout<< "motionSmoother::movePoints : testing sync of newPoints."
            << endl;
        testSyncPositions(newPoints, 1e-6*mesh_.bounds().mag());
    }

    // Move actual mesh points. Make sure to delete tetBasePtIs so it
    // gets rebuilt.
    mesh_.clearAdditionalGeom();
    tmp<scalarField> tsweptVol = mesh_.movePoints(newPoints);

    pp_.movePoints(mesh_.points());

    return tsweptVol;
}


Foam::scalar Foam::motionSmoother::setErrorReduction
(
    const scalar errorReduction
)
{
    scalar oldErrorReduction = readScalar(paramDict_.lookup("errorReduction"));

    paramDict_.remove("errorReduction");
    paramDict_.add("errorReduction", errorReduction);

    return oldErrorReduction;
}


bool Foam::motionSmoother::scaleMesh
(
    labelList& checkFaces,
    const bool smoothMesh,
    const label nAllowableErrors
)
{
    List<labelPair> emptyBaffles;
    return scaleMesh
    (
        checkFaces,
        emptyBaffles,
        smoothMesh,
        nAllowableErrors
    );
}


bool Foam::motionSmoother::scaleMesh
(
    labelList& checkFaces,
    const List<labelPair>& baffles,
    const bool smoothMesh,
    const label nAllowableErrors
)
{
    return scaleMesh
    (
        checkFaces,
        baffles,
        paramDict_,
        paramDict_,
        smoothMesh,
        nAllowableErrors
    );
}


bool Foam::motionSmoother::scaleMesh
(
    labelList& checkFaces,
    const List<labelPair>& baffles,
    const dictionary& paramDict,
    const dictionary& meshQualityDict,
    const bool smoothMesh,
    const label nAllowableErrors
)
{
    if (!smoothMesh && adaptPatchIDs_.empty())
    {
        FatalErrorIn
        (
            "motionSmoother::scaleMesh"
            "("
                "labelList&, "
                "const List<labelPair>&, "
                "const dictionary&, "
                "const dictionary&, "
                "const bool, "
                "const label"
            ")"
        )
            << "You specified both no movement on the internal mesh points"
            << " (smoothMesh = false)" << nl
            << "and no movement on the patch (adaptPatchIDs is empty)" << nl
            << "Hence nothing to adapt."
            << exit(FatalError);
    }

    if (debug)
    {
        // Had a problem with patches moved non-synced. Check transformations.
        const polyBoundaryMesh& patches = mesh_.boundaryMesh();

        Pout<< "Entering scaleMesh : coupled patches:" << endl;
        forAll(patches, patchI)
        {
            if (patches[patchI].coupled())
            {
                const coupledPolyPatch& pp =
                    refCast<const coupledPolyPatch>(patches[patchI]);

                Pout<< '\t' << patchI << '\t' << pp.name()
                    << " parallel:" << pp.parallel()
                    << " separated:" << pp.separated()
                    << " forwardT:" << pp.forwardT().size()
                    << endl;
            }
        }
    }

    const scalar errorReduction =
        readScalar(paramDict.lookup("errorReduction"));
    const label nSmoothScale =
        readLabel(paramDict.lookup("nSmoothScale"));


    // Note: displacement_ should already be synced already from setDisplacement
    // but just to make sure.
    syncTools::syncPointList
    (
        mesh_,
        displacement_,
        maxMagEqOp(),
        vector::zero    // null value
    );

    // Set newPoints as old + scale*displacement
    pointField newPoints;
    {
        // Create overall displacement with same b.c.s as displacement_
        wordList actualPatchTypes;
        {
            const pointBoundaryMesh& pbm = displacement_.mesh().boundary();
            actualPatchTypes.setSize(pbm.size());
            forAll(pbm, patchI)
            {
                actualPatchTypes[patchI] = pbm[patchI].type();
            }
        }

        pointVectorField totalDisplacement
        (
            IOobject
            (
                "totalDisplacement",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            scale_*displacement_,
            displacement_.boundaryField().types(),
            actualPatchTypes
        );
        correctBoundaryConditions(totalDisplacement);

        if (debug)
        {
            Pout<< "scaleMesh : testing sync of totalDisplacement" << endl;
            testSyncField
            (
                totalDisplacement,
                maxMagEqOp(),
                vector::zero,   // null value
                1e-6*mesh_.bounds().mag()
            );
        }

        newPoints = oldPoints_ + totalDisplacement.internalField();
    }

    Info<< "Moving mesh using displacement scaling :"
        << " min:" << gMin(scale_.internalField())
        << "  max:" << gMax(scale_.internalField())
        << endl;


    // Move
    movePoints(newPoints);

    // Check. Returns parallel number of incorrect faces.
    faceSet wrongFaces(mesh_, "wrongFaces", mesh_.nFaces()/100+100);
    checkMesh(false, mesh_, meshQualityDict, checkFaces, baffles, wrongFaces);

    if (returnReduce(wrongFaces.size(), sumOp<label>()) <= nAllowableErrors)
    {
        return true;
    }
    else
    {
        // Sync across coupled faces by extending the set.
        wrongFaces.sync(mesh_);

        // Special case:
        // if errorReduction is set to zero, extend wrongFaces
        // to face-Cell-faces to ensure quick return to previously valid mesh

        if (mag(errorReduction) < SMALL)
        {
            labelHashSet newWrongFaces(wrongFaces);
            forAllConstIter(labelHashSet, wrongFaces, iter)
            {
                label own = mesh_.faceOwner()[iter.key()];
                const cell& ownFaces = mesh_.cells()[own];

                forAll(ownFaces, cfI)
                {
                    newWrongFaces.insert(ownFaces[cfI]);
                }

                if (iter.key() < mesh_.nInternalFaces())
                {
                    label nei = mesh_.faceNeighbour()[iter.key()];
                    const cell& neiFaces = mesh_.cells()[nei];

                    forAll(neiFaces, cfI)
                    {
                        newWrongFaces.insert(neiFaces[cfI]);
                    }
                }
            }
            wrongFaces.transfer(newWrongFaces);
            wrongFaces.sync(mesh_);
        }


        // Find out points used by wrong faces and scale displacement.
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        pointSet usedPoints(mesh_, "usedPoints", getPoints(wrongFaces));
        usedPoints.sync(mesh_);



        // Grow a few layers to determine
        // - points to be smoothed
        // - faces to be checked in next iteration
        PackedBoolList isAffectedPoint(mesh_.nPoints());
        getAffectedFacesAndPoints
        (
            nSmoothScale,       // smoothing iterations
            wrongFaces,         // error faces
            checkFaces,
            isAffectedPoint
        );

        if (debug)
        {
            Pout<< "Faces in error:" << wrongFaces.size()
                << "  with points:" << usedPoints.size()
                << endl;
        }

        if (adaptPatchIDs_.size())
        {
            // Scale conflicting patch points
            scaleField(pp_.meshPoints(), usedPoints, errorReduction, scale_);
            //subtractField(pp_.meshPoints(), usedPoints, 0.1, scale_);
        }
        if (smoothMesh)
        {
            // Scale conflicting internal points
            scaleField(usedPoints, errorReduction, scale_);
            //subtractField(usedPoints, 0.1, scale_);
        }

        scalarField eWeights(calcEdgeWeights(oldPoints_));

        for (label i = 0; i < nSmoothScale; i++)
        {
            if (adaptPatchIDs_.size())
            {
                // Smooth patch values
                pointScalarField oldScale(scale_);
                minSmooth
                (
                    eWeights,
                    isAffectedPoint,
                    pp_.meshPoints(),
                    oldScale,
                    scale_
                );
                checkFld(scale_);
            }
            if (smoothMesh)
            {
                // Smooth internal values
                pointScalarField oldScale(scale_);
                minSmooth(eWeights, isAffectedPoint, oldScale, scale_);
                checkFld(scale_);
            }
        }

        syncTools::syncPointList
        (
            mesh_,
            scale_,
            maxEqOp<scalar>(),
            -GREAT              // null value
        );


        if (debug)
        {
            Pout<< "scale_ after smoothing :"
                << " min:" << Foam::gMin(scale_)
                << " max:" << Foam::gMax(scale_)
                << endl;
        }

        return false;
    }
}


void Foam::motionSmoother::updateMesh()
{
    const pointBoundaryMesh& patches = pMesh_.boundary();

    // Check whether displacement has fixed value b.c. on adaptPatchID
    forAll(adaptPatchIDs_, i)
    {
        label patchI = adaptPatchIDs_[i];

        if
        (
           !isA<fixedValuePointPatchVectorField>
            (
                displacement_.boundaryField()[patchI]
            )
        )
        {
            FatalErrorIn
            (
                "motionSmoother::updateMesh"
            )   << "Patch " << patches[patchI].name()
                << " has wrong boundary condition "
                << displacement_.boundaryField()[patchI].type()
                << " on field " << displacement_.name() << nl
                << "Only type allowed is "
                << fixedValuePointPatchVectorField::typeName
                << exit(FatalError);
        }
    }


    // Determine internal points. Note that for twoD there are no internal
    // points so we use the points of adaptPatchIDs instead

    twoDCorrector_.updateMesh();

    const labelList& meshPoints = pp_.meshPoints();

    forAll(meshPoints, i)
    {
        isInternalPoint_.unset(meshPoints[i]);
    }

    // Calculate master edge addressing
    isMasterEdge_ = syncTools::getMasterEdges(mesh_);

    makePatchPatchAddressing();
}


// Specialisation of applyCornerConstraints for scalars because
// no constraint need be applied
template<>
void Foam::motionSmoother::applyCornerConstraints<Foam::scalar>
(
    GeometricField<scalar, pointPatchField, pointMesh>& pf
) const
{}


// ************************************************************************* //
