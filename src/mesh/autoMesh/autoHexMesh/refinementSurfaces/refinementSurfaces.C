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

#include "refinementSurfaces.H"
#include "Time.H"
#include "searchableSurfaces.H"
#include "shellSurfaces.H"
#include "triSurfaceMesh.H"
#include "labelPair.H"
#include "searchableSurfacesQueries.H"
#include "UPtrList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::refinementSurfaces::areaSelectionAlgo,
        4
    >::names[] =
    {
        "inside",
        "outside",
        "insidePoint",
        "none"
    };
}
const Foam::NamedEnum<Foam::refinementSurfaces::areaSelectionAlgo, 4>
    Foam::refinementSurfaces::areaSelectionAlgoNames;


namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::refinementSurfaces::faceZoneType,
        3
    >::names[] =
    {
        "internal",
        "baffle",
        "boundary"
    };
}
const Foam::NamedEnum<Foam::refinementSurfaces::faceZoneType, 3>
    Foam::refinementSurfaces::faceZoneTypeNames;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinementSurfaces::refinementSurfaces
(
    const searchableSurfaces& allGeometry,
    const dictionary& surfacesDict
)
:
    allGeometry_(allGeometry),
    surfaces_(surfacesDict.size()),
    names_(surfacesDict.size()),
    faceZoneNames_(surfacesDict.size()),
    cellZoneNames_(surfacesDict.size()),
    zoneInside_(surfacesDict.size(), NONE),
    zoneInsidePoints_(surfacesDict.size()),
    faceType_(surfacesDict.size(), INTERNAL),
    regionOffset_(surfacesDict.size())
{
    // Wilcard specification : loop over all surface, all regions
    // and try to find a match.

    // Count number of surfaces.
    label surfI = 0;
    forAll(allGeometry.names(), geomI)
    {
        const word& geomName = allGeometry_.names()[geomI];

        if (surfacesDict.found(geomName))
        {
            surfI++;
        }
    }

    // Size lists
    surfaces_.setSize(surfI);
    names_.setSize(surfI);
    faceZoneNames_.setSize(surfI);
    cellZoneNames_.setSize(surfI);
    zoneInside_.setSize(surfI, NONE);
    faceType_.setSize(surfI, INTERNAL),
    regionOffset_.setSize(surfI);

    labelList globalMinLevel(surfI, 0);
    labelList globalMaxLevel(surfI, 0);
    scalarField globalAngle(surfI, -GREAT);
    PtrList<dictionary> globalPatchInfo(surfI);
    List<Map<label> > regionMinLevel(surfI);
    List<Map<label> > regionMaxLevel(surfI);
    List<Map<scalar> > regionAngle(surfI);
    List<Map<autoPtr<dictionary> > > regionPatchInfo(surfI);


    HashSet<word> unmatchedKeys(surfacesDict.toc());

    surfI = 0;
    forAll(allGeometry.names(), geomI)
    {
        const word& geomName = allGeometry_.names()[geomI];

        const entry* ePtr = surfacesDict.lookupEntryPtr(geomName, false, true);

        if (ePtr)
        {
            const dictionary& dict = ePtr->dict();
            unmatchedKeys.erase(ePtr->keyword());

            names_[surfI] = geomName;
            surfaces_[surfI] = geomI;

            const labelPair refLevel(dict.lookup("level"));
            globalMinLevel[surfI] = refLevel[0];
            globalMaxLevel[surfI] = refLevel[1];

            // Global zone names per surface
            if (dict.readIfPresent("faceZone", faceZoneNames_[surfI]))
            {
                // Read optional entry to determine inside of faceZone

                word method;
                bool hasSide = dict.readIfPresent("cellZoneInside", method);
                if (hasSide)
                {
                    zoneInside_[surfI] = areaSelectionAlgoNames[method];
                    if (zoneInside_[surfI] == INSIDEPOINT)
                    {
                        dict.lookup("insidePoint") >> zoneInsidePoints_[surfI];
                    }

                }
                else
                {
                    // Check old syntax
                    bool inside;
                    if (dict.readIfPresent("zoneInside", inside))
                    {
                        hasSide = true;
                        zoneInside_[surfI] = (inside ? INSIDE : OUTSIDE);
                    }
                }

                // Read optional cellZone name

                if (dict.readIfPresent("cellZone", cellZoneNames_[surfI]))
                {
                    if
                    (
                        (
                            zoneInside_[surfI] == INSIDE
                         || zoneInside_[surfI] == OUTSIDE
                        )
                    && !allGeometry_[surfaces_[surfI]].hasVolumeType()
                    )
                    {
                        IOWarningIn
                        (
                            "refinementSurfaces::refinementSurfaces(..)",
                            dict
                        )   << "Illegal entry zoneInside "
                            << areaSelectionAlgoNames[zoneInside_[surfI]]
                            << " for faceZone "
                            << faceZoneNames_[surfI]
                            << " since surface " << names_[surfI]
                            << " is not closed." << endl;
                    }
                }
                else if (hasSide)
                {
                    IOWarningIn
                    (
                        "refinementSurfaces::refinementSurfaces(..)",
                        dict
                    )   << "Unused entry zoneInside for faceZone "
                        << faceZoneNames_[surfI]
                        << " since no cellZone specified."
                        << endl;
                }

                // How to handle faces on faceZone
                word faceTypeMethod;
                if (dict.readIfPresent("faceType", faceTypeMethod))
                {
                    faceType_[surfI] = faceZoneTypeNames[faceTypeMethod];
                }
            }



            // Global perpendicular angle
            if (dict.found("patchInfo"))
            {
                globalPatchInfo.set
                (
                    surfI,
                    dict.subDict("patchInfo").clone()
                );
            }
            dict.readIfPresent("perpendicularAngle", globalAngle[surfI]);

            if (dict.found("regions"))
            {
                const dictionary& regionsDict = dict.subDict("regions");
                const wordList& regionNames =
                    allGeometry_[surfaces_[surfI]].regions();

                forAll(regionNames, regionI)
                {
                    if (regionsDict.found(regionNames[regionI]))
                    {
                        // Get the dictionary for region
                        const dictionary& regionDict = regionsDict.subDict
                        (
                            regionNames[regionI]
                        );

                        const labelPair refLevel(regionDict.lookup("level"));

                        regionMinLevel[surfI].insert(regionI, refLevel[0]);
                        regionMaxLevel[surfI].insert(regionI, refLevel[1]);

                        if (regionDict.found("perpendicularAngle"))
                        {
                            regionAngle[surfI].insert
                            (
                                regionI,
                                readScalar
                                (
                                    regionDict.lookup("perpendicularAngle")
                                )
                            );
                        }

                        if (regionDict.found("patchInfo"))
                        {
                            regionPatchInfo[surfI].insert
                            (
                                regionI,
                                regionDict.subDict("patchInfo").clone()
                            );
                        }
                    }
                }
            }
            surfI++;
        }
    }

    if (unmatchedKeys.size() > 0)
    {
        IOWarningIn
        (
            "refinementSurfaces::refinementSurfaces(..)",
            surfacesDict
        )   << "Not all entries in refinementSurfaces dictionary were used."
            << " The following entries were not used : "
            << unmatchedKeys.sortedToc()
            << endl;
    }


    // Calculate local to global region offset
    label nRegions = 0;

    forAll(surfaces_, surfI)
    {
        regionOffset_[surfI] = nRegions;
        nRegions += allGeometry_[surfaces_[surfI]].regions().size();
    }

    // Rework surface specific information into information per global region
    minLevel_.setSize(nRegions);
    minLevel_ = 0;
    maxLevel_.setSize(nRegions);
    maxLevel_ = 0;
    perpendicularAngle_.setSize(nRegions);
    perpendicularAngle_ = -GREAT;
    patchInfo_.setSize(nRegions);


    forAll(globalMinLevel, surfI)
    {
        label nRegions = allGeometry_[surfaces_[surfI]].regions().size();

        // Initialise to global (i.e. per surface)
        for (label i = 0; i < nRegions; i++)
        {
            label globalRegionI = regionOffset_[surfI] + i;
            minLevel_[globalRegionI] = globalMinLevel[surfI];
            maxLevel_[globalRegionI] = globalMaxLevel[surfI];
            perpendicularAngle_[globalRegionI] = globalAngle[surfI];
            if (globalPatchInfo.set(surfI))
            {
                patchInfo_.set
                (
                    globalRegionI,
                    globalPatchInfo[surfI].clone()
                );
            }
        }

        // Overwrite with region specific information
        forAllConstIter(Map<label>, regionMinLevel[surfI], iter)
        {
            label globalRegionI = regionOffset_[surfI] + iter.key();

            minLevel_[globalRegionI] = iter();
            maxLevel_[globalRegionI] = regionMaxLevel[surfI][iter.key()];

            // Check validity
            if
            (
                minLevel_[globalRegionI] < 0
             || maxLevel_[globalRegionI] < minLevel_[globalRegionI]
            )
            {
                FatalErrorIn
                (
                    "refinementSurfaces::refinementSurfaces"
                    "(const searchableSurfaces&, const dictionary>&"
                )   << "Illegal level or layer specification for surface "
                    << names_[surfI]
                    << " : minLevel:" << minLevel_[globalRegionI]
                    << " maxLevel:" << maxLevel_[globalRegionI]
                    << exit(FatalError);
            }
        }
        forAllConstIter(Map<scalar>, regionAngle[surfI], iter)
        {
            label globalRegionI = regionOffset_[surfI] + iter.key();

            perpendicularAngle_[globalRegionI] = regionAngle[surfI][iter.key()];
        }

        const Map<autoPtr<dictionary> >& localInfo = regionPatchInfo[surfI];
        forAllConstIter(Map<autoPtr<dictionary> >, localInfo, iter)
        {
            label globalRegionI = regionOffset_[surfI] + iter.key();

            patchInfo_.set(globalRegionI, iter()().clone());
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Get indices of unnamed surfaces (surfaces without faceZoneName)
Foam::labelList Foam::refinementSurfaces::getUnnamedSurfaces() const
{
    labelList anonymousSurfaces(faceZoneNames_.size());

    label i = 0;
    forAll(faceZoneNames_, surfI)
    {
        if (faceZoneNames_[surfI].empty())
        {
            anonymousSurfaces[i++] = surfI;
        }
    }
    anonymousSurfaces.setSize(i);

    return anonymousSurfaces;
}


// Get indices of named surfaces (surfaces with faceZoneName)
Foam::labelList Foam::refinementSurfaces::getNamedSurfaces() const
{
   labelList namedSurfaces(faceZoneNames_.size());

    label namedI = 0;
    forAll(faceZoneNames_, surfI)
    {
        if (faceZoneNames_[surfI].size())
        {
            namedSurfaces[namedI++] = surfI;
        }
    }
    namedSurfaces.setSize(namedI);

    return namedSurfaces;
}


// Get indices of closed named surfaces
Foam::labelList Foam::refinementSurfaces::getClosedNamedSurfaces() const
{
    labelList closed(cellZoneNames_.size());

    label closedI = 0;
    forAll(cellZoneNames_, surfI)
    {
        if
        (
            cellZoneNames_[surfI].size()
         && (
                zoneInside_[surfI] == INSIDE
             || zoneInside_[surfI] == OUTSIDE
            )
         && allGeometry_[surfaces_[surfI]].hasVolumeType()
        )
        {
            closed[closedI++] = surfI;
        }
    }
    closed.setSize(closedI);

    return closed;
}


// Get indices of named surfaces with a
Foam::labelList Foam::refinementSurfaces::getInsidePointNamedSurfaces() const
{
    labelList closed(cellZoneNames_.size());

    label closedI = 0;
    forAll(cellZoneNames_, surfI)
    {
        if (cellZoneNames_[surfI].size() && zoneInside_[surfI] == INSIDEPOINT)
        {
            closed[closedI++] = surfI;
        }
    }
    closed.setSize(closedI);

    return closed;
}


// // Count number of triangles per surface region
// Foam::labelList Foam::refinementSurfaces::countRegions(const triSurface& s)
// {
//     const geometricSurfacePatchList& regions = s.patches();
//
//     labelList nTris(regions.size(), 0);
//
//     forAll(s, triI)
//     {
//         nTris[s[triI].region()]++;
//     }
//     return nTris;
// }


// // Pre-calculate the refinement level for every element
// void Foam::refinementSurfaces::wantedRefinementLevel
// (
//     const shellSurfaces& shells,
//     const label surfI,
//     const List<pointIndexHit>& info,    // Indices
//     const pointField& ctrs,             // Representative coordinate
//     labelList& minLevelField
// ) const
// {
//     const searchableSurface& geom = allGeometry_[surfaces_[surfI]];
//
//     // Get per element the region
//     labelList region;
//     geom.getRegion(info, region);
//
//     // Initialise fields to region wise minLevel
//     minLevelField.setSize(ctrs.size());
//     minLevelField = -1;
//
//     forAll(minLevelField, i)
//     {
//         if (info[i].hit())
//         {
//             minLevelField[i] = minLevel(surfI, region[i]);
//         }
//     }
//
//     // Find out if triangle inside shell with higher level
//     // What level does shell want to refine fc to?
//     labelList shellLevel;
//     shells.findHigherLevel(ctrs, minLevelField, shellLevel);
//
//     forAll(minLevelField, i)
//     {
//         minLevelField[i] = max(minLevelField[i], shellLevel[i]);
//     }
// }


// Precalculate the refinement level for every element of the searchable
// surface.
void Foam::refinementSurfaces::setMinLevelFields
(
    const shellSurfaces& shells
)
{
    forAll(surfaces_, surfI)
    {
        const searchableSurface& geom = allGeometry_[surfaces_[surfI]];

        // Precalculation only makes sense if there are different regions
        // (so different refinement levels possible) and there are some
        // elements. Possibly should have 'enough' elements to have fine
        // enough resolution but for now just make sure we don't catch e.g.
        // searchableBox (size=6)
        if (geom.regions().size() > 1 && geom.globalSize() > 10)
        {
            // Representative local coordinates
            const pointField ctrs(geom.coordinates());

            labelList minLevelField(ctrs.size(), -1);
            {
                // Get the element index in a roundabout way. Problem is e.g.
                // distributed surface where local indices differ from global
                // ones (needed for getRegion call)
                List<pointIndexHit> info;
                geom.findNearest
                (
                    ctrs,
                    scalarField(ctrs.size(), sqr(GREAT)),
                    info
                );

                // Get per element the region
                labelList region;
                geom.getRegion(info, region);

                // From the region get the surface-wise refinement level
                forAll(minLevelField, i)
                {
                    if (info[i].hit())
                    {
                        minLevelField[i] = minLevel(surfI, region[i]);
                    }
                }
            }

            // Find out if triangle inside shell with higher level
            // What level does shell want to refine fc to?
            labelList shellLevel;
            shells.findHigherLevel(ctrs, minLevelField, shellLevel);

            forAll(minLevelField, i)
            {
                minLevelField[i] = max(minLevelField[i], shellLevel[i]);
            }

            // Store minLevelField on surface
            const_cast<searchableSurface&>(geom).setField(minLevelField);
        }
    }
}


// Find intersections of edge. Return -1 or first surface with higher minLevel
// number.
void Foam::refinementSurfaces::findHigherIntersection
(
    const pointField& start,
    const pointField& end,
    const labelList& currentLevel,   // current cell refinement level

    labelList& surfaces,
    labelList& surfaceLevel
) const
{
    surfaces.setSize(start.size());
    surfaces = -1;
    surfaceLevel.setSize(start.size());
    surfaceLevel = -1;

    if (surfaces_.empty())
    {
        return;
    }

    if (surfaces_.size() == 1)
    {
        // Optimisation: single segmented surface. No need to duplicate
        // point storage.

        label surfI = 0;

        const searchableSurface& geom = allGeometry_[surfaces_[surfI]];

        // Do intersection test
        List<pointIndexHit> intersectionInfo(start.size());
        geom.findLineAny(start, end, intersectionInfo);

        // See if a cached level field available
        labelList minLevelField;
        geom.getField(intersectionInfo, minLevelField);
        bool haveLevelField =
        (
            returnReduce(minLevelField.size(), sumOp<label>())
          > 0
        );

        if (!haveLevelField && geom.regions().size() == 1)
        {
            minLevelField = labelList
            (
                intersectionInfo.size(),
                minLevel(surfI, 0)
            );
            haveLevelField = true;
        }

        if (haveLevelField)
        {
            forAll(intersectionInfo, i)
            {
                if
                (
                    intersectionInfo[i].hit()
                 && minLevelField[i] > currentLevel[i]
                )
                {
                    surfaces[i] = surfI;    // index of surface
                    surfaceLevel[i] = minLevelField[i];
                }
            }
            return;
        }
    }



    // Work arrays
    pointField p0(start);
    pointField p1(end);
    labelList intersectionToPoint(identity(start.size()));
    List<pointIndexHit> intersectionInfo(start.size());

    forAll(surfaces_, surfI)
    {
        const searchableSurface& geom = allGeometry_[surfaces_[surfI]];

        // Do intersection test
        geom.findLineAny(p0, p1, intersectionInfo);

        // See if a cached level field available
        labelList minLevelField;
        geom.getField(intersectionInfo, minLevelField);

        // Copy all hits into arguments, In-place compact misses.
        label missI = 0;
        forAll(intersectionInfo, i)
        {
            // Get the minLevel for the point
            label minLocalLevel = -1;

            if (intersectionInfo[i].hit())
            {
                // Check if minLevelField for this surface.
                if (minLevelField.size())
                {
                    minLocalLevel = minLevelField[i];
                }
                else
                {
                    // Use the min level for the surface instead. Assume
                    // single region 0.
                    minLocalLevel = minLevel(surfI, 0);
                }
            }


            label pointI = intersectionToPoint[i];

            if (minLocalLevel > currentLevel[pointI])
            {
                // Mark point for refinement
                surfaces[pointI] = surfI;
                surfaceLevel[pointI] = minLocalLevel;
            }
            else
            {
                p0[missI] = start[pointI];
                p1[missI] = end[pointI];
                intersectionToPoint[missI] = pointI;
                missI++;
            }
        }

        // All done? Note that this decision should be synchronised
        if (returnReduce(missI, sumOp<label>()) == 0)
        {
            break;
        }

        // Trim misses
        p0.setSize(missI);
        p1.setSize(missI);
        intersectionToPoint.setSize(missI);
        intersectionInfo.setSize(missI);
    }
}


void Foam::refinementSurfaces::findAllHigherIntersections
(
    const pointField& start,
    const pointField& end,
    const labelList& currentLevel,   // current cell refinement level

    List<vectorList>& surfaceNormal,
    labelListList& surfaceLevel
) const
{
    surfaceLevel.setSize(start.size());
    surfaceNormal.setSize(start.size());

    if (surfaces_.empty())
    {
        return;
    }

    // Work arrays
    List<List<pointIndexHit> > hitInfo;
    labelList pRegions;
    vectorField pNormals;

    forAll(surfaces_, surfI)
    {
        allGeometry_[surfaces_[surfI]].findLineAll(start, end, hitInfo);

        // Repack hits for surface into flat list
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // To avoid overhead of calling getRegion for every point

        label n = 0;
        forAll(hitInfo, pointI)
        {
            n += hitInfo[pointI].size();
        }

        List<pointIndexHit> surfInfo(n);
        labelList pointMap(n);
        n = 0;

        forAll(hitInfo, pointI)
        {
            const List<pointIndexHit>& pHits = hitInfo[pointI];

            forAll(pHits, i)
            {
                surfInfo[n] = pHits[i];
                pointMap[n] = pointI;
                n++;
            }
        }

        labelList surfRegion(n);
        vectorField surfNormal(n);
        allGeometry_[surfaces_[surfI]].getRegion(surfInfo, surfRegion);
        allGeometry_[surfaces_[surfI]].getNormal(surfInfo, surfNormal);

        surfInfo.clear();


        // Extract back into pointwise
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

        forAll(surfRegion, i)
        {
            label region = globalRegion(surfI, surfRegion[i]);
            label pointI = pointMap[i];

            if (maxLevel_[region] > currentLevel[pointI])
            {
                // Append to pointI info
                label sz = surfaceNormal[pointI].size();
                surfaceNormal[pointI].setSize(sz+1);
                surfaceNormal[pointI][sz] = surfNormal[i];

                surfaceLevel[pointI].setSize(sz+1);
                surfaceLevel[pointI][sz] = maxLevel_[region];
            }
        }
    }
}


void Foam::refinementSurfaces::findNearestIntersection
(
    const labelList& surfacesToTest,
    const pointField& start,
    const pointField& end,

    labelList& surface1,
    List<pointIndexHit>& hit1,
    labelList& region1,
    labelList& surface2,
    List<pointIndexHit>& hit2,
    labelList& region2
) const
{
    // 1. intersection from start to end
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Initialize arguments
    surface1.setSize(start.size());
    surface1 = -1;
    hit1.setSize(start.size());
    region1.setSize(start.size());

    // Current end of segment to test.
    pointField nearest(end);
    // Work array
    List<pointIndexHit> nearestInfo(start.size());
    labelList region;

    forAll(surfacesToTest, testI)
    {
        label surfI = surfacesToTest[testI];

        // See if any intersection between start and current nearest
        allGeometry_[surfaces_[surfI]].findLine
        (
            start,
            nearest,
            nearestInfo
        );
        allGeometry_[surfaces_[surfI]].getRegion
        (
            nearestInfo,
            region
        );

        forAll(nearestInfo, pointI)
        {
            if (nearestInfo[pointI].hit())
            {
                hit1[pointI] = nearestInfo[pointI];
                surface1[pointI] = surfI;
                region1[pointI] = region[pointI];
                nearest[pointI] = hit1[pointI].hitPoint();
            }
        }
    }


    // 2. intersection from end to last intersection
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Find the nearest intersection from end to start. Note that we initialize
    // to the first intersection (if any).
    surface2 = surface1;
    hit2 = hit1;
    region2 = region1;

    // Set current end of segment to test.
    forAll(nearest, pointI)
    {
        if (hit1[pointI].hit())
        {
            nearest[pointI] = hit1[pointI].hitPoint();
        }
        else
        {
            // Disable testing by setting to end.
            nearest[pointI] = end[pointI];
        }
    }

    forAll(surfacesToTest, testI)
    {
        label surfI = surfacesToTest[testI];

        // See if any intersection between end and current nearest
        allGeometry_[surfaces_[surfI]].findLine
        (
            end,
            nearest,
            nearestInfo
        );
        allGeometry_[surfaces_[surfI]].getRegion
        (
            nearestInfo,
            region
        );

        forAll(nearestInfo, pointI)
        {
            if (nearestInfo[pointI].hit())
            {
                hit2[pointI] = nearestInfo[pointI];
                surface2[pointI] = surfI;
                region2[pointI] = region[pointI];
                nearest[pointI] = hit2[pointI].hitPoint();
            }
        }
    }


    // Make sure that if hit1 has hit something, hit2 will have at least the
    // same point (due to tolerances it might miss its end point)
    forAll(hit1, pointI)
    {
        if (hit1[pointI].hit() && !hit2[pointI].hit())
        {
            hit2[pointI] = hit1[pointI];
            surface2[pointI] = surface1[pointI];
            region2[pointI] = region1[pointI];
        }
    }
}


void Foam::refinementSurfaces::findNearestIntersection
(
    const labelList& surfacesToTest,
    const pointField& start,
    const pointField& end,

    labelList& surface1,
    List<pointIndexHit>& hit1,
    labelList& region1,
    vectorField& normal1,

    labelList& surface2,
    List<pointIndexHit>& hit2,
    labelList& region2,
    vectorField& normal2
) const
{
    // 1. intersection from start to end
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Initialize arguments
    surface1.setSize(start.size());
    surface1 = -1;
    hit1.setSize(start.size());
    region1.setSize(start.size());
    normal1.setSize(start.size());

    // Current end of segment to test.
    pointField nearest(end);
    // Work array
    List<pointIndexHit> nearestInfo(start.size());
    labelList region;
    vectorField normal;

    forAll(surfacesToTest, testI)
    {
        label surfI = surfacesToTest[testI];
        const searchableSurface& geom = allGeometry_[surfaces_[surfI]];

        // See if any intersection between start and current nearest
        geom.findLine(start, nearest, nearestInfo);
        geom.getRegion(nearestInfo, region);
        geom.getNormal(nearestInfo, normal);

        forAll(nearestInfo, pointI)
        {
            if (nearestInfo[pointI].hit())
            {
                hit1[pointI] = nearestInfo[pointI];
                surface1[pointI] = surfI;
                region1[pointI] = region[pointI];
                normal1[pointI] = normal[pointI];
                nearest[pointI] = hit1[pointI].hitPoint();
            }
        }
    }


    // 2. intersection from end to last intersection
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Find the nearest intersection from end to start. Note that we initialize
    // to the first intersection (if any).
    surface2 = surface1;
    hit2 = hit1;
    region2 = region1;
    normal2 = normal1;

    // Set current end of segment to test.
    forAll(nearest, pointI)
    {
        if (hit1[pointI].hit())
        {
            nearest[pointI] = hit1[pointI].hitPoint();
        }
        else
        {
            // Disable testing by setting to end.
            nearest[pointI] = end[pointI];
        }
    }

    forAll(surfacesToTest, testI)
    {
        label surfI = surfacesToTest[testI];
        const searchableSurface& geom = allGeometry_[surfaces_[surfI]];

        // See if any intersection between end and current nearest
        geom.findLine(end, nearest, nearestInfo);
        geom.getRegion(nearestInfo, region);
        geom.getNormal(nearestInfo, normal);

        forAll(nearestInfo, pointI)
        {
            if (nearestInfo[pointI].hit())
            {
                hit2[pointI] = nearestInfo[pointI];
                surface2[pointI] = surfI;
                region2[pointI] = region[pointI];
                normal2[pointI] = normal[pointI];
                nearest[pointI] = hit2[pointI].hitPoint();
            }
        }
    }


    // Make sure that if hit1 has hit something, hit2 will have at least the
    // same point (due to tolerances it might miss its end point)
    forAll(hit1, pointI)
    {
        if (hit1[pointI].hit() && !hit2[pointI].hit())
        {
            hit2[pointI] = hit1[pointI];
            surface2[pointI] = surface1[pointI];
            region2[pointI] = region1[pointI];
            normal2[pointI] = normal1[pointI];
        }
    }
}


void Foam::refinementSurfaces::findAnyIntersection
(
    const pointField& start,
    const pointField& end,

    labelList& hitSurface,
    List<pointIndexHit>& hitInfo
) const
{
    searchableSurfacesQueries::findAnyIntersection
    (
        allGeometry_,
        surfaces_,
        start,
        end,
        hitSurface,
        hitInfo
    );
}


void Foam::refinementSurfaces::findNearest
(
    const labelList& surfacesToTest,
    const pointField& samples,
    const  scalarField& nearestDistSqr,
    labelList& hitSurface,
    List<pointIndexHit>& hitInfo
) const
{
    labelList geometries(UIndirectList<label>(surfaces_, surfacesToTest));

    // Do the tests. Note that findNearest returns index in geometries.
    searchableSurfacesQueries::findNearest
    (
        allGeometry_,
        geometries,
        samples,
        nearestDistSqr,
        hitSurface,
        hitInfo
    );

    // Rework the hitSurface to be surface (i.e. index into surfaces_)
    forAll(hitSurface, pointI)
    {
        if (hitSurface[pointI] != -1)
        {
            hitSurface[pointI] = surfacesToTest[hitSurface[pointI]];
        }
    }
}


void Foam::refinementSurfaces::findNearestRegion
(
    const labelList& surfacesToTest,
    const pointField& samples,
    const scalarField& nearestDistSqr,
    labelList& hitSurface,
    labelList& hitRegion
) const
{
    labelList geometries(UIndirectList<label>(surfaces_, surfacesToTest));

    // Do the tests. Note that findNearest returns index in geometries.
    List<pointIndexHit> hitInfo;
    searchableSurfacesQueries::findNearest
    (
        allGeometry_,
        geometries,
        samples,
        nearestDistSqr,
        hitSurface,
        hitInfo
    );

    // Rework the hitSurface to be surface (i.e. index into surfaces_)
    forAll(hitSurface, pointI)
    {
        if (hitSurface[pointI] != -1)
        {
            hitSurface[pointI] = surfacesToTest[hitSurface[pointI]];
        }
    }

    // Collect the region
    hitRegion.setSize(hitSurface.size());
    hitRegion = -1;

    forAll(surfacesToTest, i)
    {
        label surfI = surfacesToTest[i];

        // Collect hits for surfI
        const labelList localIndices(findIndices(hitSurface, surfI));

        List<pointIndexHit> localHits
        (
            UIndirectList<pointIndexHit>
            (
                hitInfo,
                localIndices
            )
        );

        labelList localRegion;
        allGeometry_[surfaces_[surfI]].getRegion(localHits, localRegion);

        forAll(localIndices, i)
        {
            hitRegion[localIndices[i]] = localRegion[i];
        }
    }
}


void Foam::refinementSurfaces::findNearestRegion
(
    const labelList& surfacesToTest,
    const pointField& samples,
    const scalarField& nearestDistSqr,
    labelList& hitSurface,
    List<pointIndexHit>& hitInfo,
    labelList& hitRegion,
    vectorField& hitNormal
) const
{
    labelList geometries(UIndirectList<label>(surfaces_, surfacesToTest));

    // Do the tests. Note that findNearest returns index in geometries.
    searchableSurfacesQueries::findNearest
    (
        allGeometry_,
        geometries,
        samples,
        nearestDistSqr,
        hitSurface,
        hitInfo
    );

    // Rework the hitSurface to be surface (i.e. index into surfaces_)
    forAll(hitSurface, pointI)
    {
        if (hitSurface[pointI] != -1)
        {
            hitSurface[pointI] = surfacesToTest[hitSurface[pointI]];
        }
    }

    // Collect the region
    hitRegion.setSize(hitSurface.size());
    hitRegion = -1;
    hitNormal.setSize(hitSurface.size());
    hitNormal = vector::zero;

    forAll(surfacesToTest, i)
    {
        label surfI = surfacesToTest[i];

        // Collect hits for surfI
        const labelList localIndices(findIndices(hitSurface, surfI));

        List<pointIndexHit> localHits
        (
            UIndirectList<pointIndexHit>
            (
                hitInfo,
                localIndices
            )
        );

        // Region
        labelList localRegion;
        allGeometry_[surfaces_[surfI]].getRegion(localHits, localRegion);

        forAll(localIndices, i)
        {
            hitRegion[localIndices[i]] = localRegion[i];
        }

        // Normal
        vectorField localNormal;
        allGeometry_[surfaces_[surfI]].getNormal(localHits, localNormal);

        forAll(localIndices, i)
        {
            hitNormal[localIndices[i]] = localNormal[i];
        }
    }
}


//// Find intersection with max of edge. Return -1 or the surface
//// with the highest maxLevel above currentLevel
//Foam::label Foam::refinementSurfaces::findHighestIntersection
//(
//    const point& start,
//    const point& end,
//    const label currentLevel,   // current cell refinement level
//
//    pointIndexHit& maxHit
//) const
//{
//    // surface with highest maxlevel
//    label maxSurface = -1;
//    // maxLevel of maxSurface
//    label maxLevel = currentLevel;
//
//    forAll(*this, surfI)
//    {
//        pointIndexHit hit = operator[](surfI).findLineAny(start, end);
//
//        if (hit.hit())
//        {
//            const triSurface& s = operator[](surfI);
//
//            label region = globalRegion(surfI, s[hit.index()].region());
//
//            if (maxLevel_[region] > maxLevel)
//            {
//                maxSurface = surfI;
//                maxLevel = maxLevel_[region];
//                maxHit = hit;
//            }
//        }
//    }
//
//    if (maxSurface == -1)
//    {
//        // maxLevel unchanged. No interesting surface hit.
//        maxHit.setMiss();
//    }
//
//    return maxSurface;
//}


void Foam::refinementSurfaces::findInside
(
    const labelList& testSurfaces,
    const pointField& pt,
    labelList& insideSurfaces
) const
{
    insideSurfaces.setSize(pt.size());
    insideSurfaces = -1;

    forAll(testSurfaces, i)
    {
        label surfI = testSurfaces[i];

        if (zoneInside_[surfI] != INSIDE && zoneInside_[surfI] != OUTSIDE)
        {
            FatalErrorIn("refinementSurfaces::findInside(..)")
                << "Trying to use surface "
                << allGeometry_[surfaces_[surfI]].name()
                << " which has non-geometric inside selection method "
                << areaSelectionAlgoNames[zoneInside_[surfI]]
                << exit(FatalError);
        }

        if (allGeometry_[surfaces_[surfI]].hasVolumeType())
        {
            List<searchableSurface::volumeType> volType;
            allGeometry_[surfaces_[surfI]].getVolumeType(pt, volType);

            forAll(volType, pointI)
            {
                if (insideSurfaces[pointI] == -1)
                {
                    if
                    (
                        (
                            volType[pointI] == triSurfaceMesh::INSIDE
                         && zoneInside_[surfI] == INSIDE
                        )
                     || (
                            volType[pointI] == triSurfaceMesh::OUTSIDE
                         && zoneInside_[surfI] == OUTSIDE
                        )
                    )
                    {
                        insideSurfaces[pointI] = surfI;
                    }
                }
            }
        }
    }
}


// ************************************************************************* //
