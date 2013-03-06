
IOdictionary planeToPatchDict
(
    IOobject
    (
        "planeToPatchDict",
        runTime.system(),
        runTime,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Switch planeMapping
(
    planeToPatchDict.lookup("mappingPlane")
);

Info << "Constructing planePatch" << endl;
cuttingPlane cut(plane(planeToPatchDict), mesh);

primitiveFacePatch planePatch
(
    cut.faces(),
    cut.points()
);

Info << "Finding target patch" << endl;

word toPatchName
(
    planeToPatchDict.lookup("targetPatchName")
);

const fvPatchList& patches = mesh.boundary();
bool targetPatchNamePresent(false);
label targetPatchNumber = -1;

forAll(patches, patchI)
{
    if (patches[patchI].name() == toPatchName)
    {
        targetPatchNamePresent = true;
        targetPatchNumber = patchI;
    }
}

const labelList& cutCells(cut.cells());

if (!targetPatchNamePresent)
{
    targetPatchNumber = 1;
    Serr<< "Target patch not present. " << endl;
}
else
{
    Serr<< "Target patch name: " << patches[targetPatchNumber].name() << endl;
}

patchToPatchInterpolation interPatch
(
    planePatch,
    patches[targetPatchNumber].patch()
);


PrimitivePatchInterpolation
<
    PrimitivePatch<face, List, const pointField&>
> planePatchInter
(
    planePatch
);


PrimitivePatchInterpolation
<
    PrimitivePatch<face, List, const pointField&>
> toPatchInter
(
    patches[targetPatchNumber].patch()
);
