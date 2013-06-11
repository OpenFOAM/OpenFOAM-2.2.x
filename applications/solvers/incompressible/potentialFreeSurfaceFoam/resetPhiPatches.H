FieldField<fvsPatchField, scalar>& phibf = phi.boundaryField();
const FieldField<fvPatchField, vector>& Ubf = U.boundaryField();
const FieldField<fvsPatchField, vector>& Sfbf = mesh.Sf().boundaryField();

forAll(phibf, patchI)
{
    phibf[patchI] = (Ubf[patchI] & Sfbf[patchI]);
}
