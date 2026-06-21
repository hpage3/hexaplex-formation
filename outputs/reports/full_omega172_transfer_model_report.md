# Full Omega172 Transfer Model Report

## Purpose

This run attempted a conservative transfer of Asem's omega172 local-coordinate adjustment into the full 7,146-atom ideal antiparallel 30-degree model.

## Strict Stop Condition

The workflow is allowed to build a derived full PDB only if residue/atom correspondence is clear and validation can prove that atom count, rise, twist, and unmoved coordinates are preserved.

## Result

Transfer stopped. No derived full PDB, XYZ, diffraction profile, metrics, or plots were generated.

## Mapping Method Tested

The script compared the local omega167 and omega172 files by atom serial/name/residue identity, then attempted to map each local CYP/MEP residue to the best full-model CYP/MEP residue of the same residue name using shared unmoved atom names. Direct absolute-coordinate RMSD was used as the conservative transfer validation criterion.

## Stop Reason

No high-confidence direct coordinate correspondence was established. Although the full model contains 90 CYP/MEP residues and the local struct2 model also contains 90 CYP/MEP residues, the best absolute-coordinate matches of shared unmoved atoms have multi-angstrom RMSD (range 4.383-50.078 A), far above the 0.10 A transfer tolerance. Copying raw omega172 displacement vectors into the full model would therefore perturb the wrong coordinate frame.

## Validation Headline

- Local moved atom names: H'':84;N'':84;O':84
- Full atom count inspected: 7146
- Derived full atom count: not applicable; no derived model was written.
- Rise/twist/global coordinates: unchanged because no transfer was performed.
- Corrected diffraction: not run.

## Recommendation

Reject this transfer attempt as ambiguous and ask Asem to build the official full periodic omega172 model directly. The derived coordinate-transfer model should not be used for discussion because it was not generated.

## Limitations

- Coordinate-transfer is not energy minimization.
- Atom-name mapping is not a chemically exact torsion definition.
- Asem/Nick visual chemistry review is still needed.
- No diffraction result exists from this stopped transfer attempt.
