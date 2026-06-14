# Base-Length Variant Generation Report

## Scope and cautions

- This is a geometry sensitivity study, not a structural determination.
- CYP/MEP candidate arm atoms were axis-facing in the fitted-axis inspection; variants scale local anchor-to-atom vectors, not global outward radial vectors.
- The diffraction sweep is not run by this workflow.
- Native PDB residue names, chain IDs, atom order, and fixed atom coordinates are preserved by construction before PDB formatting.

## Baseline

- Input PDB: outputs/intermediates/ai_candidate_inputs/full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain.pdb

## Operational definition

Base/hexad-arm length is defined as the local distance from each CYP/MEP residue backbone anchor to its selected non-backbone candidate arm atoms.

## Atom-selection rule

- Fixed atoms: all GLU atoms and CYP/MEP N, CA, C, O, OXT atoms.
- Transformable atoms from geometry inspection: CYP:C2,C4,C6,CB,CC,HB2,HB3,HC2,HC3,HN3,HN5,N1,N3,N5,OC2,OC4,OC6; MEP:C1,C3,C5,CB,CC,HB2,HB3,HC2,HC3,HX11,HX21,HX22,HX31,HX32,N2,N4,N6,NX1,NX2,NX3

## Anchor and transformation rules

- Anchor rule: CA if present; otherwise centroid of available backbone-like N/CA/C/O/OXT atoms; otherwise residue centroid.
- Transformation rule: `new_xyz = anchor_xyz + scale_factor * (old_xyz - anchor_xyz)`.

## Generated scale factors

- 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20

## Geometry sanity-check summary

- Variants generated: 8
- Transformed atom count(s): 1665
- Scale 1.00 matched baseline coordinates before formatting: scale 1.00 baseline coordinate match
- Total suspicious overlap counts across variants: 135

| variant_id | scale_factor | transformed_atom_count | min_heavy_atom_distance_A | suspicious_overlap_count | mean_transformed_atom_displacement_A | max_transformed_atom_displacement_A | notes | warnings |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| hexaplex_base_length_scale_0p85 | 0.85 | 1665 | 1.081819 | 0 | 0.689436 | 1.313214 |  |  |
| hexaplex_base_length_scale_0p90 | 0.90 | 1665 | 1.145455 | 0 | 0.459624 | 0.875476 |  |  |
| hexaplex_base_length_scale_0p95 | 0.95 | 1665 | 1.209092 | 0 | 0.229812 | 0.437738 |  |  |
| hexaplex_base_length_scale_1p00 | 1.00 | 1665 | 1.227463 | 0 | 0.000000 | 0.000000 | scale 1.00 baseline coordinate match |  |
| hexaplex_base_length_scale_1p05 | 1.05 | 1665 | 1.227463 | 0 | 0.229812 | 0.437738 |  |  |
| hexaplex_base_length_scale_1p10 | 1.10 | 1665 | 1.227463 | 0 | 0.459624 | 0.875476 |  |  |
| hexaplex_base_length_scale_1p15 | 1.15 | 1665 | 1.151122 | 0 | 0.689436 | 1.313214 |  |  |
| hexaplex_base_length_scale_1p20 | 1.20 | 1665 | 0.643109 | 135 | 0.919248 | 1.750952 |  | 135 heavy-atom pair(s) below 1.00 A |
