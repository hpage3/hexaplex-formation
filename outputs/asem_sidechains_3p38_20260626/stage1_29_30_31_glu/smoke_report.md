# Asem 3.38 A Side-Chain 29/30/31 Smoke Report

## Scope

This smoke validation converted the planned 29/30/31 complete Glu side-chain structures to heavy-atom XYZ files, then ran corrected diffraction and radial averaging for exactly one representative model: `angle30_cand0`.

No raw PDB files were modified. The remaining 49 corrected-diffraction jobs were not run.

## XYZ Conversion

- Planned rows converted: 50
- Converted status rows: 50
- Hydrogens were excluded using the existing `reference/asem_corrected_diffraction_engine/pdb_to_xyz.py` converter with `dedupe_exact=True`.
- Each converted XYZ has 2310 heavy atoms and 0 hydrogens.
- Element counts in XYZ use the converter element policy: `C:1152; N:624; O:534`. This differs from the inventory fallback that reported alpha-carbon atom names as `Ca`; the XYZ files correctly contain no calcium atoms.

## Corrected Diffraction Smoke

- Candidate: `angle30_cand0`
- Input XYZ: `outputs/asem_sidechains_3p38_20260626/stage1_29_30_31_glu/xyz/angle30_cand0.xyz`
- Grid size: 65
- Grid limit: 100 mm
- Theta count: 4
- Phi count: 8
- Psi count: 1
- Theta max: 180 deg
- Orientation count: 32
- Mode: serial, workers 1
- Status: success
- Runtime: 1.755729400 seconds
- Atom count in diffraction metadata: 2310
- Image shape: [65, 65]
- Image intensity range: 258089.10249976328 to 7739670527.9875965

## Radial Smoke

- Radial bins requested: 120
- Non-empty radial profile rows written: 113
- d-spacing range in written rows: 1.977706 to 445.013056 A
- q range in written rows: 0.014119 to 3.177007 A^-1
- Radial q-space and d-space plots were written for diagnostic review.

## Output Files

- `outputs/asem_sidechains_3p38_20260626/stage1_29_30_31_glu/xyz/angle30_cand0.xyz`
- `outputs/asem_sidechains_3p38_20260626/stage1_29_30_31_glu/sidechain_3p38_xyz_manifest.csv`
- `outputs/asem_sidechains_3p38_20260626/stage1_29_30_31_glu/smoke_corrected_diffraction/analysis_manifest.csv`
- `outputs/asem_sidechains_3p38_20260626/stage1_29_30_31_glu/smoke_corrected_diffraction/angle30_cand0.npy`
- `outputs/asem_sidechains_3p38_20260626/stage1_29_30_31_glu/smoke_corrected_diffraction/angle30_cand0.png`
- `outputs/asem_sidechains_3p38_20260626/stage1_29_30_31_glu/smoke_corrected_diffraction/angle30_cand0.json`
- `outputs/asem_sidechains_3p38_20260626/stage1_29_30_31_glu/smoke_radial_profiles/angle30_cand0_radial.csv`
- `outputs/asem_sidechains_3p38_20260626/stage1_29_30_31_glu/smoke_radial_profiles/angle30_cand0_radial_q.png`
- `outputs/asem_sidechains_3p38_20260626/stage1_29_30_31_glu/smoke_radial_profiles/angle30_cand0_radial_d.png`

## Consistency Check

The smoke output is structurally consistent with the prior 3.40 A corrected-diffraction workflow: heavy-atom XYZ input, corrected non-accumulating orientation averaging, 65 x 65 grid, 32 orientations, and 120-bin radial averaging. This smoke pass does not evaluate scientific fit quality by itself.

## Recommendation

If the smoke outputs are visually and procedurally acceptable, the next step is to run the remaining planned 29/30/31 corrected-diffraction jobs using the existing 50-row XYZ manifest. Keep those generated diffraction/radial outputs uncommitted until reviewed.
