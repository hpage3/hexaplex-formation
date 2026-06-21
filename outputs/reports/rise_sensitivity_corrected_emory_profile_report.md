# Rise Sensitivity Against Emory-Corrected Profile

## Purpose

This controlled test evaluates Nick's hypothesis that the clarified ideal antiparallel 30-degree model may fit the Emory-corrected powder profile better at a helical rise near 3.38 A than at the current ideal 3.40 A.

## Inputs

- Source ideal PDB: `inputs/nick_ideal_models/Hexaplex_AntiParallel_30deg_Ideal.pdb`
- Corrected experimental profile: `inputs\experimental\nick_powder_profile_corrected_emory.csv`
- Generated rise variants: `inputs\nick_ideal_models\rise_variants`

## Rise Variant Generation

The source PDB does not store explicit layer IDs. The generator inferred 15 six-residue CYP/MEP base planes from residue z-centroids; each plane contains three CYP and three MEP residues. Non-base GLU residues were assigned to the nearest inferred base plane and translated rigidly with that plane. The central inferred base plane was held fixed as the z anchor.

This perturbation preserves atom order, x/y coordinates, twist geometry, and intralayer geometry. It only changes rigid layer z positions. It does not adjust peptide omega, side-chain removal logic, twist angle, or perform minimization.

## Validation

| Target rise (A) | Atom count preserved | Measured rise (A) | Max rise error (A) | Max layer z delta (A) |
| ---: | --- | ---: | ---: | ---: |
| 3.36 | True | 3.3600 | 1.33e-15 | 0.2892 |
| 3.37 | True | 3.3700 | 4.44e-15 | 0.2192 |
| 3.38 | True | 3.3800 | 4.44e-15 | 0.1492 |
| 3.39 | True | 3.3900 | 1.33e-15 | 0.0792 |
| 3.40 | True | 3.4000 | 2.22e-15 | 0.0092 |

## Corrected Diffraction Settings

- Engine: `reference/asem_corrected_diffraction_engine/`
- Asem correction: azimuthal rotations are applied to independent tilted coordinate stacks.
- Tilts: `[0]`
- Rotations: `range(0, 181, 5)`
- Grid: 129 x 129 over +/-100.0 mm
- Detector distance: 338.4 mm
- Wavelength: 0.7749 A
- Radial bins: 420
- Hydrogens excluded and exact heavy-atom deduplication applied in XYZ inputs.

## Ranking

| Rise (A) | Mean abs primary offset (A) | Max abs primary offset (A) | Base-stack offset (A) | 3.77 offset (A) | 4.4 offset (A) | 5.6 offset (A) | 7.3 offset (A) | Rank |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 3.40 | 0.0634 | 0.1374 | 0.0348 | -0.0003 | 0.1374 | 0.1264 | 0.0181 | 1 |
| 3.38 | 0.0648 | 0.1264 | 0.0348 | -0.0864 | -0.0582 | 0.1264 | 0.0181 | 2 |
| 3.39 | 0.0648 | 0.1264 | 0.0348 | -0.0864 | -0.0582 | 0.1264 | 0.0181 | 3 |
| 3.36 | 0.0740 | 0.1264 | 0.0348 | -0.0864 | -0.1044 | 0.1264 | 0.0181 | 4 |
| 3.37 | 0.0740 | 0.1264 | 0.0348 | -0.0864 | -0.1044 | 0.1264 | 0.0181 | 5 |

## Interpretation

- Best rise by mean primary-window peak offset: 3.40 A.
- The 3.38 A variant does not improve relative to 3.40 A by mean primary-window peak offset.
- 3.38 A mean abs primary offset: 0.0648 A.
- 3.40 A mean abs primary offset: 0.0634 A.

Feature-specific offsets are in `outputs/metrics/rise_sensitivity_feature_summary_corrected_emory_profile.csv`. The comparison is a controlled sensitivity test; a favorable 3.38 A result would support asking Asem for a larger/refined 3.38 A model, not replacing a chemically relaxed model.

## Plots

- `outputs/rise_sensitivity_corrected_emory_profile/plots/rise_variant_profile_overlay.png`
- `outputs/rise_sensitivity_corrected_emory_profile/plots/rise_variant_peak_offsets.png`
- `outputs/rise_sensitivity_corrected_emory_profile/plots/rise_variant_mean_error.png`
- `outputs/rise_sensitivity_corrected_emory_profile/plots/rise_3p38_vs_3p40_nick_style.png`

The Nick-style plot uses Gaussian smoothing in d-space with sigma 0.10 A for visualization only. Scoring uses unsmoothed radial profiles.

## Limitations

- Rigid-layer z-translation does not relax backbone geometry.
- This does not address peptide omega.
- This does not alter twist angle.
- Powder/radial matching is a sensitivity test, not full refinement.
- Layer inference is based on CYP/MEP base-plane z-centroids because explicit layer IDs are not present in the source PDB.
