# Parallel-Sheet Control Against Emory-Corrected Profile

## Purpose

Run the idealized parallel-sheet model as a control/falsification comparison against Nick's clarified ideal antiparallel 30-degree hexaplex model.

## Inputs

- Parallel control PDB: `inputs\nick_ideal_models\parallel_control\TwoBetaSheetBackbones180_180.pdb`
- Parallel control XYZ: `inputs\nick_ideal_models\parallel_control\TwoBetaSheetBackbones180_180.xyz`
- Antiparallel baseline radial profile: `outputs\nick_ideal_16mer_corrected_emory_profile\radial_profiles\ideal_16mer_antiparallel_30deg_radial_profile.csv`
- Corrected experimental profile: `inputs\experimental\nick_powder_profile_corrected_emory.csv`
- Comparison plot: `outputs\parallel_control_corrected_emory_profile\plots\parallel_vs_antiparallel_corrected_emory_profile.png`

## Conversion

- PDB ATOM/HETATM count: 120
- PDB record counts: `{"ATOM": 120}`
- Source hydrogen count: 20
- Heavy atoms before deduplication: 100
- Duplicate heavy records removed: 0
- Final XYZ atom count: 100
- XYZ element counts: `{"C": 60, "N": 20, "O": 20}`

## Corrected Diffraction Settings

- Engine: `reference/asem_corrected_diffraction_engine/`.
- Asem correction: non-accumulating/vectorized orientation path.
- Tilts: `[0]`.
- Rotations: `range(0, 181, 5)`, 37 rotations.
- Grid size: 129 x 129; detector half-width 100 mm; detector distance 338.4 mm; wavelength 0.7749 A; radial bins 420.

## Primary Feature Windows

| feature_window | experimental_peak_d_A | antiparallel_peak_d_A | antiparallel_offset_d_A | parallel_peak_d_A | parallel_offset_d_A | closer_model_by_abs_offset |
| --- | --- | --- | --- | --- | --- | --- |
| 3.4 A | 3.381 | 3.41583473048 | 0.03483473048 | 3.488064464 | 0.107064464 | antiparallel |
| 3.7 A | 3.795 | 3.79471639446 | -0.00028360554 | 3.7771670273 | -0.0178329727 | antiparallel |
| 4.4 A | 4.401 | 4.34277828279 | -0.05822171721 | 4.46291994427 | 0.06191994427 | antiparallel |
| 5.6 A | 5.577 | 5.58327265748 | 0.00627265748 | 5.62272176729 | 0.04572176729 | antiparallel |
| 7.25 A | 7.325 | 7.34311339524 | 0.01811339524 | 7.34311339524 | 0.01811339524 | tie |

## Interpretation

The feature-window comparison supports the antiparallel 30-degree ideal model over the idealized parallel-sheet control under the corrected workflow.

Across the five primary windows, antiparallel is closer in 4 window(s), while parallel is closer in 0 window(s). Avoid reading this as a full structural refinement; it is a falsification-style control against one idealized wrong-geometry model.

## Limitations

- Powder profiles are orientationally/radially averaged.
- The result depends on the idealized alanine parallel-sheet control geometry.
- This is a negative-control/falsification test, not a global phase or structure refinement.
- The plot uses d-space Gaussian smoothing with sigma 0.10 A for visualization only; metrics use the unsmoothed radial/profile CSVs.
