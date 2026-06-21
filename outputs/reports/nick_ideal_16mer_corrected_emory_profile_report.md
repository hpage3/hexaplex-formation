# Nick Ideal 16-mer Corrected Emory-Profile Comparison

## Purpose

Nick clarified that `Hexaplex_AntiParallel_30deg_Ideal.pdb` is the coordinate file to treat as the 16-mer simulation. This run establishes that ideal antiparallel 30-degree full-hexaplex baseline against the Emory-corrected experimental powder profile.

This commit focuses only on the full ideal 16-mer baseline. No-side-chain/no-COO and side-chain-subtraction comparisons require a matched derived ideal no-side-chain model and should remain a separate follow-up.

## Inputs

- Ideal PDB: `inputs\nick_ideal_models\Hexaplex_AntiParallel_30deg_Ideal.pdb`
- Ideal XYZ: `inputs\nick_ideal_models\Hexaplex_AntiParallel_30deg_Ideal.xyz`
- Corrected experimental profile: `inputs\experimental\nick_powder_profile_corrected_emory.csv`
- Experimental rows: 442
- Experimental d-spacing range: 2.222 A to 9.828 A

## PDB-to-XYZ Conversion

- PDB ATOM/HETATM count: 7146
- PDB record counts: `{"ATOM": 7146}`
- PDB element counts: `{"C": 2160, "H": 2814, "N": 1170, "O": 1002}`
- Heavy atoms before exact deduplication: 4332
- XYZ atom count: 2166
- XYZ element counts: `{"C": 1080, "N": 585, "O": 501}`
- Hydrogens included in diffraction XYZ: no; source PDB hydrogens: 2814
- Exact deduplication applied: yes; removed 2166 duplicate heavy-atom records.

## Corrected Diffraction Path

- Engine: `reference/asem_corrected_diffraction_engine/`.
- Asem correction: azimuthal rotations are applied to independent tilted coordinate stacks rather than accumulated through the rotation loop.
- Vectorized path: `make_oriented_coords` plus `generate_fiber_diffraction_series`.
- Nick-style tilts: `[0]`.
- Nick-style rotations: `range(0, 181, 5)`, producing 37 rotations from 0 to 180 degrees.
- The ambiguous `range(0, 5, 180)` form was not used.

## Detector and Radial Settings

- Grid size: 129 x 129
- Detector half-width: 100 mm
- Detector distance: 338.4 mm
- Wavelength: 0.7749 A
- Radial bins: 420
- Profile normalization: max mean intensity at q >= 0.15 A^-1.
- Plotting choice: both traces are independently normalized; the theoretical trace is vertically offset by +1.15 for visual comparison.

## Primary Feature Windows

| feature_window | simulated_peak_d_A | simulated_peak_intensity_norm | experimental_peak_d_A | experimental_peak_intensity_norm | peak_offset_d_A | window_area_simulated | window_area_experimental |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 3.4 A | 3.41583473048 | 0.1672543003 | 3.381 | 0.939561294337 | 0.0348347304828 | 0.00924415682505 | 0.145579534848 |
| 3.7 A | 3.79471639446 | 0.0171786178889 | 3.795 | 0.553671437461 | -0.000283605539219 | 0.00160209779862 | 0.0660399424393 |
| 4.4 A | 4.34277828279 | 0.0168672582233 | 4.401 | 1 | -0.0582217172089 | 0.00243247005096 | 0.173910547604 |
| 5.6 A | 5.58327265748 | 0.0205122756678 | 5.577 | 0.506922837586 | 0.00627265747814 | 0.00259644821034 | 0.0921098319851 |
| 7.25 A | 7.34311339524 | 0.0132520342606 | 7.325 | 0.40969197262 | 0.0181133952385 | 0.00337070319866 | 0.175837429994 |

## Full Feature Summary

| feature_window | feature_group | simulated_peak_d_A | simulated_peak_intensity_norm | experimental_peak_d_A | experimental_peak_intensity_norm | peak_offset_d_A | window_area_simulated | window_area_experimental |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 3.4 A | primary | 3.41583473048 | 0.1672543003 | 3.381 | 0.939561294337 | 0.0348347304828 | 0.00924415682505 | 0.145579534848 |
| 3.7 A | primary | 3.79471639446 | 0.0171786178889 | 3.795 | 0.553671437461 | -0.000283605539219 | 0.00160209779862 | 0.0660399424393 |
| 4.4 A | primary | 4.34277828279 | 0.0168672582233 | 4.401 | 1 | -0.0582217172089 | 0.00243247005096 | 0.173910547604 |
| 5.6 A | primary | 5.58327265748 | 0.0205122756678 | 5.577 | 0.506922837586 | 0.00627265747814 | 0.00259644821034 | 0.0921098319851 |
| 7.25 A | primary | 7.34311339524 | 0.0132520342606 | 7.325 | 0.40969197262 | 0.0181133952385 | 0.00337070319866 | 0.175837429994 |
| 3.0 A | optional | 3.07477148552 | 0.0110191747821 | 3.095 | 0.0515712507778 | -0.020228514475 | 0.000884107505373 | 0.00229873210952 |
| 4.1 A | optional | 4.16388921441 | 0.010003539525 | 4.198 | 0.315494710641 | -0.0341107855909 | 0.00128514759222 | 0.0162584007467 |
| 4.5-5.0 A | optional | 4.53836680174 | 0.0179278266289 | 4.5 | 0.660003111388 | 0.0383668017393 | 0.00335376040779 | 0.0629276991288 |
| 5.5 A | optional | 5.58327265748 | 0.0205122756678 | 5.577 | 0.506922837586 | 0.00627265747814 | 0.00252901864423 | 0.0755768123833 |
| 7.0 A | optional | 7.07955639782 | 0.00723436120388 | 7.166 | 0.360843186061 | -0.086443602184 | 0.00182812855232 | 0.0801477520224 |
| 8.4 A | optional | 8.4458409689 | 0.0105532512644 | 8.206 | 0.106798382078 | 0.239840968896 | 0.0015448337912 | 0.0281979231487 |

## Headline Interpretation

- 3.4 A: the simulated ideal 16-mer peak at 3.41583473048 A is close to the corrected experimental base-stacking peak at 3.381 A, with offset 0.0348347304828 A.
- Treat this as the corrected full ideal baseline only, not a side-chain attribution result.
