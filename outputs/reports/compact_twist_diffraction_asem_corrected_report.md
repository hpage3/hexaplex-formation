# Compact Twist Diffraction Report: Asem-Corrected

## Inputs

Compact twist PDBs were copied into `inputs/compact_twist_variants` and the original source PDBs were not modified.

## Corrected Workflow

- Diffraction code path: `../fiber-diffraction/scripts.py` and `../fiber-diffraction/orientation_average.py`.
- Correction: azimuthal rotations are applied to the same tilted coordinate frame rather than to already-rotated coordinates.
- Vectorization: `numexpr` is installed in `.venv` and `scripts._HAVE_NUMEXPR` is true.
- PDB conversion: heavy atoms only with exact duplicate removal.
- Diffraction settings: grid size 41, grid limit 100 mm, theta-count 4, phi-count 8, psi-count 3, theta max 180 degrees, workers 1.
- Radial normalization: maximum mean intensity at q >= 0.4 A^-1 per twist.

## Outputs

- Corrected output root: `outputs\compact_twist_diffraction_asem_corrected`
- `outputs/metrics/compact_twist_feature_summary_asem_corrected.csv`
- `outputs/metrics/compact_twist_feature_rankings_asem_corrected.csv`
- `outputs/reports/compact_twist_diffraction_asem_corrected_report.md`

## Corrected Feature Summary

| twist | d_3p4_max_norm_intensity | d_3p0_max_norm_intensity | d_4p5_5p0_max_norm_intensity | d_4p1_max_norm_intensity | d_5p5_max_norm_intensity | d_7p0_max_norm_intensity | d_8p4_max_norm_intensity |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 24.0 | 0.484621613095 | 0.0882802164232 | 0.227380723545 | 0.11200844793 | 0.191131433999 | 0.242704353075 | 0.650807057975 |
| 26.0 | 0.310529844568 | 0.0448880235664 | 0.149993070238 | 0.0816967017879 | 0.0991201456123 | 0.153835517233 | 0.190165008932 |
| 28.0 | 0.413180266706 | 0.0551322297491 | 0.149196617654 | 0.103496916764 | 0.171921463687 | 0.187494173726 | 0.180101618064 |
| 30.0 | 1 | 0.163802245368 | 0.290542680173 | 0.229061774559 | 0.512411856945 | 0.391869083071 | 0.340482271437 |
| 32.0 | 0.402296101395 | 0.0552232747844 | 0.145663271183 | 0.102821969032 | 0.171331748982 | 0.179403018141 | 0.180788984166 |
| 34.0 | 0.305709051216 | 0.0423627263951 | 0.138140290427 | 0.0820704404507 | 0.0965434397843 | 0.141058978593 | 0.176868616134 |
| 36.0 | 0.475783419483 | 0.0890746490716 | 0.212579872924 | 0.102454288828 | 0.182821831871 | 0.239029256436 | 0.63040862087 |

## Corrected Top Feature Responses

| feature_window | rank | twist_deg | normalized_intensity | d_A_at_max |
| --- | --- | --- | --- | --- |
| 3.4 A | 1 | 30 | 1 | 3.45606988144 |
| 3.0 A | 1 | 30 | 0.163802245368 | 2.94348441188 |
| 4.5-5.0 A | 1 | 30 | 0.290542680173 | 4.57871619776 |
| 4.1 A | 1 | 30 | 0.229061774559 | 4.03819213836 |
| 5.5 A | 1 | 30 | 0.512411856945 | 5.58572165085 |
| 7.0 A | 1 | 30 | 0.391869083071 | 7.18019485812 |
| 8.4 A | 1 | 24 | 0.650807057975 | 8.38392522748 |

## Pre/Post Asem-Fix Comparison

| feature_window | pre_Asem_fix_top_twist | corrected_top_twist | changed |
| --- | --- | --- | --- |
| 3.4 A | 30 | 30 | no |
| 3.0 A | 30 | 30 | no |
| 4.5-5.0 A | 30 | 30 | no |
| 4.1 A | 30 | 30 | no |
| 5.5 A | 30 | 30 | no |
| 7.0 A | 30 | 30 | no |
| 8.4 A | 24 | 24 | no |

## Interpretation

- Strongest twist per feature window changed after the correction: no.
- 30 degrees remains strongest in the primary windows (3.4 A, 3.0 A, 4.5-5.0 A): yes.
- The 8.4 A exception remains: yes.
- New corrected full rankings:
  - 3.4 A: 30 deg (1.000000), 24 deg (0.484622), 36 deg (0.475783), 28 deg (0.413180), 32 deg (0.402296), 26 deg (0.310530), 34 deg (0.305709).
  - 3.0 A: 30 deg (0.163802), 36 deg (0.089075), 24 deg (0.088280), 32 deg (0.055223), 28 deg (0.055132), 26 deg (0.044888), 34 deg (0.042363).
  - 4.5-5.0 A: 30 deg (0.290543), 24 deg (0.227381), 36 deg (0.212580), 26 deg (0.149993), 28 deg (0.149197), 32 deg (0.145663), 34 deg (0.138140).
  - 4.1 A: 30 deg (0.229062), 24 deg (0.112008), 28 deg (0.103497), 32 deg (0.102822), 36 deg (0.102454), 34 deg (0.082070), 26 deg (0.081697).
  - 5.5 A: 30 deg (0.512412), 24 deg (0.191131), 36 deg (0.182822), 28 deg (0.171921), 32 deg (0.171332), 26 deg (0.099120), 34 deg (0.096543).
  - 7.0 A: 30 deg (0.391869), 24 deg (0.242704), 36 deg (0.239029), 28 deg (0.187494), 32 deg (0.179403), 26 deg (0.153836), 34 deg (0.141059).
  - 8.4 A: 24 deg (0.650807), 36 deg (0.630409), 30 deg (0.340482), 26 deg (0.190165), 32 deg (0.180789), 28 deg (0.180102), 34 deg (0.176869).

## Cautions

- These are comparative powder-style radial summaries, not direct fits to oriented/fiber experimental images.
- The corrected run uses the old pre-fix Track B report's 41 x 41 / 96-orientation / 160-bin settings for direct feature-window comparison.
- The input structures are compact rigid twist variants, not pNAB-generated conformer ensembles.
