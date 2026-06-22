# Asem Twist Series Falsification Update

## Purpose

Apply the prior falsification-style framework to Asem's new 29, 30, and 31 degree twist-series models.

## Inputs

- Imported Asem source folders: `inputs/asem_twist_series_29_30_31/raw/29`, `raw/30`, `raw/31`.
- Corrected experimental profile: `inputs/experimental/nick_powder_profile_corrected_emory.csv`.
- Reused benchmark metrics from `outputs/metrics/asem_twist_series_*`.
- Reused current ideal baseline, parallel control, rise-sensitivity, atom-contribution, omega audit, and omega172 stop-context outputs.

## Methods

No diffraction was rerun. The script reuses existing corrected-diffraction CSVs produced with the Asem-corrected non-accumulating/vectorized path, `tilts = [0]`, `rotations = range(0, 181, 5)`, hydrogen exclusion, exact heavy-atom deduplication, and the current baseline grid/radial settings. Primary windows are 3.38/3.4 A, 3.77 A, 4.4 A, 5.6 A, and 7.3 A.

The combined ranking score is documented in the CSV as: mean absolute primary-window offset + 0.25 times max offset + structural warning penalty + H-bond plausibility penalty + role/context penalty. This is a transparent screen, not a fitted physical score.

## Main Findings

- 29 degree and 30 degree tie by the current mean primary-window peak offset.
- 31 degree is worse by the same metric and is comparatively disfavored.
- 30 degree remains supported, but it is not uniquely selected by this metric.
- The result supports a near-30 degree twist family rather than exact unique 30 degree from this run alone.
- All Asem candidates show plausible H-bond proxies, but all also carry short-contact flags.

## Falsification Summary

| model_label | control_role | mean_abs_primary_peak_offset_A | max_abs_primary_peak_offset_A | passes_primary_peak_position_filter | passes_multi_window_consistency_filter | structural_warning_flag | falsification_status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| twist_29_candidate_01_29_pdb_txt | new_asem_candidate | 0.0235452211895 | 0.0582217172089 | yes | yes | yes | survives_current_filters |
| twist_30_candidate_01_30_pdb_txt | new_asem_candidate | 0.0235452211895 | 0.0582217172089 | yes | yes | yes | survives_current_filters |
| twist_31_candidate_01_31_pdb_txt | new_asem_candidate | 0.0407607423262 | 0.0863612112229 | no | yes | yes | disfavored_by_peak_offsets |
| ideal_16mer_antiparallel_30deg | current_positive_baseline | 0.0235452211895 | 0.0582217172089 | yes | yes | no | survives_current_filters |
| parallel_sheet_control | parallel_control | 0.0501305087 | 0.107064464 | no | yes | no | disfavored_by_peak_offsets |
| rise_3.36_control | rise_variant_control | 0.0740174359674 | 0.126353636346 | no | no | no | disfavored_by_peak_offsets |
| rise_3.37_control | rise_variant_control | 0.0740174359674 | 0.126353636346 | no | no | no | disfavored_by_peak_offsets |
| rise_3.38_control | rise_variant_control | 0.0647769380998 | 0.126353636346 | no | no | no | disfavored_by_peak_offsets |
| rise_3.39_control | rise_variant_control | 0.0647769380998 | 0.126353636346 | no | no | no | disfavored_by_peak_offsets |
| rise_3.40_control | rise_variant_control | 0.0633904338691 | 0.137366801739 | no | no | no | disfavored_by_peak_offsets |
| atom_contribution_eight_hexads | atom_contribution_control | 0.0940565474156 | 0.16843775124 | no | no | no | disfavored_by_peak_offsets |
| atom_contribution_only_bases | atom_contribution_control | 0.106700976625 | 0.308365456918 | no | no | no | disfavored_by_peak_offsets |
| atom_contribution_with_coo | atom_contribution_control | 0.0452401261598 | 0.0871388102752 | yes | yes | no | survives_but_not_unique |
| omega172_structural_context | structural_only_not_diffraction_ready |  |  | no | no | no | context_only |
| asem_struct2_omega_series_context | structural_only_not_diffraction_ready |  |  | no | no | no | context_only |

## Ranking

| rank | model_label | control_role | combined_screen_score | mean_abs_primary_peak_offset_A | max_abs_primary_peak_offset_A | falsification_status |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | ideal_16mer_antiparallel_30deg | current_positive_baseline | 0.0431006504917 | 0.0235452211895 | 0.0582217172089 | survives_current_filters |
| 2 | twist_29_candidate_01_29_pdb_txt | new_asem_candidate | 0.0681006504917 | 0.0235452211895 | 0.0582217172089 | survives_current_filters |
| 3 | twist_30_candidate_01_30_pdb_txt | new_asem_candidate | 0.0681006504917 | 0.0235452211895 | 0.0582217172089 | survives_current_filters |
| 4 | twist_31_candidate_01_31_pdb_txt | new_asem_candidate | 0.0923510451319 | 0.0407607423262 | 0.0863612112229 | disfavored_by_peak_offsets |
| 5 | rise_3.38_control | rise_variant_control | 0.116365347186 | 0.0647769380998 | 0.126353636346 | disfavored_by_peak_offsets |
| 6 | rise_3.39_control | rise_variant_control | 0.116365347186 | 0.0647769380998 | 0.126353636346 | disfavored_by_peak_offsets |
| 7 | parallel_sheet_control | parallel_control | 0.1168966247 | 0.0501305087 | 0.107064464 | disfavored_by_peak_offsets |
| 8 | atom_contribution_with_coo | atom_contribution_control | 0.117024828729 | 0.0452401261598 | 0.0871388102752 | survives_but_not_unique |
| 9 | rise_3.40_control | rise_variant_control | 0.117732134304 | 0.0633904338691 | 0.137366801739 | disfavored_by_peak_offsets |
| 10 | rise_3.36_control | rise_variant_control | 0.125605845054 | 0.0740174359674 | 0.126353636346 | disfavored_by_peak_offsets |
| 11 | rise_3.37_control | rise_variant_control | 0.125605845054 | 0.0740174359674 | 0.126353636346 | disfavored_by_peak_offsets |
| 12 | atom_contribution_eight_hexads | atom_contribution_control | 0.186165985226 | 0.0940565474156 | 0.16843775124 | disfavored_by_peak_offsets |
| 13 | atom_contribution_only_bases | atom_contribution_control | 0.233792340855 | 0.106700976625 | 0.308365456918 | disfavored_by_peak_offsets |

## Interpretation

These results strengthen the case that the corrected experimental profile is selecting a narrow helical/twist family rather than arbitrary models. The current screen should not be read as proof of exact unique 30 degree twist because 29 and 30 tie on peak positions and only one full 30 degree candidate is available.

Prior control framing remains intact: the parallel control is comparatively disfavored by peak offsets, rise 3.38 did not beat rise 3.40 overall in the previous controlled test, and omega172 remains a local structural recommendation because no full diffraction-ready omega172 model exists yet.

## Recommended Next Step

Ask Asem to send several completed full 30 degree PDB candidates with carboxylates built. Then rerun this same falsification screen as an ensemble-stability test.

## Limitations

- No AmberTools/tleap generation was run; AmberTools/tleap is not available here.
- No minimization was run.
- No notebooks were executed.
- Powder/radial scores are screening metrics, not full structural refinement.
- H-bond and steric/contact metrics are proxies.
- The 29/30 tie may depend on the current radial-window metric and should not be overclaimed.

## Appendix: Short Nick/Asem Update Draft

Asem's provided 29, 30, and 31 degree full PDBs all benchmark cleanly through the corrected Emory-profile screen. The 29 and 30 degree models tie by the current mean primary-window peak-offset metric, while 31 degree is worse, so the data support a near-30 degree twist family rather than uniquely selecting exact 30 degree from this single-candidate run. All three show plausible backbone H-bond proxy distances, but they also retain short-contact flags, so I would treat them as screening candidates. The most decisive next test would be several independent completed 30 degree full PDB candidates with carboxylates built, then rerunning this same falsification screen for ensemble stability.
