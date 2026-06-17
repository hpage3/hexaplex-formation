# HXC590 S1 Powder Falsification Screen

## Purpose

This analysis is a falsification-style screening analysis, not a definitive phase assignment.

A candidate is more persuasive if it reproduces the diagnostic windows while plausible alternatives fail or match materially worse.

## Experimental input and limitations

The input peak list is John Bacsa's approximate HXC590 S1 powder list for TM_TC without salt. The experimental background is broad, peak positions are approximate, and relative intensities were not treated as reliable constraints.

Uniform rings are expected for an unoriented powder sample, whereas oriented fibers give arcs. The ring-versus-arc difference is not by itself evidence for a different molecular phase.

## Candidate families tested

| candidate_id | family | role | profile_path |
|---|---|---|---|
| central6_units_30deg | existing_hxc590_scored | candidate | outputs\mini_hexaplex\radial_profiles\central6_units_radial.csv |
| central7_units_30deg | existing_hxc590_scored | candidate | outputs\mini_hexaplex\radial_profiles\central7_units_radial.csv |
| central8_units_30deg | existing_hxc590_scored | candidate | outputs\mini_hexaplex\radial_profiles\central8_units_radial.csv |
| central12_units_30deg | existing_hxc590_scored | candidate | outputs\mini_hexaplex\radial_profiles\central12_units_radial.csv |
| full_length_twist_30 | full_length_twist_variant | candidate | outputs\mini_hexaplex\radial_profiles\full_length_baseline_radial.csv |
| base_length_scale_0p85 | base_length_variant | candidate | outputs\base_length_sweep\radial_profiles\hexaplex_base_length_scale_0p85_radial.csv |
| base_length_scale_0p90 | base_length_variant | candidate | outputs\base_length_sweep\radial_profiles\hexaplex_base_length_scale_0p90_radial.csv |
| base_length_scale_0p95 | base_length_variant | candidate | outputs\base_length_sweep\radial_profiles\hexaplex_base_length_scale_0p95_radial.csv |
| base_length_scale_1p00 | base_length_variant | candidate | outputs\base_length_sweep\radial_profiles\hexaplex_base_length_scale_1p00_radial.csv |
| base_length_scale_1p05 | base_length_variant | candidate | outputs\base_length_sweep\radial_profiles\hexaplex_base_length_scale_1p05_radial.csv |
| base_length_scale_1p10 | base_length_variant | candidate | outputs\base_length_sweep\radial_profiles\hexaplex_base_length_scale_1p10_radial.csv |
| base_length_scale_1p15 | base_length_variant | candidate | outputs\base_length_sweep\radial_profiles\hexaplex_base_length_scale_1p15_radial.csv |
| negative_hexads_only | wrong_geometry_control | negative_control | outputs\metrics\ladder_diffraction\profiles\reference_hexads_only_heavy_deduped_debye_profile.csv |
| negative_scaffold_only | wrong_geometry_control | negative_control | outputs\metrics\ladder_diffraction\profiles\reference_scaffold_only_complement_heavy_deduped_debye_profile.csv |

## Candidate families unavailable

| candidate_id | family | reason |
|---|---|---|
| negative_alanine_beta_sheet | wrong_geometry_control | No existing radial profile was found; not generating a new control profile in this pass. |
| full_length_twist_24 | full_length_twist_variant | Manifest marks this twist family as pending locally; coordinate/profile files are not present. |
| full_length_twist_26 | full_length_twist_variant | Manifest marks this twist family as pending locally; coordinate/profile files are not present. |
| full_length_twist_28 | full_length_twist_variant | Manifest marks this twist family as pending locally; coordinate/profile files are not present. |
| full_length_twist_32 | full_length_twist_variant | Manifest marks this twist family as pending locally; coordinate/profile files are not present. |
| full_length_twist_34 | full_length_twist_variant | Manifest marks this twist family as pending locally; coordinate/profile files are not present. |
| full_length_twist_36 | full_length_twist_variant | Manifest marks this twist family as pending locally; coordinate/profile files are not present. |
| rise_3p20_synthetic_control | rise_variant | No existing rise-variant generation/profile workflow was found; generation would require model-generation support. |
| rise_3p30_synthetic_control | rise_variant | No existing rise-variant generation/profile workflow was found; generation would require model-generation support. |
| rise_3p35_synthetic_control | rise_variant | No existing rise-variant generation/profile workflow was found; generation would require model-generation support. |
| rise_3p40_synthetic_control | rise_variant | No existing rise-variant generation/profile workflow was found; generation would require model-generation support. |
| rise_3p50_synthetic_control | rise_variant | No existing rise-variant generation/profile workflow was found; generation would require model-generation support. |
| rise_3p60_synthetic_control | rise_variant | No existing rise-variant generation/profile workflow was found; generation would require model-generation support. |

## Scoring method

Observed-peak recovery uses the same local-profile-maximum d-window method as the HXC590 S1 peak comparison. The diagnostic windows are 3.35 A, 4.33 A, 3.90 A, and 3.71 A.

Predicted unmatched peaks were counted from the top local maxima within the experimental d-range. A predicted peak is unmatched when it does not fall inside any current experimental d-window. Unmatched predicted peaks are treated as screening diagnostics only, because weak or broad experimental features may be obscured by background.

Survival means matching all 4 diagnostic windows and at least 6 of 8 total peak windows. Strict survival means matching all 4 diagnostic windows and all 8 total windows.

## Observed-peak recovery

| candidate_id | match_count | diagnostic_match_count | screen_survives | strict_survives | mean_abs_d_error_angstrom | discrimination_score |
|---|---|---|---|---|---|---|
| central8_units_30deg | 8 | 4 | yes | yes | 0.052881 | 8.823967 |
| central12_units_30deg | 8 | 4 | yes | yes | 0.054269 | 8.823611 |
| base_length_scale_0p85 | 7 | 4 | yes | no | 0.067834 | 7.821615 |
| central6_units_30deg | 7 | 4 | yes | no | 0.049721 | 7.532241 |
| central7_units_30deg | 7 | 4 | yes | no | 0.049721 | 7.532241 |
| full_length_twist_30 | 7 | 4 | yes | no | 0.049721 | 7.532241 |
| base_length_scale_0p90 | 7 | 4 | yes | no | 0.049721 | 7.532241 |
| base_length_scale_0p95 | 7 | 4 | yes | no | 0.049721 | 7.532241 |
| base_length_scale_1p00 | 7 | 4 | yes | no | 0.049721 | 7.532241 |
| base_length_scale_1p05 | 7 | 4 | yes | no | 0.049721 | 7.532241 |
| base_length_scale_1p10 | 7 | 4 | yes | no | 0.049721 | 7.532241 |
| base_length_scale_1p15 | 7 | 4 | yes | no | 0.049721 | 7.532241 |
| negative_hexads_only | 5 | 3 | no | no | 0.047204 | 5.847599 |
| negative_scaffold_only | 3 | 3 | no | no | 0.021327 | 3.794421 |

## Diagnostic-window recovery

`central8_units_30deg` recovered 4 of 4 diagnostic windows under current tolerances.

Current-tolerance survivors:
- `central8_units_30deg`
- `central12_units_30deg`
- `base_length_scale_0p85`
- `central6_units_30deg`
- `central7_units_30deg`
- `full_length_twist_30`
- `base_length_scale_0p90`
- `base_length_scale_0p95`
- `base_length_scale_1p00`
- `base_length_scale_1p05`
- `base_length_scale_1p10`
- `base_length_scale_1p15`

Current-tolerance failures:
- `negative_hexads_only` missed survival criteria
- `negative_scaffold_only` missed survival criteria

## Predicted unmatched peaks

| candidate_id | predicted_peak_rank | predicted_d_angstrom | relative_intensity |
|---|---|---|---|
| central6_units_30deg | 1 | 3.449745 | 1.000000 |
| central6_units_30deg | 3 | 3.581075 | 0.403869 |
| central6_units_30deg | 8 | 4.613187 | 0.257817 |
| central6_units_30deg | 11 | 6.931946 | 0.222575 |
| central6_units_30deg | 12 | 4.193540 | 0.198542 |
| central7_units_30deg | 1 | 3.449745 | 1.000000 |
| central7_units_30deg | 3 | 3.581075 | 0.361197 |
| central7_units_30deg | 7 | 6.931946 | 0.226204 |
| central7_units_30deg | 9 | 4.613187 | 0.200237 |
| central7_units_30deg | 12 | 5.188157 | 0.173936 |
| central8_units_30deg | 1 | 3.449745 | 1.000000 |
| central8_units_30deg | 4 | 3.581075 | 0.290006 |
| central8_units_30deg | 9 | 4.613187 | 0.211070 |
| central8_units_30deg | 12 | 4.046998 | 0.165612 |
| central12_units_30deg | 1 | 3.449745 | 1.000000 |
| central12_units_30deg | 9 | 4.613187 | 0.143182 |
| central12_units_30deg | 10 | 4.012012 | 0.139841 |
| central12_units_30deg | 11 | 3.581075 | 0.129852 |

## Tolerance / identifiability audit

| tolerance_setting | surviving_candidate_count | strict_surviving_candidate_count | best_candidate | central8_units_30deg_survives | central8_units_30deg_uniquely_best |
|---|---|---|---|---|---|
| narrow | 2 | 0 | central8_units_30deg | yes | yes |
| current | 12 | 2 | central8_units_30deg | yes | yes |
| broad | 12 | 4 | central8_units_30deg | yes | yes |

## Interpretation

Under the current tolerance setting, central8_units_30deg is the best-scoring available candidate in this screen. The result is more persuasive when considered against failed or weaker alternatives, but it remains a screening result rather than a phase assignment.

Best current-tolerance candidate: `central8_units_30deg`. central8_units_30deg is best under the current scoring.

## What would be needed for stronger falsification

Stronger falsification would require calibrated q values, uncertainty estimates, background-subtracted experimental profiles, and generated profiles for full-length twist variants, rise variants, and chemically meaningful wrong-register compact controls.

## Limitations

This analysis is a falsification-style screening analysis, not a definitive phase assignment.

The experimental background is broad, peak positions are approximate, and relative intensities were not treated as reliable constraints.

Unmatched predicted peaks are treated as screening diagnostics only, because weak or broad experimental features may be obscured by background.

The available alternative set is incomplete: missing full-length twist variants and rise variants could not be generated safely from current local assets.

## Outputs

- `outputs/metrics/hxc590_s1_falsification_candidate_manifest.csv`
- `outputs/metrics/hxc590_s1_falsification_scores.csv`
- `outputs/metrics/hxc590_s1_predicted_unmatched_peaks.csv`
- `outputs/metrics/hxc590_s1_tolerance_survival_summary.csv`
- `outputs\plots\hxc590_s1_falsification\candidate_discrimination_scores.svg`
- `outputs\plots\hxc590_s1_falsification\diagnostic_window_survival.svg`
- `outputs\plots\hxc590_s1_falsification\tolerance_survival_counts.svg`
- `outputs\plots\hxc590_s1_falsification\unmatched_predicted_peaks.svg`
