# HXC590 S1 Powder Falsification Screen

## Purpose

This analysis is a falsification-style screening analysis, not a definitive phase assignment.

This remains a falsification-style screen, not a definitive phase assignment.

A candidate is more persuasive if it reproduces the diagnostic windows while plausible alternatives fail or match materially worse.

## Experimental input and limitations

The input peak list is John Bacsa / Nick's corrected HXC590 S1 powder target list for TM_TC without salt. The corrected distance scale shifts target positions modestly and places the base-stacking feature close to the expected 3.4 A region. Relative intensities were preserved in the corrected target table but were not treated as reliable constraints.

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
| nick_hexaplex_8hexads | nick_provided_8hexad | candidate | outputs\metrics\hxc590_s1_nick_hexaplex_8hexads_profile.csv |
| nick_16mer_antiparallel_30deg_ideal | nick_confirmed_16mer_simulation_input | candidate | outputs\metrics\hxc590_s1_nick_16mer_antiparallel_30deg_profile.csv |
| negative_hexads_only | wrong_geometry_control | negative_control | outputs\metrics\ladder_diffraction\profiles\reference_hexads_only_heavy_deduped_debye_profile.csv |
| negative_scaffold_only | wrong_geometry_control | negative_control | outputs\metrics\ladder_diffraction\profiles\reference_scaffold_only_complement_heavy_deduped_debye_profile.csv |

## Nick 8-Hexad Candidate

The `nick_hexaplex_8hexads` row is Nick's included `Hexaplex_8Hexads.xyz` candidate. It is treated as an 8-hexad candidate, not a 16-mer.

Tolerance survival for Nick's 8-hexad candidate:

| tolerance_setting | match_count | diagnostic_match_count | screen_survives | strict_survives | discrimination_score |
|---|---|---|---|---|---|
| narrow | 1 | 1 | no | no | 1.499130 |
| current | 3 | 2 | no | no | 3.987464 |
| broad | 3 | 2 | no | no | 3.987464 |

## Nick 16-Mer Simulation Input

The `nick_16mer_antiparallel_30deg_ideal` row is Nick's confirmed `Hexaplex_AntiParallel_30deg_Ideal.pdb` 16-mer simulation input.

Tolerance survival for Nick's 16-mer simulation input:

| tolerance_setting | match_count | diagnostic_match_count | screen_survives | strict_survives | discrimination_score |
|---|---|---|---|---|---|
| narrow | 1 | 0 | no | no | -1.374514 |
| current | 3 | 2 | no | no | 1.620709 |
| broad | 4 | 3 | no | no | 3.119145 |

## Candidate families unavailable

| candidate_id | family | reason |
|---|---|---|
| negative_alanine_beta_sheet | wrong_geometry_control | No existing radial profile was found; not generating a new control profile in this pass. |
| full_length_twist_24 | full_length_twist_variant | No coordinate/profile files were found for this twist. Existing pNAB helper support requires the current-model baseline YAML and pNAB runtime before generation is safe. |
| full_length_twist_26 | full_length_twist_variant | No coordinate/profile files were found for this twist. Existing pNAB helper support requires the current-model baseline YAML and pNAB runtime before generation is safe. |
| full_length_twist_28 | full_length_twist_variant | No coordinate/profile files were found for this twist. Existing pNAB helper support requires the current-model baseline YAML and pNAB runtime before generation is safe. |
| full_length_twist_32 | full_length_twist_variant | No coordinate/profile files were found for this twist. Existing pNAB helper support requires the current-model baseline YAML and pNAB runtime before generation is safe. |
| full_length_twist_34 | full_length_twist_variant | No coordinate/profile files were found for this twist. Existing pNAB helper support requires the current-model baseline YAML and pNAB runtime before generation is safe. |
| full_length_twist_36 | full_length_twist_variant | No coordinate/profile files were found for this twist. Existing pNAB helper support requires the current-model baseline YAML and pNAB runtime before generation is safe. |
| rise_3p20_synthetic_control | rise_variant | No coordinate/profile files were found for this rise. No safe existing rise-generation workflow or audited stack-axis transform was found in this repo. |
| rise_3p30_synthetic_control | rise_variant | No coordinate/profile files were found for this rise. No safe existing rise-generation workflow or audited stack-axis transform was found in this repo. |
| rise_3p35_synthetic_control | rise_variant | No coordinate/profile files were found for this rise. No safe existing rise-generation workflow or audited stack-axis transform was found in this repo. |
| rise_3p40_synthetic_control | rise_variant | No coordinate/profile files were found for this rise. No safe existing rise-generation workflow or audited stack-axis transform was found in this repo. |
| rise_3p50_synthetic_control | rise_variant | No coordinate/profile files were found for this rise. No safe existing rise-generation workflow or audited stack-axis transform was found in this repo. |
| rise_3p60_synthetic_control | rise_variant | No coordinate/profile files were found for this rise. No safe existing rise-generation workflow or audited stack-axis transform was found in this repo. |

## Scoring method

Observed-peak recovery uses the same local-profile-maximum d-window method as the HXC590 S1 peak comparison. This screen uses 5 target windows, including 3 diagnostic windows.

Predicted unmatched peaks were counted from the top local maxima within the experimental d-range. A predicted peak is unmatched when it does not fall inside any current experimental d-window. Unmatched predicted peaks are treated as screening diagnostics only, because weak or broad experimental features may be obscured by background.

Survival means matching all 3 diagnostic windows and at least 5 of 5 total peak windows. Strict survival means matching all 3 diagnostic windows and all 5 total windows.

## Observed-peak recovery

| candidate_id | match_count | diagnostic_match_count | screen_survives | strict_survives | mean_abs_d_error_angstrom | discrimination_score |
|---|---|---|---|---|---|---|
| central12_units_30deg | 5 | 3 | yes | yes | 0.061268 | 4.445075 |
| nick_hexaplex_8hexads | 3 | 2 | no | no | 0.064363 | 3.987464 |
| base_length_scale_0p85 | 4 | 3 | no | no | 0.061309 | 3.417331 |
| negative_scaffold_only | 3 | 2 | no | no | 0.095241 | 3.285702 |
| negative_hexads_only | 3 | 3 | no | no | 0.041907 | 3.204031 |
| central8_units_30deg | 4 | 2 | no | no | 0.066163 | 2.652840 |
| full_length_twist_30 | 3 | 2 | no | no | 0.067850 | 1.620709 |
| base_length_scale_0p90 | 3 | 2 | no | no | 0.067850 | 1.620709 |
| base_length_scale_0p95 | 3 | 2 | no | no | 0.067850 | 1.620709 |
| base_length_scale_1p00 | 3 | 2 | no | no | 0.067850 | 1.620709 |
| nick_16mer_antiparallel_30deg_ideal | 3 | 2 | no | no | 0.067850 | 1.620709 |
| central6_units_30deg | 3 | 2 | no | no | 0.067850 | 1.359346 |
| central7_units_30deg | 3 | 2 | no | no | 0.067850 | 1.359346 |
| base_length_scale_1p05 | 3 | 2 | no | no | 0.067850 | 1.359346 |
| base_length_scale_1p10 | 3 | 2 | no | no | 0.067850 | 1.359346 |
| base_length_scale_1p15 | 3 | 2 | no | no | 0.067850 | 1.359346 |

## Diagnostic-window recovery

`central8_units_30deg` recovered 2 of 3 diagnostic windows under current tolerances.

Current-tolerance survivors:
- `central12_units_30deg`

Current-tolerance failures:
- `nick_hexaplex_8hexads` missed survival criteria
- `base_length_scale_0p85` missed survival criteria
- `negative_scaffold_only` missed survival criteria
- `negative_hexads_only` missed survival criteria
- `central8_units_30deg` missed survival criteria
- `full_length_twist_30` missed survival criteria
- `base_length_scale_0p90` missed survival criteria
- `base_length_scale_0p95` missed survival criteria
- `base_length_scale_1p00` missed survival criteria
- `nick_16mer_antiparallel_30deg_ideal` missed survival criteria
- `central6_units_30deg` missed survival criteria
- `central7_units_30deg` missed survival criteria
- `base_length_scale_1p05` missed survival criteria
- `base_length_scale_1p10` missed survival criteria
- `base_length_scale_1p15` missed survival criteria

## Predicted unmatched peaks

| candidate_id | predicted_peak_rank | predicted_d_angstrom | relative_intensity |
|---|---|---|---|
| central6_units_30deg | 1 | 3.449745 | 1.000000 |
| central6_units_30deg | 2 | 3.693994 | 0.510103 |
| central6_units_30deg | 3 | 3.581075 | 0.403869 |
| central6_units_30deg | 5 | 6.626807 | 0.281742 |
| central6_units_30deg | 6 | 3.878156 | 0.277024 |
| central6_units_30deg | 7 | 4.613187 | 0.257817 |
| central6_units_30deg | 10 | 6.931946 | 0.222575 |
| central6_units_30deg | 11 | 4.193540 | 0.198542 |
| central6_units_30deg | 12 | 5.188157 | 0.160758 |
| central7_units_30deg | 1 | 3.449745 | 1.000000 |
| central7_units_30deg | 2 | 3.693994 | 0.407006 |
| central7_units_30deg | 3 | 3.581075 | 0.361197 |
| central7_units_30deg | 5 | 6.626807 | 0.227322 |
| central7_units_30deg | 6 | 6.931946 | 0.226204 |
| central7_units_30deg | 8 | 4.613187 | 0.200237 |
| central7_units_30deg | 10 | 3.878156 | 0.185330 |
| central7_units_30deg | 11 | 5.188157 | 0.173936 |
| central7_units_30deg | 12 | 4.046998 | 0.163657 |

## Tolerance / identifiability audit

| tolerance_setting | surviving_candidate_count | strict_surviving_candidate_count | best_candidate | central8_units_30deg_survives | central8_units_30deg_uniquely_best |
|---|---|---|---|---|---|
| narrow | 0 | 0 |  | no | no |
| current | 1 | 1 | central12_units_30deg | no | no |
| broad | 2 | 2 | central12_units_30deg | yes | no |

## Twist/rise sensitivity status

Missing full-length non-30-degree twist variants were not generated in this pass. The repo contains pNAB twist scaffolding, but no current-model baseline pNAB YAML and matching runtime inputs were found, and the local Python environment cannot currently import pNAB.

Rise variants were not generated in this pass. No safe existing rise-generation workflow or audited stack-axis transform was found for the current candidate model.

Synthetic twist/rise variants are controls for diffraction sensitivity, not chemically optimized structures.

If nearby twist or rise variants survive under current tolerances, the current powder peak list supports the conformation family but does not uniquely determine those parameters.

Because nearby twist and rise variants are unavailable here, this screen can compare the current conformation family against available length and negative-control alternatives, but it cannot uniquely determine twist or rise parameters.

Non-30-degree twist rows:

| candidate_id | parameter_value | status | reason |
|---|---|---|---|
| full_length_twist_24 | 24 deg | unavailable | No coordinate/profile files were found for this twist. The repo contains a pNAB twist helper, but no current-model baseline YAML with HelicalParameters.h_twist and CYP content was found. Also, python cannot import pnab. |
| full_length_twist_26 | 26 deg | unavailable | No coordinate/profile files were found for this twist. The repo contains a pNAB twist helper, but no current-model baseline YAML with HelicalParameters.h_twist and CYP content was found. Also, python cannot import pnab. |
| full_length_twist_28 | 28 deg | unavailable | No coordinate/profile files were found for this twist. The repo contains a pNAB twist helper, but no current-model baseline YAML with HelicalParameters.h_twist and CYP content was found. Also, python cannot import pnab. |
| full_length_twist_32 | 32 deg | unavailable | No coordinate/profile files were found for this twist. The repo contains a pNAB twist helper, but no current-model baseline YAML with HelicalParameters.h_twist and CYP content was found. Also, python cannot import pnab. |
| full_length_twist_34 | 34 deg | unavailable | No coordinate/profile files were found for this twist. The repo contains a pNAB twist helper, but no current-model baseline YAML with HelicalParameters.h_twist and CYP content was found. Also, python cannot import pnab. |
| full_length_twist_36 | 36 deg | unavailable | No coordinate/profile files were found for this twist. The repo contains a pNAB twist helper, but no current-model baseline YAML with HelicalParameters.h_twist and CYP content was found. Also, python cannot import pnab. |

Rise rows:

| candidate_id | parameter_value | status | reason |
|---|---|---|---|
| rise_3p20_synthetic_control | 3.20 A | unavailable | No coordinate/profile files were found for this rise. No safe existing rise-generation workflow or audited stack-axis transform was found in this repo. |
| rise_3p30_synthetic_control | 3.30 A | unavailable | No coordinate/profile files were found for this rise. No safe existing rise-generation workflow or audited stack-axis transform was found in this repo. |
| rise_3p35_synthetic_control | 3.35 A | unavailable | No coordinate/profile files were found for this rise. No safe existing rise-generation workflow or audited stack-axis transform was found in this repo. |
| rise_3p40_synthetic_control | 3.40 A | unavailable | No coordinate/profile files were found for this rise. No safe existing rise-generation workflow or audited stack-axis transform was found in this repo. |
| rise_3p50_synthetic_control | 3.50 A | unavailable | No coordinate/profile files were found for this rise. No safe existing rise-generation workflow or audited stack-axis transform was found in this repo. |
| rise_3p60_synthetic_control | 3.60 A | unavailable | No coordinate/profile files were found for this rise. No safe existing rise-generation workflow or audited stack-axis transform was found in this repo. |

## Interpretation

If many alternatives survive under current or broad tolerances, the current powder peak list is compatible with the proposed conformation but does not distinguish it uniquely.

Best current-tolerance candidate: `central12_units_30deg`. central8_units_30deg is not best under the current scoring.

The corrected current-tolerance screen ranks `central12_units_30deg` ahead of `central8_units_30deg` within the available screened profiles. This is a corrected-target rank change, not a unique refined phase assignment.

Nearby non-30-degree twists and requested rise variants neither survive nor fail in this pass because their coordinate/profile files are unavailable.

The current powder peak list helps screen the available conformation family against available alternatives, but it does not distinguish the exact twist or rise without generated nearby controls.

## What would be needed for stronger falsification

Stronger falsification would require calibrated q values, uncertainty estimates, background-subtracted experimental profiles, and generated profiles for full-length twist variants, rise variants, and chemically meaningful wrong-register compact controls.

## Limitations

This analysis is a falsification-style screening analysis, not a definitive phase assignment.

The experimental background is broad, peak positions are approximate, and relative intensities were not treated as reliable constraints.

Unmatched predicted peaks are treated as screening diagnostics only, because weak or broad experimental features may be obscured by background.

The available alternative set is incomplete: missing full-length twist variants and rise variants could not be generated safely from current local assets.

## Outputs

- `outputs/metrics/hxc590_s1_corrected_falsification_candidate_manifest.csv`
- `outputs/metrics/hxc590_s1_corrected_falsification_scores.csv`
- `outputs/metrics/hxc590_s1_corrected_predicted_unmatched_peaks.csv`
- `outputs/metrics/hxc590_s1_corrected_tolerance_survival_summary.csv`
- `outputs/metrics/hxc590_s1_twist_rise_sensitivity.csv` when the twist/rise audit has been run
- `outputs\plots\hxc590_s1_corrected_with_nick_16mer_falsification\candidate_discrimination_scores.svg`
- `outputs\plots\hxc590_s1_corrected_with_nick_16mer_falsification\diagnostic_window_survival.svg`
- `outputs\plots\hxc590_s1_corrected_with_nick_16mer_falsification\tolerance_survival_counts.svg`
- `outputs\plots\hxc590_s1_corrected_with_nick_16mer_falsification\unmatched_predicted_peaks.svg`
