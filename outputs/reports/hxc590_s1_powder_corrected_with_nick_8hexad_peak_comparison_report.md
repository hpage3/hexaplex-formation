# HXC590 S1 Powder Peak Comparison

## Purpose

This report compares John Bacsa / Nick's corrected HXC590 S1 powder distance scale against existing simulated radial profiles for hexaflex / stacked-hexad candidate coordinate models. It is a peak-window diagnostic-spacing comparison, not a full Rietveld refinement or phase-refinement workflow.

## Experimental Context

Sample: `TM_TC without salt`, `HXC590 S1 powder`.

John Bacsa provided a corrected distance scale for the HXC590 S1 powder data. The correction shifts the target peak positions modestly and places the base-stacking feature close to the expected 3.4 A region. Relative intensities are preserved in the target table but are treated as approximate rather than reliable scoring constraints.

Uniform rings are expected for an unoriented powder sample, whereas oriented fibers give arcs. Therefore, the ring-versus-arc difference is not by itself evidence for a different molecular phase.

## Experimental Peaks and q Conversion

| d_A | q_Ainv | window_A | confidence |
|---|---|---|---|
| 7.33 | 0.857773 | +/-0.25 | lower |
| 5.58 | 1.126625 | +/-0.15 | lower |
| 4.40 | 1.427672 | +/-0.10 | medium_high |
| 3.79 | 1.655648 | +/-0.08 | medium_high |
| 3.38 | 1.858381 | +/-0.06 | medium_high |

The q conversion used `q = 2*pi/d`. Window matching was performed in d-space, with q-space bounds recorded in `outputs/metrics/hxc590_s1_powder_peak_targets.csv`.

## Candidate Models

| candidate_id | category | profile_path |
|---|---|---|
| central6_units_30deg | central_formed_endpoint | outputs\mini_hexaplex\radial_profiles\central6_units_radial.csv |
| central7_units_30deg | central_formed_endpoint | outputs\mini_hexaplex\radial_profiles\central7_units_radial.csv |
| central8_units_30deg | central_formed_endpoint | outputs\mini_hexaplex\radial_profiles\central8_units_radial.csv |
| central12_units_30deg | central_formed_endpoint | outputs\mini_hexaplex\radial_profiles\central12_units_radial.csv |
| full_length_twist_30 | full_length_30deg | outputs\mini_hexaplex\radial_profiles\full_length_baseline_radial.csv |
| base_length_scale_0p85 | base_length_variant | outputs\base_length_sweep\radial_profiles\hexaplex_base_length_scale_0p85_radial.csv |
| base_length_scale_0p90 | base_length_variant | outputs\base_length_sweep\radial_profiles\hexaplex_base_length_scale_0p90_radial.csv |
| base_length_scale_0p95 | base_length_variant | outputs\base_length_sweep\radial_profiles\hexaplex_base_length_scale_0p95_radial.csv |
| base_length_scale_1p00 | base_length_variant | outputs\base_length_sweep\radial_profiles\hexaplex_base_length_scale_1p00_radial.csv |
| base_length_scale_1p05 | base_length_variant | outputs\base_length_sweep\radial_profiles\hexaplex_base_length_scale_1p05_radial.csv |
| base_length_scale_1p10 | base_length_variant | outputs\base_length_sweep\radial_profiles\hexaplex_base_length_scale_1p10_radial.csv |
| base_length_scale_1p15 | base_length_variant | outputs\base_length_sweep\radial_profiles\hexaplex_base_length_scale_1p15_radial.csv |
| nick_hexaplex_8hexads | nick_provided_8hexad | outputs\metrics\hxc590_s1_nick_hexaplex_8hexads_profile.csv |

## Nick 8-Hexad Candidate

The `nick_hexaplex_8hexads` row is Nick's included `Hexaplex_8Hexads.xyz` candidate. It is treated as an 8-hexad candidate, not a 16-mer.
It scores with 3 of 5 corrected windows and 2 diagnostic windows matched.
By the existing rank score, Nick's 8-hexad candidate is below `central12_units_30deg` in this corrected screen.
By the existing rank score, Nick's 8-hexad candidate is below `central8_units_30deg` in this corrected screen.

Unavailable twist variants:

- Length/twist manifest rows for 24, 26, 28, 32, 34, and 36 degree full-length twists are marked pending locally and were not scored because their coordinate/profile files are not present.
- No local 15 or 26.7 degree coordinate/profile outputs were found in this branch.

## Scoring Method

Each experimental d-spacing was assigned a conservative d-window tolerance from the prompt. A candidate was counted as matching that peak when its radial profile contained at least one local profile maximum inside the d-window. The reported matched d-spacing is the closest local feature in that window.

Ranking uses peak-window coverage, diagnostic-window coverage, and spacing error. Relative intensity is not used in the ranking because the experimental background is broad and relative intensities were not treated as reliable constraints.

## Best-Matching Candidates

| candidate_id | candidate_label | model_path | profile_path | match_count | diagnostic_match_count | mean_abs_d_error_angstrom | mean_fractional_d_error | rank_score |
|---|---|---|---|---|---|---|---|---|
| central12_units_30deg | 12-unit central 30-degree model | outputs\mini_hexaplex\structures\mini_hexaplex_central12_units.pdb | outputs\mini_hexaplex\radial_profiles\central12_units_radial.csv | 5 | 3 | 0.061268 | 0.013258 | 6.486742 |
| base_length_scale_0p85 | Base-length variant scale 0.85 | outputs\base_length_sweep\structures\hexaplex_base_length_scale_0p85.pdb | outputs\base_length_sweep\radial_profiles\hexaplex_base_length_scale_0p85_radial.csv | 4 | 3 | 0.061309 | 0.014487 | 5.485513 |
| central8_units_30deg | 8-unit formed endpoint, central 30-degree model | outputs\mini_hexaplex\structures\mini_hexaplex_central8_units.pdb | outputs\mini_hexaplex\radial_profiles\central8_units_radial.csv | 4 | 2 | 0.066163 | 0.013826 | 4.986174 |
| nick_hexaplex_8hexads | Nick-provided Hexaplex_8Hexads.xyz 8-hexad candidate (not a 16-mer) | inputs\candidates\nick_hexaplex_8hexads.xyz | outputs\metrics\hxc590_s1_nick_hexaplex_8hexads_profile.csv | 3 | 2 | 0.064363 | 0.012536 | 3.987464 |
| base_length_scale_0p90 | Base-length variant scale 0.90 | outputs\base_length_sweep\structures\hexaplex_base_length_scale_0p90.pdb | outputs\base_length_sweep\radial_profiles\hexaplex_base_length_scale_0p90_radial.csv | 3 | 2 | 0.067850 | 0.015654 | 3.984346 |
| base_length_scale_0p95 | Base-length variant scale 0.95 | outputs\base_length_sweep\structures\hexaplex_base_length_scale_0p95.pdb | outputs\base_length_sweep\radial_profiles\hexaplex_base_length_scale_0p95_radial.csv | 3 | 2 | 0.067850 | 0.015654 | 3.984346 |
| base_length_scale_1p00 | Base-length variant scale 1.00 | outputs\base_length_sweep\structures\hexaplex_base_length_scale_1p00.pdb | outputs\base_length_sweep\radial_profiles\hexaplex_base_length_scale_1p00_radial.csv | 3 | 2 | 0.067850 | 0.015654 | 3.984346 |
| base_length_scale_1p05 | Base-length variant scale 1.05 | outputs\base_length_sweep\structures\hexaplex_base_length_scale_1p05.pdb | outputs\base_length_sweep\radial_profiles\hexaplex_base_length_scale_1p05_radial.csv | 3 | 2 | 0.067850 | 0.015654 | 3.984346 |

Compared with the earlier approximate-target run, the corrected five-target comparison changes the top-ranked available candidate from `central8_units_30deg` to `central12_units_30deg` under the existing scoring logic. This is a rank change within the available screened profiles, not a unique refined phase assignment.

The corrected near-3.4 A stacking target is included as a diagnostic window. `central8_units_30deg` remains scored in the corrected comparison with 4 of 5 total windows and 2 diagnostic windows matched.

## Diagnostic Peak-Window Match Table

| candidate_id | target_d_angstrom | matched | matched_d_angstrom | abs_d_error_angstrom | window_point_count |
|---|---|---|---|---|---|
| central12_units_30deg | 4.40 | yes | 4.310904 | 0.090096 | 1 |
| central12_units_30deg | 3.79 | yes | 3.753313 | 0.041687 | 1 |
| central12_units_30deg | 3.38 | yes | 3.328208 | 0.052792 | 1 |

## Conservative Same-Phase Interpretation

The corrected HXC590 S1 powder targets strengthen compatibility with a related stacked-hexad/hexaflex-like model family if the corrected diagnostic windows, including the near-3.4 A stacking feature, are reproduced by the simulated model profile.

For the current top-ranked candidate, `central12_units_30deg`, the diagnostic windows are reproduced as follows: 4.40 A = yes, 3.79 A = yes, 3.38 A = yes.

A cautious answer to John's same-phase question is that the available powder peak list is compatible with a related stacked-hexad/hexaflex-like model family when these diagnostic spacing regions are matched, but the current comparison is not sufficient by itself to assign the phase.

## Limitations

This comparison does not establish a definitive phase assignment. The experimental background is broad, peak positions are approximate, and relative intensities were not treated as reliable constraints.

The model profiles are simplified radial-profile diagnostics. They are useful for comparing spacing windows, but they are not a full powder diffraction refinement, not a full fiber diffraction simulator, and not a calibrated experimental fit.

Improved comparison would benefit from calibrated q values, detector radii or pixel positions, beam center, wavelength, sample-to-detector distance, background subtraction details, and uncertainty estimates for each peak.

## Outputs

- `outputs/metrics/hxc590_s1_powder_peak_targets.csv`
- `outputs/metrics/hxc590_s1_powder_peak_match_scores.csv`
- `outputs\plots\hxc590_s1_powder_corrected_with_nick_8hexad_peak_comparison\peak_window_coverage.svg`
- `outputs\plots\hxc590_s1_powder_corrected_with_nick_8hexad_peak_comparison\best_candidate_d_errors.svg`
