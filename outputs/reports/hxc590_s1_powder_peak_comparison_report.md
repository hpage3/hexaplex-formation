# HXC590 S1 Powder Peak Comparison

## Purpose

This report compares John Bacsa's approximate HXC590 S1 powder peak list against existing simulated radial profiles for hexaflex / stacked-hexad candidate coordinate models. It is a peak-window diagnostic-spacing comparison, not a full Rietveld refinement or phase-refinement workflow.

## Experimental Context

Sample: `TM_TC without salt`, `HXC590 S1 powder`.

John reported that merged raw frames show uniform powder rings rather than oriented fiber arcs, with large broad background under peaks. Peak positions are approximate and relative intensities are treated as approximate rather than reliable constraints.

Uniform rings are expected for an unoriented powder sample, whereas oriented fibers give arcs. Therefore, the ring-versus-arc difference is not by itself evidence for a different molecular phase.

## Experimental Peaks and q Conversion

| d_A | q_Ainv | window_A | confidence |
|---|---|---|---|
| 7.90 | 0.795340 | +/-0.25 | lower |
| 7.30 | 0.860710 | +/-0.25 | lower |
| 6.50 | 0.966644 | +/-0.20 | lower |
| 5.50 | 1.142397 | +/-0.15 | lower |
| 4.33 | 1.451082 | +/-0.10 | medium_high |
| 3.90 | 1.611073 | +/-0.08 | medium_high |
| 3.71 | 1.693581 | +/-0.08 | medium_high |
| 3.35 | 1.875578 | +/-0.06 | medium_high |

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

Unavailable twist variants:

- Length/twist manifest rows for 24, 26, 28, 32, 34, and 36 degree full-length twists are marked pending locally and were not scored because their coordinate/profile files are not present.
- No local 15 or 26.7 degree coordinate/profile outputs were found in this branch.

## Scoring Method

Each experimental d-spacing was assigned a conservative d-window tolerance from the prompt. A candidate was counted as matching that peak when its radial profile contained at least one local profile maximum inside the d-window. The reported matched d-spacing is the closest local feature in that window.

Ranking uses peak-window coverage, diagnostic-window coverage, and spacing error. Relative intensity is not used in the ranking because the experimental background is broad and relative intensities were not treated as reliable constraints.

## Best-Matching Candidates

| candidate_id | candidate_label | model_path | profile_path | match_count | diagnostic_match_count | mean_abs_d_error_angstrom | mean_fractional_d_error | rank_score |
|---|---|---|---|---|---|---|---|---|
| central8_units_30deg | 8-unit formed endpoint, central 30-degree model | outputs\mini_hexaplex\structures\mini_hexaplex_central8_units.pdb | outputs\mini_hexaplex\radial_profiles\central8_units_radial.csv | 8 | 4 | 0.052881 | 0.009366 | 9.990634 |
| central12_units_30deg | 12-unit central 30-degree model | outputs\mini_hexaplex\structures\mini_hexaplex_central12_units.pdb | outputs\mini_hexaplex\radial_profiles\central12_units_radial.csv | 8 | 4 | 0.054269 | 0.009722 | 9.990278 |
| base_length_scale_0p90 | Base-length variant scale 0.90 | outputs\base_length_sweep\structures\hexaplex_base_length_scale_0p90.pdb | outputs\base_length_sweep\radial_profiles\hexaplex_base_length_scale_0p90_radial.csv | 7 | 4 | 0.049721 | 0.009426 | 8.990574 |
| base_length_scale_0p95 | Base-length variant scale 0.95 | outputs\base_length_sweep\structures\hexaplex_base_length_scale_0p95.pdb | outputs\base_length_sweep\radial_profiles\hexaplex_base_length_scale_0p95_radial.csv | 7 | 4 | 0.049721 | 0.009426 | 8.990574 |
| base_length_scale_1p00 | Base-length variant scale 1.00 | outputs\base_length_sweep\structures\hexaplex_base_length_scale_1p00.pdb | outputs\base_length_sweep\radial_profiles\hexaplex_base_length_scale_1p00_radial.csv | 7 | 4 | 0.049721 | 0.009426 | 8.990574 |
| base_length_scale_1p05 | Base-length variant scale 1.05 | outputs\base_length_sweep\structures\hexaplex_base_length_scale_1p05.pdb | outputs\base_length_sweep\radial_profiles\hexaplex_base_length_scale_1p05_radial.csv | 7 | 4 | 0.049721 | 0.009426 | 8.990574 |
| base_length_scale_1p10 | Base-length variant scale 1.10 | outputs\base_length_sweep\structures\hexaplex_base_length_scale_1p10.pdb | outputs\base_length_sweep\radial_profiles\hexaplex_base_length_scale_1p10_radial.csv | 7 | 4 | 0.049721 | 0.009426 | 8.990574 |
| base_length_scale_1p15 | Base-length variant scale 1.15 | outputs\base_length_sweep\structures\hexaplex_base_length_scale_1p15.pdb | outputs\base_length_sweep\radial_profiles\hexaplex_base_length_scale_1p15_radial.csv | 7 | 4 | 0.049721 | 0.009426 | 8.990574 |

## Diagnostic Peak-Window Match Table

| candidate_id | target_d_angstrom | matched | matched_d_angstrom | abs_d_error_angstrom | window_point_count |
|---|---|---|---|---|---|
| central8_units_30deg | 4.33 | yes | 4.310904 | 0.019096 | 1 |
| central8_units_30deg | 3.90 | yes | 3.910739 | 0.010739 | 1 |
| central8_units_30deg | 3.71 | yes | 3.693994 | 0.016006 | 1 |
| central8_units_30deg | 3.35 | yes | 3.328208 | 0.021792 | 1 |

## Conservative Same-Phase Interpretation

The HXC590 S1 powder peak list is consistent with a related stacked-hexad/hexaflex phase if the diagnostic 3.35 A, 4.33 A, and 3.7-3.9 A windows are reproduced by the simulated model profile.

For the current top-ranked candidate, `central8_units_30deg`, the diagnostic windows are reproduced as follows: 3.35 A = yes, 4.33 A = yes, 3.90 A = yes, and 3.71 A = yes.

A cautious answer to John's same-phase question is that the available powder peak list is compatible with a related stacked-hexad/hexaflex phase when these diagnostic spacing regions are matched, but the current comparison is not sufficient by itself to assign the phase.

## Limitations

This comparison does not establish a definitive phase assignment. The experimental background is broad, peak positions are approximate, and relative intensities were not treated as reliable constraints.

The model profiles are simplified radial-profile diagnostics. They are useful for comparing spacing windows, but they are not a full powder diffraction refinement, not a full fiber diffraction simulator, and not a calibrated experimental fit.

Improved comparison would benefit from calibrated q values, detector radii or pixel positions, beam center, wavelength, sample-to-detector distance, background subtraction details, and uncertainty estimates for each peak.

## Outputs

- `outputs/metrics/hxc590_s1_powder_peak_targets.csv`
- `outputs/metrics/hxc590_s1_powder_peak_match_scores.csv`
- `outputs\plots\hxc590_s1_powder_peak_comparison\peak_window_coverage.svg`
- `outputs\plots\hxc590_s1_powder_peak_comparison\best_candidate_d_errors.svg`
