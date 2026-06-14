# Base-Length Diffraction Sweep Report

## Purpose

This controlled sweep asks whether changing local CYP/MEP base/hexad-arm length in a six-strand hexaplex model generates or strengthens a simulated diffraction feature near d ~= 4.5-5.0 A.

This is a computational sensitivity study only. It does not determine the structure.

## Scientific cautions

- The d ~= 4.5 A feature is reciprocal-space-like, not a literal atom-distance assignment.
- The unit convention is q = 2*pi/d.
- The simulation is powder-averaged and simplified relative to oriented/fiber-like experimental arcs.
- The variants are idealized and depend on the atom-selection and local-anchor transformation rule.
- CYP/MEP candidate arm atoms are axis-facing in the fitted-axis inspection; the transform scales local anchor-to-atom vectors, not global outward radial vectors.
- Raw experimental radial data are not included; comparisons are limited to q or d-spacing regions.

## Inputs and variants

- Baseline structure: outputs/intermediates/ai_candidate_inputs/full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain.pdb
- Variant manifest: outputs/base_length_sweep/structures/base_length_variant_manifest.csv
- Geometry sanity table: outputs/metrics/base_length_variant_geometry.csv
- Scale factors analyzed: 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15
- Transformed atom count(s): 1665
- Fixed atom count(s): 1908
- Scale 1.20 is excluded from default interpretation because geometry checks found suspicious heavy-atom overlaps below 1.00 A. It should be treated only as a stress-test variant if explicitly included.

Operationally, base/hexad-arm length is the local distance from each CYP/MEP residue anchor to selected non-backbone candidate arm atoms. GLU atoms and CYP/MEP backbone-like atoms remain fixed.

## Feature windows

- d_4p5_5p0: d = 4.5-5 A, q ~= 1.257-1.396 A^-1
- d_4p1: d = 3.95-4.25 A, q ~= 1.478-1.591 A^-1
- d_5p5: d = 5.35-5.75 A, q ~= 1.093-1.174 A^-1
- d_7p0: d = 6.75-7.25 A, q ~= 0.867-0.931 A^-1
- d_8p4: d = 8.1-8.7 A, q ~= 0.722-0.776 A^-1
- d_3p4: d = 3.25-3.55 A, q ~= 1.770-1.933 A^-1
- d_3p0: d = 2.9-3.1 A, q ~= 2.027-2.167 A^-1

## Summary table

| variant_id | scale_factor | geometry_warning | d_A_at_max_in_4p5_5A_window | integrated_intensity_4p5_5A | has_local_maximum_4p5_5A | integrated_intensity_3p4A | integrated_intensity_3p0A |
| --- | --- | --- | --- | --- | --- | --- | --- |
| hexaplex_base_length_scale_0p85 | 0.85 |  | 4.806469 | 24510728.941474 | yes | 54120810.839480 | 7850979.371604 |
| hexaplex_base_length_scale_0p90 | 0.90 |  | 4.806469 | 18966772.424768 | yes | 56381958.493204 | 8242381.094191 |
| hexaplex_base_length_scale_0p95 | 0.95 |  | 4.613187 | 14982802.828695 | no | 58483985.757161 | 8579717.420188 |
| hexaplex_base_length_scale_1p00 | 1.00 |  | 4.613187 | 12329206.192042 | no | 60652814.836791 | 8540328.003391 |
| hexaplex_base_length_scale_1p05 | 1.05 |  | 4.613187 | 10660883.743836 | no | 63412324.598186 | 8091219.878109 |
| hexaplex_base_length_scale_1p10 | 1.10 |  | 4.613187 | 9663438.904723 | no | 67115801.211582 | 7487561.946458 |
| hexaplex_base_length_scale_1p15 | 1.15 |  | 4.613187 | 9123194.391505 | no | 71756525.712569 | 7280454.905806 |

## Plots

- outputs/base_length_sweep/plots/base_length_sweep_q_profile_overlay.png
- outputs/base_length_sweep/plots/base_length_sweep_d_profile_overlay.png
- outputs/base_length_sweep/plots/base_length_sweep_d_4p1_8p4_zoom.png
- outputs/base_length_sweep/plots/base_length_sweep_d_4p5_5p0_zoom.png
- outputs/base_length_sweep/plots/base_length_scale_vs_integrated_4p5_5p0.png
- outputs/base_length_sweep/plots/base_length_scale_vs_peak_d_4p5_5p0.png

## Conservative interpretation

- 4.5-5.0 A local maximum call: at least one analyzed variant has a conservative local maximum.
- The integrated 4.5-5.0 A-window intensity decreases across the analyzed scale range.
- The 3.4 A reference window has nonzero intensity in all analyzed variants.
- The 3.0 A reference window has nonzero intensity in all analyzed variants.

In this simplified comparative dataset, a change in the 4.5-5.0 A-window score with base-length scale is consistent with sensitivity to local CYP/MEP arm geometry. It does not establish causality or structural identity.

## Geometry warnings

- No geometry warnings were reported for the analyzed variants.

## Limitations

- Powder-averaged simulation is not equivalent to Emory oriented/fiber-like arcs.
- Structural variants are idealized local-anchor perturbations.
- Results depend on the selected CYP/MEP transformable atoms.
- Raw experimental radial profiles are unavailable, so this is not a direct fit.
- Scale 1.20 has a known overlap warning and is not part of the default main interpretation.
