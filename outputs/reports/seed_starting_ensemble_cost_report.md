# Seed Starting-Ensemble Cost Audit

## Purpose

This audit compares existing synthetic loose-start endpoint classes against formed-perturbed endpoints in existing order-parameter space. It is intended to identify which cost groups dominate the distance from each loose-start class to the formed endpoint distribution.

This analysis is hypothesis-generating. It is not an atomistic Schrodinger bridge, not molecular dynamics, and not evidence of a physical nucleation pathway.

## Inputs Used

- Order-parameter CSV: `outputs\metrics\seed_formation_order_parameters.csv`
- Reference class: `formed_perturbed`

## Cost Definition

For each unit count, the formed-perturbed rows define a reference distribution. For each available feature, the audit computes the formed-reference mean and standard deviation, compares the loose-start class mean against that formed mean, and records signed and absolute standardized deviations.

Group costs are root-mean-square absolute standardized deviations across available features. Lower cost means closer to the formed-perturbed endpoint distribution in the selected exploratory order parameters.

A small standard-deviation floor is used when the formed-reference standard deviation is zero or tiny, so constant synthetic reference features do not crash the audit.

## Feature Groups Used

- `geometric`: `RMSD_to_formed_seed_A`, `axial_extent_A`, `compactness_score`, `helical_axis_alignment_score`, `radial_extent_A`, `radius_of_gyration_A`
- `contact_recovery`: `CYP_MEP_contact_fraction_vs_target`, `backbone_contact_fraction_vs_target`, `contact_fraction_vs_target`
- `register_phase`: `angular_phase_order_score`, `axial_register_score`, `refined_angular_phase_score`

Skipped or unavailable feature notes:
- `register_phase`: refined_angular_phase_score: no loose_initial values for unit 4; refined_angular_phase_score: no loose_initial values for unit 5; refined_angular_phase_score: no loose_initial values for unit 6; refined_angular_phase_score: no loose_initial values for unit 7

## 6-vs-7 Focused Summary

| unit_count | start_class | geometric_cost | contact_recovery_cost | register_phase_cost | overall_cost | dominant_cost_group |
|---|---|---|---|---|---|---|
| 6 | angular_randomized_loose_initial | 102.142571 | 13.338226 | 55.382855 | 77.639195 | geometric |
| 6 | loose_initial | 94.012577 | 12.756948 | 9.969406 | 69.881188 | geometric |
| 7 | angular_randomized_loose_initial | 106.974284 | 14.327610 | 17.460214 | 76.480609 | geometric |
| 7 | loose_initial | 97.128347 | 14.216781 | 9.917234 | 72.241104 | geometric |

## Ranking By Overall Cost

| unit_count | start_class | overall_cost | dominant_cost_group |
|---|---|---|---|
| 4 | loose_initial | 61.627019 | geometric |
| 5 | loose_initial | 63.891927 | geometric |
| 6 | loose_initial | 69.881188 | geometric |
| 5 | angular_randomized_loose_initial | 71.871530 | geometric |
| 7 | loose_initial | 72.241104 | geometric |
| 7 | angular_randomized_loose_initial | 76.480609 | geometric |
| 6 | angular_randomized_loose_initial | 77.639195 | geometric |
| 4 | angular_randomized_loose_initial | 90.926871 | register_phase |

## Cost Decomposition

- Unit 4 `loose_initial`: overall cost 61.627019; dominant group `geometric`.
- Unit 5 `loose_initial`: overall cost 63.891927; dominant group `geometric`.
- Unit 6 `loose_initial`: overall cost 69.881188; dominant group `geometric`.
- Unit 5 `angular_randomized_loose_initial`: overall cost 71.871530; dominant group `geometric`.
- Unit 7 `loose_initial`: overall cost 72.241104; dominant group `geometric`.
- Unit 7 `angular_randomized_loose_initial`: overall cost 76.480609; dominant group `geometric`.
- Unit 6 `angular_randomized_loose_initial`: overall cost 77.639195; dominant group `geometric`.
- Unit 4 `angular_randomized_loose_initial`: overall cost 90.926871; dominant group `register_phase`.

## Conservative Interpretation

The cost values are order-parameter distance proxies. They indicate how far existing synthetic loose-start classes are from formed-perturbed endpoints in selected summary coordinates. They do not establish a physical mechanism or a validated assembly route.

A lower cost class is closer to the formed-perturbed endpoint distribution in this feature set. A dominant cost group suggests which type of coordinate difference is largest: geometry, contact recovery, or register/phase.

For the current 6-vs-7 focus, this audit should be read alongside the known endpoint construction: mini-hexaplex endpoints are coordinate-derived fragments cut from already-formed models, so formed geometry is present by construction.

## Recommended Next Step

After this audit is validated, add `axially_misregistered` as the next single start class and rerun the same cost-decomposition report before adding any other start classes.

## Plots

- `outputs\plots\seed_starting_ensemble_costs\seed_starting_ensemble_overall_costs.svg`
