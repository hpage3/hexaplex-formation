# Seed Starting-Ensemble Cost Audit

## Purpose

This audit compares synthetic loose-start endpoint classes against formed-perturbed endpoints in existing order-parameter space. It is intended to identify which cost groups dominate the distance from each loose-start class to the formed endpoint distribution.

This analysis is hypothesis-generating. It is not an atomistic Schrodinger bridge, not molecular dynamics, and not evidence of a physical nucleation pathway.

The `radially_separated` class is a dry-down/concentration-compaction proxy: it keeps rough angular phase and small axial offsets while moving chains farther apart radially. It does not model solvent evaporation, crystallization, molecular dynamics, or a full atomistic Schrodinger bridge.

The `axially_misregistered` class is a synthetic register-perturbation proxy. It starts from the loose baseline and applies symmetric chain-specific offsets along the fitted stack axis, asking whether axial register disruption gives a distinct cost signature from radial separation.

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

## 6-vs-7 Focused Summary

| unit_count | start_class | geometric_cost | contact_recovery_cost | register_phase_cost | overall_cost | dominant_cost_group |
|---|---|---|---|---|---|---|
| 6 | angular_randomized_loose_initial | 99.323359 | 12.780915 | 71.255052 | 79.010907 | geometric |
| 6 | axially_misregistered | 98.356960 | 12.671418 | 11.797784 | 70.085547 | geometric |
| 6 | loose_initial | 88.004264 | 12.150853 | 8.674412 | 62.674536 | geometric |
| 6 | radially_separated | 153.400504 | 12.631154 | 2.391708 | 108.660821 | geometric |
| 7 | angular_randomized_loose_initial | 118.147811 | 14.527616 | 26.983085 | 84.936666 | geometric |
| 7 | axially_misregistered | 120.668655 | 14.523890 | 11.314873 | 85.820770 | geometric |
| 7 | loose_initial | 102.988181 | 14.327182 | 8.039076 | 73.285445 | geometric |
| 7 | radially_separated | 184.855881 | 14.616580 | 2.025232 | 130.920911 | geometric |

## Ranking By Overall Cost

| unit_count | start_class | overall_cost | dominant_cost_group |
|---|---|---|---|
| 4 | loose_initial | 58.812113 | geometric |
| 6 | loose_initial | 62.674536 | geometric |
| 4 | axially_misregistered | 66.146135 | geometric |
| 6 | axially_misregistered | 70.085547 | geometric |
| 5 | loose_initial | 71.133251 | geometric |
| 7 | loose_initial | 73.285445 | geometric |
| 5 | axially_misregistered | 78.704382 | geometric |
| 6 | angular_randomized_loose_initial | 79.010907 | geometric |
| 5 | angular_randomized_loose_initial | 84.550127 | geometric |
| 7 | angular_randomized_loose_initial | 84.936666 | geometric |
| 7 | axially_misregistered | 85.820770 | geometric |
| 8 | loose_initial | 89.027799 | geometric |
| 4 | angular_randomized_loose_initial | 96.936553 | register_phase |
| 8 | angular_randomized_loose_initial | 97.405568 | geometric |
| 8 | axially_misregistered | 101.391323 | geometric |
| 4 | radially_separated | 105.824290 | geometric |
| 6 | radially_separated | 108.660821 | geometric |
| 5 | radially_separated | 122.111838 | geometric |
| 7 | radially_separated | 130.920911 | geometric |
| 8 | radially_separated | 168.150484 | geometric |

## Cost Decomposition

- Unit 4 `loose_initial`: overall cost 58.812113; dominant group `geometric`.
- Unit 6 `loose_initial`: overall cost 62.674536; dominant group `geometric`.
- Unit 4 `axially_misregistered`: overall cost 66.146135; dominant group `geometric`.
- Unit 6 `axially_misregistered`: overall cost 70.085547; dominant group `geometric`.
- Unit 5 `loose_initial`: overall cost 71.133251; dominant group `geometric`.
- Unit 7 `loose_initial`: overall cost 73.285445; dominant group `geometric`.
- Unit 5 `axially_misregistered`: overall cost 78.704382; dominant group `geometric`.
- Unit 6 `angular_randomized_loose_initial`: overall cost 79.010907; dominant group `geometric`.
- Unit 5 `angular_randomized_loose_initial`: overall cost 84.550127; dominant group `geometric`.
- Unit 7 `angular_randomized_loose_initial`: overall cost 84.936666; dominant group `geometric`.
- Unit 7 `axially_misregistered`: overall cost 85.820770; dominant group `geometric`.
- Unit 8 `loose_initial`: overall cost 89.027799; dominant group `geometric`.
- Unit 4 `angular_randomized_loose_initial`: overall cost 96.936553; dominant group `register_phase`.
- Unit 8 `angular_randomized_loose_initial`: overall cost 97.405568; dominant group `geometric`.
- Unit 8 `axially_misregistered`: overall cost 101.391323; dominant group `geometric`.
- Unit 4 `radially_separated`: overall cost 105.824290; dominant group `geometric`.
- Unit 6 `radially_separated`: overall cost 108.660821; dominant group `geometric`.
- Unit 5 `radially_separated`: overall cost 122.111838; dominant group `geometric`.
- Unit 7 `radially_separated`: overall cost 130.920911; dominant group `geometric`.
- Unit 8 `radially_separated`: overall cost 168.150484; dominant group `geometric`.

## Radial-Separation Readout

- Unit 4 `radially_separated`: overall cost 105.824290, higher than `loose_initial` (58.812); dominant group `geometric`.
- Unit 5 `radially_separated`: overall cost 122.111838, higher than `loose_initial` (71.133); dominant group `geometric`.
- Unit 6 `radially_separated`: overall cost 108.660821, higher than `loose_initial` (62.675); dominant group `geometric`.
- Unit 7 `radially_separated`: overall cost 130.920911, higher than `loose_initial` (73.285); dominant group `geometric`.
- Unit 8 `radially_separated`: overall cost 168.150484, higher than `loose_initial` (89.028); dominant group `geometric`.

This readout asks whether a radial compaction proxy changes the cost relative to the older loose baseline. It should be interpreted as a geometric endpoint comparison, not as a simulation of the lab dry-down process.

## Axial-Misregistration Readout

- Unit 6 `axially_misregistered`: overall cost 70.085547, higher than `loose_initial` (62.675); lower than `radially_separated` (108.661); lower than `angular_randomized_loose_initial` (79.011); dominant group `geometric`.
- Unit 7 `axially_misregistered`: overall cost 85.820770, higher than `loose_initial` (73.285); lower than `radially_separated` (130.921); higher than `angular_randomized_loose_initial` (84.937); dominant group `geometric`.

This readout asks whether stack-axis offsets create a distinct register or phase cost signature. It should be interpreted as a synthetic endpoint comparison, not as a physical trajectory.

## Conservative Interpretation

The cost values are order-parameter distance proxies. They indicate how far synthetic loose-start classes are from formed-perturbed endpoints in selected summary coordinates. They do not establish a physical mechanism or a validated assembly route.

A lower cost class is closer to the formed-perturbed endpoint distribution in this feature set. A dominant cost group suggests which type of coordinate difference is largest: geometry, contact recovery, or register/phase.

For the current 6-vs-7 focus, this audit should be read alongside the known endpoint construction: mini-hexaplex endpoints are coordinate-derived fragments cut from already-formed models, so formed geometry is present by construction.

## Recommended Next Step

Keep the start-class set limited to `loose_initial`, `angular_randomized_loose_initial`, `radially_separated`, and `axially_misregistered` until the axial-register readout is reviewed. Additional start classes should be added only one at a time with the same cost-decomposition checks.

## Plots

- `outputs\plots\seed_starting_ensemble_costs\seed_starting_ensemble_overall_costs.svg`
