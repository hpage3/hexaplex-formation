# Seed Starting-Ensemble Cost Audit

## Purpose

This audit compares synthetic loose-start endpoint classes against formed-perturbed endpoints in existing order-parameter space. It is intended to identify which cost groups dominate the distance from each loose-start class to the formed endpoint distribution.

This analysis is hypothesis-generating. It is not an atomistic Schrodinger bridge, not molecular dynamics, and not evidence of a physical nucleation pathway.

The `radially_separated` class is a dry-down/concentration-compaction proxy: it keeps rough angular phase and small axial offsets while moving chains farther apart radially. It does not model solvent evaporation, crystallization, molecular dynamics, or a full atomistic Schrodinger bridge.

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
| 6 | angular_randomized_loose_initial | 100.014094 | 14.078397 | 65.992123 | 78.356237 | geometric |
| 6 | loose_initial | 98.512946 | 13.564798 | 7.135016 | 70.079443 | geometric |
| 6 | radially_separated | 168.327884 | 13.831249 | 1.868881 | 119.230186 | geometric |
| 7 | angular_randomized_loose_initial | 126.097471 | 13.221978 | 34.150204 | 91.024999 | geometric |
| 7 | loose_initial | 111.093408 | 13.080670 | 8.119697 | 78.931179 | geometric |
| 7 | radially_separated | 201.243490 | 13.332773 | 2.222928 | 142.461037 | geometric |

## Ranking By Overall Cost

| unit_count | start_class | overall_cost | dominant_cost_group |
|---|---|---|---|
| 4 | loose_initial | 58.812113 | geometric |
| 5 | loose_initial | 60.246335 | geometric |
| 6 | loose_initial | 70.079443 | geometric |
| 5 | angular_randomized_loose_initial | 73.337982 | geometric |
| 6 | angular_randomized_loose_initial | 78.356237 | geometric |
| 7 | loose_initial | 78.931179 | geometric |
| 8 | loose_initial | 80.380937 | geometric |
| 8 | angular_randomized_loose_initial | 85.560038 | geometric |
| 7 | angular_randomized_loose_initial | 91.024999 | geometric |
| 4 | angular_randomized_loose_initial | 96.936553 | register_phase |
| 5 | radially_separated | 99.982092 | geometric |
| 4 | radially_separated | 105.824290 | geometric |
| 6 | radially_separated | 119.230186 | geometric |
| 7 | radially_separated | 142.461037 | geometric |
| 8 | radially_separated | 145.442953 | geometric |

## Cost Decomposition

- Unit 4 `loose_initial`: overall cost 58.812113; dominant group `geometric`.
- Unit 5 `loose_initial`: overall cost 60.246335; dominant group `geometric`.
- Unit 6 `loose_initial`: overall cost 70.079443; dominant group `geometric`.
- Unit 5 `angular_randomized_loose_initial`: overall cost 73.337982; dominant group `geometric`.
- Unit 6 `angular_randomized_loose_initial`: overall cost 78.356237; dominant group `geometric`.
- Unit 7 `loose_initial`: overall cost 78.931179; dominant group `geometric`.
- Unit 8 `loose_initial`: overall cost 80.380937; dominant group `geometric`.
- Unit 8 `angular_randomized_loose_initial`: overall cost 85.560038; dominant group `geometric`.
- Unit 7 `angular_randomized_loose_initial`: overall cost 91.024999; dominant group `geometric`.
- Unit 4 `angular_randomized_loose_initial`: overall cost 96.936553; dominant group `register_phase`.
- Unit 5 `radially_separated`: overall cost 99.982092; dominant group `geometric`.
- Unit 4 `radially_separated`: overall cost 105.824290; dominant group `geometric`.
- Unit 6 `radially_separated`: overall cost 119.230186; dominant group `geometric`.
- Unit 7 `radially_separated`: overall cost 142.461037; dominant group `geometric`.
- Unit 8 `radially_separated`: overall cost 145.442953; dominant group `geometric`.

## Radial-Separation Readout

- Unit 4 `radially_separated`: overall cost 105.824290, higher than `loose_initial` (58.812); dominant group `geometric`.
- Unit 5 `radially_separated`: overall cost 99.982092, higher than `loose_initial` (60.246); dominant group `geometric`.
- Unit 6 `radially_separated`: overall cost 119.230186, higher than `loose_initial` (70.079); dominant group `geometric`.
- Unit 7 `radially_separated`: overall cost 142.461037, higher than `loose_initial` (78.931); dominant group `geometric`.
- Unit 8 `radially_separated`: overall cost 145.442953, higher than `loose_initial` (80.381); dominant group `geometric`.

This readout asks whether a radial compaction proxy changes the cost relative to the older loose baseline. It should be interpreted as a geometric endpoint comparison, not as a simulation of the lab dry-down process.

## Conservative Interpretation

The cost values are order-parameter distance proxies. They indicate how far synthetic loose-start classes are from formed-perturbed endpoints in selected summary coordinates. They do not establish a physical mechanism or a validated assembly route.

A lower cost class is closer to the formed-perturbed endpoint distribution in this feature set. A dominant cost group suggests which type of coordinate difference is largest: geometry, contact recovery, or register/phase.

For the current 6-vs-7 focus, this audit should be read alongside the known endpoint construction: mini-hexaplex endpoints are coordinate-derived fragments cut from already-formed models, so formed geometry is present by construction.

## Recommended Next Step

Keep the start-class set limited to `loose_initial`, `angular_randomized_loose_initial`, and `radially_separated` until the radial-compaction readout is reviewed. Additional start classes should be added only one at a time with the same cost-decomposition checks.

## Plots

- `outputs\plots\seed_starting_ensemble_costs\seed_starting_ensemble_overall_costs.svg`
