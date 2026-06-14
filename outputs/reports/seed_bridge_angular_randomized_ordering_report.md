# Seed Bridge Feature-Ordering Analysis

## Purpose

This is a first-pass, SB-inspired order-parameter bridge analysis for the loose-to-formed mini-hexaplex seed transition. It estimates which order parameters become formed-like earliest when `angular_randomized_loose_initial` endpoints are paired to `formed_perturbed` endpoints by a low-cost entropic transport coupling in feature space.

This is not a full atomistic Schrödinger bridge. It does not include molecular dynamics, a force field, solvent, counterions, or physical transition probabilities.

## Inputs

- `outputs/metrics/seed_formation_order_parameters.csv`
- Existing `angular_randomized_loose_initial` and `formed_perturbed` ensembles summarized in that CSV.
- Contact-network outputs are used as prior context only, not as path constraints in this first bridge script.

## Endpoint Pairing

Endpoint pairs were built per unit count using methods: {'sinkhorn_greedy_unique': 400}. The primary method is Sinkhorn entropic optimal transport in standardized oriented-feature space, followed by unique pair extraction from highest coupling weights. Greedy minimum-cost matching is recorded if Sinkhorn fails.

## Feature Transformations

All selected features are oriented so larger means more formed-like:
- `compactness_score` from `compactness_score`: larger compactness_score is more formed-like.
- `contact_fraction_vs_target` from `contact_fraction_vs_target`: larger target-contact recovery is more formed-like.
- `CYP_MEP_contact_fraction_vs_target` from `CYP_MEP_contact_fraction_vs_target`: larger CYP/MEP-involving target-contact recovery is more formed-like.
- `axial_register_score` from `axial_register_score`: larger axial register score is more formed-like.
- `rmsd_formedness_score` from `RMSD_to_formed_seed_A`: RMSD converted to 1 / (1 + RMSD_to_formed_seed_A / 5).
- `seed_formation_score` from `seed_formation_score`: larger exploratory seed formation score is more formed-like.
- `refined_angular_phase_score` from `refined_angular_phase_score`: chain-label-aware angular agreement with the formed seed after optimizing one global rotation around the formed fitted axis.

Angular phase mode: `refined`. The legacy `angular_phase_order_score` is retained only as a comparison coordinate and should not be interpreted as a physical ordering claim. The refined coordinate measures chain-label-aware angular agreement to the formed seed after optimizing one global rotation around the formed fitted axis.

Coordinate-backed sample availability:
- unit 4: refined_angular_phase_score available for 203 coordinate-backed samples; 97 rows lacked saved ensemble PDBs and were excluded when refined phase was selected
- unit 5: refined_angular_phase_score available for 203 coordinate-backed samples; 97 rows lacked saved ensemble PDBs and were excluded when refined phase was selected
- unit 6: refined_angular_phase_score available for 203 coordinate-backed samples; 97 rows lacked saved ensemble PDBs and were excluded when refined phase was selected
- unit 7: refined_angular_phase_score available for 203 coordinate-backed samples; 97 rows lacked saved ensemble PDBs and were excluded when refined phase was selected

## Thresholds And Activation Times

For each feature, the formed-like threshold is `loose_mean + fraction * (formed_mean - loose_mean)`. Activation time is the first bridge time where the paired linear order-parameter path reaches that threshold. Threshold fractions evaluated: 0.50, 0.75, 0.80.

## Feature Ordering

| unit_count | first | second | third | confidence | summary |
|---:|---|---|---|---|---|
| 4 | contact_fraction_vs_target | rmsd_formedness_score | seed_formation_score | ambiguous | contact_fraction_vs_target -> rmsd_formedness_score -> seed_formation_score -> CYP_MEP_contact_fraction_vs_target -> axial_register_score |
| 5 | refined_angular_phase_score | rmsd_formedness_score | seed_formation_score | ambiguous | refined_angular_phase_score -> rmsd_formedness_score -> seed_formation_score -> contact_fraction_vs_target -> CYP_MEP_contact_fraction_vs_target |
| 6 | rmsd_formedness_score | CYP_MEP_contact_fraction_vs_target | axial_register_score | ambiguous | rmsd_formedness_score -> CYP_MEP_contact_fraction_vs_target -> axial_register_score -> compactness_score -> contact_fraction_vs_target |
| 7 | CYP_MEP_contact_fraction_vs_target | compactness_score | contact_fraction_vs_target | ambiguous | CYP_MEP_contact_fraction_vs_target -> compactness_score -> contact_fraction_vs_target -> refined_angular_phase_score -> rmsd_formedness_score |

Ordering table uses the primary threshold fraction `0.75`.

## Plots

- `outputs/seed_formation/plots/seed_bridge_mean_trajectories_unit4.png`
- `outputs/seed_formation/plots/seed_bridge_activation_times_unit4.png`
- `outputs/seed_formation/plots/seed_bridge_mean_trajectories_unit5.png`
- `outputs/seed_formation/plots/seed_bridge_activation_times_unit5.png`
- `outputs/seed_formation/plots/seed_bridge_mean_trajectories_unit6.png`
- `outputs/seed_formation/plots/seed_bridge_activation_times_unit6.png`
- `outputs/seed_formation/plots/seed_bridge_mean_trajectories_unit7.png`
- `outputs/seed_formation/plots/seed_bridge_activation_times_unit7.png`
- `outputs/seed_formation/plots/seed_bridge_feature_ordering_heatmap.png`
- `outputs/seed_formation/plots/seed_bridge_first_feature_frequency.png`
- `outputs/seed_formation/plots/seed_bridge_mean_trajectories_4_vs_7.png`

## Conservative Interpretation

Earliest activation in this workflow means earliest crossing along a paired linear interpolation in order-parameter space. It should be read as a candidate ordering for later SB collective-variable design, not as proof of physical time ordering.

If compactness or RMSD formedness activates early, closure-like coordinates may be good bridge progress variables. If CYP/MEP contact recovery activates early or frequently, base/contact specificity may need to be represented explicitly. Axial register and angular phase order are especially useful if they activate later or distinguish 7-unit behavior from 4-6 units, because they encode register and six-strand phasing beyond compactness.

For `angular_randomized_loose_initial`, refined angular phase should be interpreted as a test coordinate for deliberately broken chain-label angular placement. If it becomes early, angular phasing is informative under randomized loose conditions; if contacts, compactness, register, or RMSD formedness remain first, those coordinates are still the stronger bridge-ordering candidates.

The known contact-network result that 4 units is already six-chain connected means this bridge analysis should emphasize ordering among redundancy, contact recovery, register, and phase-order metrics rather than binary connectivity.

## Limitations

- Order-parameter interpolation only.
- No atomistic dynamics.
- No solvent or counterions.
- No energy model.
- Endpoint ensembles are synthetic.
- Activation thresholds are exploratory.
- Does not prove spontaneous formation mechanism.
