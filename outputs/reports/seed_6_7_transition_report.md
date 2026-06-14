# Focused 6-vs-7 Seed Transition Analysis

This report reuses existing mini-hexaplex, seed-formation, contact-network, and bridge-ordering outputs. No ensembles or structures are regenerated.

The framing is cautious: this is an SB-inspired/order-parameter comparison, not a physical proof of nucleation and not an atomistic Schrodinger bridge. The defensible question is whether the existing 6- and 7-unit data show stronger, earlier, or more perturbation-resistant structural organization in the 7-unit seed than the 6-unit seed.

## Inputs

- Seed order parameters: `outputs\metrics\seed_formation_order_parameters.csv`
- Contact-network summary: `outputs\metrics\seed_contact_network_summary.csv`
- Bridge activation summary: `outputs\metrics\seed_bridge_activation_summary.csv`
- Angular-randomized bridge activation summary: `outputs\metrics\seed_bridge_activation_summary_angular_randomized.csv`
- Mini-hexaplex geometry/helicity summaries: `outputs\metrics\mini_hexaplex_geometry_summary.csv`, `outputs\metrics\mini_hexaplex_helicity_summary.csv`

## Primary Formed-Perturbed 7-vs-6 Comparison

| Feature | 6-unit mean | 7-unit mean | 7-minus-6 | Direction |
|---|---:|---:|---:|---|
| CYP_MEP_contact_fraction_vs_target | 0.807726 | 0.809685 | 0.001960 | 7-unit stronger |
| contact_fraction_vs_target | 0.805178 | 0.805585 | 0.000407 | 7-unit stronger |
| compactness_score | 0.997380 | 0.996952 | -0.000428 | 6-unit stronger |
| axial_register_score | 0.894824 | 0.885768 | -0.009056 | 6-unit stronger |
| RMSD_to_formed_seed_A | 0.618511 | 0.638805 | 0.020294 | 6-unit stronger |
| seed_formation_score | 0.909636 | 0.904958 | -0.004678 | 6-unit stronger |

## Contact-Network And Retention Summary

| Metric | 6-unit | 7-unit | 7-minus-6 | Direction |
|---|---:|---:|---:|---|
| CYP/MEP-involving contacts | 1605.000000 | 2067.000000 | 462.000000 | 7-unit stronger |
| total interchain contacts | 1740.000000 | 2256.000000 | 516.000000 | 7-unit stronger |
| contact retention under perturbation | 0.801724 | 0.784427 | -0.017297 | 6-unit stronger |

## Bridge Ordering

| Source | Feature | 6-unit median t | 7-unit median t | Direction |
|---|---|---:|---:|---|
| loose-initial bridge median activation time | CYP_MEP_contact_fraction_vs_target | 0.750000 | 0.800000 | 6-unit earlier |
| loose-initial bridge median activation time | contact_fraction_vs_target | 0.750000 | 0.800000 | 6-unit earlier |
| loose-initial bridge median activation time | compactness_score | 0.800000 | 0.800000 | similar |
| loose-initial bridge median activation time | axial_register_score | 0.750000 | 0.800000 | 6-unit earlier |
| loose-initial bridge median activation time | rmsd_formedness_score | 0.750000 | 0.800000 | 6-unit earlier |
| loose-initial bridge median activation time | seed_formation_score | 0.750000 | 0.800000 | 6-unit earlier |
| angular-randomized bridge median activation time | CYP_MEP_contact_fraction_vs_target | 0.775000 | 0.750000 | 7-unit earlier |
| angular-randomized bridge median activation time | contact_fraction_vs_target | 0.800000 | 0.750000 | 7-unit earlier |
| angular-randomized bridge median activation time | compactness_score | 0.800000 | 0.750000 | 7-unit earlier |
| angular-randomized bridge median activation time | axial_register_score | 0.800000 | 0.775000 | 7-unit earlier |
| angular-randomized bridge median activation time | rmsd_formedness_score | 0.750000 | 0.750000 | similar |
| angular-randomized bridge median activation time | seed_formation_score | 0.800000 | 0.750000 | 7-unit earlier |

## Interpretation

The primary formed-perturbed ensemble means are mixed. The existing data do not support a simple universal 7-unit strengthening across every feature.
The contact-network summary separates CYP/MEP-involving, GLU-involving, backbone-like, interchain, per-unit, and perturbation-retention measures where those columns are available. It does not provide a full CYP/MEP-CYP/MEP versus CYP/MEP-GLU pair-type decomposition, so that finer split should not be overclaimed here.
Activation times are order-parameter diagnostics along previously generated paired paths. Earlier activation should be read as earlier threshold crossing in this constructed feature space, not as a mechanistic assembly step.

## Output Files

- Summary CSV: `outputs\metrics\seed_6_7_transition_summary.csv`
- Feature comparison CSV: `outputs\metrics\seed_6_7_feature_comparison.csv`
- Plot directory: `outputs\plots\seed_6_7_transition`
- Plot: `outputs\plots\seed_6_7_transition\CYP_MEP_contact_fraction_vs_target_6_vs_7.svg`
- Plot: `outputs\plots\seed_6_7_transition\RMSD_to_formed_seed_A_6_vs_7.svg`
- Plot: `outputs\plots\seed_6_7_transition\axial_register_score_6_vs_7.svg`
- Plot: `outputs\plots\seed_6_7_transition\compactness_score_6_vs_7.svg`
- Plot: `outputs\plots\seed_6_7_transition\contact_fraction_vs_target_6_vs_7.svg`
- Plot: `outputs\plots\seed_6_7_transition\seed_formation_score_6_vs_7.svg`

## Limitations

- Existing ensembles have finite sample counts and are coordinate-derived perturbations.
- The analysis compares order parameters and static contact summaries; it does not establish a physical nucleation threshold.
- Contact classes are limited by the columns available in the prior contact-network CSV.
- Bridge ordering is SB-inspired feature interpolation/ordering, not atomistic Schrodinger-bridge dynamics.
