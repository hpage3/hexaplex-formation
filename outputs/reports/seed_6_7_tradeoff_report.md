# 6-vs-7 Seed Tradeoff Analysis

This follow-up reuses existing CSV outputs only. It does not regenerate ensembles, mini-hexaplex structures, or bridge paths.

The goal is explanatory rather than mechanistic: compare absolute contact-network growth, normalized formedness metrics, and perturbation-retention metrics in coordinate-derived 6- and 7-unit mini-seeds. This report does not claim a physical nucleation threshold or a true assembly mechanism.

## Category Summary

| Reasoning category | Metrics favoring 6 | Metrics favoring 7 | Net 7-minus-6 | Overall classification |
|---|---:|---:|---:|---|
| absolute_network_growth | 0 | 8 | 8 | favors_7 |
| geometric_register | 5 | 3 | -2 | favors_6 |
| normalized_fraction | 5 | 9 | 4 | favors_7 |
| perturbation_retention | 2 | 0 | -2 | favors_6 |
| possible_edge_effect | 1 | 4 | 3 | favors_7 |

## Absolute Contact-Network Growth

| Metric | 6-unit | 7-unit | 7-minus-6 | Classification |
|---|---:|---:|---:|---|
| total interchain contacts | 1740.000000 | 2256.000000 | 516.000000 | favors_7 |
| CYP/MEP-involving contacts | 1605.000000 | 2067.000000 | 462.000000 | favors_7 |
| GLU-involving contacts | 150.000000 | 210.000000 | 60.000000 | favors_7 |
| backbone-like contacts | 132.000000 | 180.000000 | 48.000000 | favors_7 |
| chain-summed total contacts | 1740.000000 | 2256.000000 | 516.000000 | favors_7 |
| CYP/MEP vs CYP/MEP contacts | 1590.000000 | 2046.000000 | 456.000000 | favors_7 |
| CYP/MEP vs GLU contacts | 15.000000 | 21.000000 | 6.000000 | favors_7 |
| GLU vs GLU contacts | 135.000000 | 189.000000 | 54.000000 | favors_7 |

## Normalized, Retention, And Register Metrics

| Metric | 6-unit | 7-unit | 7-minus-6 | Classification | Category |
|---|---:|---:|---:|---|---|
| perturbation contact retention | 0.801724 | 0.784427 | -0.017297 | favors_6 | perturbation_retention |
| perturbation contact retention variability | 0.023336 | 0.039299 | 0.015963 | favors_6 | perturbation_retention |
| formed_perturbed compactness | 0.997380 | 0.996952 | -0.000428 | mixed_or_negligible | geometric_register |
| formed_perturbed axial register | 0.894824 | 0.885768 | -0.009056 | favors_6 | geometric_register |
| formed_perturbed RMSD to formed seed | 0.618511 | 0.638805 | 0.020294 | favors_6 | geometric_register |
| loose_initial compactness | 0.575031 | 0.600787 | 0.025756 | favors_7 | geometric_register |
| loose_initial axial register | 0.437444 | 0.427521 | -0.009923 | favors_6 | geometric_register |
| loose_initial RMSD to formed seed | 14.793849 | 15.105454 | 0.311605 | favors_6 | geometric_register |
| angular_randomized_loose_initial compactness | 0.597233 | 0.620286 | 0.023053 | favors_7 | geometric_register |
| angular_randomized_loose_initial axial register | 0.428850 | 0.440355 | 0.011505 | favors_7 | geometric_register |
| angular_randomized_loose_initial RMSD to formed seed | 19.035960 | 19.521899 | 0.485939 | favors_6 | geometric_register |

## Edge-Effect Proxies

| Metric | 6-unit | 7-unit | 7-minus-6 | Classification |
|---|---:|---:|---:|---|
| contacts per unit | 290.000000 | 322.285714 | 32.285714 | favors_7 |
| contacts per chain pair mean | 116.000000 | 150.400000 | 34.400000 | favors_7 |
| contacts per chain pair max | 304.000000 | 442.000000 | 138.000000 | favors_7 |
| contacts involving terminal/edge units | 708.000000 | 708.000000 | 0.000000 | mixed_or_negligible |
| contacts involving only internal units | 1032.000000 | 1548.000000 | 516.000000 | favors_7 |
| fraction of unit-edge contacts involving terminal/edge units | 0.406897 | 0.313830 | -0.093067 | favors_6 |

## Interpretation

The 7-unit seed increases contact-network opportunity and redundancy: raw total, CYP/MEP-involving, and chemistry-resolved chain-edge contact counts all rise in the existing static contact network.
The 6-unit seed can still score slightly better on normalized formedness, geometric register, or endpoint-retention metrics. Those quantities are fractions, fitted scores, or retention ratios rather than absolute contact opportunities.
This pattern supports a transition/tradeoff regime rather than a single best seed length. Added length creates more absolute contacts, but it also adds more degrees of freedom and edge-sensitive contacts that can reduce normalized retention or register scores.
The unit-level edge proxy can flag whether terminal/edge contacts change between 6 and 7 units, but it does not prove that edge effects cause the normalized-score differences.

## Outputs

- Summary CSV: `outputs\metrics\seed_6_7_tradeoff_summary.csv`
- Feature table CSV: `outputs\metrics\seed_6_7_tradeoff_feature_table.csv`
- Plot directory: `outputs\plots\seed_6_7_tradeoff`
- Plot: `outputs\plots\seed_6_7_tradeoff\seed_6_7_tradeoff_classification_counts.svg`

## Limitations

- The contact-type split is available from `seed_contact_network_edges.csv`; any absent pair classes are reported only where existing columns support them.
- Endpoint retention is based on existing perturbation summaries and bridge diagnostics, not new sampling.
- The analysis is order-parameter and coordinate-derived; it should not be read as an atomistic mechanism.
