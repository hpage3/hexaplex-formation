# Exploratory Backbone Distance-Window Attribution

This analysis uses provisional distance windows chosen without detector calibration. The results should be interpreted as a sensitivity screen for Nick's backbone hypothesis, not as an assignment of the circled experimental features. Calibrated d-spacings, q values, detector radii, beam center, wavelength, and sample-to-detector distance would allow these windows to be replaced with measured experimental windows.

This is a geometry-only heavy-atom pair analysis. Hydrogens, exact duplicate atom records, and same-residue pairs are excluded. Directly bonded or near-bonded local pairs in different residues are not excluded because no bond/topology logic is used here.

## Inputs

- `full`: `outputs\intermediates\ai_candidate_inputs\full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain.pdb` (2166 heavy atoms)
- `central6`: `outputs\mini_hexaplex\structures\mini_hexaplex_central6_units.pdb` (864 heavy atoms)
- `central7`: `outputs\mini_hexaplex\structures\mini_hexaplex_central7_units.pdb` (1008 heavy atoms)

## Backbone-Involving Fractions

| Window | Type | central6 | central7 | full |
|---|---|---:|---:|---:|
| 3p0 | known | 1.000 | 1.000 | 1.000 |
| 3p4 | known | 0.455 | 0.463 | 0.458 |
| 4p5_5p0 | known | 0.387 | 0.388 | 0.387 |
| guessed_5p0_6p0 | exploratory_guessed | 0.499 | 0.499 | 0.501 |
| guessed_6p5_7p5 | exploratory_guessed | 0.524 | 0.516 | 0.513 |
| guessed_7p8_9p0 | exploratory_guessed | 0.536 | 0.534 | 0.535 |

## Dominant Atom-Pair Classes

| Window | Structure | Dominant atom-pair class | Pair count |
|---|---|---|---:|
| 3p0 | central6 | backbone_like_vs_base_like | 18 |
| 3p0 | central7 | backbone_like_vs_base_like | 21 |
| 3p0 | full | backbone_like_vs_base_like | 45 |
| 3p4 | central6 | base_like_vs_base_like | 90 |
| 3p4 | central7 | base_like_vs_base_like | 108 |
| 3p4 | full | base_like_vs_base_like | 252 |
| 4p5_5p0 | central6 | base_like_vs_base_like | 1374 |
| 4p5_5p0 | central7 | base_like_vs_base_like | 1689 |
| 4p5_5p0 | full | base_like_vs_base_like | 3873 |
| guessed_5p0_6p0 | central6 | base_like_vs_base_like | 2613 |
| guessed_5p0_6p0 | central7 | base_like_vs_base_like | 3273 |
| guessed_5p0_6p0 | full | base_like_vs_base_like | 7473 |
| guessed_6p5_7p5 | central6 | base_like_vs_base_like | 4350 |
| guessed_6p5_7p5 | central7 | base_like_vs_base_like | 5424 |
| guessed_6p5_7p5 | full | base_like_vs_base_like | 13119 |
| guessed_7p8_9p0 | central6 | base_like_vs_base_like | 5451 |
| guessed_7p8_9p0 | central7 | base_like_vs_base_like | 6705 |
| guessed_7p8_9p0 | full | base_like_vs_base_like | 16104 |

## Backbone-Hypothesis Screen

| Structure | Guessed window with highest backbone-involving fraction | Fraction |
|---|---|---:|
| central6 | guessed_7p8_9p0 | 0.536 |
| central7 | guessed_7p8_9p0 | 0.534 |
| full | guessed_7p8_9p0 | 0.535 |

## Interpretation

- Backbone-involving pairs are present across both known and guessed windows, but this screen does not assign any experimental feature.
- The guessed longer windows are more backbone-involving than the known 3.4 A and 4.5-5.0 A windows by fraction, except that the 3.0 A local window is entirely backbone-involving under this atom-name heuristic.
- Among the guessed windows, guessed_7p8_9p0 is most consistent with Nick's backbone hypothesis by backbone-involving fraction in central6, central7, and full.
- Even in the guessed windows, the single largest atom-pair class is base_like_vs_base_like, so the backbone hypothesis is supported only as an enrichment/sensitivity signal, not as a clean class assignment.
- The full model has backbone-involving fractions similar to central6/central7 in the guessed windows, suggesting scaling of geometric opportunities rather than a qualitatively new full-length backbone-only contribution.
- Experimental calibration would make this stronger by replacing guessed windows with measured d/q windows from detector geometry.

## Outputs

- Summary CSV: `outputs\metrics\backbone_distance_window_attribution_summary.csv`
- Atom-class CSV: `outputs\metrics\backbone_distance_window_attribution_by_atom_class.csv`
- Chain-relation CSV: `outputs\metrics\backbone_distance_window_attribution_by_chain_relation.csv`
- Pair sample CSV: `outputs\metrics\backbone_distance_window_attribution_pairs_sample.csv`
- Plot directory: `outputs\plots\backbone_distance_window_attribution`
- Plot: `outputs\plots\backbone_distance_window_attribution\3p0_atom_pair_class_by_structure.svg`
- Plot: `outputs\plots\backbone_distance_window_attribution\3p4_atom_pair_class_by_structure.svg`
- Plot: `outputs\plots\backbone_distance_window_attribution\4p5_5p0_atom_pair_class_by_structure.svg`
- Plot: `outputs\plots\backbone_distance_window_attribution\guessed_5p0_6p0_atom_pair_class_by_structure.svg`
- Plot: `outputs\plots\backbone_distance_window_attribution\guessed_6p5_7p5_atom_pair_class_by_structure.svg`
- Plot: `outputs\plots\backbone_distance_window_attribution\guessed_7p8_9p0_atom_pair_class_by_structure.svg`
- Plot: `outputs\plots\backbone_distance_window_attribution\guessed_5p0_6p0_backbone_involving_by_structure.svg`
- Plot: `outputs\plots\backbone_distance_window_attribution\guessed_5p0_6p0_chain_relation_by_structure.svg`
- Plot: `outputs\plots\backbone_distance_window_attribution\guessed_6p5_7p5_backbone_involving_by_structure.svg`
- Plot: `outputs\plots\backbone_distance_window_attribution\guessed_6p5_7p5_chain_relation_by_structure.svg`
- Plot: `outputs\plots\backbone_distance_window_attribution\guessed_7p8_9p0_backbone_involving_by_structure.svg`
- Plot: `outputs\plots\backbone_distance_window_attribution\guessed_7p8_9p0_chain_relation_by_structure.svg`
- Plot: `outputs\plots\backbone_distance_window_attribution\full_atom_pair_class_across_windows.svg`
