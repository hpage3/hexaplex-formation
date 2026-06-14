# Distance-Window Attribution Analysis

This is a structural heavy-atom distance attribution analysis. It does not compute diffraction intensities directly. Pair counts identify enriched geometric contributors within X-ray-relevant distance windows in the molecular models.

## Inputs

- `full`: `outputs\intermediates\ai_candidate_inputs\full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain.pdb` (2166 heavy atoms; 3573 raw atom records)
- `central6`: `outputs\mini_hexaplex\structures\mini_hexaplex_central6_units.pdb` (864 heavy atoms; 1422 raw atom records)
- `central7`: `outputs\mini_hexaplex\structures\mini_hexaplex_central7_units.pdb` (1008 heavy atoms; 1659 raw atom records)

Selected full model: `outputs/intermediates/ai_candidate_inputs/full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain.pdb`, because it is the six-chain deduped full model used by prior workflows.

Hydrogens and exact duplicate atom records are excluded. Same-residue atom pairs are excluded. No bond/topology information is used, so directly bonded or near-bonded pairs in different residues are not removed; this is geometry-only.

Distance windows: 3p0 = 2.9-3.1 A, 3p4 = 3.3-3.5 A, 4p5_5p0 = 4.5-5.0 A.

## Dominant Residue-Pair Classes

| Window | Structure | Dominant class | Pair count |
|---|---|---|---:|
| 3p0 | central6 | CYP_MEP_vs_GLU | 18 |
| 3p0 | central7 | CYP_MEP_vs_GLU | 21 |
| 3p0 | full | CYP_MEP_vs_GLU | 45 |
| 3p4 | central6 | CYP_MEP_vs_CYP_MEP | 90 |
| 3p4 | central7 | CYP_MEP_vs_CYP_MEP | 108 |
| 3p4 | full | CYP_MEP_vs_CYP_MEP | 252 |
| 4p5_5p0 | central6 | CYP_MEP_vs_CYP_MEP | 1572 |
| 4p5_5p0 | central7 | CYP_MEP_vs_CYP_MEP | 1935 |
| 4p5_5p0 | full | CYP_MEP_vs_CYP_MEP | 4455 |

## Growth Across 6, 7, And Full Models

| Window | Structure | Total pairs | CYP/MEP-CYP/MEP | CYP/MEP-GLU | GLU-GLU | Dominant sequence relation |
|---|---|---:|---:|---:|---:|---|
| 3p0 | central6 | 18 | 0 | 18 | 0 | same_unit |
| 3p0 | central7 | 21 | 0 | 21 | 0 | same_unit |
| 3p0 | full | 45 | 0 | 45 | 0 | same_unit |
| 3p4 | central6 | 165 | 90 | 60 | 15 | adjacent_unit |
| 3p4 | central7 | 201 | 108 | 72 | 21 | adjacent_unit |
| 3p4 | full | 465 | 252 | 168 | 45 | adjacent_unit |
| 4p5_5p0 | central6 | 2286 | 1572 | 612 | 102 | adjacent_unit |
| 4p5_5p0 | central7 | 2816 | 1935 | 749 | 132 | adjacent_unit |
| 4p5_5p0 | full | 6447 | 4455 | 1692 | 300 | adjacent_unit |

## 4.5-5.0 A Chain Relation

| Structure | Same-chain pairs | Interchain pairs | Interpretation |
|---|---:|---:|---|
| central6 | 1248 | 1038 | mostly same-chain |
| central7 | 1490 | 1326 | mostly same-chain |
| full | 3438 | 3009 | mostly same-chain |

## Interpretation For Nick

- These pair counts are geometric contributors inside distance windows, not direct diffraction intensities and not proof that a class causes an experimental feature.
- In these models, the 3.0 A window is dominated by CYP/MEP-GLU pairs, while the 3.4 A window is dominated by CYP/MEP-CYP/MEP pairs.
- The 4.5-5.0 A window is dominated by CYP/MEP-CYP/MEP pairs in central6, central7, and full; CYP/MEP-GLU and GLU-GLU pairs are present but smaller contributors.
- The 4.5-5.0 A window is mixed but slightly same-chain enriched in all three structures under this geometry-only counting scheme.
- The 4.5-5.0 A count grows from central6 to central7 to full primarily because CYP/MEP-involving geometric opportunities grow with model length, with accompanying GLU/scaffold contributions.
- Growth from 6 to 7 to full indicates increasing geometric opportunities for those pair classes in the model, not necessarily proportional intensity growth.

## Outputs

- Summary CSV: `outputs\metrics\distance_window_attribution_summary.csv`
- Pair sample CSV: `outputs\metrics\distance_window_attribution_pairs_sample.csv`
- Chain-pair CSV: `outputs\metrics\distance_window_attribution_by_chain_pair.csv`
- Residue-class CSV: `outputs\metrics\distance_window_attribution_by_residue_class.csv`
- Plot directory: `outputs\plots\distance_window_attribution`
- Plot: `outputs\plots\distance_window_attribution\3p0_residue_class_contributions.svg`
- Plot: `outputs\plots\distance_window_attribution\3p4_residue_class_contributions.svg`
- Plot: `outputs\plots\distance_window_attribution\4p5_5p0_residue_class_contributions.svg`
- Plot: `outputs\plots\distance_window_attribution\4p5_5p0_chain_relation_contributions.svg`
- Plot: `outputs\plots\distance_window_attribution\4p5_5p0_chain_pair_contributions.svg`
