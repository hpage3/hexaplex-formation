# Asem Twist Registry-Shift Hypothesis Report

## Purpose

Test whether the 31 degree Asem model is disfavored because over-twist shifts medium-range inter-hexad registry, while 29 and 30 degrees remain in the same radial-distance basin.

## Inputs

- `inputs/asem_twist_series_29_30_31/raw/29/29.pdb.txt`
- `inputs/asem_twist_series_29_30_31/raw/30/30.pdb.txt`
- `inputs/asem_twist_series_29_30_31/raw/31/31.pdb.txt`
- Existing corrected diffraction metrics: `outputs/metrics/asem_twist_series_feature_summary_corrected_emory_profile.csv`.

## Methods

No diffraction was rerun. The analysis uses existing feature metrics plus direct structural pair-distance audits of the imported PDBs. Layers are inferred as consecutive groups of six residues, giving 30 layers for 29/30 degree models and 32 layers for the 31 degree model. Atom classes are heuristic: backbone (`N`, `CA`, `C`, `O`), base core for CYP/MEP ring-like atoms, side-chain/carboxylate for GLU and carboxylate-like atoms, hydrogen, and unknown.

Distance bands: 3.2-3.6 A, 3.6-3.95 A, 4.15-4.65 A, 5.35-5.85 A, and 7.0-7.6 A. Pair-distance summaries focus on adjacent layer heavy-atom pairs.

## Feature Failure

31 degree is worse mainly where its absolute peak offset exceeds the 29/30 tied/best basin:

| feature | offset 29 | offset 30 | offset 31 | 31 penalty | drives? |
| --- | --- | --- | --- | --- | --- |
| 3.38/3.4 A | 0.0348347304828 | 0.0348347304828 | 0.0348347304828 | 0 | no |
| 3.77 A | -0.000283605539219 | -0.000283605539219 | -0.0863612112229 | 0.0860776056837 | yes |
| 4.4 A | -0.0582217172089 | -0.0582217172089 | -0.0582217172089 | 0 | no |
| 5.6 A | 0.00627265747814 | 0.00627265747814 | 0.00627265747814 | 0 | no |
| 7.3 A | 0.0181133952385 | 0.0181133952385 | 0.0181133952385 | 0 | no |

## Pair-Distance Findings

Flagged 31-specific shifted families: 12 supported and 5 contextual. Top supported examples:

| band | atom-class pair | metric | 29 | 30 | 31 | delta 31 vs 30 |
| --- | --- | --- | --- | --- | --- | --- |
| base_stack_3p2_3p6 | backbone|backbone | hist_peak_distance_A | 3.29 | 3.29 | 3.39 | 0.1 |
| base_stack_3p2_3p6 | backbone|side_chain_or_carboxylate | hist_peak_distance_A | 3.57 | 3.57 | 3.43 | -0.14 |
| feature_3p77_3p6_3p95 | backbone|side_chain_or_carboxylate | hist_peak_distance_A | 3.67 | 3.67 | 3.91 | 0.24 |
| feature_3p77_3p6_3p95 | backbone|side_chain_or_carboxylate | pair_count | 48 | 48 | 64 | 16 |
| feature_3p77_3p6_3p95 | base_core|side_chain_or_carboxylate | median_distance_A | 3.86600543197 | 3.8934063492 | 3.83083391704 | -0.06257243216 |
| feature_3p77_3p6_3p95 | base_core|side_chain_or_carboxylate | pair_count | 91 | 97 | 68 | -29 |
| feature_4p4_4p15_4p65 | backbone|side_chain_or_carboxylate | median_distance_A | 4.38005821879 | 4.37548122992 | 4.48825907207 | 0.11277784215 |
| feature_5p6_5p35_5p85 | backbone|side_chain_or_carboxylate | hist_peak_distance_A | 5.6 | 5.6 | 5.36 | -0.24 |

## Short-Contact Findings

Non-same-residue contacts below 2.2 A: 29 degree = 174, 30 degree = 174, 31 degree = 186. These follow the existing structural-warning convention by excluding same-residue covalent-neighbor pairs.

## Hypothesis Verdict

Verdict: `partially_supported`.

The feature deltas and pair-distance audit support a cautious interpretation that 31 degree moves one or more medium-range adjacent-layer pair-distance families out of the 29/30 basin. This is consistent with over-rotation disrupting inter-hexad registry, but it does not prove a literal leading/trailing end effect.

## Interpretation

The results should be interpreted as a proxy-level registry analysis. If supported or partially supported, they suggest asymmetric tolerance: under-twist to 29 degree is tolerated, but over-twist to 31 degree begins to degrade registry. This does not replace Asem/Nick visual chemistry review.

## Limitations

- Atom-class assignment is heuristic.
- Pair-distance counts are proxies for diffraction contributions.
- Powder/radial features aggregate many atom pairs.
- This is not an energy calculation.
- No AmberTools/tleap, minimization, notebooks, pNAB, diffraction rerun, or model generation was performed.
