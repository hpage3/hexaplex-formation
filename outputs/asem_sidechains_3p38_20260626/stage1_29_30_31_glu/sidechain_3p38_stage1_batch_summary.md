# Asem 3.38 A Side-Chain 29/30/31 Stage 1 Batch Summary

## Scope

This batch used only complete Glu side-chain structures at 29, 30, and 31 degrees from the committed Asem 3.38 A input set. No 28/32 structures and no half-degree twists were run.

## Inputs And Settings

- Input XYZ count: 50
- Rise: 3.38 A
- Corrected-diffraction engine: `scripts/run_asem_corrected_diffraction.py` with `reference/asem_corrected_diffraction_engine/`
- Grid size: 65
- Grid limit: 100 mm
- Theta count: 4
- Phi count: 8
- Psi count: 1
- Theta max: 180 deg
- Orientation count: 32
- Radial bins: 120
- Scoring targets: `inputs/experimental_peak_windows/hxc590_twist_rise_scan_targets.csv`

## Output Counts

- Corrected-diffraction NPY files: 50
- Corrected-diffraction PNG files: 50
- Corrected-diffraction JSON files: 50
- Radial CSV files: 50
- Observed peak rows: 250
- Wide observed peak rows: 50
- Scored manifest rows: 50
- Ranked manifest rows: 50

## Runtime

- Corrected-diffraction runtime sum from manifest: 76.046 seconds
- Mean per-model corrected-diffraction runtime: 1.521 seconds
- Min/max per-model runtime: 1.477 / 1.739 seconds

## Top Ranked Candidates

| rank | model_id | twist_deg | candidate_id | base_rmsd | helical_rmsd | combined_rmsd | completeness |
| ---: | --- | ---: | --- | ---: | ---: | ---: | ---: |
| 1 | angle30_cand9 | 30 | cand9 | 0.016294 | 0.068086 | 0.061333 | 1.000000 |
| 2 | angle31_cand19 | 31 | cand19 | 0.016294 | 0.068086 | 0.061333 | 1.000000 |
| 3 | angle31_cand21 | 31 | cand21 | 0.016294 | 0.068086 | 0.061333 | 1.000000 |
| 4 | angle31_cand23 | 31 | cand23 | 0.016294 | 0.068086 | 0.061333 | 1.000000 |
| 5 | angle31_cand24 | 31 | cand24 | 0.016294 | 0.068086 | 0.061333 | 1.000000 |
| 6 | angle31_cand25 | 31 | cand25 | 0.016294 | 0.068086 | 0.061333 | 1.000000 |
| 7 | angle31_cand16 | 31 | cand16 | 0.016294 | 0.086137 | 0.077387 | 1.000000 |
| 8 | angle31_cand18 | 31 | cand18 | 0.016294 | 0.086137 | 0.077387 | 1.000000 |
| 9 | angle31_cand22 | 31 | cand22 | 0.016294 | 0.086137 | 0.077387 | 1.000000 |
| 10 | angle29_cand6 | 29 | cand6 | 0.016294 | 0.098398 | 0.088311 | 1.000000 |

## Best Candidate Per Twist

| twist_deg | rank | model_id | candidate_id | base_rmsd | helical_rmsd | combined_rmsd | completeness |
| ---: | ---: | --- | --- | ---: | ---: | ---: | ---: |
| 29 | 10 | angle29_cand6 | cand6 | 0.016294 | 0.098398 | 0.088311 | 1.000000 |
| 30 | 1 | angle30_cand9 | cand9 | 0.016294 | 0.068086 | 0.061333 | 1.000000 |
| 31 | 2 | angle31_cand19 | cand19 | 0.016294 | 0.068086 | 0.061333 | 1.000000 |

## Cautious Interpretation

The 29/30/31 set does not fully resolve to a single twist at this stage. The best 30-degree row and several 31-degree rows have identical current scoring metrics, while the best 29-degree row is worse by combined RMSD. This supports treating 29 degrees as less favorable in this screen, but 30 vs 31 remains ambiguous under the present radial-window metric.

The 2D corrected-diffraction images are diagnostic outputs. The ranking above is based on the target-window radial profile score.

## Recommended Next Step

Compare these 3.38 A 29/30/31 Stage 1 scores against the prior 3.40 A Stage 2 refined results, then produce comparison plots before deciding whether expanding to 28/32 or testing half-degree twists is justified.
