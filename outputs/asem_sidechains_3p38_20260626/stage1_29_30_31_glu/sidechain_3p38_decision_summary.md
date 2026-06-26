# 3.38 A Side-Chain Decision Summary

## What Was Tested

Complete Glu side-chain structures at 3.38 A rise for twists 29, 30, and 31 degrees. The screen used 50 existing heavy-atom XYZ files and existing corrected-diffraction/radial outputs. No 28/32 structures, half-degree twists, or new diffraction jobs were run during this comparison step.

## Validation Passed

- 50 corrected-diffraction NPY/PNG/JSON outputs were present.
- 50 radial CSV profiles were present.
- 50 scored manifest rows and 50 ranked rows were present.
- All scored rows were complete across the five target windows.

## Best Candidates

| twist | best model | rank | combined RMSD | helical RMSD | base RMSD |
| ---: | --- | ---: | ---: | ---: | ---: |
| 29 | angle29_cand6 | 10 | 0.088311 | 0.098398 | 0.016294 |
| 30 | angle30_cand9 | 1 | 0.061333 | 0.068086 | 0.016294 |
| 31 | angle31_cand19 | 2 | 0.061333 | 0.068086 | 0.016294 |

## Current Read

29 degrees is disfavored relative to 30 and 31 in the current 3.38 A radial-window score because its best combined/helical RMSD is worse than the tied 30/31 best rows. 30 and 31 remain unresolved because the best 30-degree row and multiple 31-degree rows share the same current combined RMSD and completeness.

## Expansion Decision

- 28/32 expansion is not the most urgent next step from this 29/30/31-only comparison. It could become useful if plot review suggests edge behavior remains important, but the current 3.38 screen points first to resolving the 30/31 ambiguity.
- A focused half-degree test around 30-31, especially 30.5 degrees, may be more informative than broad 28/32 expansion if the overlay plots support pursuing the 30/31 boundary.
- This remains a screening result, not final structural proof.

## Recommended Next Scientific Step

Review the generated radial overlays and 3.38-vs-3.40 comparison plots, then decide whether a focused 30.5-degree 3.38 A side-chain generation/scoring pass is warranted before broadening to 28/32.