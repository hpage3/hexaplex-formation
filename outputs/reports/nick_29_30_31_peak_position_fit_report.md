# Nick 29-31 Peak-Position Fit Report

## Inputs

Scoring command:

`.\.venv\Scripts\python.exe scripts\score_peak_position_fit.py --targets inputs\experimental_peak_windows\hxc590_rise_twist_targets.csv --assignments inputs\theoretical_peak_assignments\nick_example_29_30_31_peak_assignments.csv --summary-output outputs\metrics\nick_29_30_31_peak_position_fit_summary.csv --per-peak-output outputs\metrics\nick_29_30_31_peak_position_fit_per_peak_errors.csv`

Input files:

- `inputs/experimental_peak_windows/hxc590_rise_twist_targets.csv`
- `inputs/theoretical_peak_assignments/nick_example_29_30_31_peak_assignments.csv`

Experimental targets:

| target_label | target_d_A | target_group |
| --- | ---: | --- |
| base | 3.38 | base_stacking |
| A | 3.8 | backbone_associated |
| B | 4.4 | backbone_associated |
| C | 5.65 | backbone_associated |
| D | 7.3 | backbone_associated |

Nick measured theoretical peak assignments:

| model_id | twist_deg | base | A | B | C | D |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| nick_29 | 29 | 3.396 | 3.895 | 4.46 | 5.79 | 7.465 |
| nick_30 | 30 | 3.407 | 3.853 | 4.465 | 5.73 | 7.63 |
| nick_31 | 31 | 3.407 | 3.77 | 4.39 | 5.55 | 7.865 |

## Base-Stacking Score

| rank | model_id | root_sum_relative_error | rms_relative_error | weighted_rms_relative_error | missing_peak_count |
| ---: | --- | ---: | ---: | ---: | ---: |
| 1 | nick_29 | 0.004734 | 0.004734 | 0.004734 | 0 |
| 2 | nick_30 | 0.007988 | 0.007988 | 0.007988 | 0 |
| 3 | nick_31 | 0.007988 | 0.007988 | 0.007988 | 0 |

## Backbone-Associated Score

| rank | model_id | root_sum_relative_error | rms_relative_error | weighted_rms_relative_error | missing_peak_count |
| ---: | --- | ---: | ---: | ---: | ---: |
| 1 | nick_29 | 0.043998 | 0.021999 | 0.021999 | 0 |
| 2 | nick_30 | 0.051544 | 0.025772 | 0.025772 | 0 |
| 3 | nick_31 | 0.079819 | 0.039910 | 0.039910 | 0 |

## Combined Score

| rank | model_id | root_sum_relative_error | rms_relative_error | weighted_rms_relative_error | missing_peak_count |
| ---: | --- | ---: | ---: | ---: | ---: |
| 1 | nick_29 | 0.044252 | 0.019790 | 0.019790 | 0 |
| 2 | nick_30 | 0.052159 | 0.023326 | 0.023326 | 0 |
| 3 | nick_31 | 0.080218 | 0.035874 | 0.035874 | 0 |

## Interpretation

The calculated scores reproduce Nick's reported ordering and approximate values. By the backbone-associated A-D peak-position metric, 29 degrees ranks best, 30 degrees is close to 29 degrees, and 31 degrees is worse. The 31 degree model is especially penalized by the D peak: Nick's measured D = 7.865 A is farther from the experimental D = 7.3 A than the 29 and 30 degree D assignments.

The base-stacking score contributes little to twist discrimination here. The base peak mainly checks rise consistency; after rise is fixed, the A-D backbone-associated peaks are the more useful twist/pitch discriminator.

## Limitations

This is a peak-position comparison, not a full intensity fit. The assignments are Nick-measured assignments, not automated de novo assignments. Powder/radial data alone may restrict a twist range, but it does not establish a unique structure.
