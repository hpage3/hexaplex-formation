# HXC590 Peak-Position Fit Report

## Target peaks

The scoring target set uses four backbone-associated experimental d-spacings:

| target_label | experimental_d_A | default_weight | structural_note |
| --- | ---: | ---: | --- |
| D | 7.3 | 1.0 | longer-distance backbone-associated feature |
| C | 5.65 | 1.0 | intermediate backbone-associated feature |
| B | 4.4 | 1.0 | intermediate/short backbone-associated feature |
| A | 3.8 | 1.0 | short backbone-associated feature |

The target table is stored at `inputs/experimental_peak_windows/hxc590_backbone_peak_targets.csv`.

## Metric definitions

For each assigned target, the signed relative error is:

`(experimental_d_A - theoretical_d_A) / experimental_d_A`

The primary root-sum score is:

`sqrt(sum(relative_error^2))`

The RMS score is:

`sqrt(mean(relative_error^2))`

The weighted RMS score is:

`sqrt(sum(weight * relative_error^2) / sum(weight))`

The summary output also reports maximum absolute relative error, mean absolute relative error, and missing peak count. Missing target assignments receive the configured relative-error penalty, which defaults to 0.25.

## Example scores

The example input file at `inputs/theoretical_peak_assignments/nick_example_29_30_31_peak_assignments.csv` is header-only because the per-target theoretical A/B/C/D peak positions needed to reproduce the aggregate values were not available in this repository. Those peak positions should be added when the original assignments are available.

Reported aggregate example scores:

| twist_deg | root_sum_relative_error |
| ---: | ---: |
| 29 | 0.044 |
| 30 | 0.052 |
| 31 | 0.080 |

This pattern indicates that the 29 degree and 30 degree examples are close by this peak-position metric, while the 31 degree example is worse. The metric is useful for detecting poor peak-position fits above or below the preferred twist range, but it is not a full intensity fit.

## Ranked output

No ranked table is shown here because the included example assignment CSV has no per-target theoretical peak positions. Running `scripts/score_peak_position_fit.py` after adding assignments writes:

- `outputs/metrics/hxc590_peak_position_fit_summary.csv`
- `outputs/metrics/hxc590_peak_position_fit_per_peak_errors.csv`

The summary CSV is sorted best-to-worst by `weighted_rms_relative_error`.

## Limitations

Peak assignment choices matter. Intensity is not included in the score. Broad or overlapping peaks may shift with background subtraction. Full 2D fiber-pattern information is not captured by this one-dimensional peak-position metric.
