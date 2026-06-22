# Theoretical peak assignments

`nick_example_29_30_31_peak_assignments.csv` is intentionally header-only because the per-target A/B/C/D theoretical peak positions needed to reproduce the reported aggregate scores were not available in this repository. Do not infer or back-calculate those positions from the aggregate scores.

The expected assignment columns are:

`model_id,twist_deg,target_label,theoretical_d_A,assignment_method,notes`

Each model should provide one row per target label (`A`, `B`, `C`, `D`) when a peak can be assigned. Missing targets are scored with the configured missing-peak penalty by `scripts/score_peak_position_fit.py`.

Reported aggregate example scores to carry as context:

| twist_deg | root_sum_relative_error |
| ---: | ---: |
| 29 | 0.044 |
| 30 | 0.052 |
| 31 | 0.080 |
