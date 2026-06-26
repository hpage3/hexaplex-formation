# 3.38 A vs 3.40 A Side-Chain Comparison Summary

This comparison uses the 3.38 A Stage 1 29/30/31 screen and the prior 3.40 A Stage 2 top-five-per-angle subset. Candidate IDs are not assumed to identify the same physical structures across rises.

## Candidate Counts

| twist_deg | 3.38 Stage 1 count | 3.40 source count | 3.40 Stage 2 subset count |
| ---: | ---: | ---: | ---: |
| 29 | 14 | 15 | 5 |
| 30 | 10 | 12 | 5 |
| 31 | 26 | 22 | 5 |

## Best Per Twist

| twist | 3.38 best | 3.38 combined | 3.40 best | 3.40 combined | delta 3.38-3.40 |
| ---: | --- | ---: | --- | ---: | ---: |
| 29 | angle29_cand6 | 0.088311 | angle29_cand11 | 0.101369 | -0.013058 |
| 30 | angle30_cand9 | 0.061333 | angle30_cand10 | 0.088311 | -0.026978 |
| 31 | angle31_cand19 | 0.061333 | angle31_cand10 | 0.088311 | -0.026978 |

## Top-10 Rank Distribution

- 3.38 top-10 twist counts: 29=1, 30=1, 31=8
- 3.40 top-10 among 29/30/31 twist counts: 29=4, 30=1, 31=5

## Observed Peak Position Notes

All 3.38 Stage 1 models scored complete across the five target windows: base, A, B, C, and D. The base-window observed peak is identical across top rows under the current binning, so separation is driven by backbone-associated windows. The observed peak error plot shows the spread by target and twist.

## Cautious Interpretation

- Relative to the prior 3.40 Stage 2 subset, the 3.38 best-per-twist combined RMSD is lower for 29, 30, and 31 in this scoring table.
- The comparison is not candidate-ID matched and should be interpreted as a screening comparison by twist/rise, not as a one-to-one structural provenance match.
- Within the 3.38 Stage 1 screen, 29 degrees remains weaker than the best 30/31 rows.
- 30 and 31 remain tied or near-tied under the current radial-window score; this does not prove 30 degrees and does not exclude 31 degrees.

## Plot Outputs

- `plots/sidechain_3p38_combined_rmsd_by_twist.png`
- `plots/sidechain_3p38_helical_rmsd_by_twist.png`
- `plots/sidechain_3p38_best_per_twist_bar.png`
- `plots/sidechain_3p38_vs_3p40_best_per_twist_combined_rmsd.png`
- `plots/sidechain_3p38_observed_peak_errors_by_target.png`
- `plots/sidechain_3p38_best_per_twist_radial_overlays.png`
- `plots/sidechain_3p38_best_per_twist_2d_contact_sheet.png`
