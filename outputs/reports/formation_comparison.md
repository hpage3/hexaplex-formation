# Formation Comparison Report

## Scientific caution

The d ~= 4.5 A diffraction feature is a reciprocal-space feature and should not be read as a literal 4.5 A atom contact. Contact maps and pair-distance summaries are real-space summaries. Hydrogen-bond outputs are rough geometric candidates, not definitive hydrogen bonds. Diffraction intensities from hexads-only and scaffold-only structures are comparative controls, not additive decompositions of the full intensity.

## Comparison table

| structure_base | residue_count | heavy_deduped_atom_count | contact_count_4p5A | motif_GLU_GLU | motif_GLU_any | hbond_candidate_count | GLU_involved_hbond_candidate_count | mean_radius_xy | z_span | angular_coverage_rad |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| alanine_beta_sheet_central5_res3_to_7 | 5 | 50 | 8 | 0 | 0 | 4 | 0 | 15.254623 | 0.035200 | 0.007353 |
| alanine_beta_sheet_central5_res4_to_8 | 5 | 50 | 8 | 0 | 0 | 4 | 0 | 15.254970 | 0.035200 | 0.007353 |
| alanine_beta_sheet_full | 10 | 100 | 21 | 0 | 0 | 10 | 0 | 15.254893 | 0.050600 | 0.010628 |
| full_hexaplex_anti_parallel_30deg_ideal | 180 | 2166 | 645 | 129 | 348 | 270 | 0 | 10.643939 | 52.434200 | 6.163350 |
| hexaplex_hexads_only | 90 | 810 | 258 | 0 | 0 | 0 | 0 | 5.565043 | 47.619111 | 6.279224 |
| hexaplex_scaffold_only_complement | 180 | 1356 | 513 | 129 | 348 | 0 | 0 | 12.097015 | 52.434200 | 6.079052 |

## Automatically generated notes

- The largest GLU-GLU motif count is in `full_hexaplex_anti_parallel_30deg_ideal` (129). This is consistent with a GLU-rich real-space contact pattern, but does not prove a specific stabilizing mechanism.
- The alanine controls have zero GLU motif counts, as expected for GLU-free controls.
- The largest angular coverage is in `hexaplex_hexads_only` (6.279224 rad), supporting the working hypothesis only as a coarse z-axis ordering screen.
- The Hexaplex scaffold-only complement has higher GLU-any and GLU-GLU counts than the alanine controls. This supports the working hypothesis that GLU-rich scaffold contacts distinguish the organized scaffold, but requires follow-up and does not prove formation pathway.
