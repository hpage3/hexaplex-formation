# Ladder Diffraction Approximation Report

## Scientific cautions

This report uses a simplified Debye-style isotropic scattering approximation for comparative radial profiles. It is not a replacement for full fiber-diffraction simulations and does not model orientation, lattice effects, form factors, hydration, or experimental instrument response.

The d ~= 4.5 A feature remains a reciprocal-space scaffold signature in the working hypothesis, not a literal 4.5 A atom-contact assignment. Component structure intensities are comparative controls, not additive decompositions of the full intensity.

The candidate block map remains unvalidated against PyMOL colored strand paths, and ladder models do not prove temporal assembly order.

## Runtime metadata

- Profile methods: histogram.
- Distance bin widths: 0.05.
- 15 of 15 profile(s) used full post-filter atom sets.

## Window fraction table

| model_name | included_blocks | includes_hexads | d_8p4_fraction | d_5p5_6p0_fraction | d_4p5_fraction | d_4p1_fraction | d_3p4_fraction | d_3p0_fraction |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| hexads_plus_scaffold_blocks_1_2_3_4_5_6_heavy_deduped | 1,2,3,4,5,6 | yes | 0.018272 | 0.082189 | 0.021626 | 0.024174 | 0.077104 | 0.020522 |
| hexads_plus_scaffold_blocks_1_2_3_4_5_heavy_deduped | 1,2,3,4,5 | yes | 0.019885 | 0.060936 | 0.018306 | 0.020706 | 0.072219 | 0.017077 |
| hexads_plus_scaffold_blocks_1_2_3_4_heavy_deduped | 1,2,3,4 | yes | 0.016292 | 0.040289 | 0.014856 | 0.017679 | 0.064457 | 0.013320 |
| hexads_plus_scaffold_blocks_1_2_3_heavy_deduped | 1,2,3 | yes | 0.014722 | 0.028167 | 0.011554 | 0.014713 | 0.061876 | 0.011453 |
| hexads_plus_scaffold_blocks_1_2_heavy_deduped | 1,2 | yes | 0.014623 | 0.016514 | 0.009785 | 0.013054 | 0.059538 | 0.009314 |
| hexads_plus_scaffold_blocks_1_heavy_deduped | 1 | yes | 0.016600 | 0.009404 | 0.007370 | 0.010637 | 0.061783 | 0.008418 |
| reference_full_hexaplex_heavy_deduped |  | yes | 0.018272 | 0.082189 | 0.021626 | 0.024174 | 0.077104 | 0.020522 |
| reference_hexads_only_heavy_deduped |  | yes | 0.012565 | 0.002742 | 0.004412 | 0.008499 | 0.061171 | 0.007143 |
| reference_scaffold_only_complement_heavy_deduped | 1,2,3,4,5,6 | no | 0.015260 | 0.086576 | 0.019043 | 0.024993 | 0.015259 | 0.021512 |
| scaffold_blocks_1_2_3_4_5_6_heavy_deduped | 1,2,3,4,5,6 | no | 0.015260 | 0.086576 | 0.019043 | 0.024993 | 0.015259 | 0.021512 |
| scaffold_blocks_1_2_3_4_5_heavy_deduped | 1,2,3,4,5 | no | 0.018675 | 0.073344 | 0.017841 | 0.023641 | 0.015272 | 0.019299 |
| scaffold_blocks_1_2_3_4_heavy_deduped | 1,2,3,4 | no | 0.015140 | 0.054661 | 0.016290 | 0.021407 | 0.013753 | 0.015870 |
| scaffold_blocks_1_2_3_heavy_deduped | 1,2,3 | no | 0.013742 | 0.045187 | 0.013679 | 0.019481 | 0.013829 | 0.014327 |
| scaffold_blocks_1_2_heavy_deduped | 1,2 | no | 0.015062 | 0.032544 | 0.014381 | 0.018924 | 0.012721 | 0.011961 |
| scaffold_blocks_1_heavy_deduped | 1 | no | 0.030075 | 0.030418 | 0.014388 | 0.019305 | 0.014846 | 0.012157 |

## Automatically generated notes

- d_4p5_fraction is present in scaffold block 1 (0.014388), so this approximate profile does not require complete multi-block assembly to show a 4.5 A-window score.
- Across scaffold-only ladder rows, d_4p5_fraction increases from block 1 to all blocks.
- Adding hexads changes d_3p4_fraction by a mean paired difference of 0.051883.
- Pearson correlation between d_4p5_fraction and GLU-GLU motif count: 0.904.
- Pearson correlation between d_4p5_fraction and real-space contact count: 0.609.
- Pearson correlation between d_3p4_fraction and hexad inclusion flag: 0.982.
- The d_4p5 window should be interpreted cautiously as a reciprocal-space-like comparative score, not as a literal 4.5 A contact count.
