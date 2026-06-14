# Fitted-Axis Helical Order Report

## Scientific cautions

The fitted-axis helical metrics are geometric descriptors computed from residue centroids. They do not prove a physical helix, formation pathway, or molecular dynamics. The d ~= 4.5 A diffraction feature remains a reciprocal-space scaffold signature in the working hypothesis, not a literal atom-contact assignment.

The current block map is a candidate contiguous-residue map. Fitted-axis results help test whether candidate blocks behave like folded or twisted paths, but the mapping still requires validation against PyMOL colored strand paths.

Unwrapped angle span, approximate turns, and pitch-like values depend on residue order, so they should be interpreted as order-aware descriptors rather than invariant helix parameters.

## Comparison table

| model | residue_count | mean_radius_fitted | axial_span | angular_coverage_rad | approximate_turns | approximate_pitch_per_turn | z_axis_angle_degrees |
| --- | --- | --- | --- | --- | --- | --- | --- |
| alanine_beta_sheet_central5_res3_to_7_heavy_deduped | 5 | 0.000074 | 0.119109 | 3.250530 | 0.957390 | 0.124410 | 72.856160 |
| alanine_beta_sheet_central5_res4_to_8_heavy_deduped | 5 | 0.000065 | 0.119109 | 3.623654 | 0.897286 | 0.132744 | 72.857751 |
| alanine_beta_sheet_full_heavy_deduped | 10 | 0.000087 | 0.172090 | 4.688804 | 1.445676 | 0.119038 | 72.874865 |
| full_hexaplex_anti_parallel_30deg_ideal_heavy_deduped | 180 | 10.643928 | 52.382828 | 6.156728 | 0.758371 | 69.072804 | 0.109918 |
| hexads_plus_scaffold_blocks_1_2_3_4_5_6_heavy_deduped | 180 | 10.643928 | 52.382828 | 6.156728 | 1.758371 | 29.790540 | 0.109918 |
| hexads_plus_scaffold_blocks_1_2_3_4_5_heavy_deduped | 165 | 10.172044 | 52.858258 | 6.256888 | 2.879391 | 18.357442 | 1.918655 |
| hexads_plus_scaffold_blocks_1_2_3_4_heavy_deduped | 150 | 9.584398 | 53.335628 | 6.241806 | 1.406890 | 37.910297 | 3.974386 |
| hexads_plus_scaffold_blocks_1_2_3_heavy_deduped | 135 | 8.916554 | 51.676227 | 6.184021 | 2.569137 | 20.114237 | 3.694505 |
| hexads_plus_scaffold_blocks_1_2_heavy_deduped | 120 | 8.055510 | 51.403309 | 6.106985 | 1.096362 | 46.885350 | 4.369672 |
| hexads_plus_scaffold_blocks_1_heavy_deduped | 105 | 7.001921 | 49.898590 | 6.217255 | 2.256204 | 22.116166 | 2.418610 |
| hexaplex_hexads_only_heavy_deduped | 90 | 5.565012 | 47.600912 | 5.767437 | 0.834589 | 57.035158 | 0.109911 |
| hexaplex_scaffold_only_complement_heavy_deduped | 180 | 12.097009 | 52.382856 | 6.071303 | 0.750925 | 69.757728 | 0.109866 |
| reference_full_hexaplex_heavy_deduped | 180 | 10.643928 | 52.382828 | 6.156728 | 0.758371 | 69.072804 | 0.109918 |
| reference_hexads_only_heavy_deduped | 90 | 5.565012 | 47.600912 | 5.767437 | 0.834589 | 57.035158 | 0.109911 |
| reference_scaffold_only_complement_heavy_deduped | 180 | 12.097009 | 52.382856 | 6.071303 | 0.750925 | 69.757728 | 0.109866 |
| scaffold_blocks_1_2_3_4_5_6_heavy_deduped | 180 | 12.097009 | 52.382856 | 6.071303 | 0.750925 | 69.757728 | 0.109866 |
| scaffold_blocks_1_2_3_4_5_heavy_deduped | 150 | 12.045789 | 53.215110 | 6.174892 | 1.861070 | 28.593819 | 3.468185 |
| scaffold_blocks_1_2_3_4_heavy_deduped | 120 | 11.806947 | 54.115822 | 6.214849 | 0.376982 | 143.550255 | 8.305521 |
| scaffold_blocks_1_2_3_heavy_deduped | 90 | 11.654513 | 49.947035 | 6.237544 | 1.566860 | 31.877143 | 10.018321 |
| scaffold_blocks_1_2_heavy_deduped | 60 | 11.225443 | 48.093436 | 6.075752 | 0.070355 | 683.579556 | 14.209179 |
| scaffold_blocks_1_heavy_deduped | 30 | 11.186964 | 44.149981 | 6.233895 | 1.346288 | 32.793854 | 14.350630 |

## Automatically generated notes

- For the scaffold-only complement, the fitted principal axis is 0.110 degrees from the global z-axis.
- Scaffold block 1 still shows near-full angular coverage under the fitted-axis metric (6.234 rad), supporting the current individual folded/twisted path interpretation.
- These fitted-axis results are geometric descriptors, not proof of a physical helix or temporal assembly order.
- The candidate block mapping still requires validation against PyMOL colored strand paths.
