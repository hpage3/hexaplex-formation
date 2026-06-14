# Mini-Hexaplex Length Response Report

## Purpose

This pass asks whether fewer than 8 base/GLU units per strand can still preserve a six-strand mini-hexaplex geometry and whether the 4.5-5.0 A diffraction-window signal appears below the 8-unit models.

This is a sensitivity and compatibility analysis only. It does not determine the true structure.

## Baseline And Convention

- Baseline full-length structure: outputs/intermediates/ai_candidate_inputs/full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain.pdb
- Variant manifest: outputs/mini_hexaplex/mini_hexaplex_variant_manifest.csv
- Geometry sanity table: outputs/metrics/mini_hexaplex_geometry_summary.csv
- Analysis mode: fiber
- Unit convention reminder: q = 2*pi/d.
- Full cleaned baseline length: 15 base/GLU units per chain, 30 residues per chain, 180 residues total.
- These variants are coordinate truncations, not independently relaxed or minimized mini-structures.
- CentralN and lower_end_firstN are the physically meaningful truncations; literal_firstN is a sequence-order control and can be misleading in anti-parallel geometry.

## Geometry Summary

- Coherent counts among the primary non-literal variants: 7, 8, 12.
- Borderline counts among the primary non-literal variants: 6.
- Control-only counts: 4, 5.
- Shortest primary count flagged coherent by the conservative geometry heuristic: 7.

## Feature Windows

- 4p5_5A: d = 4.5-5 A, q ~= 1.257-1.396 A^-1
- 3p4A: d = 3.25-3.55 A, q ~= 1.770-1.933 A^-1
- 3p0A: d = 2.9-3.1 A, q ~= 2.027-2.167 A^-1
- 4p1_8p4A: d = 4.1-8.4 A, q ~= 0.748-1.532 A^-1

## Feature Summary

| variant_id | units_per_chain | structural_coherence_flag | d_A_at_max_in_4p5_5A_window | integrated_intensity_4p5_5A | ratio_to_full_length_4p5_5A | ratio_to_matching_8unit_model_4p5_5A | integrated_intensity_3p4A | integrated_intensity_3p0A | comparison_to_8unit_models |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| central4_units | 4 | not_compact/control_only | 4.613187 | 3444006.874750 | 0.279337 | 0.559749 | 7339796.804455 | 3302018.344523 | weaker than matching 8-unit model; 4.5-5.0 A integrated intensity ratio vs 8-unit = 0.560 |
| lower_end_first4_units | 4 | not_compact/control_only | 4.613187 | 3037422.968126 | 0.246360 | 0.507938 | 7023346.989272 | 3269048.322854 | weaker than matching 8-unit model; 4.5-5.0 A integrated intensity ratio vs 8-unit = 0.508 |
| central5_units | 5 | not_compact/control_only | 4.613187 | 4338648.731845 | 0.351900 | 0.705154 | 10048019.188233 | 3734306.002264 | weaker than matching 8-unit model; 4.5-5.0 A integrated intensity ratio vs 8-unit = 0.705 |
| lower_end_first5_units | 5 | not_compact/control_only | 4.613187 | 4341638.908187 | 0.352143 | 0.726037 | 10548836.733031 | 3936485.500936 | weaker than matching 8-unit model; 4.5-5.0 A integrated intensity ratio vs 8-unit = 0.726 |
| central6_units | 6 | borderline | 4.613187 | 4571997.760169 | 0.370827 | 0.743080 | 14680431.815599 | 4422994.042244 | weaker than matching 8-unit model; 4.5-5.0 A integrated intensity ratio vs 8-unit = 0.743 |
| lower_end_first6_units | 6 | borderline | 4.613187 | 4829934.172087 | 0.391747 | 0.807693 | 13995806.281239 | 4258544.181266 | weaker than matching 8-unit model; 4.5-5.0 A integrated intensity ratio vs 8-unit = 0.808 |
| central7_units | 7 | coherent | 4.613187 | 5112151.537978 | 0.414638 | 0.830870 | 18586021.083765 | 4524187.480741 | weaker than matching 8-unit model; 4.5-5.0 A integrated intensity ratio vs 8-unit = 0.831 |
| lower_end_first7_units | 7 | coherent | 4.613187 | 5132924.530402 | 0.416322 | 0.858361 | 18545808.663313 | 4568137.980928 | weaker than matching 8-unit model; 4.5-5.0 A integrated intensity ratio vs 8-unit = 0.858 |
| central8_units | 8 | coherent | 4.613187 | 6152765.982527 | 0.499040 |  | 22997849.857152 | 4807186.342551 | 8-unit reference |
| lower_end_first8_units | 8 | coherent | 4.613187 | 5979913.762520 | 0.485020 |  | 22690277.909479 | 4884927.972859 | 8-unit reference |
| central12_units | 12 | coherent | 4.613187 | 8907574.386349 | 0.722478 | 1.447735 | 43524524.035499 | 6736996.903022 | stronger than matching 8-unit model; 4.5-5.0 A integrated intensity ratio vs 8-unit = 1.448 |
| lower_end_first12_units | 12 | coherent | 4.613187 | 9034690.426107 | 0.732788 | 1.510840 | 43421239.834538 | 6897059.111317 | stronger than matching 8-unit model; 4.5-5.0 A integrated intensity ratio vs 8-unit = 1.511 |
| full_length_baseline | 15 | coherent | 4.613187 | 12329206.192042 | 1.000000 |  | 60652814.836791 | 8540328.003391 | full-length reference |

## Length Response

- 4.5-5.0 A signal below 8 units: present in the primary central/lower-end truncations.
- Central-family trend from 4 to 12 units: monotonic nondecreasing across the sampled counts.
- Lower-end family trend from 4 to 12 units: monotonic nondecreasing across the sampled counts.
- 3.4 A persistence: The 3.4 A reference window has nonzero intensity in all analyzed mini variants.
- 3.0 A persistence: The 3.0 A reference window has nonzero intensity in all analyzed mini variants.

## Plots

- outputs/mini_hexaplex/plots/mini_hexaplex_q_profile_overlay.png
- outputs/mini_hexaplex/plots/mini_hexaplex_d_profile_overlay.png
- outputs/mini_hexaplex/plots/mini_hexaplex_central_d_profile_overlay.png
- outputs/mini_hexaplex/plots/mini_hexaplex_lower_end_d_profile_overlay.png
- outputs/mini_hexaplex/plots/mini_hexaplex_central_d_4p1_8p4_zoom.png
- outputs/mini_hexaplex/plots/mini_hexaplex_lower_end_d_4p1_8p4_zoom.png
- outputs/mini_hexaplex/plots/mini_hexaplex_central_d_4p5_5p0_zoom.png
- outputs/mini_hexaplex/plots/mini_hexaplex_lower_end_d_4p5_5p0_zoom.png
- outputs/mini_hexaplex/plots/mini_hexaplex_integrated_4p5_5p0_by_variant.png
- outputs/mini_hexaplex/plots/mini_hexaplex_units_vs_integrated_4p5_5p0.png
- outputs/mini_hexaplex/plots/mini_hexaplex_units_vs_ratio_to_full_4p5_5p0.png
- outputs/mini_hexaplex/plots/mini_hexaplex_units_vs_3p4A_intensity.png
- outputs/mini_hexaplex/plots/mini_hexaplex_units_vs_3p0A_intensity.png
- outputs/mini_hexaplex/plots/mini_hexaplex_reference_window_intensity_by_variant.png
- outputs/mini_hexaplex/plots/mini_hexaplex_reference_window_intensity_by_unit_count.png

## Conservative Interpretation

- Fewer than 8 units: yes, the central/lower-end variants retain six-strand geometry at 6-7 units; 4-5 units fall into the not-compact/control-only category.
- The 4.5-5.0 A signal appears below 8 units: yes; the signal rises from 4 to 12 units rather than appearing only at 8.
- The 4.5-5.0 A window does not receive a conservative local-maximum call in any of the analyzed short variants.
- The signal grows with length in a broadly smooth way across the central and lower-end families, with 12-unit models recovering more of the full-length signal than 8-unit models.
- 3.4 A and 3.0 A signals persist across the sampled lengths.
- A minimum coherent length suggested by this conservative geometry flag: 7.

## Limitations

- No minimization or molecular dynamics was performed.
- End effects may matter for the shortest truncations.
- The diffraction calculation is idealized and comparative, not a full experiment model.
- Emory data are oriented/fiber-like, not random powder.
- Raw experimental radial data are needed for stronger conclusions.
