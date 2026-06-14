# Mini-Hexaplex Helicity Report

## Purpose

This geometry-only analysis asks how helix-like each coordinate-truncated mini-hexaplex is as a function of base/GLU unit count, without rerunning diffraction.

## Metric Definition

- The common axis is fitted once from the full cleaned six-chain baseline using the same principal-axis logic as the mini truncation workflow.
- Each repeat unit is represented by the base-like CYP/MEP residue CA atom when present, otherwise by the CYP/MEP heavy-atom centroid.
- Representative points are converted to cylindrical coordinates around the full-baseline fitted axis.
- For each chain, theta is unwrapped after sorting by axial coordinate, then theta = omega*z + phi is fit by least squares.
- The raw helical score is mean_helical_r2 across chains, clipped to 0-1 as helical_coherence_score; it is retained as a local helix-fit diagnostic, not the main length-response metric.
- Because these mini-hexaplexes are coordinate truncations from an ideal helix, theta-vs-z linearity saturates at 1.0 and cannot distinguish structural length across the sampled truncations.
- The main structural-length metric is coherent_helical_turns = axial_extent_A / mean_pitch_A, with normalized_coherent_helical_turns_vs_full comparing each truncation against the full-length baseline.
- normalized_axial_extent_vs_full uses the existing geometry axial extent where available and the full-baseline PDB axial extent as the denominator.
- Six-strand phase coherence is estimated from angular spacing consistency at comparable axial ranks; it is secondary to the raw chain R2 score.

## Helicity Summary

| variant_id | units_per_chain | structural_coherence_flag | axial_extent_A | normalized_axial_extent_vs_full | coherent_helical_turns | normalized_coherent_helical_turns_vs_full | helical_coherence_score | mean_helical_r2 | mean_angular_residual_deg | std_twist_per_unit_deg | mean_pitch_A | ratio_to_full_length_4p5_5A |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| central4_units | 4 | not_compact/control_only | 15.451115 | 0.277828 | 0.378710 | 0.277832 | 1.000000 | 1.000000 | 0.004285 | 0.003994 | 40.799323 | 0.279337 |
| lower_end_first4_units | 4 | not_compact/control_only | 18.214105 | 0.327509 | 0.446397 | 0.327489 | 1.000000 | 1.000000 | 0.004680 | 0.002891 | 40.802468 | 0.246360 |
| central5_units | 5 | not_compact/control_only | 21.612096 | 0.388609 | 0.529710 | 0.388610 | 1.000000 | 1.000000 | 0.004864 | 0.003724 | 40.799841 | 0.351900 |
| lower_end_first5_units | 5 | not_compact/control_only | 21.612887 | 0.388623 | 0.529708 | 0.388608 | 1.000000 | 1.000000 | 0.005083 | 0.003803 | 40.801519 | 0.352143 |
| central6_units | 6 | borderline | 22.250214 | 0.400083 | 0.545353 | 0.400086 | 1.000000 | 1.000000 | 0.005256 | 0.004632 | 40.799634 | 0.370827 |
| lower_end_first6_units | 6 | borderline | 25.013991 | 0.449779 | 0.613078 | 0.449770 | 1.000000 | 1.000000 | 0.005093 | 0.004200 | 40.800695 | 0.391747 |
| central7_units | 7 | coherent | 28.414575 | 0.510925 | 0.696437 | 0.510926 | 1.000000 | 1.000000 | 0.005249 | 0.004607 | 40.799893 | 0.414638 |
| lower_end_first7_units | 7 | coherent | 28.414011 | 0.510915 | 0.696414 | 0.510908 | 1.000000 | 1.000000 | 0.005136 | 0.004301 | 40.800466 | 0.416322 |
| central8_units | 8 | coherent | 29.051919 | 0.522385 | 0.712064 | 0.522390 | 1.000000 | 1.000000 | 0.005108 | 0.004477 | 40.799598 | 0.499040 |
| lower_end_first8_units | 8 | coherent | 31.812792 | 0.572029 | 0.779720 | 0.572024 | 1.000000 | 1.000000 | 0.005190 | 0.004246 | 40.800265 | 0.485020 |
| central12_units | 12 | coherent | 42.651285 | 0.766916 | 1.045381 | 0.766920 | 1.000000 | 1.000000 | 0.005694 | 0.004359 | 40.799733 | 0.722478 |
| lower_end_first12_units | 12 | coherent | 45.413816 | 0.816590 | 1.113088 | 0.816592 | 1.000000 | 1.000000 | 0.005643 | 0.004297 | 40.799827 | 0.732788 |
| full_length_baseline | 15 | coherent | 55.613995 | 1.000000 | 1.363090 | 1.000000 | 1.000000 | 1.000000 | 0.005923 | 0.004253 | 40.799950 | 1.000000 |

## Plots

- outputs/mini_hexaplex/plots/mini_hexaplex_units_vs_helical_coherence_score.png
- outputs/mini_hexaplex/plots/mini_hexaplex_units_vs_mean_angular_residual_deg.png
- outputs/mini_hexaplex/plots/mini_hexaplex_units_vs_std_twist_per_unit_deg.png
- outputs/mini_hexaplex/plots/mini_hexaplex_units_vs_coherent_helical_turns.png
- outputs/mini_hexaplex/plots/mini_hexaplex_units_vs_normalized_axial_extent_vs_full.png
- outputs/mini_hexaplex/plots/mini_hexaplex_units_vs_ratio_to_full_4p5_5p0.png
- outputs/mini_hexaplex/plots/mini_hexaplex_helicity_vs_4p5_5p0_response.png
- outputs/mini_hexaplex/plots/mini_hexaplex_coherent_turns_and_4p5_5p0_response.png

## Conservative Interpretation

- Raw helix-linearity trend: central: saturated or effectively flat across sampled counts; lower-end: saturated or effectively flat across sampled counts.
- Coherent helical turns trend: central: monotonic nondecreasing; lower-end: monotonic nondecreasing.
- Normalized axial extent trend: central: monotonic nondecreasing; lower-end: monotonic nondecreasing.
- Shortest unit count still flagged geometrically coherent by the existing structural summary: 7.
- In this coordinate-truncation set, the raw theta-vs-z R2 score is saturated at 1.0 for all sampled lengths, so it does not by itself define the minimum meaningful mini-hexaplex length.
- The useful structural metric is coherent helical extent: axial length expressed as fitted helical turns, interpreted alongside normalized axial extent and the existing structural_coherence_flag.
- The 4.5-5.0 A intensity response rises with length while the raw helix-linearity score is already saturated; the diffraction response therefore appears to track structural extent more than local helix fit quality.
- The metric is most useful for comparing central and lower-end truncation families across length, not for proving independent stability.

## Limitations

- These are coordinate truncations only; no relaxation or minimization was performed.
- The metric depends on representative atom choice.
- Anti-parallel chains require careful angle unwrapping, and short chains can have inflated linear-fit scores.
- This analysis does not prove independent stability or a formation mechanism.
