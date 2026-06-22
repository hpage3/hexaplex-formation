# HXC590 3.38 A Rise Sensitivity Report

## Scope

This run tests whether changing the ideal antiparallel 30 degree model rise from 3.40 A to 3.38 A changes the match to Nick's powder/radial peak positions. The models are rigid-layer geometric variants, not relaxed physical structures.

The baseline model is `inputs/nick_ideal_models/Hexaplex_AntiParallel_30deg_Ideal.pdb`. The focused generator writes only:

| model_id | twist_deg | rise_A | atom_count_preserved |
| --- | ---: | ---: | --- |
| twist30_rise3p40 | 30 | 3.40 | true |
| twist30_rise3p38 | 30 | 3.38 | true |

The generator preserves atom count, atom identity, residue identity, x/y coordinates, and intralayer geometry. It changes inferred layer z positions only.

## Inputs and scoring

Target peaks are in `inputs/experimental_peak_windows/hxc590_rise_twist_targets.csv`.

| target_label | experimental_d_A | target_group |
| --- | ---: | --- |
| base_stack | 3.35 | base_stacking |
| D | 7.3 | backbone_associated |
| C | 5.65 | backbone_associated |
| B | 4.4 | backbone_associated |
| A | 3.8 | backbone_associated |

Profiles were generated with the basic Debye radial-profile script and scored with the basic window scorer plus grouped peak-position scoring. Local maxima from each radial profile were matched to targets with a 0.40 A tolerance because the coarse q grid places the 7.3 A nearest local maximum at 7.662421 A.

## Peak-position fit summary

| model_id | target_group | root_sum_relative_error | rms_relative_error | weighted_rms_relative_error | missing_peak_count |
| --- | --- | ---: | ---: | ---: | ---: |
| twist30_rise3p38 | combined | 0.053861 | 0.024087 | 0.024087 | 0 |
| twist30_rise3p40 | combined | 0.054699 | 0.024462 | 0.024462 | 0 |
| twist30_rise3p38 | backbone_associated | 0.053779 | 0.026889 | 0.026889 | 0 |
| twist30_rise3p40 | backbone_associated | 0.054054 | 0.027027 | 0.027027 | 0 |
| twist30_rise3p38 | base_stacking | 0.002983 | 0.002983 | 0.002983 | 0 |
| twist30_rise3p40 | base_stacking | 0.008375 | 0.008375 | 0.008375 | 0 |

Per-peak matched positions:

| model_id | target_label | target_group | experimental_d_A | theoretical_d_A | abs_relative_error |
| --- | --- | --- | ---: | ---: | ---: |
| twist30_rise3p40 | base_stack | base_stacking | 3.35 | 3.378057 | 0.008375 |
| twist30_rise3p38 | base_stack | base_stacking | 3.35 | 3.359992 | 0.002983 |
| twist30_rise3p40 | D | backbone_associated | 7.3 | 7.662421 | 0.049647 |
| twist30_rise3p38 | D | backbone_associated | 7.3 | 7.662421 | 0.049647 |
| twist30_rise3p40 | C | backbone_associated | 5.65 | 5.764390 | 0.020246 |
| twist30_rise3p38 | C | backbone_associated | 5.65 | 5.764390 | 0.020246 |
| twist30_rise3p40 | B | backbone_associated | 4.4 | 4.424778 | 0.005631 |
| twist30_rise3p38 | B | backbone_associated | 4.4 | 4.393836 | 0.001401 |
| twist30_rise3p40 | A | backbone_associated | 3.8 | 3.785051 | 0.003934 |
| twist30_rise3p38 | A | backbone_associated | 3.8 | 3.785051 | 0.003934 |

## Interpretation

Changing rise from 3.40 A to 3.38 A improves the base-stacking peak-position match to the 3.35 A target. The matched peak shifts from 3.378057 A to 3.359992 A, reducing the base-stacking relative error from 0.008375 to 0.002983.

The backbone-associated score changes only slightly. The 7.3 A, 5.65 A, and 3.8 A matched positions are unchanged in this coarse basic-Debye profile. The 4.4 A matched peak improves from 4.424778 A to 4.393836 A, which slightly improves the backbone-associated RMS score.

The combined score improves modestly, but most of that clear, expected improvement is in the base-stacking region. This supports using 3.38 A as the fixed rise value for the next twist scan, while keeping base-stacking and backbone-associated scores separate.

## Limitations

These are geometric models, not MD-relaxed structures. Improving the base-stacking peak by changing rise is expected and does not by itself determine twist. The longer backbone-associated peaks are likely more informative for twist/pitch discrimination. The basic Debye radial profile is a simplified comparison, not a full fiber-pattern or intensity fit.

## Next twist scan plan

Use fixed rise 3.38 A for the next scan unless a chemically relaxed model set suggests otherwise. The next scan should cover twists from 25 to 35 degrees at 1 degree increments if model generation is computationally reasonable; otherwise use 2 degree increments for a first pass.

Exact inputs needed:

- one coordinate model per twist value at fixed 3.38 A rise
- preserved residue and atom naming across the model set
- a manifest with `model_id`, `twist_deg`, `rise_A`, and source file path
- radial profiles and peak assignments scored with the grouped target file used here

The assignment template is `inputs/theoretical_peak_assignments/hxc590_rise_twist_scan_assignments.csv`.
