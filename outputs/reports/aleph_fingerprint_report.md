# Aleph Geometric Fingerprint Prototype

Aleph is an exploratory one-dimensional geometric fingerprint along the fitted hexaplex axis. It is not a diffraction simulator and does not prove formation, stability, or experimental correctness.

Aleph asks whether each model has an ordered repeating geometric signature along its axis that can be plotted as per-unit signals and summarized with a discrete FFT.

## Inputs

- `full`: `outputs\intermediates\ai_candidate_inputs\full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain.pdb` (6 chains, 15 units, 2166 heavy atoms)
- `central6`: `outputs\mini_hexaplex\structures\mini_hexaplex_central6_units.pdb` (6 chains, 6 units, 864 heavy atoms)
- `central7`: `outputs\mini_hexaplex\structures\mini_hexaplex_central7_units.pdb` (6 chains, 7 units, 1008 heavy atoms)
- Optional compact twist variants: optional compact 24/30/36 variants were not found in this repo

## Aleph Summary

| Structure | Units | Mean signed twist | Mean abs twist | Twist std | Mean rise | Rise std | Regular score |
|---|---:|---:|---:|---:|---:|---:|---:|
| full | 15 | -14.638537 | 16.251148 | 15.025284 | 0.167200 | 2.389845 | 0.399653 |
| central6 | 6 | -18.623268 | 18.623268 | 20.542040 | 4.065739 | 0.549846 | 0.359578 |
| central7 | 7 | -30.000705 | 30.000705 | 34.422308 | 3.857624 | 0.721071 | 0.301553 |

## Geometry QC

| Structure | Axis flipped | Phase unwrap ok | Expected units | Observed units | Missing-chain units | Negative rises | Twist warning | Rise warning | Plane warning |
|---|---|---|---:|---:|---:|---:|---|---|---|
| full | true | true | 15 | 15 | 0 | 7 | mean absolute local twist differs from nominal 30 deg by >10 deg | unit 1 to 2 has negative rise -2.35 A; unit 2 to 3 has negative rise -2.48 A; unit 3 to 4 has negative rise -1.77 A; unit 4 to 5 has negative rise -0.44 A; unit 12 to 13 has negative rise -1.44 A; unit 13 to 14 has negative rise -2.35 A; unit 14 to 15 has negative rise -2.48 A; 7 negative local rise values after axis orientation |  |
| central6 | true | true | 6 | 6 | 0 | 0 | mean absolute local twist differs from nominal 30 deg by >10 deg |  |  |
| central7 | true | true | 7 | 7 | 0 | 0 |  |  |  |

## Series-style Aleph fingerprint

The Aleph series fingerprint is the primary visual fingerprint in this prototype. It is closer to a true structural fingerprint because it represents the hexaplex as an ordered one-dimensional trace rather than a multi-feature dashboard.

The primary Aleph series is signed local twist between adjacent units, in degrees. The x-axis labels show unit transitions along the hexaplex, so `U1->U2` is the step from unit 1 to unit 2. The y-axis is signed local twist after phase unwrapping. The dashed 30 deg line is a geometric reference, not an experimental target, and markers with QC warnings are highlighted.

Richer unit labels are stored in the per-unit CSV. `unit_label` gives the concise unit name, `transition_label` gives the plotted transition label where a local transition exists, and `residue_label_summary` records the chain/residue labels contributing to each unit.

Under the current Aleph definitions, central7 currently looks like the cleanest 30 deg-like Aleph series fingerprint because its mean absolute local twist is near 30 deg, its rise is positive, and it has no QC warnings. central6 is shorter and has positive rise, but its signed twist trace and mean absolute twist deviate from the nominal 30 deg value. The full model remains a geometry-definition diagnostic case because its current warnings indicate that layer assignment and antiparallel ordering require further inspection.

This representation provides a more natural ordered signal for future DFT/FFT exploration than the feature comparison panel. Whether spectral analysis adds value remains an open question.

Series outputs:

- `outputs\plots\aleph_fingerprint\aleph_series_fingerprint_full.svg`
- `outputs\plots\aleph_fingerprint\aleph_series_fingerprint_central6.svg`
- `outputs\plots\aleph_fingerprint\aleph_series_fingerprint_central7.svg`
- `outputs\plots\aleph_fingerprint\aleph_series_fingerprint_comparison.svg`
- `outputs\plots\aleph_fingerprint\aleph_series_companion_traces_full.svg`
- `outputs\plots\aleph_fingerprint\aleph_series_companion_traces_central7.svg`

## Visual Aleph feature comparison panel

Aleph converts the hexaplex into ordered per-unit geometric traces. The feature comparison panel places unit index on the x-axis and stacks Aleph features as rows, so changes in twist, rise, radial spread, phase, plane bend, chain coherence, and QC flags can be scanned compactly.

Rows labeled `abs twist` and `signed twist` show the local angular step after phase unwrapping; `rise` shows the axis-oriented local axial step; `radial spread` tracks base-like radial variability; `phase` shows the unwrapped angular progression; `base bend` and `scaffold bend` compare adjacent fitted plane normals; `chain coherence` reports the circular mean resultant length; and the QC row marks per-unit warnings.

The visualization is a structural fingerprint, not a diffraction simulation. It reveals axial ordering, local geometric irregularity, phase progression, and QC behavior that pair-distance counts do not show directly.

Feature comparison panel outputs:

- `outputs\plots\aleph_fingerprint\aleph_fingerprint_strip_full.svg`
- `outputs\plots\aleph_fingerprint\aleph_fingerprint_strip_central6.svg`
- `outputs\plots\aleph_fingerprint\aleph_fingerprint_strip_central7.svg`
- `outputs\plots\aleph_fingerprint\aleph_fingerprint_comparison.svg`
- `outputs\plots\aleph_fingerprint\aleph_fingerprint_qc_summary.svg`

## FFT Diagnostic Summary

| Structure | Signal | Samples | Dominant index | Normalized amplitude | Warning |
|---|---|---:|---:|---:|---|
| full | local_twist_deg | 14 | 1 | 0.681473 |  |
| full | local_rise_A | 14 | 1 | 0.816357 |  |
| full | radial_spread_A | 15 | 1 | 0.384917 |  |
| central6 | local_twist_deg | 5 | 1 | 0.703753 | short signal; spectral interpretation is limited |
| central6 | local_rise_A | 5 | 1 | 0.834914 | short signal; spectral interpretation is limited |
| central6 | radial_spread_A | 6 | 1 | 0.631820 | short signal; spectral interpretation is limited |
| central7 | local_twist_deg | 6 | 1 | 0.655957 | short signal; spectral interpretation is limited |
| central7 | local_rise_A | 6 | 1 | 0.806969 | short signal; spectral interpretation is limited |
| central7 | radial_spread_A | 7 | 1 | 0.681759 | short signal; spectral interpretation is limited |

## Interpretation

The full model remains the best first Aleph diagnostic target because it has the longest ordered per-unit signal; shorter central6 and central7 fragments are useful for local comparison but have limited spectral resolution.
central6 and central7 can be compared to the full model by local twist, local rise, radial spread, and phase progression rather than by pair-distance counts alone.
Backbone/scaffold Aleph bend and base-like Aleph bend are reported separately where per-layer plane normals are computable; disagreement between them would indicate different scaffold and base-like geometric signals.
This QC pass should be read before expanding FFT interpretation: stable axis orientation, phase unwrapping, positive rise convention, and bounded circular dispersion are prerequisites for treating Aleph signals as repeat fingerprints.

## Assumptions And Cautions

- Unit assignment uses base-like CYP/MEP residues as repeat anchors by residue order within each chain.
- The Aleph phase and axial position use the first available chain-specific base-like centroid as an angular/axial anchor for each unit, while the layer centroid is still reported separately. This avoids symmetry cancellation of the six-strand layer centroid.
- The fitted axis is flipped when needed so the chain-specific anchor axial coordinate generally increases with unit index.
- Chain angular dispersion is computed with bounded circular statistics and paired with mean resultant length.
- FFT summaries for 6- and 7-unit models are explicitly marked as short-signal diagnostics.
- Aleph is a geometric fingerprint, not a diffraction simulator or structural mechanism.

## Outputs

- Per-unit CSV: `outputs\metrics\aleph_fingerprint_per_unit.csv`
- Summary CSV: `outputs\metrics\aleph_fingerprint_summary.csv`
- FFT CSV: `outputs\metrics\aleph_fingerprint_fft_summary.csv`
- QC CSV: `outputs\metrics\aleph_fingerprint_qc.csv`
- Plot directory: `outputs\plots\aleph_fingerprint`
- Plot: `outputs\plots\aleph_fingerprint\aleph_local_twist_vs_unit.svg`
- Plot: `outputs\plots\aleph_fingerprint\aleph_local_rise_vs_unit.svg`
- Plot: `outputs\plots\aleph_fingerprint\aleph_radial_spread_vs_unit.svg`
- Plot: `outputs\plots\aleph_fingerprint\aleph_phase_raw_vs_unit.svg`
- Plot: `outputs\plots\aleph_fingerprint\aleph_phase_progression_vs_unit.svg`
- Plot: `outputs\plots\aleph_fingerprint\aleph_chain_resultant_vs_unit.svg`
- Plot: `outputs\plots\aleph_fingerprint\aleph_series_fingerprint_full.svg`
- Plot: `outputs\plots\aleph_fingerprint\aleph_series_fingerprint_central6.svg`
- Plot: `outputs\plots\aleph_fingerprint\aleph_series_fingerprint_central7.svg`
- Plot: `outputs\plots\aleph_fingerprint\aleph_series_fingerprint_comparison.svg`
- Plot: `outputs\plots\aleph_fingerprint\aleph_series_companion_traces_full.svg`
- Plot: `outputs\plots\aleph_fingerprint\aleph_series_companion_traces_central7.svg`
- Plot: `outputs\plots\aleph_fingerprint\aleph_fingerprint_strip_full.svg`
- Plot: `outputs\plots\aleph_fingerprint\aleph_fingerprint_strip_central6.svg`
- Plot: `outputs\plots\aleph_fingerprint\aleph_fingerprint_strip_central7.svg`
- Plot: `outputs\plots\aleph_fingerprint\aleph_fingerprint_comparison.svg`
- Plot: `outputs\plots\aleph_fingerprint\aleph_fingerprint_qc_summary.svg`
- Plot: `outputs\plots\aleph_fingerprint\aleph_fft_dominant_amplitudes.svg`
