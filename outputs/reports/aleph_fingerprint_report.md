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

## FFT Diagnostic Summary

| Structure | Signal | Samples | Dominant index | Normalized amplitude | Warning |
|---|---|---:|---:|---:|---|
| full | local_twist_deg | 14 | 1 | 0.681473 | numpy unavailable; used standard-library DFT fallback |
| full | local_rise_A | 14 | 1 | 0.816357 | numpy unavailable; used standard-library DFT fallback |
| full | radial_spread_A | 15 | 1 | 0.384917 | numpy unavailable; used standard-library DFT fallback |
| central6 | local_twist_deg | 5 | 1 | 0.703753 | numpy unavailable; used standard-library DFT fallback; short signal; spectral interpretation is limited |
| central6 | local_rise_A | 5 | 1 | 0.834914 | numpy unavailable; used standard-library DFT fallback; short signal; spectral interpretation is limited |
| central6 | radial_spread_A | 6 | 1 | 0.631820 | numpy unavailable; used standard-library DFT fallback; short signal; spectral interpretation is limited |
| central7 | local_twist_deg | 6 | 1 | 0.655957 | numpy unavailable; used standard-library DFT fallback; short signal; spectral interpretation is limited |
| central7 | local_rise_A | 6 | 1 | 0.806969 | numpy unavailable; used standard-library DFT fallback; short signal; spectral interpretation is limited |
| central7 | radial_spread_A | 7 | 1 | 0.681759 | numpy unavailable; used standard-library DFT fallback; short signal; spectral interpretation is limited |

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
- Plot: `outputs\plots\aleph_fingerprint\aleph_fft_dominant_amplitudes.svg`
