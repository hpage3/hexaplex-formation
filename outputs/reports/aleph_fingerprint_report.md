# Aleph Geometric Fingerprint Prototype

Aleph is an exploratory one-dimensional geometric fingerprint along the fitted hexaplex axis. It is not a diffraction simulator and does not prove formation, stability, or experimental correctness.

Aleph asks whether each model has an ordered repeating geometric signature along its axis that can be plotted as per-unit signals and summarized with a discrete FFT.

## Inputs

- `full`: `outputs\intermediates\ai_candidate_inputs\full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain.pdb` (6 chains, 15 units, 2166 heavy atoms)
- `central6`: `outputs\mini_hexaplex\structures\mini_hexaplex_central6_units.pdb` (6 chains, 6 units, 864 heavy atoms)
- `central7`: `outputs\mini_hexaplex\structures\mini_hexaplex_central7_units.pdb` (6 chains, 7 units, 1008 heavy atoms)
- Optional compact twist variants: optional compact 24/30/36 variants were not found in this repo

## Aleph Summary

| Structure | Units | Mean local twist | Twist std | Mean rise | Rise std | Regular score |
|---|---:|---:|---:|---:|---:|---:|
| full | 15 | 14.638537 | 15.025284 | -0.167200 | 2.389845 | 0.399653 |
| central6 | 6 | 18.623268 | 20.542040 | -4.065739 | 0.549846 | 0.359578 |
| central7 | 7 | 30.000705 | 34.422308 | -3.857624 | 0.721071 | 0.301553 |

## FFT Summary

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

The full model is the best first Aleph FFT target because it has the longest ordered per-unit signal; shorter central6 and central7 fragments are useful for local comparison but have limited spectral resolution.
central6 and central7 can be compared to the full model by local twist, local rise, radial spread, and phase progression rather than by pair-distance counts alone.
Backbone/scaffold Aleph bend and base-like Aleph bend are reported separately where per-layer plane normals are computable; disagreement between them would indicate different scaffold and base-like geometric signals.
This first prototype looks promising only if the plots reveal repeat regularity or variant differences that were not already obvious from distance-window attribution.

## Assumptions And Cautions

- Unit assignment uses base-like CYP/MEP residues as repeat anchors by residue order within each chain.
- The Aleph phase and axial position use the first available chain-specific base-like centroid as an angular/axial anchor for each unit, while the layer centroid is still reported separately. This avoids symmetry cancellation of the six-strand layer centroid.
- FFT summaries for 6- and 7-unit models are explicitly marked as short-signal diagnostics.
- Aleph is a geometric fingerprint, not a diffraction simulator or structural mechanism.

## Outputs

- Per-unit CSV: `outputs\metrics\aleph_fingerprint_per_unit.csv`
- Summary CSV: `outputs\metrics\aleph_fingerprint_summary.csv`
- FFT CSV: `outputs\metrics\aleph_fingerprint_fft_summary.csv`
- Plot directory: `outputs\plots\aleph_fingerprint`
- Plot: `outputs\plots\aleph_fingerprint\aleph_local_twist_vs_unit.svg`
- Plot: `outputs\plots\aleph_fingerprint\aleph_local_rise_vs_unit.svg`
- Plot: `outputs\plots\aleph_fingerprint\aleph_radial_spread_vs_unit.svg`
- Plot: `outputs\plots\aleph_fingerprint\aleph_phase_progression_vs_unit.svg`
- Plot: `outputs\plots\aleph_fingerprint\aleph_fft_dominant_amplitudes.svg`
