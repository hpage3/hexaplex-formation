# pNAB Fixed-Rise Twist-Series Analysis

- Source workflow: `inputs/pnab_asem_30deg`
- Original source: `C:\Users\hpage3\OneDrive - Georgia Institute of Technology\Documents\GitHub\research\30_A`
- Fixed rise: 3.38 A
- Attempted twist grid: 28.0 to 32.0 degrees in 0.5-degree increments
- Generated structure directory: `outputs/pnab_twist_series_rise3p38/structures`
- Analysis mode: parallel
- Worker count: 16

## Analysis Status

- Successful computations: 0
- Existing outputs reused: 5
- Available analyzed model outputs: 5 (28.5, 29.5, 30.0, 30.5, 31.0)
- Missing generation outputs skipped: 4
- Generation timeouts without promoted structures: 28.0, 29.0, 31.5, 32.0
- Analysis failures: 0
- Analysis timeouts: 0

## Comparative Readout

- Best available penalized combined local-window score: 31.0 degrees (0.0260574470638).
- The available 31.0-degree model ranks 1 of 5 by combined score.
- This does not establish that 31.0 degrees is preferred or refute the earlier disfavored assessment: it is the only available model with an accepted interior B-window maximum, while the other available models receive a missing-peak penalty.
- The 29.5-degree combined score is lower than the 30.0-degree score, but both lack an accepted B-window assignment and the numerical separation is small.
- The 30.5-degree combined score is lower than the 30.0-degree score.

These rankings are comparative and are not proof of a unique structure. Missing pNAB models are not excluded scientifically; they were unavailable because generation timed out.

## Caveats

- pNAB structures are not MD- or QM-relaxed in this workflow.
- The radial profiles use a simplified isotropic Debye approximation.
- Powder/radial profiles discard directional two-dimensional information.
- Peak assignments require an interior local maximum in a predefined target window; ambiguous windows are left unassigned.
- Parallel execution distributes independent model jobs. Each model calculation remains serial, avoiding nested worker pools.
- Generated PDB structures and per-twist work directories are retained locally but are not intended for this commit; the generator and manifests reproduce their locations.
