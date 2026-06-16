# Aleph Exploratory Status

Aleph is a geometric fingerprint for visualizing and quantifying ordered axial organization in the hexaplex model. It converts a hexaplex model into ordered per-unit structural traces along the helix axis.

Aleph is not a diffraction simulator. It does not prove formation, stability, or experimental correctness.

## Implementation Status

The Aleph workflow was implemented and validated on branch `aleph-fingerprint`, including the handedness-normalized twist update in commit `a23f4c4 Add handedness-normalized Aleph twist fingerprints`.

The workflow computes signed and absolute local twist, deviation from nominal 30-degree twist, rise, radial/phase traces, bend measures, chain coherence, QC flags, and SVG visualizations.

Validation passed with:

- analyzer compile
- analyzer run
- test compile
- pytest

The latest validation reported 28 per-unit rows, 23 plots, and 18 passing tests.

## Current Interpretation

The `central7` fragment recovers the intended approximately 30-degree local repeat by absolute local twist. This is useful as a geometric sanity check, but it mostly reflects the coordinate-derived construction of the current endpoint fragment.

The `central6` and `full` cases expose edge sensitivity and possible layer-assignment / antiparallel-ordering issues in the current Aleph definition. These are useful diagnostics, but they do not yet provide a clear independent structural-discrimination signal.

The current mini-hexaplex fragments are probably too small, too edge-sensitive, and too construction-derived for Aleph to provide strong additional scientific insight at this stage.

## Conclusion

Aleph should be parked as a validated exploratory side branch rather than merged into the main project narrative.

Aleph may become useful as a model-comparison tool if longer models, better layer assignment, full-chain directional handling, simulated assembly intermediates, or multiple independent candidate structures become available.

Spectral analysis remains optional and exploratory; it should only be revisited after the primary Aleph trace is geometrically trustworthy.
