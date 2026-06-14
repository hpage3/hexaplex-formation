# Research Roadmap

## Phase 1: Organize existing structures and outputs

Collect:

- full Hexaplex no-H deduped PDB
- hexads-only no-H deduped PDB
- scaffold-only no-H deduped PDB
- alanine beta-sheet controls
- central 5-alanine controls
- radial profile CSVs
- pair-attribution outputs
- experimental radial profiles or X-ray data
- PyMOL selections / strand mapping

## Phase 2: Build assembly-intermediate ladder

Generate structural subsets:

- isolated local peptide-pair controls
- native two-strand scaffold fragments
- three-strand scaffold subset
- four-strand scaffold subset
- six-strand scaffold-only model
- hexads-only model
- hexads + partial scaffold models
- full Hexaplex model

## Phase 3: Compute features for each model

For every model:

- contact map
- residue-pair distance histograms
- GLU–GLU motif recurrence
- backbone/side-chain oxygen motif recurrence
- helical axis
- radius and phase around axis
- twist / pitch estimate
- strand curvature
- scaffold closure metric
- powder diffraction radial profile
- intensity summaries in selected d-spacing windows

Key windows:

- d ≈ 8.4 Å
- d ≈ 5.5–6 Å
- d ≈ 4–5 Å
- d ≈ 3.3–3.4 Å
- d ≈ 3.0 Å

## Phase 4: Compare against controls and experimental data

Compare:

- alanine pair controls
- native scaffold fragments
- scaffold-only complement
- hexads-only
- full Hexaplex
- possible assembly intermediates

Primary question:

At what structural rung does the Hexaplex-like 4.5 Å scaffold signature appear?

## Phase 5: Schrödinger bridge / stochastic optimal transport

Use only after the intermediate ladder gives useful structural features.

Recommended first abstraction:

- contact-map bridge
- pair-distance-matrix bridge
- low-dimensional order-parameter bridge

Do not start with full atomistic coordinate bridge.

The bridge should generate and rank plausible low-action assembly pathways, not claim to prove the real physical assembly mechanism.
