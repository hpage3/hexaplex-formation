# pNAB Twist Diffraction Sensitivity Report

## Purpose

This report stages the official pNAB path for helical-twist variants around the current 30 degree hexaplex baseline. It is a forward-model sensitivity workflow, not structure determination.

## pNAB Source And Availability

- Upstream checkout: external/pnab
- Upstream commit: e82a3297cc3a744efec6134f8ef8f6b959194be1
- Python import status: not available (python cannot import pnab)
- `which pnab`: not on PATH
- `conda env list`: conda not installed/on PATH

Install command documented by pNAB README:

```bash
conda create -n pnab -c conda-forge pnab
conda activate pnab
```

No environment changes were made by this report.

## pNAB YAML Format

- pNAB input files are YAML mappings with `Backbone`, `HelicalParameters`, and `RuntimeParameters` sections.
- `HelicalParameters.h_twist` is in degrees.
- A single fixed value can be represented as `[30, 30, 1]`.
- A range can be represented as `[28, 32, 5]`, which pNAB expands into uniformly spaced values including endpoints.
- pNAB writes `results.csv`, `prefix.yaml`, and accepted conformer PDBs named `<Prefix>_<Conformer Index>.pdb`.

## Current Baseline Input Status

- Status: no baseline YAML supplied
- Current baseline PDB: outputs/intermediates/ai_candidate_inputs/full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain.pdb
- Current baseline model has CYP/MEP/GLU residues, six chains, and 15 base/GLU units per chain.
- No current-project baseline pNAB YAML was found or supplied.

## Official pNAB Hexad Example

- Status: baseline YAML supplied
- YAML: external/pnab/pnab/data/Hexad_Antiparallel.yaml
- Backbone file: threoninol_backbone.pdb
- Strand sequence: OOOOOOOOOO
- is_hexad: True
- h_twist: [24, 24, 1]
- h_rise: [3.4, 3.4, 1]
- x_displacement: [0, 0, 1]
- y_displacement: [0, 0, 1]
- inclination: [0, 0, 1]
- tip: [0, 0, 1]
- This official example is schema-relevant, but it is not automatically assumed to be the current 30 degree baseline because its sequence/components differ from the current CYP/MEP/GLU model.

## Missing Inputs Before Generating Current Twist Variants

- Current 30 degree pNAB YAML or equivalent builder parameter file.
- Backbone component file used for the current CYP/MEP/GLU model.
- Base/nucleobase component definitions for CYP/MEP or their pNAB one-letter codes.
- Sequence corresponding to 15 units per strand.
- Runtime parameters, energy filters, search algorithm, random seed, strand orientation, and build_strand flags used for the current model.
- Confirmation that only `h_twist` should vary while h_rise, x/y displacement, inclination, tip, sequence, and component files stay fixed.

## Generated Twist Manifest

- No pNAB twist manifest was found yet.

## Diffraction Workflow

- Once selected pNAB PDBs exist, run the existing normalization/heavy-atom conversion, detector simulation, and `radial_average.py` path.
- Native `.npy` detector/image-plate outputs remain unchanged.
- Comparisons to John Bacsa's experimental features must use radial `q_Ainv` or `d_A`, with q = 2*pi/d.

## Conservative Interpretation

- No twist-dependent structural or diffraction conclusion is made until pNAB-generated structures are available and validated.
- Width metrics should remain simulated width proxies, not Scherrer/domain-size estimates.
- The Emory data are fiber-like/oriented while current simulations are powder-averaged, limiting direct interpretation.

## Next Steps

1. Provide or reconstruct the exact baseline pNAB YAML and component files for the current 30 degree model.
2. Run `python3 scripts/generate_pnab_twist_variants.py --baseline-yaml <baseline.yaml>` in a pNAB-enabled environment.
3. Validate the regenerated 30 degree pNAB model against the current baseline PDB before interpreting twist variants.
4. Run diffraction/radial averaging and feature extraction on accepted pNAB structures.
5. Only after the single-parameter twist sweep works, run a small length-by-twist grid.
