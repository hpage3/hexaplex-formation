# Hexaplex Formation

This project investigates how the proposed Hexaplex structure forms, with emphasis on:

1. When the helical Hexaplex scaffold emerges during assembly.
2. What interactions stabilize the final helical scaffold.
3. Whether Schrödinger bridge / stochastic optimal transport methods can help infer plausible formation pathways.

This is a standalone project. It may reuse outputs or ideas from the sibling `fiber-diffraction` project, but it should not depend on that repository's internal directory structure.

## Current working interpretation

The d ≈ 4.5 Å diffraction feature is treated as a useful signature of the Hexaplex scaffold architecture because it differentiates the Hexaplex scaffold from isolated beta-strand / alanine-pair controls.

Important caution: the d ≈ 4.5 Å feature is a reciprocal-space feature. It should not be interpreted as meaning that the responsible atoms are separated by a literal 4.5 Å PyMOL distance.

## Current scientific summary

The current state-of-the-project mechanistic interpretation is summarized in `outputs/reports/current_mechanistic_summary.md`.

## Initial research focus

The first computational goal is to build an assembly-intermediate ladder and determine when the following features emerge:

- helical scaffold order
- repeated GLU-rich scaffold geometry
- structured diffraction features in the 4–6 Å region
- the hexad-associated d ≈ 3.3–3.4 Å feature

## Relationship to fiber-diffraction

The sibling project `../fiber-diffraction` contains diffraction simulation tools and prior radial-profile work.

This project may consume exported PDBs, CSVs, JSON metadata, and figures from that work, but this repository should remain logically independent.

## Setup

Use the existing virtual environment for this repo, then install the local package in editable mode:

```bash
source .venv/bin/activate
python -m pip install -e .
```

The first scripts intentionally use only the Python standard library. `pytest` is needed for the test suite.

## Tests

```bash
pytest
```

## Initial Metrics

Place raw `.pdb` files in `inputs/structures/`, then run:

```bash
scripts/run_initial_metrics.sh
```

Raw input structures may include hydrogens and exact duplicate atom records. This project intentionally preserves hydrogen-containing structures for formation and stabilization questions, including hydrogen bonding, protonation hypotheses, scaffold stabilization, and emergence of helical architecture.

The runner first writes normalized working structures to:

- `outputs/intermediates/normalized_structures/<basename>_allatom_deduped.pdb`
- `outputs/intermediates/normalized_structures/<basename>_heavy.pdb`
- `outputs/intermediates/normalized_structures/<basename>_heavy_deduped.pdb`

Main fast structure metrics use the `_heavy_deduped.pdb` files for speed and comparison to prior no-H diffraction calculations. Hydrogen-aware GLU-only local geometry metrics use the `_allatom_deduped.pdb` files.

The runner creates output directories if needed and writes:

- `outputs/metrics/structure_inventory.csv`
- `outputs/metrics/structure_normalization_summary.csv`
- `outputs/metrics/normalized_structure_inventory.csv`
- `outputs/metrics/<basename>_pair_distances.csv`
- `outputs/metrics/<basename>_GLU_pair_distances.csv`
- `outputs/metrics/<basename>_GLU_allatom_pair_distances.csv`

Baseline pair-distance histograms are computed from normalized `_heavy_deduped.pdb` structures. GLU all-atom pair-distance histograms are computed from normalized `_allatom_deduped.pdb` structures that contain GLU residues.

You can also run the scripts directly:

```bash
python3 scripts/normalize_structures.py
python3 scripts/inventory_structures.py
python3 scripts/inventory_structures.py --input-dir outputs/intermediates/normalized_structures --out outputs/metrics/normalized_structure_inventory.csv
python3 scripts/pair_distance_histogram.py --pdb outputs/intermediates/normalized_structures/example_heavy_deduped.pdb --out outputs/metrics/example_pair_distances.csv
python3 scripts/pair_distance_histogram.py --pdb outputs/intermediates/normalized_structures/example_allatom_deduped.pdb --out outputs/metrics/example_GLU_allatom_pair_distances.csv --all-atoms --residue-filter GLU
```

Pair-distance histograms are real-space atom-pair summaries. They should not be confused with reciprocal-space d-spacing features such as the d ≈ 4.5 Å scaffold signature.

For the current Hexaplex working set, useful heavy-deduped sanity-check counts are approximately:

- full Hexaplex: 2166 atoms
- hexads-only: 810 atoms
- scaffold-only complement: 1356 atoms

These real-space normalized-structure counts do not imply that hexads-only and scaffold-only diffraction intensities are additive decompositions of the full structure. Scattering amplitudes add before intensities are squared.

## Formation Metrics

Run the first formation-oriented structural summaries with:

```bash
scripts/run_formation_metrics.sh
```

This refreshes normalized structures, then writes CSV outputs under `outputs/metrics/formation/`.

Contact maps are real-space residue proximity summaries using minimum inter-residue atom distances. They are useful for comparing scaffold closure and contact recurrence, but they are not reciprocal-space d-spacing assignments.

GLU motif summaries count GLU-containing contact patterns from those contact maps. These are intended to test whether GLU-rich contacts distinguish the organized Hexaplex scaffold from alanine beta-strand controls.

Hydrogen-bond candidate files use all-atom deduped structures and simple donor-H-acceptor distance cutoffs. They are rough geometric candidate lists only; definitive hydrogen-bond assignment still requires angle, protonation, and chemical validation.

Helical order summaries currently use a simple global z-axis approximation from residue centroids. This is not a fitted helical axis and does not yet include strand mapping; it is an early screen for scaffold organization before more specific formation models are added.

### Formation comparison report

The formation runner also writes a cross-structure comparison table to `outputs/metrics/formation_comparison.csv` and a human-readable summary to `outputs/reports/formation_comparison.md`.

These outputs aggregate structural metrics across the normalized models to compare the Hexaplex scaffold structure against alanine beta-sheet controls. The report highlights scaffold-distinguishing features such as GLU-rich contact motifs and coarse helical order metrics, but it does not infer a formation pathway by itself.

## Intermediate Ladder

The intermediate ladder workflow creates candidate assembly models by dividing the scaffold-only complement into contiguous residue-order blocks. The initial map is written to `inputs/metadata/strand_map_candidate.csv` and is a candidate block map, not a validated PyMOL strand map. It should be replaced or revised when the true colored strand paths are available.

Run the ladder workflow with:

```bash
bash scripts/run_intermediate_ladder_metrics.sh
```

This writes candidate scaffold-only and hexads-plus-scaffold block structures under `outputs/intermediates/ladder_structures/`, per-model metrics under `outputs/metrics/intermediate_ladder/`, a ladder summary at `outputs/metrics/intermediate_ladder_summary.csv`, and a human-readable report at `outputs/reports/intermediate_ladder_report.md`.

These candidate intermediates are for testing when Hexaplex-like scaffold order appears in structural subsets. They do not prove temporal assembly order by themselves.

## Block Contact Decomposition

The candidate ladder showed that one contiguous 30-residue scaffold block already spans nearly the full circumference and axial length of the scaffold. This means the candidate blocks are not simple angular wedge units; they may correspond to long folded or twisted scaffold paths.

The block contact decomposition asks what multi-block assembly adds beyond that individual-block helical geometry. It separates contacts into within-block scaffold contacts, between-block scaffold contacts, hexad/other internal contacts, scaffold-hexad/other contacts, and GLU-rich inter-block contact patterns.

In full structures where scaffold and hexad atoms share residue labels, `hexad_or_other` is inferred from base-like atom names such as `N1`, `C2`, `N3`, `C4`, `N5`, `C6`, `OC2`, `OC4`, and `OC6`. This is an operational annotation for comparison, not a final chemical assignment.

Run the analysis with:

```bash
bash scripts/run_block_contact_analysis.sh
```

Outputs are written under `outputs/metrics/block_contacts/`, with a human-readable report at `outputs/reports/block_contact_analysis.md`. The block map is still a candidate contiguous-residue map and requires validation against PyMOL colored strand paths.

## Fitted-Axis Helical Metrics

The original `scripts/helical_order.py` uses a simple global z-axis approximation. The fitted-axis workflow in `scripts/helical_order_fitted_axis.py` estimates a principal axis from residue centroids, then computes residue radii, angles, axial span, approximate turns, and pitch-like descriptors around that fitted axis.

This is a better geometric test of whether candidate blocks behave like folded or twisted scaffold paths, while still remaining a descriptor rather than proof of a physical helix or assembly pathway.

Run fitted-axis metrics with:

```bash
bash scripts/run_fitted_helical_metrics.sh
```

Outputs are written under `outputs/metrics/fitted_helical_order/`, with a human-readable report at `outputs/reports/fitted_helical_order_report.md`.

## Ladder Diffraction Approximation

The ladder diffraction workflow computes simplified Debye-style isotropic scattering profiles for heavy-deduped ladder structures. This is a comparative approximation, not a replacement for full fiber-diffraction simulations.

The purpose is to ask whether 4-6 A window scores emerge with individual folded/twisted scaffold paths or with multi-block contact-network assembly. The d ~= 4.5 A feature remains a reciprocal-space feature and should not be interpreted as a literal atom-contact distance.

Run the workflow with:

```bash
bash scripts/run_ladder_diffraction.sh
```

The default profile method is a histogram-based Debye approximation. It bins weighted pair distances once, then evaluates the q grid from those binned counts. This keeps full heavy-deduped runs feasible more often. The direct Debye method remains available in `scripts/debye_radial_profile.py --method direct` for small structures and tests.

Full heavy-deduped runs should be preferred when feasible. For runtime-limited previews, set `SAMPLE_ATOMS` to use deterministic stratified sampling across residue/order groups:

```bash
SAMPLE_ATOMS=800 bash scripts/run_ladder_diffraction.sh
```

`--sample-atoms`/`SAMPLE_ATOMS` avoids simply taking the first atoms and preserves approximate coverage across the full atom list, but sampled outputs remain runtime-limited approximations. `MAX_ATOMS` is retained only as a legacy fallback; it uses the first atoms in PDB order after filtering. For hexads-plus-scaffold structures, this can preferentially sample hexads before scaffold atoms and should not be used for scientific interpretation.

Outputs are written under `outputs/metrics/ladder_diffraction/`, with a human-readable report at `outputs/reports/ladder_diffraction_report.md`.

## Combined Evidence Table

The combined evidence workflow joins fitted helical geometry, intermediate ladder metrics, block contact decomposition, GLU motif recurrence, and available diffraction-window scores into one CSV and Markdown report.

Run the workflow with:

```bash
bash scripts/run_combined_evidence.sh
```

Main outputs:

- `outputs/metrics/combined_evidence_table.csv`
- `outputs/reports/combined_evidence_report.md`

This table is intended to summarize evidence for the current formation mechanism: individual folded/twisted scaffold paths, multi-block GLU-rich scaffold assembly, and hexad-associated diffraction-window signals. It is not proof of temporal assembly order, and the candidate block map still requires validation against PyMOL colored strand paths.

## Seed-Formation Order Parameters

The seed-formation workflow asks a different question from the mini-hexaplex truncation analyses. The central mini-hexaplex structures show that short segments cut from the already-formed structure can remain locally helical; the seed-formation workflow instead builds synthetic rigid-body ensembles to test whether order parameters can distinguish formed-like six-strand seeds from loose arrangements of the same strand fragments.

Run the workflow with:

```bash
MPLCONFIGDIR=/tmp/hexaplex_matplotlib .venv/bin/python scripts/analyze_seed_formation_order_parameters.py --unit-counts 4,5,6,7,8 --samples-per-ensemble 25
```

Outputs are written to `outputs/seed_formation/ensembles/`, `outputs/seed_formation/plots/`, `outputs/metrics/seed_formation_order_parameters.csv`, and `outputs/reports/seed_formation_order_parameters_report.md`.

The generated ensembles treat each chain fragment as a rigid body. They are synthetic preparation for later Schrödinger bridge or formation-pathway modeling, not molecular dynamics and not evidence of spontaneous assembly by themselves.

## Seed Starting-Ensemble Cost Audit

The seed starting-ensemble cost audit compares synthetic loose-start endpoint classes against `formed_perturbed` endpoints in the current order-parameter space. It decomposes the standardized distance to the formed endpoint distribution into geometric, contact-recovery, and register/phase cost groups. The current start classes are `loose_initial`, `angular_randomized_loose_initial`, and `radially_separated`; the radial class is a dry-down/concentration-compaction proxy, not a solvent-evaporation simulation.

Run the audit with:

```bash
.venv/bin/python scripts/analyze_seed_starting_ensemble_costs.py
```

Outputs are written to `outputs/metrics/seed_starting_ensemble_cost_components.csv`, `outputs/metrics/seed_starting_ensemble_cost_summary.csv`, `outputs/plots/seed_starting_ensemble_costs/`, and `outputs/reports/seed_starting_ensemble_cost_report.md`.

This audit is an exploratory order-parameter cost proxy. It is not an atomistic Schrodinger bridge, not molecular dynamics, and not evidence of a physical nucleation pathway. Its immediate use is to identify which loose-start class is closest to the formed endpoint distribution and whether geometry, contact recovery, or register/phase terms dominate that distance.

## Seed Contact Networks

The seed contact-network workflow analyzes the formed mini-hexaplex targets directly. It asks whether the 4-, 5-, 6-, 7-, and 8-unit seeds already contain a connected six-chain interchain contact graph, and whether longer seeds add contact redundancy, CYP/MEP contact count, unit-level connectivity, or axial contact span.

Run the workflow with:

```bash
MPLCONFIGDIR=/tmp/hexaplex_matplotlib .venv/bin/python scripts/analyze_seed_contact_networks.py --unit-counts 4,5,6,7,8
```

Outputs are written to `outputs/metrics/seed_contact_network_summary.csv`, `outputs/metrics/seed_contact_network_edges.csv`, `outputs/seed_formation/plots/`, and `outputs/reports/seed_contact_network_report.md`.

The contact-network score is exploratory. It summarizes static heavy-atom interchain contact topology and should not be interpreted as proof of spontaneous nucleation, stability, or a kinetic threshold.

## Seed Bridge Feature Ordering

The seed bridge feature-ordering workflow is a Schrödinger-bridge-facing order-parameter analysis. It pairs `loose_initial` and `formed_perturbed` endpoints in standardized formedness-feature space, then interpolates order parameters to ask which features cross loose-to-formed thresholds first.

Run the workflow with:

```bash
MPLCONFIGDIR=/tmp/hexaplex_matplotlib .venv/bin/python scripts/analyze_seed_bridge_ordering.py --unit-counts 4,5,6,7 --n-time-points 21 --formed-fractions 0.50,0.75,0.80
```

Outputs are written to `outputs/metrics/seed_bridge_endpoint_pairs.csv`, `outputs/metrics/seed_bridge_order_parameter_paths.csv`, `outputs/metrics/seed_bridge_activation_summary.csv`, `outputs/metrics/seed_bridge_feature_ordering.csv`, `outputs/seed_formation/plots/`, and `outputs/reports/seed_bridge_ordering_report.md`.

The default angular phase coordinate is `refined_angular_phase_score`, which compares chain-label-aware angular positions to the formed seed after allowing one global rotation around the formed fitted axis. The bridge script can recompute this coordinate from saved sample PDBs when available, and newer seed-formation runs also write it directly into the order-parameter CSV. A fair legacy comparison can be run with `--angular-phase-mode legacy --coordinate-backed-only` and separate output paths.

This is not a full atomistic Schrödinger bridge. It uses synthetic endpoint ensembles and order-parameter interpolation only, so feature ordering should be treated as a guide for later collective-variable design rather than proof of physical assembly order. The legacy `angular_phase_order_score` should not be used as a physical ordering claim.

To test whether refined angular phase becomes informative when loose chain-label phasing is intentionally broken, generate an angular-randomized loose ensemble and rerun the refined bridge:

```bash
MPLCONFIGDIR=/tmp/hexaplex_matplotlib .venv/bin/python scripts/analyze_seed_formation_order_parameters.py --unit-counts 4,5,6,7 --samples-per-ensemble 100 --loose-mode angular_randomized
MPLCONFIGDIR=/tmp/hexaplex_matplotlib .venv/bin/python scripts/analyze_seed_bridge_ordering.py --unit-counts 4,5,6,7 --n-time-points 21 --formed-fractions 0.50,0.75,0.80 --angular-phase-mode refined --loose-ensemble-type angular_randomized_loose_initial
```

The angular-randomized bridge writes distinct `_angular_randomized` metric outputs and `outputs/reports/seed_bridge_angular_randomized_ordering_report.md`.

## Scaffold Path Map Validation

The current contiguous-residue block map is a candidate scaffold path map. It must be validated against the PyMOL colored strand/path representation before the block-specific mechanistic interpretation should be treated as map-stable.

Run the workflow with:

```bash
bash scripts/run_scaffold_path_map_workflow.sh
```

Generated files:

- `inputs/metadata/scaffold_path_map_candidate.csv`
- `inputs/metadata/scaffold_path_map_manual_template.csv`
- `outputs/reports/pymol_scaffold_path_map_helper.pml`
- `outputs/reports/scaffold_path_map_validation.md`
- `outputs/reports/scaffold_path_map_comparison.md`

User workflow:

1. Run `bash scripts/run_scaffold_path_map_workflow.sh`.
2. Open `outputs/reports/pymol_scaffold_path_map_helper.pml` in PyMOL.
3. Compare generated colors to the known colored strand/path representation.
4. Edit `inputs/metadata/scaffold_path_map_manual_template.csv` into `inputs/metadata/scaffold_path_map_manual.csv` if needed.
5. Rerun validation/comparison.

These checks validate map consistency and compare alternatives; they do not prove temporal assembly order or biological truth by themselves.

## AlphaFold/ESM Candidate Inputs

The original blank-chain PDB can be deduplicated and split into six inferred chains for candidate-only AlphaFold/ESM exploration. The residue identities are preserved in the PDB and CSV outputs as `CYP`, `MEP`, and `GLU`, but `CYP` and `MEP` are nonstandard residues and are represented as `X` in the proxy FASTA.

Run the workflow with:

```bash
bash scripts/run_ai_candidate_inputs.sh
```

Main outputs:

- `outputs/intermediates/ai_candidate_inputs/full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain.pdb`
- `outputs/intermediates/ai_candidate_inputs/full_hexaplex_anti_parallel_30deg_ideal_chain_residue_table.csv`
- `outputs/intermediates/ai_candidate_inputs/full_hexaplex_anti_parallel_30deg_ideal_chain_pattern_summary.csv`
- `outputs/intermediates/ai_candidate_inputs/full_hexaplex_anti_parallel_30deg_ideal_alphafold_esm_proxy.fasta`
- `outputs/reports/ai_candidate_inputs_report.md`

This is a proxy export for exploratory use only. It preserves the nonstandard residue identities in the structural outputs, but it does not convert the system into a standard AlphaFold/ESM-ready amino-acid sequence without loss.

## Base-Length Variant Generation

The base-length variant workflow generates local CYP/MEP arm-length sensitivity structures without running the diffraction sweep. It uses the six-chain deduplicated baseline from the AI candidate export and the transformable atom recommendations from `outputs/metrics/hexad_arm_geometry_summary.csv`.

Run the workflow with:

```bash
python3 scripts/generate_base_length_variants.py
```

Main outputs:

- `outputs/base_length_sweep/structures/hexaplex_base_length_scale_*.pdb`
- `outputs/base_length_sweep/structures/base_length_variant_manifest.csv`
- `outputs/metrics/base_length_variant_geometry.csv`
- `outputs/reports/base_length_variant_generation_report.md`

Operationally, base/hexad-arm length is the local distance from each CYP/MEP residue anchor to selected non-backbone candidate arm atoms. The transform is local: `new_xyz = anchor_xyz + scale_factor * (old_xyz - anchor_xyz)`. This keeps GLU atoms and CYP/MEP backbone-like atoms fixed and treats the generated structures as a sensitivity study, not a structural determination.

## Base-Length Diffraction Sweep

The base-length diffraction sweep runs the generated CYP/MEP arm-length variants through the sibling `../fiber-diffraction` powder-style detector workflow, then radial-averages the native `.npy` detector maps and scores q/d-spacing feature windows. The default main sweep excludes scale `1.20` because the variant geometry sanity check found suspicious heavy-atom overlaps; use it only as an explicit stress-test.

Run a small smoke sweep with:

```bash
bash scripts/run_base_length_diffraction_sweep.sh --scales 0.95,1.00,1.05
```

Run the default main sweep with:

```bash
bash scripts/run_base_length_diffraction_sweep.sh
```

Main outputs:

- `outputs/base_length_sweep/diffraction/hexaplex_base_length_scale_*.npy`
- `outputs/base_length_sweep/radial_profiles/hexaplex_base_length_scale_*_radial.csv`
- `outputs/base_length_sweep/plots/`
- `outputs/metrics/base_length_sweep_feature_summary.csv`
- `outputs/reports/base_length_sweep_report.md`

This is a computational sensitivity study: it asks whether changing local CYP/MEP base/hexad-arm length strengthens or weakens simulated intensity near d ≈ 4.5-5.0 Å using the convention `q = 2π/d`. It does not determine the structure, and it should not be read as a direct fit to experimental detector positions.

## Mini-Hexaplex Length Response

The mini-hexaplex workflow asks whether short six-strand coordinate truncations with N repeated base/GLU units per strand still show the simulated radial features of the full-length model. The cleaned full model has 15 base/GLU units per chain, so requests above 15 cannot be generated by simple truncation. The generator preserves original coordinates, residue names, atom names, chain IDs, and residue numbering; output atom serials are regenerated by the repo PDB writer.

Because the six strands are anti-parallel or overlapping from different physical ends, `firstN` is a sequence-order control rather than a compact physical segment. The workflow prioritizes:

- `lower_end_firstN_units`
- `centralN_units`

Literal `firstN` variants are available as optional controls when needed.

Run a smoke analysis for a single short variant with:

```bash
PYTHON_BIN=.venv/bin/python bash scripts/run_mini_hexaplex_analysis.sh --variants central6_units
```

Run the short-length response set with 4, 5, 6, 7, 8, and 12 units per chain:

```bash
PYTHON_BIN=.venv/bin/python bash scripts/run_mini_hexaplex_analysis.sh --unit-counts 4,5,6,7,8,12
```

Main outputs:

- `outputs/mini_hexaplex/structures/mini_hexaplex_*.pdb`
- `outputs/mini_hexaplex/mini_hexaplex_variant_manifest.csv`
- `outputs/metrics/mini_hexaplex_geometry_summary.csv`
- `outputs/mini_hexaplex/radial_profiles/*_radial.csv`
- `outputs/mini_hexaplex/plots/`
- `outputs/metrics/mini_hexaplex_feature_summary.csv`
- `outputs/reports/mini_hexaplex_length_response_report.md`

The analysis scores d = 4.5-5.0 Å, d ≈ 3.4 Å, d ≈ 3.0 Å, and the broader d = 4.1-8.4 Å region using `q = 2π/d`. These mini models are coordinate truncations, not relaxed or minimized structures, so the result should be read as a sensitivity/compatibility check rather than a structural determination.

### Mini-hexaplex helicity analysis

The helicity workflow is geometry-only and reuses existing generated mini-hexaplex PDBs. It fits the central axis from the full cleaned baseline, maps each CYP/MEP base residue representative point into cylindrical coordinates, and scores per-chain theta-vs-z linearity as a conservative helix-like coherence metric.

Run it after generating mini-hexaplex structures:

```bash
python3 scripts/analyze_mini_hexaplex_helicity.py
```

Main outputs:

- `outputs/metrics/mini_hexaplex_helicity_summary.csv`
- `outputs/mini_hexaplex/plots/mini_hexaplex_units_vs_helical_coherence_score.png`
- `outputs/mini_hexaplex/plots/mini_hexaplex_units_vs_mean_angular_residual_deg.png`
- `outputs/mini_hexaplex/plots/mini_hexaplex_units_vs_std_twist_per_unit_deg.png`
- `outputs/mini_hexaplex/plots/mini_hexaplex_helicity_vs_4p5_5p0_response.png`
- `outputs/reports/mini_hexaplex_helicity_report.md`
