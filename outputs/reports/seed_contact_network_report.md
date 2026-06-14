# Seed Contact-Network Analysis

## Purpose

This workflow analyzes the formed mini-hexaplex seed targets as interchain contact networks. It asks whether short 4-8 unit structures already contain a connected and redundant six-chain contact graph, and whether increasing length adds a plausible nucleation threshold signal.

## Contact Definition

Contacts are heavy-atom, interchain atom pairs within the requested cutoff. The primary report uses 4.5 A. Duplicate atom-index pairs are excluded. Contacts are classified as CYP/MEP vs CYP/MEP, CYP/MEP vs GLU, GLU vs GLU, backbone-like, or non-backbone. Candidate arm/core contact classes are not assigned because no stable atom-selection map is available in this workflow.

## Graph Definitions

Chain-level graph: nodes are chains A-F. An edge exists when any interchain contact is present between a chain pair, and edge weight is the number of heavy-atom contacts.

Unit-level graph: nodes are chain plus repeat-unit index. An edge exists when residues/units from different chains contact, and edge weight is the number of heavy-atom contacts.

## Exploratory Score

`nucleation_network_score` is a transparent average of chain connectivity, largest chain component size, normalized chain average degree, normalized contacts per unit, unit-graph largest component fraction, normalized CYP/MEP contact count, and coordinate-z contact span normalized to the 8-unit central reference. It is exploratory, not a validated free-energy or kinetic coordinate.

## Summary Table

| units | chain connected | contacts | contacts/unit | avg degree | unit largest frac | CYP/MEP contacts | axial span A | network score | perturb connected |
|---:|:---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 4 | True | 1032 | 258.00 | 2.00 | 1.000 | 951 | 13.63 | 0.777 | 1.000 |
| 5 | True | 1548 | 309.60 | 2.00 | 1.000 | 1413 | 18.44 | 0.845 | 1.000 |
| 6 | True | 1740 | 290.00 | 2.00 | 1.000 | 1605 | 20.43 | 0.840 | 1.000 |
| 7 | True | 2256 | 322.29 | 2.00 | 1.000 | 2067 | 25.25 | 0.892 | 1.000 |
| 8 | True | 2448 | 306.00 | 2.00 | 1.000 | 2259 | 27.22 | 0.889 | 1.000 |

## Plots

- `outputs/seed_formation/plots/seed_contacts_units_vs_total_interchain_contacts.png`
- `outputs/seed_formation/plots/seed_contacts_units_vs_contacts_per_unit.png`
- `outputs/seed_formation/plots/seed_contacts_units_vs_chain_graph_average_degree.png`
- `outputs/seed_formation/plots/seed_contacts_units_vs_unit_graph_largest_component_fraction.png`
- `outputs/seed_formation/plots/seed_contacts_units_vs_CYP_MEP_contact_count.png`
- `outputs/seed_formation/plots/seed_contacts_units_vs_nucleation_network_score.png`
- `outputs/seed_formation/plots/seed_contacts_units_vs_perturbation_connected_probability.png`
- `outputs/seed_formation/plots/seed_contacts_category_stacked_by_unit_count.png`

## Conservative Interpretation

A connected six-chain contact network in the 4-unit target means the coordinate-truncated formed seed is already topologically closed by this contact definition. Additional units should therefore be interpreted mainly through redundancy, contact density, axial span, unit-graph component growth, and perturbation robustness rather than binary chain connectivity alone.

If the 4-8 unit curves grow smoothly, the contact-network evidence supports gradual redundancy and axial-span growth rather than a sharp length threshold. A 7-unit coherent geometry flag is supported only if the 7-unit point shows a distinct increase in redundancy, unit-graph consolidation, axial span, or perturbation robustness compared with 4-6 units.

Metrics most relevant for later Schrödinger bridge modeling are contact counts by category, chain-pair weighted edges, unit-level largest component fraction, axial contact span, and perturbation contact retention because these can define progress coordinates for closure and register.

## Limitations

- Static coordinate-truncated structures.
- Contact cutoff sensitivity.
- No energy model.
- No solvent or counterions.
- No dynamics.
- Does not prove spontaneous formation or stability.
