# Current Mechanistic Summary: Hexaplex Formation

## Executive summary

The current working model is that the proposed Hexaplex forms through cooperative organization of a GLU-rich scaffold with hexad/core coupling, rather than through a simple accumulation of local 4.5 A atom contacts. The d ~= 4.5 A diffraction feature is treated as a reciprocal-space signature of organized Hexaplex scaffold geometry because it differentiates scaffold-containing Hexaplex models from isolated beta-strand / alanine controls.

The intermediate ladder changed the interpretation of the candidate scaffold blocks. A contiguous 30-residue block already spans nearly the full angular coverage and axial length of the scaffold, so these blocks are not simple angular wedge sectors. They are better treated as candidate long folded or twisted scaffold paths until validated against PyMOL colored strand mapping.

Multi-block scaffold assembly appears to add contact density and repeated GLU-rich recurrence. In the current candidate map, the scaffold-only complement has 342 within-block contacts and 171 between-block contacts, with 84 within-block GLU-GLU contacts and 45 between-block GLU-GLU contacts. This supports a model in which individual paths carry helical geometry, while multi-path association reinforces a recurring GLU-rich contact network.

The full Hexaplex adds additional component coupling. Relative to scaffold-only, the full model contains substantial scaffold-hexad/other contacts and hexad/other internal contacts under the current atom-level annotation. These contacts are consistent with cooperative scaffold/core stabilization, but they do not by themselves prove a temporal assembly order.

## Evidence table

| comparison model | role in comparison | atom/residue scale | GLU-rich contacts | angular/helical implication | mechanistic implication |
| --- | --- | --- | --- | --- | --- |
| alanine beta-sheet controls | Negative controls for GLU-rich scaffold behavior | Small beta-sheet controls; GLU-free | Zero GLU motif counts in formation comparison | Low angular coverage relative to Hexaplex scaffold models | Supports the claim that GLU-rich scaffold metrics distinguish Hexaplex-like organization from isolated alanine beta-strand controls |
| scaffold block 1 candidate | Candidate individual scaffold path | 30 residues; 226 heavy atoms | 14 GLU-GLU motifs; 43 GLU-any motifs | Angular coverage about 6.07 rad; z-span about 49.4 A | Indicates an individual candidate block is already a long helical/twisted path, not a partial circumference sector |
| scaffold-only complement | Organized scaffold without hexads | 180 residues; 1356 heavy atoms | 129 GLU-GLU motifs; 348 GLU-any motifs; 84 within-block and 45 between-block GLU-GLU contacts | Angular coverage about 6.08 rad | Multi-path scaffold association adds dense recurring contacts and GLU-rich contact-network reinforcement |
| hexads-only | Component/control for core contribution | 90 residues; 810 heavy atoms | No GLU motifs in current formation comparison | Angular coverage about 6.28 rad from z-axis screen, but not a GLU-rich scaffold path | Comparative control for core/hexad geometry; not an additive decomposition of full intensity |
| full Hexaplex | Assembled scaffold plus hexad/core model | 180 residues; 2166 heavy atoms | 129 GLU-GLU motifs; 348 GLU-any motifs; scaffold-hexad/other contacts: 174 | Angular coverage about 6.16 rad | Adds scaffold-hexad/other coupling and core packing to the GLU-rich scaffold contact network |

## Current working model

1. Individual folded or twisted scaffold paths: candidate contiguous blocks appear to be long scaffold paths with near-full angular coverage and substantial axial span.
2. Multi-path scaffold association: adding blocks increases contact density and creates between-block scaffold contacts.
3. GLU-rich contact-network reinforcement: GLU-any and GLU-GLU motifs recur across the scaffold, with both within-block and between-block contributions.
4. Scaffold-hexad/other coupling: the full Hexaplex adds contacts between scaffold atoms and base-like hexad/other atoms, plus internal hexad/other contacts.
5. Final cooperative Hexaplex stabilization: the most plausible current model is cooperative stabilization from preorganized scaffold paths, GLU-rich contact recurrence, and core/hexad coupling.

## What the 4.5 A feature means and does not mean

The d ~= 4.5 A feature is currently interpreted as a likely reciprocal-space signature of organized Hexaplex scaffold geometry. It appears useful because it differentiates Hexaplex scaffold-containing models from isolated beta-strand / alanine controls.

It does not mean that the responsible atoms are separated by literal 4.5 A real-space distances. Contact maps and pair-distance histograms are real-space summaries and should not be treated as direct assignments of reciprocal-space d-spacings.

It also is not a clean subtraction of scaffold intensity from full intensity. Hexads-only and scaffold-only diffraction simulations are comparative controls; scattering amplitudes combine before intensities are squared, so component intensities are not additive decompositions of the full model.

## Formation hypotheses now worth testing

- H1: Individual scaffold paths are preorganized into helical trajectories.
- H2: Multi-path assembly strengthens GLU-rich contact recurrence.
- H3: Hexads stabilize the scaffold through scaffold-hexad/other contacts and/or core packing.
- H4: The 4.5 A reciprocal-space feature grows with contact-network completeness, not merely with local peptide-pair geometry.
- H5: Hydrogen bonding, protonation state, ions, or water may be needed to explain GLU-rich stabilization.

## Next computational steps

1. Validate the candidate block map against PyMOL colored strand paths.
2. Compute block-pair and contact-network graphs from the block contact decomposition.
3. Add a fitted helical axis instead of the current global z-axis approximation.
4. Compare simulated diffraction profiles for intermediate ladder structures.
5. Add hydrogen-bond geometry with angle filters.
6. Evaluate protonation, ion mediation, and water-mediated stabilization hypotheses.
7. Only then introduce Schrödinger bridge or contact-map pathway modeling.

## Schrödinger bridge positioning

Schrödinger bridge modeling is not the next immediate step. The current priority is to validate the structural state variables and ensure that block/path definitions are biologically meaningful.

The best first abstraction is likely a contact-map or block-contact-state bridge, not an atomistic dynamics claim. It should be used to rank plausible assembly paths under chosen endpoints and constraints, not to prove true dynamics.

Possible bridge state variables include:

- within-block contacts
- between-block contacts
- scaffold-hexad contacts
- GLU-rich motif counts
- helical order parameters
- diffraction-window scores

## Conservative conclusion

The current data support a cooperative scaffold/hexad formation hypothesis: individual candidate scaffold paths appear geometrically helical, multi-path scaffold assembly adds recurring GLU-rich contact-network density, and the full Hexaplex adds scaffold-hexad/other coupling. This is a working mechanistic model, not a proven formation pathway. It requires validation against PyMOL strand mapping, improved helical-axis analysis, diffraction comparisons for ladder intermediates, and more chemically specific hydrogen-bond/protonation/ion modeling.
