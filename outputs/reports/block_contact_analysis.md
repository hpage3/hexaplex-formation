# Block Contact Analysis

## Scientific cautions

The d ~= 4.5 A diffraction feature is treated as a reciprocal-space scaffold signature in the working hypothesis. It does not imply literal 4.5 A real-space atom contacts. Contact maps are real-space proximity summaries only. Diffraction intensities from component structures are comparative controls, not additive decompositions.

The block mapping is a candidate contiguous-residue mapping and must be validated against PyMOL colored strand paths. Current ladder output shows that one candidate block already spans nearly full angular coverage and axial length, so these blocks may represent long folded or twisted scaffold paths rather than angular wedge units.

`hexad_or_other` assignments use base-like atom names such as N1/C2/N3/C4/N5/C6/OC2/OC4/OC6 when scaffold and hexad atoms share residue labels in full structures. This is an atom-level operational classification for comparison, not a final chemical annotation.

## Scaffold And Full Hexaplex Summary

| model_name | scaffold_within_block_contacts | scaffold_within_block_GLU_involved | scaffold_within_block_GLU_GLU | scaffold_between_blocks_contacts | scaffold_between_blocks_GLU_involved | scaffold_between_blocks_GLU_GLU | scaffold_hexad_or_other_contacts | scaffold_hexad_or_other_GLU_involved | scaffold_hexad_or_other_GLU_GLU | hexad_or_other_internal_contacts | hexad_or_other_internal_GLU_involved | hexad_or_other_internal_GLU_GLU |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| hexaplex_scaffold_only_complement_heavy_deduped | 342 | 258 | 84 | 171 | 90 | 45 | 0 | 0 | 0 | 0 | 0 | 0 |
| full_hexaplex_anti_parallel_30deg_ideal_heavy_deduped | 300 | 258 | 84 | 129 | 90 | 45 | 174 | 0 | 0 | 42 | 0 | 0 |

## Automatically generated notes

- Full scaffold GLU-GLU contacts are mostly within-block (84) rather than between-block (45) in the candidate map.
- Full Hexaplex adds scaffold-hexad/other contacts relative to scaffold-only, consistent with component coupling in the assembled model.
- Inter-block contacts increase as scaffold blocks are added in the ladder, consistent with repeated contact-network stabilization during multi-block assembly.
- Individual blocks appear geometrically helical, while multi-block assembly appears to add repeated contact-network stabilization; this does not prove temporal order.
