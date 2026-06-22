# Asem Twist Extended Stack Build Report

## Purpose

Build synthetic layer-equivalent extensions of Asem's 29, 30, and 31 degree twist models for finite-length diffraction sensitivity testing.

## Build Method

Each source PDB was grouped into consecutive six-residue layer-equivalent units. A helical axis was estimated from layer centroids by PCA. The mean adjacent-layer rise and twist were estimated from centroid projections around that axis. Complete source layers were then repeated by the model's own mean screw transform. New residue IDs and chain `A` were assigned to avoid duplicate PDB identifiers.

These are synthetic periodic/stack extensions for diffraction sensitivity only. They are not chemically searched, independently generated, minimized, or refined structures.

## Build Summary

| twist | source layers | target | status | atoms | heavy atoms | rise mean A | twist mean deg | output |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 29 | 30 | 32 | built | 3803 | 2310 | 8.78044560747 | -9.41556748095 | inputs/asem_twist_series_29_30_31/extended_stacks/asem_29deg_extended_32layer_equiv.pdb |
| 29 | 30 | 64 | built | 7604 | 4620 | 8.78044560747 | -9.41556748095 | inputs/asem_twist_series_29_30_31/extended_stacks/asem_29deg_extended_64layer_equiv.pdb |
| 29 | 30 | 100 | built | 11910 | 7220 | 8.78044560747 | -9.41556748095 | inputs/asem_twist_series_29_30_31/extended_stacks/asem_29deg_extended_100layer_equiv.pdb |
| 30 | 30 | 32 | built | 3803 | 2310 | 8.78792981703 | -9.42581573516 | inputs/asem_twist_series_29_30_31/extended_stacks/asem_30deg_extended_32layer_equiv.pdb |
| 30 | 30 | 64 | built | 7604 | 4620 | 8.78792981703 | -9.42581573516 | inputs/asem_twist_series_29_30_31/extended_stacks/asem_30deg_extended_64layer_equiv.pdb |
| 30 | 30 | 100 | built | 11910 | 7220 | 8.78792981703 | -9.42581573516 | inputs/asem_twist_series_29_30_31/extended_stacks/asem_30deg_extended_100layer_equiv.pdb |
| 31 | 32 | 32 | built | 3810 | 2310 | 9.23803464988 | -8.82808248596 | inputs/asem_twist_series_29_30_31/extended_stacks/asem_31deg_extended_32layer_equiv.pdb |
| 31 | 32 | 64 | built | 7620 | 4620 | 9.23803464988 | -8.82808248596 | inputs/asem_twist_series_29_30_31/extended_stacks/asem_31deg_extended_64layer_equiv.pdb |
| 31 | 32 | 100 | built | 11888 | 7218 | 9.23803464988 | -8.82808248596 | inputs/asem_twist_series_29_30_31/extended_stacks/asem_31deg_extended_100layer_equiv.pdb |

## Validation Notes

Built variants: 9 of 9 requested. Layer grouping was accepted only when the residue count was divisible by six. Heavy-atom XYZ files exclude hydrogens and apply exact heavy-atom deduplication.

## Limitations

- Synthetic repeated stacks may exaggerate coherent diffraction relative to disordered real fibers.
- No AmberTools/tleap, minimization, pNAB, notebooks, or new candidate search was run.
- Layer-equivalent count is used instead of claiming independently generated hexads.
