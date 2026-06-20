# HXC590 S1 Ideal Antiparallel 30-Degree Hexaplex Inventory

## Purpose

This note inventories `Hexaplex_AntiParallel_30deg_Ideal.pdb` as the ideal antiparallel 30-degree hexaplex model for the corrected HXC590 S1 powder compatibility screen.

Nick clarified that this is the file he meant when referring to the 16-mer simulation benchmark. For scoring and falsification outputs, the repo labels it by model provenance and geometry rather than by that informal shorthand. The file was copied unchanged into the scoring repo; no pNAB generation, structural model generation, or coordinate modification was performed.

## Inventory

- Candidate label: `ideal_antiparallel_30deg_hexaplex`
- Source path: `C:\Users\Public\hexaplex-tools\repos\fiber-diffraction-original\inputs\Hexaplex_AntiParallel_30deg_Ideal.pdb`
- Copied path: `inputs\candidates\ideal_antiparallel_30deg_hexaplex.pdb`
- SHA256: `9a880544a551b3d16f9024e7897a23af0286820ed69df9e3cb805c94113b6aca`
- File size: 578837 bytes
- Total PDB atom count: 7146
- Element counts: C:2160;H:2814;N:1170;O:1002
- Residue counts: CYP:2082;GLU:2712;MEP:2352
- Chain IDs: (blank) (1 chain entries)
- Coordinate spans: x 32.885000 A, y 32.888000 A, z 55.665000 A
- Hydrogens present: yes
- Duplicate atom/coordinate records detected: yes (3573)
- Deduped atom count by exact atom identity: 3573

## Profile

- Scorer-compatible profile path: `outputs\metrics\hxc590_s1_ideal_antiparallel_30deg_hexaplex_profile.csv`
- Existing profile source: `outputs\mini_hexaplex\radial_profiles\full_length_baseline_radial.csv`
- Existing deduped reference: `outputs\intermediates\ai_candidate_inputs\full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain.pdb`

The profile is reused from the repo's existing full-length baseline radial profile for this same original PDB family after the established deduped/full-length profile workflow. This is a falsification-style compatibility input, not a unique refined phase assignment.
