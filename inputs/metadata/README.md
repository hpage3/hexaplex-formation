# Structure Metadata

Strand and scaffold block maps live here.

The initial map is a candidate map based on contiguous residue order in the current scaffold-only complement PDB. It divides the scaffold complement into evenly sized residue blocks and is intended for hypothesis-testing model generation.

This map is not yet a validated biological strand map. The true PyMOL colored strand mapping may need to replace it later.

Generated intermediate structures based on these maps are candidate assembly models, not proven assembly intermediates.

## Canonical Scaffold Path Map CSV

Use `scaffold_path_map_*.csv` files for validated, manual, or alternative scaffold path assignments.

Required columns:

- `map_name`: stable name for the path map, such as `contiguous_residue_blocks` or `pymol_manual_paths`.
- `strand_id`: path/strand identifier. This must be nonblank for validation.
- `strand_label`: readable path/strand label. This must be nonblank for validation.
- `residue_index_in_pdb_order`: one-based residue index in the source PDB order.
- `chain_id`: PDB chain ID. This may be blank for structures without chain IDs.
- `residue_name`: residue name.
- `residue_number`: PDB residue number.
- `insertion_code`: PDB insertion code. This may be blank.
- `residue_label`: readable residue label, such as `GLU2` or `A:GLU2`.
- `source`: map provenance, such as `generated_from_strand_map_candidate` or `manual_pymol_inspection`.

The legacy `strand_map_candidate.csv` may still exist for compatibility. The richer canonical candidate map is `scaffold_path_map_candidate.csv`, generated from the legacy contiguous-residue block map. Manual PyMOL-derived maps should be saved as `scaffold_path_map_manual.csv` for comparison.
