# Nick/Asem Archive Incorporation Plan

This commit incorporates the Nick/Asem powder diffraction archive as source provenance and clean inputs only. Corrected scientific reanalysis is intentionally deferred to Commit 2.

## Files Incorporated

Source provenance under `reference/nick_asem_powder_archive/`:

- `diffractionnvh1.ipynb`
- `Experimental_Powder_Diffracxtion_Pattern.rtf`
- `Hexaplex_8Hexads.xyz`
- `HexaplexOnlyBases.xyz`
- `HexaplexWithCOO.xyz`
- `scriptspowder1.py`
- `Modeling Powder Diffraction Data.docx`
- `README.md`

Clean model inputs under `inputs/nick_asem_models/`:

- `nick_hexaplex_8hexads.xyz`: 2064 atoms
- `nick_hexaplex_only_bases.xyz`: 1034 atoms
- `nick_hexaplex_with_coo.xyz`: 4620 atoms

Clean experimental input:

- `inputs/experimental/nick_powder_profile.csv`

## Experimental Profile

- Rows: 442
- d-spacing range: 2.20036382530154 A to 9.73072852773558 A
- Maximum normalized intensity: 1.0
- Source file: `reference/nick_asem_powder_archive/Experimental_Powder_Diffracxtion_Pattern.rtf`

The CSV columns are `d_A`, `intensity_raw`, and `intensity_normalized`. Numeric two-column rows were extracted from the RTF, preserving the source d-spacing order and normalizing intensities by the maximum raw intensity.

## Proposed Feature Windows

Primary windows from Nick's questions:

- 3.4 A
- 3.7 A
- 4.4 A
- 5.6 A
- 7.25 A

Additional windows already used in the compact twist corrected workflow:

- 3.0 A
- 4.1 A
- 4.5-5.0 A
- 5.5 A
- 7.0 A
- 8.4 A

## Planned Corrected Reanalysis for Commit 2

- Compare full/with-COO vs bases-only diffraction to test whether the 3.4 A feature is dominated by stacked bases.
- Compare full/with-COO vs no-side-chain/no-COO models if suitable inputs are available or can be generated, especially for the broad 7.25 A feature.
- Estimate atom-group contributions for the 3.7 A, 4.4 A, and 5.6 A feature windows, with caution that diffraction intensities are not simply additive.
- Test length response using 4, 8, 16, and longer stacked assemblies if models can be generated.
- Test twist specificity using the existing compact twist variants.
- Test antiparallel specificity if suitable parallel or alternative-orientation control models exist.
- Compare beta-sheet-like backbone constraints and nonuniform backbone spacing against the corrected diffraction features.

Experimental oriented/fiber-like features should be compared by q/d-spacing and relative feature trends, not by raw 2D detector shape.

Corrected analysis should use the Asem-corrected non-accumulating/vectorized diffraction path rather than treating old notebook/script outputs as final.
