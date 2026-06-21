# Nick/Asem Powder Diffraction Archive

Copied on: 2026-06-21

Primary source path: `C:\Users\hpage3\OneDrive - Georgia Institute of Technology\Hud Lab\fast`

Additional source path for `Modeling Powder Diffraction Data.docx`: `C:\Users\hpage3\OneDrive - Georgia Institute of Technology\Hud Lab`

This folder preserves Nick/Asem powder diffraction source files as provenance. The copied source files are intentionally preserved unchanged; cleaned input copies for analysis live outside this folder under `inputs/`.

## Copied Files

- `diffractionnvh1.ipynb`: source notebook for computing diffraction patterns from the supplied models.
- `Experimental_Powder_Diffracxtion_Pattern.rtf`: Nick's captured experimental powder diffraction numeric profile; the filename is preserved with the original spelling.
- `Hexaplex_8Hexads.xyz`: supplied 8-hexad coordinate model.
- `HexaplexOnlyBases.xyz`: supplied bases-only coordinate model.
- `HexaplexWithCOO.xyz`: supplied full/with-COO coordinate model.
- `scriptspowder1.py`: supplied powder diffraction script.
- `Modeling Powder Diffraction Data.docx`: Nick's modeling notes and scientific agenda document.

## Scientific Questions

Nick's archive frames the next corrected analysis around these questions:

- Which parts of the structure give rise to which diffraction features?
- Does the 3.4 A feature come essentially from stacked bases?
- Does subtraction of no-Glu/no-COO from the full model suggest the broad 7.25 A feature is backbone-derived?
- Do the 3.7 A, 4.4 A, and 5.6 A features have comparable backbone and Glu side-chain contributions?
- Does the pattern support beta-sheet-like backbone conformations?
- Does the pattern reveal nonuniform backbone spacing, with neighboring backbones alternating between H-bonding distance and too far apart for H-bonding?
- How specifically does the diffraction constrain antiparallel strand orientation and helical twist?
- How many stacked assemblies are required before the simulated pattern stops changing?

## Correction Context

Nick's email and the source archive predate Asem's rotation correction. Treat these files as source provenance and scientific agenda, not as final corrected simulated patterns.

Future analysis should use the Asem-corrected non-accumulating/vectorized diffraction path. In particular, azimuthal rotations should be applied independently to the tilted coordinate stack rather than repeatedly rotating an already rotated coordinate set.

The preserved notebook/code also contains a Python `range()` ambiguity: `rotations = range(0, 5, 180)` produces only `0` in Python, rather than 0-180 degrees in 5-degree steps. The archive copy is intentionally left unchanged; future corrected analysis should document the interpreted rotation sampling explicitly.
