# Nick Ideal 16-mer Derived Models

`Hexaplex_AntiParallel_30deg_Ideal_no_CH2_COOH.pdb` is derived programmatically from `../Hexaplex_AntiParallel_30deg_Ideal.pdb`.

Atom-selection rule: remove atoms named `CG`, `CD`, `OE1`, `OE2`, `HG2`, and `HG3` from every `GLU` residue. This removes the terminal GLU `CH2-COOH` group while preserving the peptide backbone and `CB` atom.

The matching XYZ file uses the same corrected ideal-baseline convention: hydrogens excluded and exact heavy-atom deduplication applied.
