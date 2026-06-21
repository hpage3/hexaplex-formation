# Ideal Antiparallel Rise Variants

Generated from `inputs\nick_ideal_models\Hexaplex_AntiParallel_30deg_Ideal.pdb`.

These are rigid-layer z-position variants for Nick's clarified ideal antiparallel 30-degree model.
The generator preserves x/y coordinates, twist geometry, atom ordering, and intralayer geometry.

Layer convention:

- The source PDB does not store explicit chain/layer IDs.
- The script infers 15 six-residue CYP/MEP base planes from residue z-centroids.
- Each base plane contains three CYP and three MEP residues.
- Non-base GLU residues are assigned to the nearest inferred base plane and moved rigidly with that plane.
- The central inferred base plane is fixed as the z anchor.

This is a controlled rise-sensitivity perturbation only. It does not adjust peptide omega, twist angle, side-chain removal logic, or perform minimization.
