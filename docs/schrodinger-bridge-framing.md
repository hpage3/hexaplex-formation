# Schrödinger Bridge Framing

## Purpose

Schrödinger bridge / stochastic optimal transport methods may help generate plausible stochastic formation pathways between an initial ensemble and the final Hexaplex-like structure.

They should be used as hypothesis-generation tools, not as proof of the actual physical assembly pathway.

## Conservative framing

Initial distribution:

- separated or weakly associated scaffold/hexad components
- loose peptide strands
- partial scaffold assemblies
- hexads-only stacks

Final distribution:

- observed Hexaplex-like scaffold geometry
- final contact map
- final pair-distance distribution
- final helical order parameters
- final diffraction features

Bridge objective:

Infer low-action stochastic paths that connect initial and final distributions while respecting structural and chemical constraints.

## Candidate representations

### 1. Contact-map bridge

State:

- which residue-residue contacts exist
- which inter-strand contacts exist
- which GLU-rich motifs are formed
- which scaffold closure contacts are present

Best first option.

### 2. Pair-distance-matrix bridge

State:

- selected residue centroid distances
- GLU side-chain distances
- backbone carbonyl / side-chain oxygen distances
- hexad center distances

Useful for connecting geometry to diffraction-like pair organization.

### 3. Torsion-angle bridge

State:

- peptide backbone torsions
- side-chain torsions
- strand-level curvature/twist variables

More chemically realistic, but harder.

### 4. Coordinate bridge

State:

- coarse-grained or atomistic coordinates

Most ambitious and highest risk. Should not be first.

## Constraints

Possible constraints:

- chain connectivity
- excluded volume
- approximate contact map
- final PDB geometry
- GLU-rich motif recurrence
- hexad stacking
- helical order parameters
- diffraction feature agreement
- electrostatic plausibility
- hydrogen-bond plausibility
- ion/water mediation hypotheses

## Avoid overclaiming

Acceptable language:

"We use Schrödinger bridge methods to generate and rank plausible low-action assembly pathways under chosen structural constraints."

Avoid:

"The Schrödinger bridge reveals the true Hexaplex formation mechanism."
