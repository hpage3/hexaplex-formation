# Asem Peptide-Dihedral Candidate Screen

## Purpose

This note audits Asem's peptide-dihedral candidate structures in light of Nick's latest guidance. The immediate question is not whether omega can be forced to 180 degrees, but whether the original H-bond-preserving candidate `1_558214.pdb` can be modestly moved by about +5 degrees toward roughly 172-173 degrees without losing the structural features that made it attractive.

No corrected diffraction was run in this step. Nick's instruction was to stop after structural and feasibility analysis unless a clean adjusted PDB could be generated and pass geometry checks.

## Nick and Literature Context

- Nick's latest note says only `1_558214.pdb` appears to preserve good backbone-backbone H-bonding.
- Nick does not require exact 180 degree omega; a modest shift toward 172-173 degrees would be easier to defend.
- The Vitagliano omega paper is relevant because peptide omega angles in real protein structures deviate from 180 degrees, those deviations correlate with local backbone geometry, and beta-sheet regions can shift omega below 180 degrees.
- Therefore, a roughly 167.5 degree omega-related geometry is strained but not automatically disqualifying.
- Helical rise refinement near 3.38 A is a separate follow-up and was not mixed into this task.

## Atom-Naming Limitation

The candidate PDB atom names expose multiple backbone-like torsions. The canonical amide `CA-C-N'-CA'` torsion is near 180 degrees in all candidate files and therefore does not reproduce Nick's cited differences. The audit reports `N-CA-C-N'` as a Nick-note omega proxy because it is the local torsion that separates `1_558214.pdb` from the higher-omega candidates in this file set. This should be confirmed with Asem before using the value as a chemically exact omega.

## Candidate Audit

| Candidate | Nick-note omega proxy (deg) | Canonical amide omega (deg) | Min heavy contact (A) | <2.2 A contacts | Best N...O (A) | H-bond plausible | Status |
| --- | ---: | ---: | ---: | ---: | ---: | --- | --- |
| 1_1241848.pdb | 148.79 | 179.98 | 1.39 | 24 | 2.82 | False | omega_improved_hbond_unclear |
| 1_4221908.pdb | 141.31 | 179.99 | 1.41 | 30 | 2.90 | False | omega_improved_hbond_unclear |
| 1_4740437.pdb | 156.58 | 179.99 | 1.41 | 24 | 2.83 | False | omega_improved_hbond_unclear |
| 1_558214.pdb | 169.80 | 179.97 | 1.47 | 12 | 3.71 | False | original_hbond_preserving_baseline |
| 1_6276960.pdb | 175.63 | 179.98 | 1.49 | 12 | 2.89 | False | omega_improved_hbond_unclear |
| 1_853635.pdb | 179.85 | 179.97 | 1.49 | 12 | 2.84 | False | omega_improved_hbond_unclear |
| 1_9754243.pdb | 164.38 | 179.99 | 1.42 | 18 | 2.83 | False | omega_improved_hbond_unclear |

## Feasibility of a +5 Degree Adjustment

- Baseline file: `1_558214.pdb`
- Current Nick-note omega proxy: 169.80 degrees
- Target range: 172-173
- Pipeline decision: `no_adjusted_pdb_generated`

The candidate atom naming supports structural auditing, but the coordinated phi/psi-like torsion controls needed to move only the peptide omega-related geometry while preserving base geometry, rise, twist, and periodicity are not defined in this repo. A constrained modeling workflow from Asem is required.

Because this repo does not have an explicit local bonding/constraint model for coordinated phi/psi-like rotations, and because the task says not to perturb base/hexad stack geometry or helical rise/twist, the script did not generate adjusted PDB variants. A naive coordinate rotation would risk inventing chemistry and changing the wrong degrees of freedom.

## Recommendation

`1_558214.pdb` should remain the structural baseline if it is the only candidate that preserves the desired backbone-backbone H-bond. The most defensible next step is to ask Asem to attempt a constrained phi/psi-only adjustment or minimization in his modeling workflow, targeting a modest +5 degree omega shift while explicitly preserving the backbone-backbone H-bond and periodic geometry.

## Outputs

- Structural audit CSV: `outputs\metrics\asem_peptide_dihedral_candidate_structural_audit.csv`
- Feasibility CSV: `outputs\metrics\asem_1_558214_omega_adjustment_feasibility.csv`

## Checks and Scope

- Candidate PDBs were copied into `inputs/asem_peptide_dihedral_candidates/`.
- No notebooks were executed.
- Preserved source archives were not modified.
- pNAB was not used or edited.
- Corrected diffraction was not rerun.
