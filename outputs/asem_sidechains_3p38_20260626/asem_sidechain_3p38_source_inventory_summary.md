# Asem 3.38 A Side-Chain Source Inventory Summary

- Source path: `C:\Users\hpage3\OneDrive - Georgia Institute of Technology\Documents\GitHub\research\tleap_structures_rise3.38`
- Import path: `inputs/asem_sidechains_3p38_20260626/raw/`
- Total PDB count: 488
- Metadata/YAML/CSV/TXT/ZIP/PDF/log files found in source: none

## Per-Angle Counts

| angle_deg | candidate_count |
| --- | ---: |
| 18 | 15 |
| 19 | 35 |
| 20 | 58 |
| 21 | 48 |
| 22 | 50 |
| 23 | 42 |
| 24 | 48 |
| 25 | 33 |
| 26 | 33 |
| 27 | 39 |
| 28 | 31 |
| 29 | 14 |
| 30 | 10 |
| 31 | 26 |
| 32 | 6 |

## Availability

- 29/30/31 present: yes (14, 10, 26)
- Initial 29/30/31 run size: 50 structures
- 28 present: yes (31 candidates)
- 32 present: yes (6 candidates)

## Representative Checks

| sample | atom_count | heavy_atom_count | element_counts | residues | has_glu_sidechain_atoms |
| --- | ---: | ---: | --- | --- | --- |
| 29/cand0/initial.pdb | 3810 | 2310 | C:960; Ca:192; H:1500; N:624; O:534 | CYP:1110; GLU:1446; MEP:1254 | True |
| 29/cand13/initial.pdb | 3810 | 2310 | C:960; Ca:192; H:1500; N:624; O:534 | CYP:1110; GLU:1446; MEP:1254 | True |
| 30/cand0/initial.pdb | 3810 | 2310 | C:960; Ca:192; H:1500; N:624; O:534 | CYP:1110; GLU:1446; MEP:1254 | True |
| 30/cand9/initial.pdb | 3810 | 2310 | C:960; Ca:192; H:1500; N:624; O:534 | CYP:1110; GLU:1446; MEP:1254 | True |
| 31/cand0/initial.pdb | 3810 | 2310 | C:960; Ca:192; H:1500; N:624; O:534 | CYP:1110; GLU:1446; MEP:1254 | True |
| 31/cand25/initial.pdb | 3810 | 2310 | C:960; Ca:192; H:1500; N:624; O:534 | CYP:1110; GLU:1446; MEP:1254 | True |

All sampled 29/30/31 structures contained GLU side-chain atom names CG, CD, OE1, and OE2. No malformed coordinate lines were found in the sampled structures.

## Comparison To Previous 3.40 A Import

| angle_deg | prior_3p40_count | new_3p38_count | delta |
| --- | ---: | ---: | ---: |
| 18 | 15 | 15 | +0 |
| 19 | 3 | 35 | +32 |
| 20 | 1 | 58 | +57 |
| 21 | 2 | 48 | +46 |
| 22 | 1 | 50 | +49 |
| 23 | 4 | 42 | +38 |
| 24 | 2 | 48 | +46 |
| 25 | 31 | 33 | +2 |
| 26 | 2 | 33 | +31 |
| 27 | 38 | 39 | +1 |
| 28 | 30 | 31 | +1 |
| 29 | 15 | 14 | -1 |
| 30 | 12 | 10 | -2 |
| 31 | 22 | 26 | +4 |
| 32 | 7 | 6 | -1 |

Specific Nick-guided starting angles:

| angle_deg | prior_3p40_count | new_3p38_count |
| --- | ---: | ---: |
| 29 | 15 | 14 |
| 30 | 12 | 10 |
| 31 | 22 | 26 |

Do not treat 3.38 candidate IDs as identical to 3.40 candidate IDs unless structure provenance later proves that mapping.

## Recommended Next Step

Prepare and review the 29/30/31 corrected-diffraction run manifest before generating XYZ files or running diffraction. The planning manifest is `outputs/asem_sidechains_3p38_20260626/sidechain_3p38_29_30_31_run_manifest.csv` and contains 50 planned rows.
