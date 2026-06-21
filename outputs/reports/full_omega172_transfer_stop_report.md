# Full Omega172 Transfer Stop Report

## Decision

Stopped without generating a derived full PDB.

## Stop Reason

No high-confidence direct coordinate correspondence was established. Although the full model contains 90 CYP/MEP residues and the local struct2 model also contains 90 CYP/MEP residues, the best absolute-coordinate matches of shared unmoved atoms have multi-angstrom RMSD (range 4.383-50.078 A), far above the 0.10 A transfer tolerance. Copying raw omega172 displacement vectors into the full model would therefore perturb the wrong coordinate frame.

## Local Omega-Series Correspondence

- omega167/omega172 atom count equality: True
- omega167/omega172 residue count equality: True
- atom serial/name/residue correspondence: True
- moved atom names: H'':84;N'':84;O':84
- moved atom count: 252
- moved names match expected `O'`, `N''`, `H''`: True
- local max displacement: 0.2772 A
- local RMS displacement: 0.0539 A

## Full Model Inspection

- full model atom count: 7146
- full model residue count: 180
- full CYP/MEP residue count: 90
- local struct2 residue count: 90
- candidate residue mappings inspected: 90
- clean absolute-coordinate matches within 0.10 A RMSD: 0

The same residue/atom-name pattern appears in repeated form, but the coordinate frames do not match closely enough for a direct displacement transfer. A local rigid-frame/Kabsch transfer might be mathematically possible, but it would require choosing among repeated symmetric motifs and rotating displacement vectors from one local frame into another. That is outside this conservative task and would risk inventing chemistry.

## Recommendation

Ask Asem to build the full periodic omega172 model directly from `1_558214_omega172.pdb`, then run corrected diffraction on that official full model. Do not use a coordinate-transfer surrogate for scientific conclusions.

## Outputs

- `outputs/metrics/full_omega172_transfer_correspondence_audit.csv`
- `outputs/reports/full_omega172_transfer_stop_report.md`
- `outputs/reports/full_omega172_transfer_model_report.md`
