# Basic Diffraction Handoff For Nick

Public folder:

`C:\Users\Public\Documents\hexaplex-diffraction-handoff`

## What Nick Should Run

From PowerShell:

```powershell
cd "C:\Users\Public\Documents\hexaplex-diffraction-handoff"
python scripts\debye_radial_profile.py --pdb input_pdbs\full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain.pdb --out outputs\basic_diffraction\full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain_profile.csv --q-min 0.2 --q-max 2.5 --q-step 0.01 --method histogram
python scripts\score_radial_windows.py --profile outputs\basic_diffraction\full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain_profile.csv --out outputs\basic_diffraction\full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain_window_scores.csv
```

## 30-Degree Model

The 30-degree antiparallel candidate model is:

`C:\Users\Public\Documents\hexaplex-diffraction-handoff\input_pdbs\full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain.pdb`

It was copied from the repo's cleaned six-chain candidate model.

## Outputs

The workflow writes:

- `C:\Users\Public\Documents\hexaplex-diffraction-handoff\outputs\basic_diffraction\full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain_profile.csv`
- `C:\Users\Public\Documents\hexaplex-diffraction-handoff\outputs\basic_diffraction\full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain_profile.metadata.json`
- `C:\Users\Public\Documents\hexaplex-diffraction-handoff\outputs\basic_diffraction\full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain_window_scores.csv`

## What This Can And Cannot Conclude

This workflow helps compare whether candidate coordinate models produce diffraction-like features consistent with experimental spacing regions. It is a simplified isotropic radial profile and d-spacing window score, not a full detector simulation and not a proof of the structure.

Molecular dynamics or minimization can help test whether the 30-degree twist and approximately 3.4 A stacking geometry are mechanically plausible under a chosen force-field model. MD or minimization would still not prove the structure.

## Experimental Information That Would Improve Comparison

Calibrated experimental data would make the comparison more meaningful, especially:

- calibrated d-spacings
- q values
- detector radii or pixel positions
- beam center
- wavelength
- sample-to-detector distance
