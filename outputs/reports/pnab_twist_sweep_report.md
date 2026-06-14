# pNAB Twist Sweep Blocker Report

## Purpose

Document why the requested 24 degree to 36 degree pNAB twist sweep cannot be run in this checkout yet, and list the exact inputs still needed.

## Local Inspection Summary

- Baseline full-length structure present: `outputs/intermediates/ai_candidate_inputs/full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain.pdb`
- Baseline source PDB present: `inputs/structures/full_hexaplex_anti_parallel_30deg_ideal.pdb`
- Existing pNAB helper scripts present:
  - `scripts/generate_pnab_twist_variants.py`
  - `scripts/run_pnab_twist_diffraction_sweep.sh`
  - `scripts/report_pnab_twist_integration.py`
- Existing analysis report already notes the missing current-project pNAB YAML: `outputs/reports/pnab_twist_diffraction_sensitivity_report.md`
- Official upstream pNAB examples are present under `external/pnab/pnab/data/`, including `Hexad_Antiparallel.yaml`

## Environment Status

- `python -c "import pnab; print(pnab.__file__)"` failed with `ModuleNotFoundError: No module named 'pnab'`
- `which pnab` returned nothing
- `conda env list` failed because `conda` is not on PATH

Install command documented by the upstream project:

```bash
conda create -n pnab -c conda-forge pnab
conda activate pnab
```

No environment changes were made here.

## What Is Missing

The exact current-model pNAB input set is not present in this repository. To generate the requested twist sweep safely, the following are still needed:

- The baseline pNAB YAML for the current 30 degree six-strand hexaplex model
- The backbone component file referenced by that YAML
- The base or hexad component definitions used for the CYP/MEP/GLU model
- The runtime parameters used by the baseline builder, including search settings and strand orientation flags
- Confirmation that only `h_twist` should vary while the rest of the baseline inputs remain fixed
- A way to run pNAB in a matching environment, since it is not installed in this workspace

## Blocking Reason

I cannot generate models, run the twist sweep, or validate the regenerated 30 degree model because the exact baseline pNAB builder inputs are missing locally and `pnab` is not available in the current environment.

## Next Step Required

Provide the exact baseline pNAB YAML and its referenced component files, or point me at the source location where they live. Once those inputs are available in a pNAB-enabled environment, the sweep can be built for twists 24, 26, 28, 30, 32, 34, and 36 degrees without changing any other model parameters.
