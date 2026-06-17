# HXC590 S1 pNAB YAML Audit

## Purpose

This audit determines whether local pNAB inputs can support reproducible twist/rise candidate generation for the HXC590 S1 falsification workflow.

A generic pNAB example YAML is not automatically the HXC590 model-generation input.

A generated 30 degree pilot must reproduce the existing 30 degree baseline before using the same workflow to generate nearby twist variants.

## pNAB location and runtime

- Explicit pNAB repo path: `C:\Users\hpage3\OneDrive - Georgia Institute of Technology\Documents\GitHub\pnab`
- `Hexad_Antiparallel.yaml` present: yes
- `Hexad.yaml` present: yes
- Active HXC590 environment pNAB import: no (python cannot import pnab)
- Local install path status: no setup.py or pyproject.toml; install.bat documents a CMake/NMake build in a conda pNAB environment with OpenBabel

The safe editable install command was checked separately with `python -m pip install -e <local-pNAB-path>`; this local checkout is not an editable Python project because it has neither `setup.py` nor `pyproject.toml`.

## YAML classification

| path | classification | backbone_exists | base_files_exist | h_rise | h_twist | is_hexad | strand_orientation | num_steps |
|---|---|---|---|---|---|---|---|---|
| C:\Users\hpage3\OneDrive - Georgia Institute of Technology\Documents\GitHub\pnab\pnab\data\Hexad.yaml | generic pNAB example YAML | yes | not_applicable | 3.4;3.4;1 | 0;60;11 | True | True;True;True;True;True;True | 10000000 |
| C:\Users\hpage3\OneDrive - Georgia Institute of Technology\Documents\GitHub\pnab\pnab\data\Hexad_Antiparallel.yaml | possible antiparallel hexad YAML | yes | yes | 3.4;3.4;1 | 24;24;1 | True | True;False;True;False;True;False | 10000000 |

## Smoke-test status

| input_yaml | status | runtime_seconds | number_of_candidates | message |
|---|---|---|---|---|
| C:\Users\hpage3\OneDrive - Georgia Institute of Technology\Documents\GitHub\pnab\pnab\data\Hexad.yaml | not_run_import_failed | 0.000006 |  | python cannot import pnab |
| C:\Users\hpage3\OneDrive - Georgia Institute of Technology\Documents\GitHub\pnab\pnab\data\Hexad_Antiparallel.yaml | not_run_import_failed | 0.000005 |  | python cannot import pnab |

## HXC590 generation decision

1. `Hexad_Antiparallel.yaml` is present locally.
2. Location: `C:\Users\hpage3\OneDrive - Georgia Institute of Technology\Documents\GitHub\pnab\pnab\data\Hexad_Antiparallel.yaml`
3. The located `Hexad_Antiparallel.yaml` is classified as a possible antiparallel hexad YAML from the public pNAB examples, not as the actual HXC590 baseline input.
4. pNAB can import/run in the active environment: no.
5. The referenced backbone/base PDB files for the located hexad YAMLs are present relative to the YAML files.
6. No generic pNAB candidate PDB was generated because pNAB is not importable in the active HXC590 virtualenv.
7. It is not safe to use the generic `Hexad_Antiparallel.yaml` to generate HXC590 full-length twist variants without provenance linking it to the current HXC590/hexaplex baseline.
8. Missing inputs/provenance: the actual HXC590 pNAB baseline YAML or equivalent builder parameters, component mapping for the current CYP/MEP/GLU model, and a pNAB runtime environment that can reproduce the existing 30 degree baseline.

Project-specific HXC590 baseline YAMLs found: 0.

No 30 degree pilot was generated, so no reproduction check against the existing full-length 30 degree baseline was performed.

Full twist/rise generation is not safe yet.

## Outputs

- `outputs\metrics\hxc590_s1_pnab_yaml_audit.csv`
- `outputs\metrics\hxc590_s1_pnab_smoke_test_results.csv`
- `outputs\reports\hxc590_s1_pnab_yaml_audit_report.md`
