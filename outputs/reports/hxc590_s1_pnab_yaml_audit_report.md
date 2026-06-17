# HXC590 S1 pNAB YAML Audit

## Purpose

This audit checks external pNAB YAML provenance for the HXC590 S1 falsification workflow.

pNAB exists as a separate local utility repository. The HXC590 branch does not depend on pNAB, and this audit does not vendor, wrap, install, or integrate pNAB.

The pNAB YAML files were inspected only for provenance. A generic pNAB example YAML is not automatically an HXC590 model-generation input.

Future pNAB work, if any, should happen outside this branch and only hand off validated candidate PDBs/profiles after a separate provenance and baseline-reproduction check.

## External pNAB location and runtime note

- Explicit pNAB repo path: `C:\Users\hpage3\OneDrive - Georgia Institute of Technology\Documents\GitHub\pnab`
- `Hexad_Antiparallel.yaml` present: yes
- `Hexad.yaml` present: yes
- pNAB import check from this Python process: no (python cannot import pnab)
- Import status is recorded only to document the local boundary; pNAB should not be installed into this branch's virtualenv as part of this audit.

No pNAB install scripts, wrappers, vendoring, or environment changes are added here.

## YAML classification

| path | classification | backbone_exists | base_files_exist | h_rise | h_twist | is_hexad | strand_orientation | num_steps |
|---|---|---|---|---|---|---|---|---|
| C:\Users\hpage3\OneDrive - Georgia Institute of Technology\Documents\GitHub\pnab\pnab\data\Hexad.yaml | generic pNAB example YAML | yes | not_applicable | 3.4;3.4;1 | 0;60;11 | True | True;True;True;True;True;True | 10000000 |
| C:\Users\hpage3\OneDrive - Georgia Institute of Technology\Documents\GitHub\pnab\pnab\data\Hexad_Antiparallel.yaml | possible antiparallel hexad YAML | yes | yes | 3.4;3.4;1 | 24;24;1 | True | True;False;True;False;True;False | 10000000 |

## Smoke-test status

No pNAB smoke test was run. Running pNAB is outside this branch-local provenance audit unless a separate external pNAB workflow first establishes provenance and reproduces the existing 30 degree baseline.

| input_yaml | status | runtime_seconds | number_of_candidates | message |
|---|---|---|---|---|
| C:\Users\hpage3\OneDrive - Georgia Institute of Technology\Documents\GitHub\pnab\pnab\data\Hexad.yaml | not_run_external_tool_not_available | 0.000001 |  | python cannot import pnab; pNAB is not required by this branch and was not installed or integrated. |
| C:\Users\hpage3\OneDrive - Georgia Institute of Technology\Documents\GitHub\pnab\pnab\data\Hexad_Antiparallel.yaml | not_run_external_tool_not_available | 0.000000 |  | python cannot import pnab; pNAB is not required by this branch and was not installed or integrated. |

## HXC590 provenance decision

1. `Hexad_Antiparallel.yaml` is present locally.
2. Location: `C:\Users\hpage3\OneDrive - Georgia Institute of Technology\Documents\GitHub\pnab\pnab\data\Hexad_Antiparallel.yaml`
3. The located `Hexad_Antiparallel.yaml` is classified as a possible antiparallel hexad YAML from the public pNAB examples, not as the actual HXC590 baseline input.
4. pNAB import status was recorded as: not available. This branch does not require pNAB.
5. The referenced backbone/base PDB files for the located hexad YAMLs are present relative to the YAML files.
6. No pNAB candidate PDB was generated because this branch is not the place to run or integrate pNAB.
7. The generic `Hexad_Antiparallel.yaml` should not be used to generate HXC590 full-length twist variants in this branch.
8. Missing provenance: a separate, provenance-clear HXC590 pNAB baseline YAML or equivalent external-builder record, component mapping for the current CYP/MEP/GLU model, and an external baseline-reproduction check.

Project-specific HXC590 baseline YAMLs found: 0.

No 30 degree pilot was generated in this branch, so no reproduction check against the existing full-length 30 degree baseline was performed.

Because no provenance-clear HXC590 pNAB input was found, no pNAB-generated twist variants should be added to this branch.

## Outputs

- `outputs\metrics\hxc590_s1_pnab_yaml_audit.csv`
- `outputs\metrics\hxc590_s1_pnab_smoke_test_results.csv`
- `outputs\reports\hxc590_s1_pnab_yaml_audit_report.md`
