# pNAB 3.38-Rise 30-Degree Handoff Status

## Summary

No suitable pNAB-generated 3.38 A rise / 30 degree Hexaplex handoff PDB was recovered or generated from the available Asem source material.

The primary source folder inspected was:

`C:\Users\hpage3\OneDrive - Georgia Institute of Technology\Documents\GitHub\research\30`

That folder contains only:

- `30.pdb.txt`

`30.pdb.txt` is a valid PDB-like coordinate file, but it appears to be output-only. It does not contain `HEADER`, `TITLE`, `REMARK`, pNAB command metadata, source YAML/JSON/script references, job IDs, generation date, rise/twist metadata, or other embedded driver information.

Because a PDB output alone is not enough evidence to regenerate the same pNAB model with a different rise, no replacement PDB was written to:

`outputs\nick_handoff\pnab_hexaplex_twist30_rise3p38.pdb`

The older geometry-adjusted file remains distinct and should not be treated as the pNAB-regenerated model:

`outputs\nick_handoff\hexaplex_16hexad_ideal_source_twist30_rise3p38.pdb`

## Source Inventory

Targeted inspection of `research\30` found:

- `30.pdb.txt`: generated coordinate structure, 239440 bytes.
- No pNAB YAML, JSON, Python driver, notebook, shell/batch script, log, `options.yaml`, `run.py`, `initial.in`, `leaprc`, `mol2`, or `frcmod` files.
- No backbone/base component files in that folder.

Targeted sibling-folder inspection found:

- `research\29\29.pdb.txt`
- `research\30\30.pdb.txt`
- `research\31\31.pdb.txt`

The repository copy under `inputs\asem_twist_series_29_30_31\raw\30` also contains only the copied `30.pdb.txt` output. Nearby reference folders contain diffraction notebooks/scripts and pNAB library examples, but not the Asem 30-degree pNAB driver/config needed to regenerate this model.

## `30.pdb.txt` Validation

- Valid PDB-like coordinate file: yes
- Metadata lines found: none among `HEADER`, `TITLE`, `REMARK`, `AUTHOR`, `COMPND`, `SOURCE`, `KEYWDS`, or `EXPDTA`
- Atom count: 3573
- Chain IDs: blank chain only
- Residue count: 180
- Residue counts: `CYP=45`, `GLU=90`, `MEP=45`
- Inferred base/hexad layers: 15
- Adjacent intervals: 14
- Mean adjacent rise: 3.4000 A
- SHA256: `87204daf5fbad8e7a6e70518777f1a941ea1aaae97f6e87e1a0445c4dbca8e44`

This looks like the 30-degree generated structure previously imported for Asem's twist series, but the file itself does not prove the original pNAB input workflow or provide enough structured input to regenerate a chemically consistent 3.38 A version.

## pNAB Availability

pNAB is available in an existing conda environment:

`C:\Users\Public\hexaplex-tools\miniforge3\envs\pnab`

Safe checks performed:

- `conda env list` shows `pnab`
- `where python` shows system Python
- `where pnab` did not find a standalone `pnab` executable
- Python import check previously confirmed `import pnab` works in the `pnab` environment

The pNAB environment should be invoked through:

```powershell
C:\Users\Public\hexaplex-tools\miniforge3\Scripts\conda.exe run -n pnab python <script-or-command>
```

Directly running the env `python.exe` can miss Open Babel plugin setup on this machine.

## Process Check

After the interrupted pNAB attempt, no live pNAB/conda generation process remained.

The only matching Python process visible was:

- PID: `4592`
- Parent: `salt-minion.exe` PID `3188`
- Command line: unavailable/blank
- Executable path: unavailable/blank
- Input/config paths visible: none
- Output/log paths visible: none
- Assessment: not identifiable as a pNAB job; left untouched

## Handoff Decision

No handoff PDB was created because the required evidence was missing:

1. No Asem pNAB input/driver workflow was found.
2. No existing pNAB-generated 3.38 A / 30 degree PDB was found.
3. The only primary source file, `30.pdb.txt`, is output-only and does not contain enough provenance to regenerate the model.

Recommendation: ask Asem or Nick for the original pNAB driver/config/notebook or exact generation command, including backbone/base input files, strand definition, orientation settings, energy filters, and any Amber/tleap post-processing steps.

## Tests

No project code was changed for this report, so no py_compile or pytest run was required.
