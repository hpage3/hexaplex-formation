# Asem 30-Degree pNAB Working Copy

Original source folder:

`C:\Users\hpage3\OneDrive - Georgia Institute of Technology\Documents\GitHub\research\30_A`

This directory is a working copy of the pNAB files needed to rerun Asem's
30-degree anti-parallel Hexaplex workflow without modifying the source folder.

The copied `options.yaml` differs from the source only in:

- `HelicalParameters.h_rise`: `[3.4, 3.4, 1]` to `[3.38, 3.38, 1]`

The following supplied settings are unchanged:

- `HelicalParameters.h_twist`: `[30, 30, 1]`
- `RuntimeParameters.strand`: `YYYYYYYYYYYYYYY` (15 units)
- `RuntimeParameters.build_strand`: all six strands enabled
- `RuntimeParameters.strand_orientation`: alternating `true/false`
- pNAB search settings and energy filters
- backbone and CYP/MEP component inputs

`run_original.py` preserves Asem's supplied script. `run.py` contains the same
workflow with a Windows multiprocessing entry-point guard so it can be run
through the installed conda environment:

```powershell
C:\Users\Public\hexaplex-tools\miniforge3\Scripts\conda.exe run -n pnab python run.py
```

The script runs pNAB, selects its accepted conformer, then reconstructs the
CYP/GLU/MEP residue layout and writes `fixed.pdb`.
