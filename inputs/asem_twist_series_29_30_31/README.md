# Asem Twist Series 29/30/31

Imported on 2026-06-22.

## Source Paths

- `C:\Users\hpage3\OneDrive - Georgia Institute of Technology\Documents\GitHub\research\29`
- `C:\Users\hpage3\OneDrive - Georgia Institute of Technology\Documents\GitHub\research\30`
- `C:\Users\hpage3\OneDrive - Georgia Institute of Technology\Documents\GitHub\research\31`

The source folders were copied verbatim into `inputs/asem_twist_series_29_30_31/raw/29`, `raw/30`, and `raw/31`. Filenames and folder structure are preserved.

## Asem's Note

Asem attached models for 29, 30, and 31 degree twist. He also attached input files for right-handed and left-handed structures and for parallel and antiparallel structures. To create alternative structures, edit `options.yaml` to adjust helical parameters and energy thresholds, execute `run.py` to generate candidate structures, and run `tleap -f initial.in` to build missing carboxylate groups. Asem notes that `tleap` is available through AmberTools. He suggests testing multiple different structures at 30 degree twist to see if the diffraction pattern stays the same. Relaxing the energy filters may allow multiple candidates with the desired backbone-backbone hydrogen bonds.

## Scope

Folders `29`, `30`, and `31` are treated as twist-angle folders. Generation inputs, when present, may include `options.yaml`, `run.py`, `initial.in`, `leaprc`, `mol2`, `frcmod`, and other Amber/tleap-related files. This commit only imports, inventories, validates, converts, and benchmarks the provided models. It does not edit generation inputs, generate new candidates, execute notebooks, run pNAB, run minimization, or run AmberTools/tleap.
