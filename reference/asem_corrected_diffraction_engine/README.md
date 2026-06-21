# Asem-Corrected Diffraction Engine Provenance

This folder preserves the corrected/vectorized diffraction engine files used for the Asem-corrected compact twist rerun. At the time of the rerun, the sibling `research\fiber-diffraction` folder was not under Git metadata, so these source copies are tracked here to keep the scientific provenance with the committed `hexaplex-formation` results.

The copied files are reference snapshots, not an active package:

- `scripts.py`
- `orientation_average.py`
- `run_powder_benchmark.py`
- `radial_average.py`
- `pdb_to_xyz.py`

The key rotation correction is that each azimuthal rotation must be applied to the same tilted coordinate frame. The old problematic form compounded rotations:

```python
current_coords = np.dot(Rz(rotation), current_coords.T).T
```

The corrected form rotates from the independent tilted stack each time:

```python
current_coords = np.dot(Rz(rotation), tilted_coords.T).T
```

The corrected/vectorized path also uses `numexpr`. `numexpr` was installed in `hexaplex-formation\.venv` for the rerun, but `.venv` is intentionally not committed.

Related rerun commit: `97997287f8c60d85e3b5973f9e07fb7dbb069a19`.

Corrected compact twist results preserve the main pre-fix trend: 30 degrees is strongest in the primary 3.4 A, 3.0 A, and 4.5-5.0 A feature windows. The 8.4 A window remains the exception, strongest at 24 degrees with 36 degrees close. Comparisons should be made by q/d-spacing feature windows and relative trends, not by raw 2D detector shape.
