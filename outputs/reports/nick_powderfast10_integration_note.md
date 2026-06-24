# Nick Powder Profile Integration

Nick's reference implementation is preserved at
`scripts\powderfast10.py`. Its `plot_powder_diffraction` function rotates and
averages a two-dimensional diffraction intensity array over 360 angles and
writes:

- `PowderPattern.png`: the powder-averaged image
- `middlerow.txt`: the middle-row intensity profile
- `plot2.png`: the input fiber pattern rendered with the same display scaling

The adapter at `scripts\run_nick_powderfast10_profile.py` calls that function
without modifying its scientific operations. It temporarily makes the requested
output directory the working directory, which preserves Nick-compatible output
names without writing files at the repository root. It also writes
`nick_powder_profile_manifest.json` with the input path, grid limits and sizes,
display scaling, 360-rotation settings, reference function, output names, run
label, and UTC timestamp.

The current integration accepts a saved two-dimensional NumPy diffraction array,
because that is the input supported directly by Nick's plotting function. It
does not reinterpret coordinate files or replace the repo's existing
diffraction-generation and scoring definitions.

Nick's routine imports SciPy for image rotation. This repository does not
currently maintain a complete scientific dependency list in `pyproject.toml`,
so enable SciPy directly in the existing project environment:

```powershell
.\.venv\Scripts\python.exe -m pip install scipy
```

The adapter reports a clear missing-package error if SciPy is unavailable. Its
manifest records `status`, `error`, and `scipy_available` in addition to the
input, output directory, grid settings, and generated files.

Example:

```powershell
python scripts\run_nick_powderfast10_profile.py `
  --diffraction-npy outputs\example\diffraction.npy `
  --output-dir outputs\example\nick_powder_profile `
  --z-min -100 --z-max 100 --x-min -100 --x-max 100 `
  --max-intensity-scaling 0.05 `
  --run-label example
```

For the committed smoke example, first create a tiny deterministic diffraction
array:

```powershell
.\.venv\Scripts\python.exe -c "import numpy as np; from pathlib import Path; p=Path(r'outputs\nick_powderfast10_smoke'); p.mkdir(parents=True, exist_ok=True); y,x=np.mgrid[-1:1:17j,-1:1:17j]; np.save(p/'smoke_diffraction.npy', np.exp(-8*(x*x+y*y)))"
```

Then run Nick's real plotting routine through the adapter:

```powershell
.\.venv\Scripts\python.exe scripts\run_nick_powderfast10_profile.py `
  --diffraction-npy outputs\nick_powderfast10_smoke\smoke_diffraction.npy `
  --output-dir outputs\nick_powderfast10_smoke `
  --z-min -1 --z-max 1 --x-min -1 --x-max 1 `
  --max-intensity-scaling 0.25 `
  --run-label nick-powderfast10-smoke
```

Expected controlled outputs are `PowderPattern.png`, `middlerow.txt`,
`plot2.png`, and `nick_powder_profile_manifest.json`. Grid limits and display
scaling are required rather than inferred, so the adapter does not introduce
new scientific defaults.

This wrapper is the intended entry point for reproducing Nick-style powder
plots from existing diffraction arrays. Parallel execution is intentionally
deferred. That future work belongs to the Asem-corrected diffraction
computation, not Nick's plotting routine.
