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

Nick's routine imports SciPy for image rotation. The current project `.venv`
does not include SciPy, so a real powder-profile run must use an environment
where the reference script's NumPy, Matplotlib, Pillow, and SciPy dependencies
are available. The adapter reports a clear missing-package error and does not
install or replace those dependencies.

Example:

```powershell
python scripts\run_nick_powderfast10_profile.py `
  --diffraction-npy outputs\example\diffraction.npy `
  --output-dir outputs\example\nick_powder_profile `
  --z-min -100 --z-max 100 --x-min -100 --x-max 100 `
  --max-intensity-scaling 0.05 `
  --run-label example
```

Run this command from a Python environment containing the reference script's
dependencies. Grid limits and display scaling are required rather than inferred,
so the adapter does not introduce new scientific defaults.

This wrapper is the intended entry point for reproducing Nick-style powder
plots from existing diffraction arrays. Parallel execution is intentionally
deferred to a later pipeline task.
