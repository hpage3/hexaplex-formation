"""Run Nick's powder plotting/profile routine in a controlled directory."""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Iterator, Sequence

import numpy as np


REFERENCE_SCRIPT = Path(__file__).with_name("powderfast10.py")
POWDER_FILENAME = "PowderPattern.png"
PROFILE_FILENAME = "middlerow.txt"
FIBER_FILENAME = "plot2.png"
MANIFEST_FILENAME = "nick_powder_profile_manifest.json"


def scipy_available() -> bool:
    return importlib.util.find_spec("scipy") is not None


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Apply Nick's powderfast10 powder plotting routine to a saved "
            "two-dimensional diffraction array."
        )
    )
    parser.add_argument(
        "--diffraction-npy",
        type=Path,
        required=True,
        help="Input two-dimensional diffraction intensity array in NumPy .npy format.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Controlled directory for Nick-compatible plots, profile, and metadata.",
    )
    parser.add_argument("--z-min", type=float, required=True)
    parser.add_argument("--z-max", type=float, required=True)
    parser.add_argument("--x-min", type=float, required=True)
    parser.add_argument("--x-max", type=float, required=True)
    parser.add_argument(
        "--max-intensity-scaling",
        type=float,
        required=True,
        help="Fraction of the input maximum used as the image display maximum.",
    )
    parser.add_argument(
        "--run-label",
        help="Optional human-readable label; defaults to the input file stem.",
    )
    return parser


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    return build_parser().parse_args(argv)


def output_paths(output_dir: Path) -> dict[str, Path]:
    output_dir = Path(output_dir)
    return {
        "powder_pattern": output_dir / POWDER_FILENAME,
        "middle_row_profile": output_dir / PROFILE_FILENAME,
        "fiber_pattern": output_dir / FIBER_FILENAME,
        "manifest": output_dir / MANIFEST_FILENAME,
    }


@contextmanager
def working_directory(path: Path) -> Iterator[None]:
    previous = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


def load_reference_plotter() -> Callable[..., None]:
    """Load Nick's committed reference function without changing its source."""

    import matplotlib

    matplotlib.use("Agg")
    spec = importlib.util.spec_from_file_location("nick_powderfast10", REFERENCE_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load Nick's reference script: {REFERENCE_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(module)
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            f"Nick's reference script requires the missing Python package {exc.name!r}"
        ) from exc
    return module.plot_powder_diffraction


def run_profile(
    diffraction_npy: Path,
    output_dir: Path,
    *,
    z_limits: tuple[float, float] = (-100.0, 100.0),
    x_limits: tuple[float, float] = (-100.0, 100.0),
    max_intensity_scaling: float = 0.05,
    run_label: str | None = None,
    plotter: Callable[..., None] | None = None,
    generated_at_utc: str | None = None,
) -> dict[str, object]:
    input_path = Path(diffraction_npy).resolve()
    destination = Path(output_dir).resolve()
    if max_intensity_scaling <= 0:
        raise ValueError("max_intensity_scaling must be positive")
    if z_limits[0] >= z_limits[1] or x_limits[0] >= x_limits[1]:
        raise ValueError("grid limits must be strictly increasing")

    diffraction = np.load(input_path, allow_pickle=False)
    if diffraction.ndim != 2:
        raise ValueError(
            f"Expected a two-dimensional diffraction array, got shape {diffraction.shape}"
        )

    destination.mkdir(parents=True, exist_ok=True)
    paths = output_paths(destination)
    timestamp = generated_at_utc or datetime.now(timezone.utc).isoformat()
    metadata: dict[str, object] = {
        "status": "pending",
        "error": "",
        "scipy_available": scipy_available(),
        "run_label": run_label or input_path.stem,
        "generated_at_utc": timestamp,
        "input_file": str(input_path),
        "output_directory": str(destination),
        "input_array": {
            "shape": list(diffraction.shape),
            "dtype": str(diffraction.dtype),
        },
        "grid_settings": {
            "z_limits": list(z_limits),
            "x_limits": list(x_limits),
            "z_grid_size": int(diffraction.shape[0]),
            "x_grid_size": int(diffraction.shape[1]),
            "max_intensity_scaling": max_intensity_scaling,
        },
        "powder_settings": {
            "rotation_count": 360,
            "rotation_angles_deg": "1 through 360 inclusive",
            "rotation_function": "scipy.ndimage.rotate",
            "reshape": False,
        },
        "reference": {
            "script": "scripts/powderfast10.py",
            "function": "plot_powder_diffraction",
            "integration": "reference routine called unchanged from a controlled working directory",
        },
        "outputs": {
            "manifest": paths["manifest"].name,
        },
    }

    try:
        selected_plotter = plotter or load_reference_plotter()
        with working_directory(destination):
            selected_plotter(
                diffraction,
                z_limits,
                x_limits,
                max_intensity_scaling,
                diffraction.shape[0],
            )

        missing = [
            str(path)
            for key, path in paths.items()
            if key not in {"manifest", "fiber_pattern"} and not path.is_file()
        ]
        if missing:
            raise RuntimeError(
                "Nick's plotting routine did not create expected output(s): "
                + ", ".join(missing)
            )
    except Exception as exc:
        metadata["status"] = "failed"
        metadata["error"] = f"{type(exc).__name__}: {exc}"
        paths["manifest"].write_text(
            json.dumps(metadata, indent=2) + "\n",
            encoding="utf-8",
        )
        raise

    metadata["status"] = "success"
    metadata["outputs"] = {
        key: path.name
        for key, path in paths.items()
        if path.is_file() or key == "manifest"
    }
    paths["manifest"].write_text(
        json.dumps(metadata, indent=2) + "\n",
        encoding="utf-8",
    )
    return metadata


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    metadata = run_profile(
        args.diffraction_npy,
        args.output_dir,
        z_limits=(args.z_min, args.z_max),
        x_limits=(args.x_min, args.x_max),
        max_intensity_scaling=args.max_intensity_scaling,
        run_label=args.run_label,
    )
    print(f"Wrote Nick powder outputs to {metadata['output_directory']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
