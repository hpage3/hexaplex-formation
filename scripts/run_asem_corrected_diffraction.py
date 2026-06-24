"""Run Asem-corrected powder diffraction serially or by orientation chunks."""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from pathlib import Path
from typing import Callable, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from reference.asem_corrected_diffraction_engine import orientation_average
from reference.asem_corrected_diffraction_engine import scripts as diffraction_engine


MANIFEST_FIELDS = [
    "model", "input_file", "status", "error", "mode", "worker_count",
    "effective_worker_count", "orientation_count", "grid_size",
    "elapsed_seconds", "npy_file", "png_file", "metadata_file",
]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--coordinate-file", type=Path, action="append", required=True,
        help="Input XYZ file. Repeat to process multiple models sequentially.",
    )
    parser.add_argument(
        "--output-dir", type=Path,
        default=Path("outputs/asem_corrected_diffraction"),
    )
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--serial", action="store_true")
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument("--max-models", type=int)
    parser.add_argument("--grid-size", type=int, default=41)
    parser.add_argument("--grid-limit-mm", type=float, default=100.0)
    parser.add_argument("--theta-count", type=int, default=4)
    parser.add_argument("--phi-count", type=int, default=8)
    parser.add_argument("--psi-count", type=int, default=3)
    parser.add_argument("--theta-max", type=float, default=180.0)
    parser.add_argument("--detector-distance-mm", type=float, default=338.4)
    parser.add_argument("--wavelength-angstrom", type=float, default=0.7749)
    parser.add_argument("--run-label", default="")
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--benchmark-csv", type=Path)
    return parser


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    args = build_parser().parse_args(argv)
    if args.workers < 1:
        raise SystemExit("--workers must be at least 1")
    if args.max_models is not None and args.max_models < 1:
        raise SystemExit("--max-models must be at least 1")
    for name in ("grid_size", "theta_count", "phi_count", "psi_count"):
        if getattr(args, name) < 1:
            raise SystemExit(f"--{name.replace('_', '-')} must be at least 1")
    return args


def analysis_mode(args: argparse.Namespace) -> tuple[str, int]:
    if args.serial or args.workers == 1:
        return "serial", 1
    return "parallel", args.workers


def configure_worker_thread_limits(mode: str) -> dict[str, str]:
    if mode != "parallel":
        return {}
    limits = {
        "NUMEXPR_NUM_THREADS": "1",
        "NUMEXPR_MAX_THREADS": "1",
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
    }
    os.environ.update(limits)
    return limits


def load_xyz(path: Path) -> tuple[np.ndarray, np.ndarray]:
    lines = path.read_text(encoding="utf-8").splitlines()
    atoms = []
    coords = []
    for line in lines[2:]:
        parts = line.split()
        if len(parts) < 4:
            continue
        atom = parts[0]
        if atom not in diffraction_engine.atomic_number:
            raise ValueError(f"Unsupported atom symbol {atom!r} in {path}")
        atoms.append(diffraction_engine.atomic_number[atom])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    if not atoms:
        raise ValueError(f"No supported atoms found in {path}")
    return np.asarray(atoms, dtype=np.float64), np.asarray(coords, dtype=np.float64)


def output_paths(output_dir: Path, model: str) -> dict[str, Path]:
    return {
        "npy": output_dir / f"{model}.npy",
        "png": output_dir / f"{model}.png",
        "metadata": output_dir / f"{model}.json",
    }


def base_row(
    path: Path,
    args: argparse.Namespace,
    mode: str,
    workers: int,
    status: str,
    error: str = "",
) -> dict[str, object]:
    paths = output_paths(args.output_dir, path.stem)
    orientation_count = args.theta_count * args.phi_count * args.psi_count
    effective_workers = 1 if mode == "serial" else min(workers, orientation_count)
    return {
        "model": path.stem,
        "input_file": str(path),
        "status": status,
        "error": error,
        "mode": mode,
        "worker_count": workers,
        "effective_worker_count": effective_workers,
        "orientation_count": orientation_count,
        "grid_size": args.grid_size,
        "elapsed_seconds": "",
        "npy_file": str(paths["npy"]),
        "png_file": str(paths["png"]),
        "metadata_file": str(paths["metadata"]),
    }


def compute_image(
    atomic_numbers: np.ndarray,
    coords_mm: np.ndarray,
    args: argparse.Namespace,
    mode: str,
    workers: int,
) -> np.ndarray:
    common = (
        atomic_numbers,
        coords_mm,
        args.wavelength_angstrom * 1e-7,
        args.detector_distance_mm,
        [-args.grid_limit_mm, args.grid_limit_mm],
        [-args.grid_limit_mm, args.grid_limit_mm],
        args.grid_size,
        args.grid_size,
        args.theta_count,
        args.phi_count,
        args.psi_count,
    )
    if mode == "serial":
        return orientation_average.average_powder_diffraction(
            *common, theta_max=args.theta_max
        )
    return orientation_average.average_powder_diffraction_parallel(
        *common, theta_max=args.theta_max, workers=workers
    )


def write_png(path: Path, image: np.ndarray, grid_limit_mm: float) -> None:
    maximum = float(np.max(image))
    fig, ax = plt.subplots(figsize=(4, 4), dpi=120)
    ax.imshow(
        image,
        extent=(-grid_limit_mm, grid_limit_mm, -grid_limit_mm, grid_limit_mm),
        cmap="gray_r",
        vmin=0.0,
        vmax=maximum * 0.01 if maximum else None,
    )
    ax.set_xlabel("x detector position, mm")
    ax.set_ylabel("z detector position, mm")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def run_model(
    path: Path,
    args: argparse.Namespace,
    mode: str,
    workers: int,
    compute_fn: Callable[..., np.ndarray] = compute_image,
) -> dict[str, object]:
    row = base_row(path, args, mode, workers, "pending")
    paths = output_paths(args.output_dir, path.stem)
    if not path.is_file():
        row["status"] = "skipped_missing_input"
        row["error"] = f"Input file does not exist: {path}"
        return row
    if args.skip_existing and all(output.is_file() for output in paths.values()):
        row["status"] = "skipped_existing"
        return row

    started = time.perf_counter()
    metadata = {
        **row,
        "run_label": args.run_label,
        "scientific_engine": "reference/asem_corrected_diffraction_engine",
        "rotation_correction": "independent non-accumulating orientations",
        "theta_count": args.theta_count,
        "phi_count": args.phi_count,
        "psi_count": args.psi_count,
        "theta_max_deg": args.theta_max,
        "grid_limit_mm": args.grid_limit_mm,
        "detector_distance_mm": args.detector_distance_mm,
        "wavelength_A": args.wavelength_angstrom,
        "worker_thread_limits": configure_worker_thread_limits(mode),
    }
    try:
        atomic_numbers, coords_angstrom = load_xyz(path)
        image = compute_fn(
            atomic_numbers, coords_angstrom * 1e-7, args, mode, workers
        )
        args.output_dir.mkdir(parents=True, exist_ok=True)
        np.save(paths["npy"], image)
        write_png(paths["png"], image, args.grid_limit_mm)
        row["status"] = "success"
        row["elapsed_seconds"] = f"{time.perf_counter() - started:.9f}"
        metadata.update(
            {
                **row,
                "atom_count": int(len(atomic_numbers)),
                "image_shape": list(image.shape),
                "image_min": float(np.min(image)),
                "image_max": float(np.max(image)),
                "image_mean": float(np.mean(image)),
            }
        )
        write_json(paths["metadata"], metadata)
    except Exception as exc:
        row["status"] = "failed"
        row["error"] = f"{type(exc).__name__}: {exc}"
        row["elapsed_seconds"] = f"{time.perf_counter() - started:.9f}"
        args.output_dir.mkdir(parents=True, exist_ok=True)
        metadata.update(row)
        write_json(paths["metadata"], metadata)
    return row


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def append_benchmark(path: Path, rows: list[dict[str, object]], run_label: str) -> None:
    fields = ["run_label", *MANIFEST_FIELDS]
    existing = []
    if path.is_file():
        with path.open("r", encoding="utf-8", newline="") as handle:
            existing = list(csv.DictReader(handle))
    write_csv(path, [*existing, *({"run_label": run_label, **row} for row in rows)], fields)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    args.output_dir = args.output_dir.resolve()
    mode, workers = analysis_mode(args)
    rows = []
    for index, path in enumerate(args.coordinate_file):
        resolved = path.resolve()
        if args.max_models is not None and index >= args.max_models:
            rows.append(base_row(resolved, args, mode, workers, "skipped_max_models"))
            continue
        row = run_model(resolved, args, mode, workers)
        rows.append(row)
        print(f"{row['model']}: {row['status']}", flush=True)

    manifest = args.manifest.resolve() if args.manifest else args.output_dir / "analysis_manifest.csv"
    write_csv(manifest, rows, MANIFEST_FIELDS)
    if args.benchmark_csv:
        append_benchmark(args.benchmark_csv.resolve(), rows, args.run_label)
    print(f"Wrote {manifest}")
    return 1 if any(row["status"] == "failed" for row in rows) else 0


if __name__ == "__main__":
    raise SystemExit(main())
