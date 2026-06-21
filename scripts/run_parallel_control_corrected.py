#!/usr/bin/env python3
"""Run the idealized parallel-sheet control against the Emory-corrected profile."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import sys
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from hexaplex_formation.pdb_utils import dedupe_exact_atoms, heavy_atoms, load_pdb_atoms  # noqa: E402


ENGINE_DIR = ROOT / "reference" / "asem_corrected_diffraction_engine"
DEFAULT_PDB = ROOT / "inputs" / "nick_ideal_models" / "parallel_control" / "TwoBetaSheetBackbones180_180.pdb"
DEFAULT_XYZ = ROOT / "inputs" / "nick_ideal_models" / "parallel_control" / "TwoBetaSheetBackbones180_180.xyz"
DEFAULT_EXPERIMENTAL = ROOT / "inputs" / "experimental" / "nick_powder_profile_corrected_emory.csv"
DEFAULT_ANTIPARALLEL_RADIAL = (
    ROOT
    / "outputs"
    / "nick_ideal_16mer_corrected_emory_profile"
    / "radial_profiles"
    / "ideal_16mer_antiparallel_30deg_radial_profile.csv"
)
OUTPUT_ROOT = ROOT / "outputs" / "parallel_control_corrected_emory_profile"
DEFAULT_METRICS = ROOT / "outputs" / "metrics" / "parallel_control_feature_comparison_corrected_emory_profile.csv"
DEFAULT_REPORT = ROOT / "outputs" / "reports" / "parallel_control_corrected_emory_profile_report.md"

FEATURE_WINDOWS = [
    ("3.4 A", 3.4, 3.3, 3.5, "primary"),
    ("3.7 A", 3.7, 3.6, 3.8, "primary"),
    ("4.4 A", 4.4, 4.3, 4.5, "primary"),
    ("5.6 A", 5.6, 5.5, 5.7, "primary"),
    ("7.25 A", 7.25, 7.0, 7.5, "primary"),
    ("3.0 A", 3.0, 2.9, 3.1, "optional"),
    ("4.1 A", 4.1, 4.0, 4.2, "optional"),
    ("4.5-5.0 A", 4.75, 4.5, 5.0, "optional"),
    ("5.5 A", 5.5, 5.4, 5.6, "optional"),
    ("7.0 A", 7.0, 6.8, 7.2, "optional"),
    ("8.4 A", 8.4, 8.2, 8.6, "optional"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parallel-pdb", type=Path, default=DEFAULT_PDB)
    parser.add_argument("--parallel-xyz", type=Path, default=DEFAULT_XYZ)
    parser.add_argument("--experimental-profile", type=Path, default=DEFAULT_EXPERIMENTAL)
    parser.add_argument("--antiparallel-radial", type=Path, default=DEFAULT_ANTIPARALLEL_RADIAL)
    parser.add_argument("--output-root", type=Path, default=OUTPUT_ROOT)
    parser.add_argument("--metrics-output", type=Path, default=DEFAULT_METRICS)
    parser.add_argument("--report-output", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--grid-size", type=int, default=129)
    parser.add_argument("--grid-limit-mm", type=float, default=100.0)
    parser.add_argument("--radial-bins", type=int, default=420)
    parser.add_argument("--detector-distance-mm", type=float, default=338.4)
    parser.add_argument("--wavelength-angstrom", type=float, default=0.7749)
    parser.add_argument("--normalize-q-min", type=float, default=0.15)
    parser.add_argument("--plot-sigma-a", type=float, default=0.10)
    return parser.parse_args()


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load {name} from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def load_engine():
    scripts_module = load_module("scripts", ENGINE_DIR / "scripts.py")
    orientation_module = load_module("orientation_average_corrected", ENGINE_DIR / "orientation_average.py")
    radial_module = load_module("radial_average_corrected", ENGINE_DIR / "radial_average.py")
    return scripts_module, orientation_module, radial_module


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path)


def counts(atoms) -> dict[str, int]:
    return dict(sorted(Counter(atom.element.upper() for atom in atoms).items()))


def convert_pdb_to_xyz(pdb_path: Path, xyz_path: Path) -> dict[str, object]:
    atoms = load_pdb_atoms(pdb_path)
    heavy = heavy_atoms(atoms)
    deduped = dedupe_exact_atoms(heavy)
    xyz_path.parent.mkdir(parents=True, exist_ok=True)
    with xyz_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f"{len(deduped)}\n")
        handle.write("Heavy atoms only, exact deduplication applied from idealized parallel-sheet control PDB.\n")
        for atom in deduped:
            handle.write(f"{atom.element.upper():<6}{atom.x:12.6f}{atom.y:12.6f}{atom.z:12.6f}\n")
    return {
        "pdb_atom_count": len(atoms),
        "pdb_record_counts": dict(sorted(Counter(atom.record_type for atom in atoms).items())),
        "source_hydrogen_count": counts(atoms).get("H", 0),
        "source_element_counts": counts(atoms),
        "heavy_atoms_before_deduplication": len(heavy),
        "duplicate_heavy_records_removed": len(heavy) - len(deduped),
        "xyz_atom_count": len(deduped),
        "xyz_element_counts": counts(deduped),
        "hydrogens_excluded": True,
        "exact_heavy_atom_deduplication_applied": True,
    }


def load_xyz(path: Path, atomic_number: dict[str, int]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    atoms = []
    coords = []
    for line in path.read_text(encoding="utf-8").splitlines()[2:]:
        parts = line.split()
        if len(parts) < 4:
            continue
        atom = parts[0]
        atoms.append(atom)
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return (
        np.asarray(atoms, dtype=str),
        np.asarray(coords, dtype=np.float64),
        np.asarray([atomic_number[atom] for atom in atoms], dtype=np.float64),
    )


def read_profile(path: Path, intensity_column: str) -> list[dict[str, float]]:
    rows = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            d_value = float(row["d_A"])
            rows.append(
                {
                    "d_A": d_value,
                    "q_Ainv": float(row.get("q_Ainv", 2.0 * np.pi / d_value)),
                    "intensity": float(row[intensity_column]),
                }
            )
    return rows


def read_experimental(path: Path) -> list[dict[str, float]]:
    return read_profile(path, "intensity_normalized")


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def normalize_radial_rows(rows: list[dict[str, float]], normalize_q_min: float) -> list[dict[str, float]]:
    finite_rows = [
        row
        for row in rows
        if np.isfinite(row["d_A"]) and row["pixel_count"] > 0 and row["q_Ainv"] >= normalize_q_min
    ]
    max_intensity = max((row["mean_intensity"] for row in finite_rows), default=0.0)
    return [
        {
            "bin_index": row["bin_index"],
            "q_Ainv": row["q_Ainv"],
            "d_A": row["d_A"],
            "mean_intensity": row["mean_intensity"],
            "intensity_norm": row["mean_intensity"] / max_intensity if max_intensity else 0.0,
            "pixel_count": row["pixel_count"],
        }
        for row in rows
        if np.isfinite(row["d_A"]) and row["pixel_count"] > 0
    ]


def radial_rows_from_image(image: np.ndarray, args: argparse.Namespace, radial_module) -> list[dict[str, float]]:
    rows = radial_module.radial_average(
        image,
        args.grid_limit_mm,
        args.detector_distance_mm,
        args.wavelength_angstrom,
        args.radial_bins,
    )
    compact_rows = [
        {
            "bin_index": int(row["bin_index"]),
            "q_Ainv": float(row["q_Ainv"]),
            "d_A": float(row["d_A"]),
            "mean_intensity": float(row["mean_intensity"]),
            "pixel_count": int(row["pixel_count"]),
        }
        for row in rows
    ]
    return normalize_radial_rows(compact_rows, args.normalize_q_min)


def write_radial(path: Path, rows: list[dict[str, float]]) -> None:
    write_csv(
        path,
        [
            {
                "bin_index": row["bin_index"],
                "q_Ainv": f"{row['q_Ainv']:.12g}",
                "d_A": f"{row['d_A']:.12g}",
                "mean_intensity": f"{row['mean_intensity']:.12g}",
                "intensity_norm": f"{row['intensity_norm']:.12g}",
                "pixel_count": row["pixel_count"],
            }
            for row in rows
        ],
        ["bin_index", "q_Ainv", "d_A", "mean_intensity", "intensity_norm", "pixel_count"],
    )


def points(rows: list[dict[str, float]], d_min: float, d_max: float) -> list[dict[str, float]]:
    return sorted([row for row in rows if d_min <= row["d_A"] <= d_max], key=lambda row: row["d_A"])


def peak_area(rows: list[dict[str, float]], d_min: float, d_max: float) -> tuple[float | None, float | None, float]:
    local = points(rows, d_min, d_max)
    if not local:
        return None, None, 0.0
    peak = max(local, key=lambda row: row["intensity"])
    x = np.asarray([row["d_A"] for row in local], dtype=np.float64)
    y = np.asarray([row["intensity"] for row in local], dtype=np.float64)
    return float(peak["d_A"]), float(peak["intensity"]), float(np.trapezoid(y, x)) if len(local) > 1 else 0.0


def build_metrics(
    experimental: list[dict[str, float]],
    antiparallel: list[dict[str, float]],
    parallel: list[dict[str, float]],
) -> list[dict[str, object]]:
    rows = []
    for feature, center, d_min, d_max, group in FEATURE_WINDOWS:
        exp_d, exp_i, exp_area = peak_area(experimental, d_min, d_max)
        anti_d, anti_i, anti_area = peak_area(antiparallel, d_min, d_max)
        par_d, par_i, par_area = peak_area(parallel, d_min, d_max)
        anti_offset = anti_d - exp_d if anti_d is not None and exp_d is not None else None
        par_offset = par_d - exp_d if par_d is not None and exp_d is not None else None
        if anti_offset is None or par_offset is None:
            closer = ""
        elif abs(anti_offset) < abs(par_offset):
            closer = "antiparallel"
        elif abs(par_offset) < abs(anti_offset):
            closer = "parallel"
        else:
            closer = "tie"
        rows.append(
            {
                "feature_window": feature,
                "feature_group": group,
                "window_center_d_A": f"{center:.12g}",
                "window_min_d_A": f"{d_min:.12g}",
                "window_max_d_A": f"{d_max:.12g}",
                "experimental_peak_d_A": f"{exp_d:.12g}" if exp_d is not None else "",
                "experimental_peak_intensity_norm": f"{exp_i:.12g}" if exp_i is not None else "",
                "antiparallel_peak_d_A": f"{anti_d:.12g}" if anti_d is not None else "",
                "antiparallel_offset_d_A": f"{anti_offset:.12g}" if anti_offset is not None else "",
                "parallel_peak_d_A": f"{par_d:.12g}" if par_d is not None else "",
                "parallel_offset_d_A": f"{par_offset:.12g}" if par_offset is not None else "",
                "closer_model_by_abs_offset": closer,
                "antiparallel_peak_intensity_norm": f"{anti_i:.12g}" if anti_i is not None else "",
                "parallel_peak_intensity_norm": f"{par_i:.12g}" if par_i is not None else "",
                "experimental_window_area": f"{exp_area:.12g}",
                "antiparallel_window_area": f"{anti_area:.12g}",
                "parallel_window_area": f"{par_area:.12g}",
            }
        )
    return rows


def interpolate(rows: list[dict[str, float]], d_grid: np.ndarray) -> np.ndarray:
    ordered = sorted(rows, key=lambda row: row["d_A"])
    x = np.asarray([row["d_A"] for row in ordered], dtype=np.float64)
    y = np.asarray([row["intensity"] for row in ordered], dtype=np.float64)
    return np.interp(d_grid, x, y)


def gaussian_kernel(sigma: float, step: float) -> np.ndarray:
    radius = max(3, int(np.ceil(4.0 * sigma / step)))
    offsets = np.arange(-radius, radius + 1, dtype=np.float64) * step
    kernel = np.exp(-0.5 * (offsets / sigma) ** 2)
    return kernel / np.sum(kernel)


def smooth(y: np.ndarray, sigma: float, step: float) -> np.ndarray:
    kernel = gaussian_kernel(sigma, step)
    padded = np.pad(y, (len(kernel) // 2, len(kernel) // 2), mode="edge")
    return np.convolve(padded, kernel, mode="valid")


def normalize(y: np.ndarray) -> np.ndarray:
    y = np.asarray(y, dtype=np.float64)
    y = y - np.nanmin(y)
    peak = np.nanmax(y)
    return y / peak if peak > 0 else y


def plot_comparison(path: Path, experimental, antiparallel, parallel, sigma: float) -> None:
    d_grid = np.linspace(3.0, 8.0, 1600)
    step = float(np.mean(np.diff(d_grid)))
    exp_y = smooth(interpolate(experimental, d_grid), sigma, step)
    anti_y = smooth(interpolate(antiparallel, d_grid), sigma, step)
    par_y = smooth(interpolate(parallel, d_grid), sigma, step)
    fig, ax = plt.subplots(figsize=(8.5, 5.2), dpi=220, facecolor="white")
    ax.set_facecolor("white")
    for _feature, _center, d_min, d_max, group in FEATURE_WINDOWS:
        if group == "primary":
            ax.axvspan(d_min, d_max, color="lightgray", alpha=0.16, linewidth=0)
    ax.plot(d_grid, normalize(exp_y) + 2.2, color="black", linewidth=1.4, label="Corrected Experimental")
    ax.plot(d_grid, normalize(anti_y) + 1.1, color="red", linewidth=1.25, label="Ideal AntiParallel 30deg")
    ax.plot(d_grid, normalize(par_y), color="#1f77b4", linewidth=1.25, label="Parallel-sheet control")
    ax.set_xlim(8.0, 3.0)
    ax.set_ylim(-0.08, 3.35)
    ax.set_xlabel("d-spacing (A)", color="black")
    ax.set_ylabel("offset normalized intensity", color="black")
    ax.set_title("Parallel-sheet control vs antiparallel ideal model", color="black")
    ax.legend(loc="upper left", frameon=True, facecolor="white", edgecolor="none", fontsize=8)
    ax.tick_params(axis="both", colors="black")
    for spine in ax.spines.values():
        spine.set_color("black")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def markdown_table(rows: list[dict[str, object]], columns: list[str]) -> str:
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join("---" for _ in columns) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(str(row.get(column, "")) for column in columns) + " |")
    return "\n".join(lines)


def write_report(path: Path, args: argparse.Namespace, conversion: dict[str, object], metrics: list[dict[str, object]], plot_path: Path) -> None:
    primary = [row for row in metrics if row["feature_group"] == "primary"]
    anti_wins = sum(1 for row in primary if row["closer_model_by_abs_offset"] == "antiparallel")
    par_wins = sum(1 for row in primary if row["closer_model_by_abs_offset"] == "parallel")
    if anti_wins > par_wins:
        interpretation = "The feature-window comparison supports the antiparallel 30-degree ideal model over the idealized parallel-sheet control under the corrected workflow."
    elif par_wins > anti_wins:
        interpretation = "The feature-window comparison does not support treating the parallel-sheet control as clearly disfavored by peak position alone."
    else:
        interpretation = "The feature-window comparison is mixed by peak position and should be treated as a limited negative-control result."
    lines = [
        "# Parallel-Sheet Control Against Emory-Corrected Profile",
        "",
        "## Purpose",
        "",
        "Run the idealized parallel-sheet model as a control/falsification comparison against Nick's clarified ideal antiparallel 30-degree hexaplex model.",
        "",
        "## Inputs",
        "",
        f"- Parallel control PDB: `{display_path(args.parallel_pdb)}`",
        f"- Parallel control XYZ: `{display_path(args.parallel_xyz)}`",
        f"- Antiparallel baseline radial profile: `{display_path(args.antiparallel_radial)}`",
        f"- Corrected experimental profile: `{display_path(args.experimental_profile)}`",
        f"- Comparison plot: `{display_path(plot_path)}`",
        "",
        "## Conversion",
        "",
        f"- PDB ATOM/HETATM count: {conversion['pdb_atom_count']}",
        f"- PDB record counts: `{json.dumps(conversion['pdb_record_counts'], sort_keys=True)}`",
        f"- Source hydrogen count: {conversion['source_hydrogen_count']}",
        f"- Heavy atoms before deduplication: {conversion['heavy_atoms_before_deduplication']}",
        f"- Duplicate heavy records removed: {conversion['duplicate_heavy_records_removed']}",
        f"- Final XYZ atom count: {conversion['xyz_atom_count']}",
        f"- XYZ element counts: `{json.dumps(conversion['xyz_element_counts'], sort_keys=True)}`",
        "",
        "## Corrected Diffraction Settings",
        "",
        "- Engine: `reference/asem_corrected_diffraction_engine/`.",
        "- Asem correction: non-accumulating/vectorized orientation path.",
        "- Tilts: `[0]`.",
        "- Rotations: `range(0, 181, 5)`, 37 rotations.",
        f"- Grid size: {args.grid_size} x {args.grid_size}; detector half-width {args.grid_limit_mm:g} mm; detector distance {args.detector_distance_mm:g} mm; wavelength {args.wavelength_angstrom:g} A; radial bins {args.radial_bins}.",
        "",
        "## Primary Feature Windows",
        "",
        markdown_table(
            primary,
            [
                "feature_window",
                "experimental_peak_d_A",
                "antiparallel_peak_d_A",
                "antiparallel_offset_d_A",
                "parallel_peak_d_A",
                "parallel_offset_d_A",
                "closer_model_by_abs_offset",
            ],
        ),
        "",
        "## Interpretation",
        "",
        interpretation,
        "",
        f"Across the five primary windows, antiparallel is closer in {anti_wins} window(s), while parallel is closer in {par_wins} window(s). Avoid reading this as a full structural refinement; it is a falsification-style control against one idealized wrong-geometry model.",
        "",
        "## Limitations",
        "",
        "- Powder profiles are orientationally/radially averaged.",
        "- The result depends on the idealized alanine parallel-sheet control geometry.",
        "- This is a negative-control/falsification test, not a global phase or structure refinement.",
        f"- The plot uses d-space Gaussian smoothing with sigma {args.plot_sigma_a:.2f} A for visualization only; metrics use the unsmoothed radial/profile CSVs.",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    scripts_module, orientation_module, radial_module = load_engine()
    radial_dir = args.output_root / "radial_profiles"
    plot_dir = args.output_root / "plots"
    table_dir = args.output_root / "tables"
    for directory in (radial_dir, plot_dir, table_dir, args.metrics_output.parent, args.report_output.parent):
        directory.mkdir(parents=True, exist_ok=True)

    conversion = convert_pdb_to_xyz(args.parallel_pdb, args.parallel_xyz)
    atom_symbols, coords_angstrom, atomic_numbers = load_xyz(args.parallel_xyz, scripts_module.atomic_number)
    tilts = [0]
    rotations = list(range(0, 181, 5))
    image = orientation_module.average_fiber_diffraction(
        atomic_numbers,
        coords_angstrom * 1e-7,
        args.wavelength_angstrom * 1e-7,
        args.detector_distance_mm,
        [-args.grid_limit_mm, args.grid_limit_mm],
        [-args.grid_limit_mm, args.grid_limit_mm],
        args.grid_size,
        args.grid_size,
        tilts,
        rotations,
    )
    radial = radial_rows_from_image(image, args, radial_module)
    radial_path = radial_dir / "parallel_control_radial_profile.csv"
    write_radial(radial_path, radial)
    metadata = {
        "parallel_pdb": display_path(args.parallel_pdb),
        "parallel_xyz": display_path(args.parallel_xyz),
        "radial_profile": display_path(radial_path),
        "conversion": conversion,
        "engine_dir": display_path(ENGINE_DIR),
        "rotation_sampling": "range(0, 181, 5)",
        "rotations_deg": rotations,
        "tilts_deg": tilts,
        "grid_size": args.grid_size,
        "grid_limit_mm": args.grid_limit_mm,
        "radial_bins": args.radial_bins,
        "detector_distance_mm": args.detector_distance_mm,
        "wavelength_A": args.wavelength_angstrom,
        "image_min": float(np.min(image)),
        "image_max": float(np.max(image)),
        "image_mean": float(np.mean(image)),
        "model_atoms": len(atom_symbols),
    }
    (table_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    experimental = read_experimental(args.experimental_profile)
    antiparallel = read_profile(args.antiparallel_radial, "intensity_norm")
    parallel = read_profile(radial_path, "intensity_norm")
    metrics = build_metrics(experimental, antiparallel, parallel)
    fields = [
        "feature_window",
        "feature_group",
        "window_center_d_A",
        "window_min_d_A",
        "window_max_d_A",
        "experimental_peak_d_A",
        "experimental_peak_intensity_norm",
        "antiparallel_peak_d_A",
        "antiparallel_offset_d_A",
        "parallel_peak_d_A",
        "parallel_offset_d_A",
        "closer_model_by_abs_offset",
        "antiparallel_peak_intensity_norm",
        "parallel_peak_intensity_norm",
        "experimental_window_area",
        "antiparallel_window_area",
        "parallel_window_area",
    ]
    write_csv(args.metrics_output, metrics, fields)
    write_csv(table_dir / "parallel_control_feature_comparison_corrected_emory_profile.csv", metrics, fields)
    plot_path = plot_dir / "parallel_vs_antiparallel_corrected_emory_profile.png"
    plot_comparison(plot_path, experimental, antiparallel, parallel, args.plot_sigma_a)
    write_report(args.report_output, args, conversion, metrics, plot_path)
    print(f"Wrote parallel radial profile: {radial_path}")
    print(f"Wrote metrics: {args.metrics_output}")
    print(f"Wrote plot: {plot_path}")
    print(f"Wrote report: {args.report_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
