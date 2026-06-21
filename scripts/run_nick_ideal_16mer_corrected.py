#!/usr/bin/env python3
"""Run Nick's clarified ideal 16-mer model against the Emory-corrected profile."""

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

from hexaplex_formation.pdb_utils import (  # noqa: E402
    dedupe_exact_atoms,
    heavy_atoms,
    load_pdb_atoms,
)


ENGINE_DIR = ROOT / "reference" / "asem_corrected_diffraction_engine"
MODEL_KEY = "ideal_16mer_antiparallel_30deg"
MODEL_LABEL = "Ideal 16-mer AntiParallel 30deg"
DEFAULT_PDB = ROOT / "inputs" / "nick_ideal_models" / "Hexaplex_AntiParallel_30deg_Ideal.pdb"
DEFAULT_XYZ = ROOT / "inputs" / "nick_ideal_models" / "Hexaplex_AntiParallel_30deg_Ideal.xyz"
DEFAULT_EXPERIMENTAL = ROOT / "inputs" / "experimental" / "nick_powder_profile_corrected_emory.csv"
DEFAULT_OUTPUT_ROOT = ROOT / "outputs" / "nick_ideal_16mer_corrected_emory_profile"
DEFAULT_FEATURE_SUMMARY = (
    ROOT / "outputs" / "metrics" / "nick_ideal_16mer_feature_summary_corrected_emory_profile.csv"
)
DEFAULT_REPORT = ROOT / "outputs" / "reports" / "nick_ideal_16mer_corrected_emory_profile_report.md"

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
    parser.add_argument("--pdb-input", type=Path, default=DEFAULT_PDB)
    parser.add_argument("--xyz-output", type=Path, default=DEFAULT_XYZ)
    parser.add_argument("--experimental-profile", type=Path, default=DEFAULT_EXPERIMENTAL)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--feature-summary-output", type=Path, default=DEFAULT_FEATURE_SUMMARY)
    parser.add_argument("--report-output", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--grid-size", type=int, default=129)
    parser.add_argument("--grid-limit-mm", type=float, default=100.0)
    parser.add_argument("--radial-bins", type=int, default=420)
    parser.add_argument("--detector-distance-mm", type=float, default=338.4)
    parser.add_argument("--wavelength-angstrom", type=float, default=0.7749)
    parser.add_argument("--normalize-q-min", type=float, default=0.15)
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


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def convert_pdb_to_heavy_deduped_xyz(pdb_path: Path, xyz_path: Path) -> dict[str, object]:
    atoms = load_pdb_atoms(pdb_path)
    heavy = heavy_atoms(atoms)
    deduped = dedupe_exact_atoms(heavy)
    element_counts_all = Counter(atom.element.upper() for atom in atoms)
    element_counts_xyz = Counter(atom.element.upper() for atom in deduped)
    xyz_path.parent.mkdir(parents=True, exist_ok=True)
    with xyz_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f"{len(deduped)}\n")
        handle.write(
            "Heavy atoms only, exact deduplication applied from "
            "Hexaplex_AntiParallel_30deg_Ideal.pdb for corrected ideal 16-mer run.\n"
        )
        for atom in deduped:
            element = atom.element.upper()
            handle.write(f"{element:<6}{atom.x:12.6f}{atom.y:12.6f}{atom.z:12.6f}\n")
    return {
        "pdb_atom_count": len(atoms),
        "pdb_record_counts": dict(sorted(Counter(atom.record_type for atom in atoms).items())),
        "pdb_element_counts": dict(sorted(element_counts_all.items())),
        "heavy_atom_count_before_dedup": len(heavy),
        "xyz_atom_count": len(deduped),
        "xyz_element_counts": dict(sorted(element_counts_xyz.items())),
        "hydrogens_in_pdb": element_counts_all.get("H", 0),
        "hydrogens_in_xyz": element_counts_xyz.get("H", 0),
        "hydrogens_included": False,
        "exact_deduplication_applied": True,
        "duplicate_heavy_atom_count_removed": len(heavy) - len(deduped),
    }


def load_xyz(path: Path, atomic_number: dict[str, int]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lines = path.read_text(encoding="utf-8").splitlines()
    atoms = []
    coords = []
    for line in lines[2:]:
        parts = line.split()
        if len(parts) < 4:
            continue
        atom = parts[0]
        if atom not in atomic_number:
            raise ValueError(f"Unsupported atom symbol {atom!r} in {path}")
        atoms.append(atom)
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return (
        np.asarray(atoms, dtype=str),
        np.asarray(coords, dtype=np.float64),
        np.asarray([atomic_number[atom] for atom in atoms], dtype=np.float64),
    )


def read_experimental_profile(path: Path) -> list[dict[str, float]]:
    rows = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            d_value = float(row["d_A"])
            rows.append(
                {
                    "d_A": d_value,
                    "q_Ainv": 2.0 * np.pi / d_value,
                    "intensity_norm": float(row["intensity_normalized"]),
                    "intensity_raw": float(row["intensity_raw"]),
                }
            )
    return rows


def normalize_radial_rows(rows: list[dict[str, float]], normalize_q_min: float) -> list[dict[str, float]]:
    finite_rows = [
        row
        for row in rows
        if np.isfinite(row["d_A"]) and row["pixel_count"] > 0 and row["q_Ainv"] >= normalize_q_min
    ]
    max_intensity = max((row["mean_intensity"] for row in finite_rows), default=0.0)
    normalized = []
    for row in rows:
        if not np.isfinite(row["d_A"]) or row["pixel_count"] <= 0:
            continue
        normalized.append(
            {
                "bin_index": row["bin_index"],
                "q_Ainv": row["q_Ainv"],
                "d_A": row["d_A"],
                "mean_intensity": row["mean_intensity"],
                "intensity_norm": row["mean_intensity"] / max_intensity if max_intensity else 0.0,
                "pixel_count": row["pixel_count"],
            }
        )
    return normalized


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


def write_profile(path: Path, rows: list[dict[str, float]]) -> None:
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


def points_in_window(rows: list[dict[str, float]], d_min: float, d_max: float) -> list[dict[str, float]]:
    return [row for row in rows if d_min <= row["d_A"] <= d_max]


def profile_peak_and_area(
    rows: list[dict[str, float]],
    d_min: float,
    d_max: float,
    intensity_key: str = "intensity_norm",
) -> tuple[float | None, float | None, float]:
    points = sorted(points_in_window(rows, d_min, d_max), key=lambda row: row["d_A"])
    if not points:
        return None, None, 0.0
    peak = max(points, key=lambda row: row[intensity_key])
    x = np.asarray([row["d_A"] for row in points], dtype=np.float64)
    y = np.asarray([row[intensity_key] for row in points], dtype=np.float64)
    area = float(np.trapezoid(y, x)) if len(points) > 1 else 0.0
    return float(peak["d_A"]), float(peak[intensity_key]), area


def build_feature_summary(
    profile: list[dict[str, float]],
    experimental: list[dict[str, float]],
) -> list[dict[str, object]]:
    rows = []
    for feature, center, d_min, d_max, group in FEATURE_WINDOWS:
        sim_peak_d, sim_peak_i, sim_area = profile_peak_and_area(profile, d_min, d_max)
        exp_peak_d, exp_peak_i, exp_area = profile_peak_and_area(experimental, d_min, d_max)
        rows.append(
            {
                "model": MODEL_KEY,
                "feature_window": feature,
                "feature_group": group,
                "window_center_d_A": f"{center:.12g}",
                "window_min_d_A": f"{d_min:.12g}",
                "window_max_d_A": f"{d_max:.12g}",
                "simulated_peak_d_A": f"{sim_peak_d:.12g}" if sim_peak_d is not None else "",
                "simulated_peak_intensity_norm": f"{sim_peak_i:.12g}" if sim_peak_i is not None else "",
                "experimental_peak_d_A": f"{exp_peak_d:.12g}" if exp_peak_d is not None else "",
                "experimental_peak_intensity_norm": f"{exp_peak_i:.12g}" if exp_peak_i is not None else "",
                "peak_offset_d_A": f"{(sim_peak_d - exp_peak_d):.12g}"
                if sim_peak_d is not None and exp_peak_d is not None
                else "",
                "window_area_simulated": f"{sim_area:.12g}",
                "window_area_experimental": f"{exp_area:.12g}",
            }
        )
    return rows


def markdown_table(rows: list[dict[str, object]], columns: list[str]) -> str:
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join("---" for _ in columns) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(str(row.get(column, "")) for column in columns) + " |")
    return "\n".join(lines)


def plot_overlay(path: Path, profile: list[dict[str, float]], experimental: list[dict[str, float]]) -> None:
    fig, ax = plt.subplots(figsize=(9.5, 5.5), dpi=180)
    exp_rows = sorted([row for row in experimental if 2.8 <= row["d_A"] <= 8.4], key=lambda row: row["d_A"])
    sim_rows = sorted([row for row in profile if 2.8 <= row["d_A"] <= 8.4], key=lambda row: row["d_A"])
    ax.plot(
        [row["d_A"] for row in exp_rows],
        [row["intensity_norm"] for row in exp_rows],
        color="black",
        linewidth=1.8,
        label="Corrected Experimental",
    )
    ax.plot(
        [row["d_A"] for row in sim_rows],
        [row["intensity_norm"] + 1.15 for row in sim_rows],
        color="#2b6cb0",
        linewidth=1.35,
        label=MODEL_LABEL,
    )
    for _feature, _center, d_min, d_max, group in FEATURE_WINDOWS:
        if group == "primary":
            ax.axvspan(d_min, d_max, color="#d8d8d8", alpha=0.13)
    ax.set_xlim(8.0, 3.0)
    ax.set_xlabel("d-spacing (A)")
    ax.set_ylabel("normalized intensity; theoretical trace offset by +1.15")
    ax.set_title("Nick ideal 16-mer vs Emory-corrected experimental profile")
    ax.grid(True, alpha=0.25, linewidth=0.6)
    ax.legend(loc="upper right")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)


def write_report(
    path: Path,
    args: argparse.Namespace,
    conversion: dict[str, object],
    experimental: list[dict[str, float]],
    feature_summary: list[dict[str, object]],
    rotations: list[int],
) -> None:
    primary_rows = [row for row in feature_summary if row["feature_group"] == "primary"]
    row_3p4 = next(row for row in feature_summary if row["feature_window"] == "3.4 A")
    offset_3p4 = float(row_3p4["peak_offset_d_A"]) if row_3p4["peak_offset_d_A"] else float("nan")
    closeness = (
        "is close to"
        if np.isfinite(offset_3p4) and abs(offset_3p4) <= 0.1
        else "is offset from"
    )
    d_values = [row["d_A"] for row in experimental]
    lines = [
        "# Nick Ideal 16-mer Corrected Emory-Profile Comparison",
        "",
        "## Purpose",
        "",
        "Nick clarified that `Hexaplex_AntiParallel_30deg_Ideal.pdb` is the coordinate file to treat as the 16-mer simulation. This run establishes that ideal antiparallel 30-degree full-hexaplex baseline against the Emory-corrected experimental powder profile.",
        "",
        "This commit focuses only on the full ideal 16-mer baseline. No-side-chain/no-COO and side-chain-subtraction comparisons require a matched derived ideal no-side-chain model and should remain a separate follow-up.",
        "",
        "## Inputs",
        "",
        f"- Ideal PDB: `{display_path(args.pdb_input)}`",
        f"- Ideal XYZ: `{display_path(args.xyz_output)}`",
        f"- Corrected experimental profile: `{display_path(args.experimental_profile)}`",
        f"- Experimental rows: {len(experimental)}",
        f"- Experimental d-spacing range: {min(d_values):.12g} A to {max(d_values):.12g} A",
        "",
        "## PDB-to-XYZ Conversion",
        "",
        f"- PDB ATOM/HETATM count: {conversion['pdb_atom_count']}",
        f"- PDB record counts: `{json.dumps(conversion['pdb_record_counts'], sort_keys=True)}`",
        f"- PDB element counts: `{json.dumps(conversion['pdb_element_counts'], sort_keys=True)}`",
        f"- Heavy atoms before exact deduplication: {conversion['heavy_atom_count_before_dedup']}",
        f"- XYZ atom count: {conversion['xyz_atom_count']}",
        f"- XYZ element counts: `{json.dumps(conversion['xyz_element_counts'], sort_keys=True)}`",
        f"- Hydrogens included in diffraction XYZ: no; source PDB hydrogens: {conversion['hydrogens_in_pdb']}",
        f"- Exact deduplication applied: yes; removed {conversion['duplicate_heavy_atom_count_removed']} duplicate heavy-atom records.",
        "",
        "## Corrected Diffraction Path",
        "",
        "- Engine: `reference/asem_corrected_diffraction_engine/`.",
        "- Asem correction: azimuthal rotations are applied to independent tilted coordinate stacks rather than accumulated through the rotation loop.",
        "- Vectorized path: `make_oriented_coords` plus `generate_fiber_diffraction_series`.",
        "- Nick-style tilts: `[0]`.",
        f"- Nick-style rotations: `range(0, 181, 5)`, producing {len(rotations)} rotations from {rotations[0]} to {rotations[-1]} degrees.",
        "- The ambiguous `range(0, 5, 180)` form was not used.",
        "",
        "## Detector and Radial Settings",
        "",
        f"- Grid size: {args.grid_size} x {args.grid_size}",
        f"- Detector half-width: {args.grid_limit_mm:g} mm",
        f"- Detector distance: {args.detector_distance_mm:g} mm",
        f"- Wavelength: {args.wavelength_angstrom:g} A",
        f"- Radial bins: {args.radial_bins}",
        f"- Profile normalization: max mean intensity at q >= {args.normalize_q_min:g} A^-1.",
        "- Plotting choice: both traces are independently normalized; the theoretical trace is vertically offset by +1.15 for visual comparison.",
        "",
        "## Primary Feature Windows",
        "",
        markdown_table(
            primary_rows,
            [
                "feature_window",
                "simulated_peak_d_A",
                "simulated_peak_intensity_norm",
                "experimental_peak_d_A",
                "experimental_peak_intensity_norm",
                "peak_offset_d_A",
                "window_area_simulated",
                "window_area_experimental",
            ],
        ),
        "",
        "## Full Feature Summary",
        "",
        markdown_table(
            feature_summary,
            [
                "feature_window",
                "feature_group",
                "simulated_peak_d_A",
                "simulated_peak_intensity_norm",
                "experimental_peak_d_A",
                "experimental_peak_intensity_norm",
                "peak_offset_d_A",
                "window_area_simulated",
                "window_area_experimental",
            ],
        ),
        "",
        "## Headline Interpretation",
        "",
        f"- 3.4 A: the simulated ideal 16-mer peak at {row_3p4['simulated_peak_d_A']} A {closeness} the corrected experimental base-stacking peak at {row_3p4['experimental_peak_d_A']} A, with offset {row_3p4['peak_offset_d_A']} A.",
        "- Treat this as the corrected full ideal baseline only, not a side-chain attribution result.",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    scripts_module, orientation_module, radial_module = load_engine()
    radial_dir = args.output_root / "radial_profiles"
    plot_dir = args.output_root / "plots"
    table_dir = args.output_root / "tables"
    for directory in (radial_dir, plot_dir, table_dir, args.feature_summary_output.parent, args.report_output.parent):
        directory.mkdir(parents=True, exist_ok=True)

    if not args.pdb_input.exists():
        raise FileNotFoundError(f"Ideal PDB not found: {args.pdb_input}")
    conversion = convert_pdb_to_heavy_deduped_xyz(args.pdb_input, args.xyz_output)

    wavelength_mm = args.wavelength_angstrom * 1e-7
    z_grid_limits = [-args.grid_limit_mm, args.grid_limit_mm]
    x_grid_limits = [-args.grid_limit_mm, args.grid_limit_mm]
    tilts = [0]
    rotations = list(range(0, 181, 5))

    print(f"Running corrected diffraction for {MODEL_KEY}: {args.xyz_output}", flush=True)
    atoms, coords_angstrom, atomic_numbers = load_xyz(args.xyz_output, scripts_module.atomic_number)
    coords_mm = coords_angstrom * 1e-7
    image = orientation_module.average_fiber_diffraction(
        atomic_numbers,
        coords_mm,
        wavelength_mm,
        args.detector_distance_mm,
        z_grid_limits,
        x_grid_limits,
        args.grid_size,
        args.grid_size,
        tilts,
        rotations,
    )
    profile = radial_rows_from_image(image, args, radial_module)
    radial_path = radial_dir / "ideal_16mer_antiparallel_30deg_radial_profile.csv"
    write_profile(radial_path, profile)

    experimental = read_experimental_profile(args.experimental_profile)
    feature_summary = build_feature_summary(profile, experimental)
    feature_fields = [
        "model",
        "feature_window",
        "feature_group",
        "window_center_d_A",
        "window_min_d_A",
        "window_max_d_A",
        "simulated_peak_d_A",
        "simulated_peak_intensity_norm",
        "experimental_peak_d_A",
        "experimental_peak_intensity_norm",
        "peak_offset_d_A",
        "window_area_simulated",
        "window_area_experimental",
    ]
    write_csv(args.feature_summary_output, feature_summary, feature_fields)
    write_csv(table_dir / "nick_ideal_16mer_feature_summary_corrected_emory_profile.csv", feature_summary, feature_fields)
    metadata = {
        "model": MODEL_KEY,
        "model_label": MODEL_LABEL,
        "pdb_input": display_path(args.pdb_input),
        "xyz_output": display_path(args.xyz_output),
        "experimental_profile": display_path(args.experimental_profile),
        "radial_profile": display_path(radial_path),
        "plot": display_path(plot_dir / "nick_ideal_16mer_vs_corrected_experimental_profile.png"),
        "engine_dir": display_path(ENGINE_DIR),
        "rotation_sampling": "range(0, 181, 5)",
        "rotations_deg": rotations,
        "tilts_deg": tilts,
        "conversion": conversion,
        "grid_size": args.grid_size,
        "grid_limit_mm": args.grid_limit_mm,
        "radial_bins": args.radial_bins,
        "detector_distance_mm": args.detector_distance_mm,
        "wavelength_A": args.wavelength_angstrom,
        "image_min": float(np.min(image)),
        "image_max": float(np.max(image)),
        "image_mean": float(np.mean(image)),
        "note": "180 degrees included; compare range(0, 180, 5) in a later endpoint sensitivity check if needed.",
    }
    (table_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    plot_overlay(plot_dir / "nick_ideal_16mer_vs_corrected_experimental_profile.png", profile, experimental)
    write_report(args.report_output, args, conversion, experimental, feature_summary, rotations)

    print(f"Wrote radial profile: {radial_path}", flush=True)
    print(f"Wrote feature summary: {args.feature_summary_output}", flush=True)
    print(f"Wrote report: {args.report_output}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
