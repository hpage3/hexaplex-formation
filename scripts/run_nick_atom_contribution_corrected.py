#!/usr/bin/env python3
"""Run Nick/Asem atom-contribution comparison with corrected diffraction code."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
ENGINE_DIR = ROOT / "reference" / "asem_corrected_diffraction_engine"
OUTPUT_ROOT = ROOT / "outputs" / "nick_atom_contribution_corrected"

MODELS = [
    ("with_coo", "With COO/full model", ROOT / "inputs" / "nick_asem_models" / "nick_hexaplex_with_coo.xyz"),
    ("only_bases", "Bases only", ROOT / "inputs" / "nick_asem_models" / "nick_hexaplex_only_bases.xyz"),
    ("eight_hexads", "8 hexads", ROOT / "inputs" / "nick_asem_models" / "nick_hexaplex_8hexads.xyz"),
]

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

DIFFERENCES = [
    ("with_coo_minus_bases_only", "with_coo", "only_bases"),
    ("with_coo_minus_8hexads", "with_coo", "eight_hexads"),
    ("8hexads_minus_bases_only", "eight_hexads", "only_bases"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid-size", type=int, default=129)
    parser.add_argument("--grid-limit-mm", type=float, default=100.0)
    parser.add_argument("--radial-bins", type=int, default=420)
    parser.add_argument("--detector-distance-mm", type=float, default=338.4)
    parser.add_argument("--wavelength-angstrom", type=float, default=0.7749)
    parser.add_argument("--normalize-q-min", type=float, default=0.15)
    parser.add_argument("--output-root", type=Path, default=OUTPUT_ROOT)
    parser.add_argument(
        "--experimental-profile",
        type=Path,
        default=ROOT / "inputs" / "experimental" / "nick_powder_profile.csv",
    )
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


def load_xyz(path: Path, atomic_number: dict[str, int]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    atoms = []
    coords = []
    with path.open("r", encoding="utf-8") as handle:
        lines = handle.read().splitlines()
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


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


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
    profiles: dict[str, list[dict[str, float]]],
    experimental: list[dict[str, float]],
) -> list[dict[str, object]]:
    rows = []
    for model_key, label, _path in MODELS:
        for feature, center, d_min, d_max, group in FEATURE_WINDOWS:
            sim_peak_d, sim_peak_i, sim_area = profile_peak_and_area(profiles[model_key], d_min, d_max)
            exp_peak_d, exp_peak_i, exp_area = profile_peak_and_area(experimental, d_min, d_max)
            rows.append(
                {
                    "model": model_key,
                    "model_label": label,
                    "feature_window": feature,
                    "feature_group": group,
                    "window_center_d_A": f"{center:.12g}",
                    "window_min_d_A": f"{d_min:.12g}",
                    "window_max_d_A": f"{d_max:.12g}",
                    "simulated_peak_d_A": f"{sim_peak_d:.12g}" if sim_peak_d is not None else "",
                    "simulated_peak_intensity_norm": f"{sim_peak_i:.12g}" if sim_peak_i is not None else "",
                    "experimental_peak_d_A": f"{exp_peak_d:.12g}" if exp_peak_d is not None else "",
                    "experimental_peak_intensity_norm": f"{exp_peak_i:.12g}" if exp_peak_i is not None else "",
                    "peak_offset_d_A": f"{(sim_peak_d - exp_peak_d):.12g}" if sim_peak_d is not None and exp_peak_d is not None else "",
                    "window_area_simulated": f"{sim_area:.12g}",
                    "window_area_experimental": f"{exp_area:.12g}",
                }
            )
    return rows


def interpolate_profile(rows: list[dict[str, float]], d_grid: np.ndarray) -> np.ndarray:
    sorted_rows = sorted(rows, key=lambda row: row["d_A"])
    x = np.asarray([row["d_A"] for row in sorted_rows], dtype=np.float64)
    y = np.asarray([row["intensity_norm"] for row in sorted_rows], dtype=np.float64)
    return np.interp(d_grid, x, y, left=np.nan, right=np.nan)


def build_difference_profiles(
    profiles: dict[str, list[dict[str, float]]],
    output_dir: Path,
) -> tuple[dict[str, list[dict[str, float]]], list[dict[str, object]]]:
    output_dir.mkdir(parents=True, exist_ok=True)
    d_grid = np.linspace(2.2, 9.8, 900)
    difference_profiles: dict[str, list[dict[str, float]]] = {}
    summary_rows = []
    for comparison, left, right in DIFFERENCES:
        left_y = interpolate_profile(profiles[left], d_grid)
        right_y = interpolate_profile(profiles[right], d_grid)
        diff_y = left_y - right_y
        valid = np.isfinite(diff_y)
        rows = [
            {
                "d_A": float(d_value),
                "q_Ainv": float(2.0 * np.pi / d_value),
                "intensity_norm": float(intensity),
            }
            for d_value, intensity in zip(d_grid[valid], diff_y[valid])
        ]
        difference_profiles[comparison] = rows
        write_csv(
            output_dir / f"{comparison}_difference_profile.csv",
            [
                {
                    "d_A": f"{row['d_A']:.12g}",
                    "q_Ainv": f"{row['q_Ainv']:.12g}",
                    "intensity_difference_norm": f"{row['intensity_norm']:.12g}",
                }
                for row in rows
            ],
            ["d_A", "q_Ainv", "intensity_difference_norm"],
        )
        for feature, center, d_min, d_max, group in FEATURE_WINDOWS:
            peak_d, peak_i, area = profile_peak_and_area(rows, d_min, d_max)
            summary_rows.append(
                {
                    "comparison": comparison,
                    "feature_window": feature,
                    "feature_group": group,
                    "window_center_d_A": f"{center:.12g}",
                    "window_min_d_A": f"{d_min:.12g}",
                    "window_max_d_A": f"{d_max:.12g}",
                    "difference_peak_d_A": f"{peak_d:.12g}" if peak_d is not None else "",
                    "difference_peak_intensity_norm": f"{peak_i:.12g}" if peak_i is not None else "",
                    "window_area_difference": f"{area:.12g}",
                }
            )
    return difference_profiles, summary_rows


def plot_overlay(path: Path, profiles: dict[str, list[dict[str, float]]], experimental: list[dict[str, float]]) -> None:
    fig, ax = plt.subplots(figsize=(9, 5.5), dpi=160)
    exp_rows = sorted([row for row in experimental if 2.2 <= row["d_A"] <= 9.8], key=lambda row: row["d_A"])
    ax.plot(
        [row["d_A"] for row in exp_rows],
        [row["intensity_norm"] for row in exp_rows],
        color="black",
        linewidth=1.8,
        label="experimental",
    )
    for model_key, label, _path in MODELS:
        rows = sorted([row for row in profiles[model_key] if 2.2 <= row["d_A"] <= 9.8], key=lambda row: row["d_A"])
        ax.plot([row["d_A"] for row in rows], [row["intensity_norm"] for row in rows], linewidth=1.2, label=label)
    for _feature, _center, d_min, d_max, group in FEATURE_WINDOWS[:5]:
        ax.axvspan(d_min, d_max, color="#d8d8d8", alpha=0.14 if group == "primary" else 0.08)
    ax.invert_xaxis()
    ax.set_xlabel("d-spacing (A)")
    ax.set_ylabel("normalized intensity")
    ax.set_title("Corrected Nick/Asem model profiles vs experimental powder profile")
    ax.grid(True, alpha=0.25, linewidth=0.6)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_feature_response(path: Path, summary_rows: list[dict[str, object]]) -> None:
    fig, ax = plt.subplots(figsize=(10, 5.8), dpi=160)
    primary = [row for row in summary_rows if row["feature_group"] == "primary"]
    labels = [row["feature_window"] for row in primary if row["model"] == "with_coo"]
    x = np.arange(len(labels))
    width = 0.24
    offsets = {"with_coo": -width, "only_bases": 0.0, "eight_hexads": width}
    for model_key, label, _path in MODELS:
        values = [
            float(row["simulated_peak_intensity_norm"]) if row["simulated_peak_intensity_norm"] else 0.0
            for row in primary
            if row["model"] == model_key
        ]
        ax.bar(x + offsets[model_key], values, width=width, label=label)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("simulated peak normalized intensity")
    ax.set_title("Primary feature-window response by model")
    ax.grid(True, axis="y", alpha=0.25, linewidth=0.6)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_differences(path: Path, difference_profiles: dict[str, list[dict[str, float]]]) -> None:
    fig, ax = plt.subplots(figsize=(9, 5.5), dpi=160)
    for comparison, rows in difference_profiles.items():
        rows = sorted([row for row in rows if 2.2 <= row["d_A"] <= 9.8], key=lambda row: row["d_A"])
        ax.plot([row["d_A"] for row in rows], [row["intensity_norm"] for row in rows], linewidth=1.2, label=comparison)
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.invert_xaxis()
    ax.set_xlabel("d-spacing (A)")
    ax.set_ylabel("difference in normalized intensity")
    ax.set_title("Corrected model difference profiles")
    ax.grid(True, alpha=0.25, linewidth=0.6)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def top_model_by_feature(summary_rows: list[dict[str, object]], feature: str) -> dict[str, object]:
    rows = [row for row in summary_rows if row["feature_window"] == feature]
    return max(rows, key=lambda row: float(row["simulated_peak_intensity_norm"] or 0.0))


def markdown_table(rows: list[dict[str, object]], columns: list[str]) -> str:
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join("---" for _ in columns) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(str(row.get(column, "")) for column in columns) + " |")
    return "\n".join(lines)


def write_report(
    path: Path,
    args: argparse.Namespace,
    atom_counts: dict[str, int],
    experimental: list[dict[str, float]],
    summary_rows: list[dict[str, object]],
    difference_rows: list[dict[str, object]],
    rotations: list[int],
) -> None:
    primary_top = [
        {
            "feature_window": feature,
            "top_simulated_model": top_model_by_feature(summary_rows, feature)["model"],
            "top_simulated_peak_norm": top_model_by_feature(summary_rows, feature)["simulated_peak_intensity_norm"],
            "experimental_peak_d_A": top_model_by_feature(summary_rows, feature)["experimental_peak_d_A"],
        }
        for feature, _center, _d_min, _d_max, group in FEATURE_WINDOWS
        if group == "primary"
    ]
    diff_primary = [row for row in difference_rows if row["feature_group"] == "primary"]
    d_values = [row["d_A"] for row in experimental]
    lines = [
        "# Nick Atom Contribution Comparison: Asem-Corrected",
        "",
        "## Purpose",
        "",
        "Redo Nick's first powder diffraction model-comparison task using the Asem-corrected non-accumulating/vectorized diffraction path and corrected rotation sampling.",
        "",
        "## Inputs",
        "",
        markdown_table(
            [
                {"model": model_key, "label": label, "atoms": atom_counts[model_key], "path": path.relative_to(ROOT)}
                for model_key, label, path in MODELS
            ],
            ["model", "label", "atoms", "path"],
        ),
        "",
        f"- Experimental profile: `{args.experimental_profile.relative_to(ROOT)}`",
        f"- Experimental rows: {len(experimental)}",
        f"- Experimental d-spacing range: {min(d_values):.12g} A to {max(d_values):.12g} A",
        "",
        "## Corrected Rotation Handling",
        "",
        "- Engine: `reference/asem_corrected_diffraction_engine/`.",
        "- Asem correction: orientation stacks are generated independently, so azimuthal rotations do not accumulate.",
        "- Vectorized path: `make_oriented_coords` plus `generate_fiber_diffraction_series`.",
        f"- Nick-style tilts: `[0]`.",
        f"- Nick-style rotations: `range(0, 181, 5)`, producing {len(rotations)} rotations from {rotations[0]} to {rotations[-1]} degrees.",
        "- 180 degrees was included to reproduce the explicit inclusive 0-180 degree interpretation. A later sensitivity check can compare `range(0, 180, 5)` if endpoint duplication after image spinning/symmetrization matters.",
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
        "",
        "## Primary Feature-Window Tops",
        "",
        markdown_table(primary_top, ["feature_window", "top_simulated_model", "top_simulated_peak_norm", "experimental_peak_d_A"]),
        "",
        "## Full Feature Summary",
        "",
        markdown_table(
            summary_rows,
            [
                "model",
                "feature_window",
                "simulated_peak_d_A",
                "simulated_peak_intensity_norm",
                "experimental_peak_d_A",
                "experimental_peak_intensity_norm",
                "peak_offset_d_A",
                "window_area_simulated",
            ],
        ),
        "",
        "## Difference Summary",
        "",
        markdown_table(
            diff_primary,
            ["comparison", "feature_window", "difference_peak_d_A", "difference_peak_intensity_norm", "window_area_difference"],
        ),
        "",
        "## Preliminary Interpretation",
        "",
        f"- 3.4 A: strongest simulated peak is in `{top_model_by_feature(summary_rows, '3.4 A')['model']}` for this corrected reproduction. Treat this as support for stacked-base contribution only if bases-only remains competitive in area and peak intensity.",
        f"- 7.25 A: strongest simulated peak is in `{top_model_by_feature(summary_rows, '7.25 A')['model']}`; compare the full-minus-bases and full-minus-8hexads differences before assigning it to backbone/full-model structure.",
        f"- 3.7 A, 4.4 A, and 5.6 A: these windows should be treated as mixed until atom-group controls separate backbone and Glu/COO contributions more directly.",
        "",
        "## Limitations",
        "",
        "- Nick's source archive predates Asem's rotation fix.",
        "- This is a corrected reproduction/sensitivity analysis, not a final structural fit.",
        "- Compare q/d-spacing positions and relative feature trends, not raw 2D detector image shape alone.",
        "- This is not yet the length-convergence or twist-refinement analysis.",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    scripts_module, orientation_module, radial_module = load_engine()
    output_root = args.output_root
    radial_dir = output_root / "radial_profiles"
    plot_dir = output_root / "plots"
    table_dir = output_root / "tables"
    for directory in (radial_dir, plot_dir, table_dir, ROOT / "outputs" / "metrics", ROOT / "outputs" / "reports"):
        directory.mkdir(parents=True, exist_ok=True)

    wavelength_mm = args.wavelength_angstrom * 1e-7
    z_grid_limits = [-args.grid_limit_mm, args.grid_limit_mm]
    x_grid_limits = [-args.grid_limit_mm, args.grid_limit_mm]
    tilts = [0]
    rotations = list(range(0, 181, 5))

    profiles: dict[str, list[dict[str, float]]] = {}
    atom_counts: dict[str, int] = {}
    run_metadata = {
        "engine_dir": str(ENGINE_DIR.relative_to(ROOT)),
        "rotation_sampling": "range(0, 181, 5)",
        "rotations_deg": rotations,
        "tilts_deg": tilts,
        "grid_size": args.grid_size,
        "grid_limit_mm": args.grid_limit_mm,
        "radial_bins": args.radial_bins,
        "detector_distance_mm": args.detector_distance_mm,
        "wavelength_A": args.wavelength_angstrom,
        "note": "180 degrees included; compare range(0, 180, 5) in a later endpoint sensitivity check if needed.",
    }

    for model_key, label, path in MODELS:
        print(f"Running corrected diffraction for {model_key}: {path}", flush=True)
        atoms, coords_angstrom, atomic_numbers = load_xyz(path, scripts_module.atomic_number)
        atom_counts[model_key] = len(atoms)
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
        profiles[model_key] = profile
        write_profile(radial_dir / f"{model_key}_radial_profile.csv", profile)
        run_metadata.setdefault("models", {})[model_key] = {
            "label": label,
            "path": str(path.relative_to(ROOT)),
            "atoms": len(atoms),
            "image_min": float(np.min(image)),
            "image_max": float(np.max(image)),
            "image_mean": float(np.mean(image)),
        }

    experimental = read_experimental_profile(args.experimental_profile)
    feature_summary = build_feature_summary(profiles, experimental)
    difference_profiles, difference_summary = build_difference_profiles(profiles, table_dir)

    feature_fields = [
        "model",
        "model_label",
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
    difference_fields = [
        "comparison",
        "feature_window",
        "feature_group",
        "window_center_d_A",
        "window_min_d_A",
        "window_max_d_A",
        "difference_peak_d_A",
        "difference_peak_intensity_norm",
        "window_area_difference",
    ]
    write_csv(ROOT / "outputs" / "metrics" / "nick_atom_contribution_feature_summary_corrected.csv", feature_summary, feature_fields)
    write_csv(table_dir / "nick_atom_contribution_feature_summary_corrected.csv", feature_summary, feature_fields)
    write_csv(
        ROOT / "outputs" / "metrics" / "nick_atom_contribution_difference_summary_corrected.csv",
        difference_summary,
        difference_fields,
    )
    write_csv(table_dir / "nick_atom_contribution_difference_summary_corrected.csv", difference_summary, difference_fields)
    (table_dir / "run_metadata.json").write_text(json.dumps(run_metadata, indent=2) + "\n", encoding="utf-8")

    plot_overlay(plot_dir / "nick_atom_contribution_profile_overlay_corrected.png", profiles, experimental)
    plot_feature_response(plot_dir / "nick_atom_contribution_feature_response_corrected.png", feature_summary)
    plot_differences(plot_dir / "nick_atom_contribution_difference_profiles_corrected.png", difference_profiles)

    write_report(
        ROOT / "outputs" / "reports" / "nick_atom_contribution_corrected_report.md",
        args,
        atom_counts,
        experimental,
        feature_summary,
        difference_summary,
        rotations,
    )

    print("Wrote corrected Nick atom-contribution outputs", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
