#!/usr/bin/env python3
"""Run corrected diffraction for ideal rise variants against the Emory profile."""

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
VARIANT_DIR = ROOT / "inputs" / "nick_ideal_models" / "rise_variants"
EXPERIMENTAL_PROFILE = ROOT / "inputs" / "experimental" / "nick_powder_profile_corrected_emory.csv"
OUTPUT_ROOT = ROOT / "outputs" / "rise_sensitivity_corrected_emory_profile"
FEATURE_SUMMARY = ROOT / "outputs" / "metrics" / "rise_sensitivity_feature_summary_corrected_emory_profile.csv"
RANKING = ROOT / "outputs" / "metrics" / "rise_sensitivity_ranking_corrected_emory_profile.csv"
REPORT = ROOT / "outputs" / "reports" / "rise_sensitivity_corrected_emory_profile_report.md"

FEATURE_WINDOWS = [
    ("3.38/3.4 A", 3.38, 3.30, 3.46, "primary"),
    ("3.77 A", 3.77, 3.60, 3.90, "primary"),
    ("4.4 A", 4.40, 4.25, 4.55, "primary"),
    ("5.6 A", 5.60, 5.45, 5.75, "primary"),
    ("7.3 A", 7.30, 7.05, 7.55, "primary"),
    ("3.0 A", 3.00, 2.90, 3.10, "optional"),
    ("4.1 A", 4.10, 4.00, 4.20, "optional"),
    ("4.5-5.0 A", 4.75, 4.50, 5.00, "optional"),
    ("5.5 A", 5.50, 5.40, 5.60, "optional"),
    ("7.0 A", 7.00, 6.80, 7.20, "optional"),
    ("8.4 A", 8.40, 8.20, 8.60, "optional"),
]
PRIMARY_FEATURES = ["3.38/3.4 A", "3.77 A", "4.4 A", "5.6 A", "7.3 A"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variant-dir", type=Path, default=VARIANT_DIR)
    parser.add_argument("--experimental-profile", type=Path, default=EXPERIMENTAL_PROFILE)
    parser.add_argument("--output-root", type=Path, default=OUTPUT_ROOT)
    parser.add_argument("--feature-summary-output", type=Path, default=FEATURE_SUMMARY)
    parser.add_argument("--ranking-output", type=Path, default=RANKING)
    parser.add_argument("--report-output", type=Path, default=REPORT)
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
    orientation_module = load_module("orientation_average_corrected_rise", ENGINE_DIR / "orientation_average.py")
    radial_module = load_module("radial_average_corrected_rise", ENGINE_DIR / "radial_average.py")
    return scripts_module, orientation_module, radial_module


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path)


def rise_from_path(path: Path) -> float:
    token = path.stem.rsplit("_rise_", 1)[1]
    return float(token.replace("p", "."))


def rise_token(value: float) -> str:
    return f"{value:.2f}".replace(".", "p")


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
    compact = [
        {
            "bin_index": int(row["bin_index"]),
            "q_Ainv": float(row["q_Ainv"]),
            "d_A": float(row["d_A"]),
            "mean_intensity": float(row["mean_intensity"]),
            "pixel_count": int(row["pixel_count"]),
        }
        for row in rows
    ]
    return normalize_radial_rows(compact, args.normalize_q_min)


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_radial(path: Path, rows: list[dict[str, float]]) -> None:
    write_csv(
        path,
        [
            {
                "bin_index": row["bin_index"],
                "q_Ainv": f"{row['q_Ainv']:.8f}",
                "d_A": f"{row['d_A']:.8f}",
                "mean_intensity": f"{row['mean_intensity']:.12g}",
                "intensity_norm": f"{row['intensity_norm']:.12g}",
                "pixel_count": row["pixel_count"],
            }
            for row in rows
        ],
    )


def window_peak(rows: list[dict[str, float]], low: float, high: float) -> tuple[float | None, float | None, float]:
    subset = [row for row in rows if low <= row["d_A"] <= high]
    if not subset:
        return None, None, 0.0
    peak = max(subset, key=lambda row: row["intensity"])
    ordered = sorted(subset, key=lambda row: row["d_A"])
    area = float(np.trapezoid([row["intensity"] for row in ordered], [row["d_A"] for row in ordered]))
    return peak["d_A"], peak["intensity"], area


def build_feature_summary(
    rise: float,
    profile: list[dict[str, float]],
    experimental: list[dict[str, float]],
) -> list[dict[str, object]]:
    rows = []
    for name, center, low, high, group in FEATURE_WINDOWS:
        exp_peak, _exp_intensity, exp_area = window_peak(experimental, low, high)
        sim_peak, sim_intensity, sim_area = window_peak(profile, low, high)
        offset = sim_peak - exp_peak if sim_peak is not None and exp_peak is not None else None
        rows.append(
            {
                "rise_A": f"{rise:.2f}",
                "feature_group": group,
                "feature_window": name,
                "window_center_A": center,
                "window_min_A": low,
                "window_max_A": high,
                "experimental_peak_d_A": exp_peak,
                "simulated_peak_d_A": sim_peak,
                "peak_offset_d_A": offset,
                "abs_peak_offset_d_A": abs(offset) if offset is not None else None,
                "simulated_peak_intensity_norm": sim_intensity,
                "window_area_simulated": sim_area,
                "window_area_experimental": exp_area,
            }
        )
    return rows


def build_ranking(summary_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    for rise in sorted({row["rise_A"] for row in summary_rows}, key=float):
        subset = [row for row in summary_rows if row["rise_A"] == rise and row["feature_group"] == "primary"]
        by_feature = {row["feature_window"]: row for row in subset}
        offsets = [float(row["abs_peak_offset_d_A"]) for row in subset if row["abs_peak_offset_d_A"] is not None]
        row = {
            "rise_A": rise,
            "mean_abs_primary_peak_offset_A": float(np.mean(offsets)) if offsets else None,
            "max_abs_primary_peak_offset_A": max(offsets) if offsets else None,
            "base_stack_peak_offset_A": by_feature["3.38/3.4 A"]["peak_offset_d_A"],
            "feature_3p77_offset_A": by_feature["3.77 A"]["peak_offset_d_A"],
            "feature_4p4_offset_A": by_feature["4.4 A"]["peak_offset_d_A"],
            "feature_5p6_offset_A": by_feature["5.6 A"]["peak_offset_d_A"],
            "feature_7p3_offset_A": by_feature["7.3 A"]["peak_offset_d_A"],
            "rank_by_mean_abs_offset": None,
            "qualitative_note": "",
        }
        rows.append(row)
    ranked = sorted(rows, key=lambda row: row["mean_abs_primary_peak_offset_A"])
    best = ranked[0]["rise_A"]
    for index, row in enumerate(ranked, start=1):
        row["rank_by_mean_abs_offset"] = index
        row["qualitative_note"] = "best mean primary-window offset" if row["rise_A"] == best else "worse than best by mean primary-window offset"
    return sorted(ranked, key=lambda row: float(row["rise_A"]))


def smooth_profile(rows: list[dict[str, float]], sigma_a: float, d_grid: np.ndarray) -> np.ndarray:
    x = np.asarray([row["d_A"] for row in rows], dtype=float)
    y = np.asarray([row["intensity"] for row in rows], dtype=float)
    order = np.argsort(x)
    x = x[order]
    y = y[order]
    interp = np.interp(d_grid, x, y, left=np.nan, right=np.nan)
    step = float(np.median(np.diff(d_grid)))
    sigma_bins = max(sigma_a / step, 0.1)
    radius = max(int(4 * sigma_bins), 1)
    kernel_x = np.arange(-radius, radius + 1)
    kernel = np.exp(-0.5 * (kernel_x / sigma_bins) ** 2)
    kernel /= np.sum(kernel)
    valid = np.isfinite(interp)
    filled = np.where(valid, interp, 0.0)
    weights = np.convolve(valid.astype(float), kernel, mode="same")
    smoothed = np.convolve(filled, kernel, mode="same")
    return np.divide(smoothed, weights, out=np.full_like(smoothed, np.nan), where=weights > 0)


def plot_overlay(path: Path, profiles: dict[float, list[dict[str, float]]], experimental: list[dict[str, float]]) -> None:
    fig, ax = plt.subplots(figsize=(10, 6), facecolor="white")
    exp = [row for row in experimental if 2.8 <= row["d_A"] <= 8.8]
    ax.plot([row["d_A"] for row in exp], [row["intensity"] for row in exp], color="black", lw=1.4, label="Corrected Experimental")
    for rise, rows in sorted(profiles.items()):
        subset = [row for row in rows if 2.8 <= row["d_A"] <= 8.8]
        ax.plot([row["d_A"] for row in subset], [row["intensity"] for row in subset], lw=1.1, label=f"rise {rise:.2f} A")
    ax.set_xlim(8.8, 2.8)
    ax.set_xlabel("d-spacing (A)")
    ax.set_ylabel("normalized intensity")
    ax.set_title("Rise variants vs Emory-corrected experimental profile")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_peak_offsets(path: Path, summary_rows: list[dict[str, object]]) -> None:
    fig, ax = plt.subplots(figsize=(9, 5.5), facecolor="white")
    for feature in PRIMARY_FEATURES:
        rows = [row for row in summary_rows if row["feature_window"] == feature]
        ax.plot(
            [float(row["rise_A"]) for row in rows],
            [float(row["peak_offset_d_A"]) for row in rows],
            marker="o",
            label=feature,
        )
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xlabel("target rise (A)")
    ax.set_ylabel("simulated - experimental peak position (A)")
    ax.set_title("Primary feature peak offsets vs rise")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_mean_error(path: Path, ranking_rows: list[dict[str, object]]) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 5), facecolor="white")
    ax.plot(
        [float(row["rise_A"]) for row in ranking_rows],
        [float(row["mean_abs_primary_peak_offset_A"]) for row in ranking_rows],
        marker="o",
        color="#2f6f9f",
    )
    ax.set_xlabel("target rise (A)")
    ax.set_ylabel("mean absolute primary peak offset (A)")
    ax.set_title("Rise sensitivity mean peak-position error")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_nick_style(path: Path, profiles: dict[float, list[dict[str, float]]], experimental: list[dict[str, float]], sigma_a: float) -> None:
    d_grid = np.linspace(3.0, 8.0, 900)
    fig, ax = plt.subplots(figsize=(9, 6), facecolor="white")
    traces = [
        ("Corrected Experimental", experimental, "black", 2.4),
        ("rise 3.38 A", profiles[3.38], "#008a2e", 1.2),
        ("rise 3.40 A", profiles[3.40], "#c62828", 0.0),
    ]
    for label, rows, color, offset in traces:
        smoothed = smooth_profile(rows, sigma_a, d_grid)
        mask = np.isfinite(smoothed)
        values = smoothed[mask]
        normed = (smoothed - np.nanmin(values)) / (np.nanmax(values) - np.nanmin(values)) if len(values) and np.nanmax(values) > np.nanmin(values) else smoothed
        ax.plot(d_grid, normed + offset, color=color, lw=1.5, label=label)
    for center in [3.38, 3.77, 4.4, 5.6, 7.3]:
        ax.axvspan(center - 0.04, center + 0.04, color="lightgray", alpha=0.18, zorder=0)
    ax.set_xlim(8.0, 3.0)
    ax.set_xlabel("d-spacing (A)")
    ax.set_ylabel("independently normalized intensity + offset")
    ax.set_title("Nick-style rise comparison: 3.38 A vs 3.40 A")
    ax.legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_report(
    path: Path,
    args: argparse.Namespace,
    ranking_rows: list[dict[str, object]],
    summary_rows: list[dict[str, object]],
    validation_rows: list[dict[str, str]],
    metadata: dict[str, object],
) -> None:
    best = min(ranking_rows, key=lambda row: float(row["mean_abs_primary_peak_offset_A"]))
    row_338 = next(row for row in ranking_rows if row["rise_A"] == "3.38")
    row_340 = next(row for row in ranking_rows if row["rise_A"] == "3.40")
    relation = "improves" if float(row_338["mean_abs_primary_peak_offset_A"]) < float(row_340["mean_abs_primary_peak_offset_A"]) else "does not improve"

    lines = [
        "# Rise Sensitivity Against Emory-Corrected Profile",
        "",
        "## Purpose",
        "",
        "This controlled test evaluates Nick's hypothesis that the clarified ideal antiparallel 30-degree model may fit the Emory-corrected powder profile better at a helical rise near 3.38 A than at the current ideal 3.40 A.",
        "",
        "## Inputs",
        "",
        "- Source ideal PDB: `inputs/nick_ideal_models/Hexaplex_AntiParallel_30deg_Ideal.pdb`",
        f"- Corrected experimental profile: `{display_path(args.experimental_profile)}`",
        f"- Generated rise variants: `{display_path(args.variant_dir)}`",
        "",
        "## Rise Variant Generation",
        "",
        "The source PDB does not store explicit layer IDs. The generator inferred 15 six-residue CYP/MEP base planes from residue z-centroids; each plane contains three CYP and three MEP residues. Non-base GLU residues were assigned to the nearest inferred base plane and translated rigidly with that plane. The central inferred base plane was held fixed as the z anchor.",
        "",
        "This perturbation preserves atom order, x/y coordinates, twist geometry, and intralayer geometry. It only changes rigid layer z positions. It does not adjust peptide omega, side-chain removal logic, twist angle, or perform minimization.",
        "",
        "## Validation",
        "",
        "| Target rise (A) | Atom count preserved | Measured rise (A) | Max rise error (A) | Max layer z delta (A) |",
        "| ---: | --- | ---: | ---: | ---: |",
    ]
    for row in validation_rows:
        lines.append(
            f"| {row['rise_A']} | {row['atom_count_preserved']} | "
            f"{float(row['measured_variant_base_layer_rise_A']):.4f} | "
            f"{float(row['max_adjacent_rise_error_A']):.3g} | "
            f"{float(row['max_abs_layer_z_delta_A']):.4f} |"
        )
    lines.extend(
        [
            "",
            "## Corrected Diffraction Settings",
            "",
            "- Engine: `reference/asem_corrected_diffraction_engine/`",
            "- Asem correction: azimuthal rotations are applied to independent tilted coordinate stacks.",
            "- Tilts: `[0]`",
            "- Rotations: `range(0, 181, 5)`",
            f"- Grid: {metadata['grid_size']} x {metadata['grid_size']} over +/-{metadata['grid_limit_mm']} mm",
            f"- Detector distance: {metadata['detector_distance_mm']} mm",
            f"- Wavelength: {metadata['wavelength_A']} A",
            f"- Radial bins: {metadata['radial_bins']}",
            "- Hydrogens excluded and exact heavy-atom deduplication applied in XYZ inputs.",
            "",
            "## Ranking",
            "",
            "| Rise (A) | Mean abs primary offset (A) | Max abs primary offset (A) | Base-stack offset (A) | 3.77 offset (A) | 4.4 offset (A) | 5.6 offset (A) | 7.3 offset (A) | Rank |",
            "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in sorted(ranking_rows, key=lambda item: int(item["rank_by_mean_abs_offset"])):
        lines.append(
            f"| {row['rise_A']} | {float(row['mean_abs_primary_peak_offset_A']):.4f} | "
            f"{float(row['max_abs_primary_peak_offset_A']):.4f} | "
            f"{float(row['base_stack_peak_offset_A']):.4f} | "
            f"{float(row['feature_3p77_offset_A']):.4f} | "
            f"{float(row['feature_4p4_offset_A']):.4f} | "
            f"{float(row['feature_5p6_offset_A']):.4f} | "
            f"{float(row['feature_7p3_offset_A']):.4f} | "
            f"{row['rank_by_mean_abs_offset']} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            f"- Best rise by mean primary-window peak offset: {best['rise_A']} A.",
            f"- The 3.38 A variant {relation} relative to 3.40 A by mean primary-window peak offset.",
            f"- 3.38 A mean abs primary offset: {float(row_338['mean_abs_primary_peak_offset_A']):.4f} A.",
            f"- 3.40 A mean abs primary offset: {float(row_340['mean_abs_primary_peak_offset_A']):.4f} A.",
            "",
            "Feature-specific offsets are in `outputs/metrics/rise_sensitivity_feature_summary_corrected_emory_profile.csv`. The comparison is a controlled sensitivity test; a favorable 3.38 A result would support asking Asem for a larger/refined 3.38 A model, not replacing a chemically relaxed model.",
            "",
            "## Plots",
            "",
            "- `outputs/rise_sensitivity_corrected_emory_profile/plots/rise_variant_profile_overlay.png`",
            "- `outputs/rise_sensitivity_corrected_emory_profile/plots/rise_variant_peak_offsets.png`",
            "- `outputs/rise_sensitivity_corrected_emory_profile/plots/rise_variant_mean_error.png`",
            "- `outputs/rise_sensitivity_corrected_emory_profile/plots/rise_3p38_vs_3p40_nick_style.png`",
            "",
            "The Nick-style plot uses Gaussian smoothing in d-space with sigma 0.10 A for visualization only. Scoring uses unsmoothed radial profiles.",
            "",
            "## Limitations",
            "",
            "- Rigid-layer z-translation does not relax backbone geometry.",
            "- This does not address peptide omega.",
            "- This does not alter twist angle.",
            "- Powder/radial matching is a sensitivity test, not full refinement.",
            "- Layer inference is based on CYP/MEP base-plane z-centroids because explicit layer IDs are not present in the source PDB.",
        ]
    )
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

    xyz_paths = sorted((args.variant_dir / "xyz").glob("Hexaplex_AntiParallel_30deg_Ideal_rise_*.xyz"), key=rise_from_path)
    if len(xyz_paths) != 5:
        raise FileNotFoundError(f"Expected five rise-variant XYZ files in {args.variant_dir / 'xyz'}, found {len(xyz_paths)}")

    experimental = read_experimental(args.experimental_profile)
    wavelength_mm = args.wavelength_angstrom * 1e-7
    z_grid_limits = [-args.grid_limit_mm, args.grid_limit_mm]
    x_grid_limits = [-args.grid_limit_mm, args.grid_limit_mm]
    tilts = [0]
    rotations = list(range(0, 181, 5))

    profiles: dict[float, list[dict[str, float]]] = {}
    metadata_rows = []
    summary_rows = []
    for xyz_path in xyz_paths:
        rise = rise_from_path(xyz_path)
        print(f"Running corrected diffraction for rise {rise:.2f} A: {xyz_path}", flush=True)
        _atoms, coords_angstrom, atomic_numbers = load_xyz(xyz_path, scripts_module.atomic_number)
        image = orientation_module.average_fiber_diffraction(
            atomic_numbers,
            coords_angstrom * 1e-7,
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
        radial_path = radial_dir / f"Hexaplex_AntiParallel_30deg_Ideal_rise_{rise_token(rise)}_radial_profile.csv"
        write_radial(radial_path, profile)
        profile_for_scoring = [
            {"d_A": row["d_A"], "q_Ainv": row["q_Ainv"], "intensity": row["intensity_norm"]}
            for row in profile
        ]
        profiles[rise] = profile_for_scoring
        summary_rows.extend(build_feature_summary(rise, profile_for_scoring, experimental))
        metadata_rows.append(
            {
                "rise_A": f"{rise:.2f}",
                "xyz_path": display_path(xyz_path),
                "radial_profile": display_path(radial_path),
                "atom_count": len(atomic_numbers),
                "image_min": float(np.min(image)),
                "image_max": float(np.max(image)),
                "image_mean": float(np.mean(image)),
            }
        )

    ranking_rows = build_ranking(summary_rows)
    write_csv(args.feature_summary_output, summary_rows)
    write_csv(args.ranking_output, ranking_rows)
    validation_rows = list(csv.DictReader((args.variant_dir / "rise_variant_validation.csv").open("r", encoding="utf-8", newline="")))
    metadata = {
        "engine_dir": display_path(ENGINE_DIR),
        "experimental_profile": display_path(args.experimental_profile),
        "rotation_sampling": "range(0, 181, 5)",
        "rotations_deg": rotations,
        "tilts_deg": tilts,
        "grid_size": args.grid_size,
        "grid_limit_mm": args.grid_limit_mm,
        "radial_bins": args.radial_bins,
        "detector_distance_mm": args.detector_distance_mm,
        "wavelength_A": args.wavelength_angstrom,
        "normalize_q_min": args.normalize_q_min,
        "runs": metadata_rows,
    }
    (table_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    plot_overlay(plot_dir / "rise_variant_profile_overlay.png", profiles, experimental)
    plot_peak_offsets(plot_dir / "rise_variant_peak_offsets.png", summary_rows)
    plot_mean_error(plot_dir / "rise_variant_mean_error.png", ranking_rows)
    plot_nick_style(plot_dir / "rise_3p38_vs_3p40_nick_style.png", profiles, experimental, args.plot_sigma_a)
    write_report(args.report_output, args, ranking_rows, summary_rows, validation_rows, metadata)

    print(f"Wrote feature summary: {args.feature_summary_output}")
    print(f"Wrote ranking: {args.ranking_output}")
    print(f"Wrote report: {args.report_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
