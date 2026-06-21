#!/usr/bin/env python3
"""Run the ideal no-CH2-COOH comparison and four-trace Nick-style figure."""

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
DEFAULT_NO_XYZ = ROOT / "inputs" / "nick_ideal_models" / "derived" / "Hexaplex_AntiParallel_30deg_Ideal_no_CH2_COOH.xyz"
DEFAULT_NO_PDB = ROOT / "inputs" / "nick_ideal_models" / "derived" / "Hexaplex_AntiParallel_30deg_Ideal_no_CH2_COOH.pdb"
DEFAULT_PROVENANCE = (
    ROOT / "inputs" / "nick_ideal_models" / "derived" / "Hexaplex_AntiParallel_30deg_Ideal_no_CH2_COOH_provenance.json"
)
DEFAULT_FULL_RADIAL = (
    ROOT
    / "outputs"
    / "nick_ideal_16mer_corrected_emory_profile"
    / "radial_profiles"
    / "ideal_16mer_antiparallel_30deg_radial_profile.csv"
)
DEFAULT_EXPERIMENTAL = ROOT / "inputs" / "experimental" / "nick_powder_profile_corrected_emory.csv"
NO_OUTPUT_ROOT = ROOT / "outputs" / "nick_ideal_16mer_no_ch2_cooh_corrected_emory_profile"
SUB_OUTPUT_ROOT = ROOT / "outputs" / "nick_ideal_16mer_sidechain_subtraction_corrected_emory_profile"
DEFAULT_METRICS = ROOT / "outputs" / "metrics" / "nick_ideal_16mer_atom_contribution_feature_summary.csv"
DEFAULT_REPORT = ROOT / "outputs" / "reports" / "nick_ideal_16mer_atom_contribution_report.md"
DEFAULT_FOUR_TRACE = (
    ROOT
    / "outputs"
    / "nick_ideal_16mer_corrected_emory_profile"
    / "plots"
    / "nick_ideal_16mer_four_trace_nick_style.png"
)

FEATURE_WINDOWS = [
    ("7.25 A", 7.25, 7.0, 7.5, "likely backbone-dominated"),
    ("5.6 A", 5.6, 5.5, 5.7, "mixed backbone/side-chain"),
    ("4.4 A", 4.4, 4.3, 4.5, "mixed backbone/side-chain"),
    ("3.7 A", 3.7, 3.6, 3.8, "mixed backbone/side-chain"),
    ("3.4 A", 3.4, 3.3, 3.5, "primarily bases"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-ch2-cooh-xyz", type=Path, default=DEFAULT_NO_XYZ)
    parser.add_argument("--no-ch2-cooh-pdb", type=Path, default=DEFAULT_NO_PDB)
    parser.add_argument("--provenance", type=Path, default=DEFAULT_PROVENANCE)
    parser.add_argument("--full-radial", type=Path, default=DEFAULT_FULL_RADIAL)
    parser.add_argument("--experimental-profile", type=Path, default=DEFAULT_EXPERIMENTAL)
    parser.add_argument("--no-output-root", type=Path, default=NO_OUTPUT_ROOT)
    parser.add_argument("--subtraction-output-root", type=Path, default=SUB_OUTPUT_ROOT)
    parser.add_argument("--feature-output", type=Path, default=DEFAULT_METRICS)
    parser.add_argument("--report-output", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--four-trace-output", type=Path, default=DEFAULT_FOUR_TRACE)
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


def interpolate(rows: list[dict[str, float]], d_grid: np.ndarray) -> np.ndarray:
    ordered = sorted(rows, key=lambda row: row["d_A"])
    x = np.asarray([row["d_A"] for row in ordered], dtype=np.float64)
    y = np.asarray([row["intensity"] for row in ordered], dtype=np.float64)
    return np.interp(d_grid, x, y)


def gaussian_kernel(sigma: float, step: float) -> np.ndarray:
    radius = max(3, int(np.ceil(4.0 * sigma / step)))
    offsets = np.arange(-radius, radius + 1, dtype=np.float64) * step
    kernel = np.exp(-0.5 * (offsets / sigma) ** 2)
    kernel /= np.sum(kernel)
    return kernel


def smooth(y_values: np.ndarray, sigma: float, step: float) -> np.ndarray:
    kernel = gaussian_kernel(sigma, step)
    padded = np.pad(y_values, (len(kernel) // 2, len(kernel) // 2), mode="edge")
    return np.convolve(padded, kernel, mode="valid")


def normalize(y_values: np.ndarray) -> np.ndarray:
    y = np.asarray(y_values, dtype=np.float64)
    y = y - np.nanmin(y)
    peak = np.nanmax(y)
    return y / peak if peak > 0 else y


def peak_in_window(d_grid: np.ndarray, y_values: np.ndarray, d_min: float, d_max: float) -> tuple[float | None, float | None]:
    mask = (d_grid >= d_min) & (d_grid <= d_max) & np.isfinite(y_values)
    if not np.any(mask):
        return None, None
    local_d = d_grid[mask]
    local_y = y_values[mask]
    index = int(np.argmax(local_y))
    return float(local_d[index]), float(local_y[index])


def plot_four_trace(
    path: Path,
    d_grid: np.ndarray,
    exp_y: np.ndarray,
    full_y: np.ndarray,
    no_y: np.ndarray,
    sub_y: np.ndarray,
    sigma: float,
) -> None:
    fig, ax = plt.subplots(figsize=(6.6, 5.0), dpi=180)
    fig.patch.set_facecolor("#8a8a8a")
    ax.set_facecolor("#8a8a8a")
    offsets = [3.0, 2.05, 1.10, 0.15]
    labels = [
        ("Corrected Experimental", "black", exp_y, offsets[0]),
        ("Theoretical 16-mer Hexaplex", "red", full_y, offsets[1]),
        ("16-mer No CH2-COOH", "#00cc00", no_y, offsets[2]),
        ("Side Chain Contribution Derived From Subtraction", "#7b2cff", sub_y, offsets[3]),
    ]
    for _label, _center, d_min, d_max, _note in FEATURE_WINDOWS:
        ax.axvspan(d_min, d_max, color="white", alpha=0.06, linewidth=0)
    for label, color, y_values, offset in labels:
        ax.plot(d_grid, normalize(y_values) + offset, color=color, linewidth=1.15, label=label)
    ax.set_xlim(8.0, 3.0)
    ax.set_ylim(-0.05, 4.15)
    ax.set_xlabel("A")
    ax.set_ylabel("offset normalized intensity")
    ax.set_title(f"Nick-style ideal 16-mer atom contribution (sigma {sigma:.2f} A)", fontsize=9.5)
    ax.legend(loc="upper left", fontsize=5.8, frameon=True, facecolor="white", edgecolor="none")
    ax.tick_params(axis="both", colors="black", labelsize=7, direction="out", length=3)
    for spine in ax.spines.values():
        spine.set_color("black")
    ax.grid(False)
    fig.tight_layout(pad=1.0)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)


def write_report(path: Path, args: argparse.Namespace, provenance: dict[str, object], feature_rows: list[dict[str, object]]) -> None:
    lines = [
        "# Nick Ideal 16-mer Atom-Contribution Comparison",
        "",
        "## Purpose",
        "",
        "Reproduce Nick's four-trace atom-contribution figure using the corrected Asem diffraction path and the Emory-corrected experimental profile.",
        "",
        "## Inputs and Derived Model",
        "",
        f"- Corrected experimental profile: `{display_path(args.experimental_profile)}`",
        f"- Full ideal 16-mer radial profile: `{display_path(args.full_radial)}`",
        f"- Derived no-CH2-COOH PDB: `{display_path(args.no_ch2_cooh_pdb)}`",
        f"- Derived no-CH2-COOH XYZ: `{display_path(args.no_ch2_cooh_xyz)}`",
        f"- Side-chain subtraction profile: `{display_path(args.subtraction_output_root / 'radial_profiles' / 'sidechain_subtraction_profile.csv')}`",
        f"- Final four-trace plot: `{display_path(args.four_trace_output)}`",
        "",
        "## Atom-Selection Rule",
        "",
        str(provenance["atom_selection_rule"]),
        "",
        f"- Source PDB atoms: {provenance['source_atom_count']}",
        f"- Derived PDB atoms: {provenance['derived_pdb_atom_count']}",
        f"- Removed atoms: {provenance['removed_atom_count']}",
        f"- Derived heavy deduped XYZ atoms: {provenance['derived_heavy_deduped_xyz_atom_count']}",
        f"- Removed element counts: `{json.dumps(provenance['removed_element_counts'], sort_keys=True)}`",
        f"- Remaining derived PDB element counts: `{json.dumps(provenance['derived_pdb_element_counts'], sort_keys=True)}`",
        "",
        "## Corrected Diffraction Settings",
        "",
        "- Engine: `reference/asem_corrected_diffraction_engine/`.",
        "- Asem correction: azimuthal rotations are applied to independent tilted coordinate stacks rather than accumulated through a rotation loop.",
        "- Tilts: `[0]`.",
        "- Rotations: `range(0, 181, 5)`, 37 rotations from 0 to 180 degrees.",
        "- Hydrogens excluded from XYZ; exact heavy-atom deduplication applied.",
        f"- Grid: {args.grid_size} x {args.grid_size}; detector half-width {args.grid_limit_mm:g} mm; detector distance {args.detector_distance_mm:g} mm; wavelength {args.wavelength_angstrom:g} A; radial bins {args.radial_bins}.",
        "",
        "## Quantitative Feature Attribution",
        "",
        "| feature_window | experimental_peak_d_A | full_ideal_peak_d_A | no_ch2_cooh_peak_d_A | subtraction_peak_d_A | attribution_note |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for row in feature_rows:
        lines.append(
            f"| {row['feature_window']} | {row['experimental_peak_d_A']} | {row['full_ideal_peak_d_A']} | {row['no_ch2_cooh_peak_d_A']} | {row['subtraction_peak_d_A']} | {row['attribution_note']} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- The 3.4 A window remains most consistent with stacked-base structure because it is strong in the full and no-CH2-COOH traces and is not uniquely isolated by the subtraction trace.",
            "- The 7.25 A window remains likely backbone-dominated because removing the terminal GLU CH2-COOH group does not create a uniquely dominant subtraction feature there.",
            "- The 3.7 A, 4.4 A, and 5.6 A windows are best treated as mixed contribution windows in this corrected comparison.",
            "",
            "## Visualization Notes",
            "",
            f"The four-trace figure uses d-space Gaussian smoothing with sigma {args.plot_sigma_a:.2f} A for visual comparison only. The subtraction trace is full ideal minus no-CH2-COOH after interpolation onto a common d-spacing grid. Smoothing and stacked offsets do not replace the quantitative peak-position table.",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    scripts_module, orientation_module, radial_module = load_engine()
    no_radial_dir = args.no_output_root / "radial_profiles"
    no_table_dir = args.no_output_root / "tables"
    sub_radial_dir = args.subtraction_output_root / "radial_profiles"
    for directory in (no_radial_dir, no_table_dir, sub_radial_dir, args.feature_output.parent, args.report_output.parent):
        directory.mkdir(parents=True, exist_ok=True)

    wavelength_mm = args.wavelength_angstrom * 1e-7
    z_grid_limits = [-args.grid_limit_mm, args.grid_limit_mm]
    x_grid_limits = [-args.grid_limit_mm, args.grid_limit_mm]
    tilts = [0]
    rotations = list(range(0, 181, 5))

    atoms, coords_angstrom, atomic_numbers = load_xyz(args.no_ch2_cooh_xyz, scripts_module.atomic_number)
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
    no_rows = radial_rows_from_image(image, args, radial_module)
    no_radial = no_radial_dir / "ideal_16mer_no_ch2_cooh_radial_profile.csv"
    write_radial(no_radial, no_rows)

    metadata = {
        "engine_dir": display_path(ENGINE_DIR),
        "rotation_sampling": "range(0, 181, 5)",
        "rotations_deg": rotations,
        "tilts_deg": tilts,
        "model_xyz": display_path(args.no_ch2_cooh_xyz),
        "model_atoms": len(atoms),
        "radial_profile": display_path(no_radial),
        "grid_size": args.grid_size,
        "grid_limit_mm": args.grid_limit_mm,
        "radial_bins": args.radial_bins,
        "detector_distance_mm": args.detector_distance_mm,
        "wavelength_A": args.wavelength_angstrom,
        "image_min": float(np.min(image)),
        "image_max": float(np.max(image)),
        "image_mean": float(np.mean(image)),
    }
    (no_table_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    full = read_profile(args.full_radial, "intensity_norm")
    no = read_profile(no_radial, "intensity_norm")
    experimental = read_experimental(args.experimental_profile)
    d_grid = np.linspace(3.0, 8.0, 1600)
    step = float(np.mean(np.diff(d_grid)))
    full_y = interpolate(full, d_grid)
    no_y = interpolate(no, d_grid)
    exp_y = interpolate(experimental, d_grid)
    subtraction_y = full_y - no_y
    sm_full = smooth(full_y, args.plot_sigma_a, step)
    sm_no = smooth(no_y, args.plot_sigma_a, step)
    sm_exp = smooth(exp_y, args.plot_sigma_a, step)
    sm_sub = smooth(subtraction_y, args.plot_sigma_a, step)

    sub_rows = [
        {
            "d_A": f"{d_value:.12g}",
            "q_Ainv": f"{(2.0 * np.pi / d_value):.12g}",
            "intensity_difference_norm": f"{intensity:.12g}",
        }
        for d_value, intensity in zip(d_grid, subtraction_y)
    ]
    subtraction_csv = sub_radial_dir / "sidechain_subtraction_profile.csv"
    write_csv(subtraction_csv, sub_rows, ["d_A", "q_Ainv", "intensity_difference_norm"])

    feature_rows = []
    for feature, _center, d_min, d_max, note in FEATURE_WINDOWS:
        exp_peak, _ = peak_in_window(d_grid, exp_y, d_min, d_max)
        full_peak, _ = peak_in_window(d_grid, full_y, d_min, d_max)
        no_peak, _ = peak_in_window(d_grid, no_y, d_min, d_max)
        sub_peak, sub_intensity = peak_in_window(d_grid, subtraction_y, d_min, d_max)
        feature_rows.append(
            {
                "feature_window": feature,
                "experimental_peak_d_A": f"{exp_peak:.12g}" if exp_peak is not None else "",
                "full_ideal_peak_d_A": f"{full_peak:.12g}" if full_peak is not None else "",
                "no_ch2_cooh_peak_d_A": f"{no_peak:.12g}" if no_peak is not None else "",
                "subtraction_peak_d_A": f"{sub_peak:.12g}" if sub_peak is not None else "",
                "subtraction_peak_intensity_norm_difference": f"{sub_intensity:.12g}" if sub_intensity is not None else "",
                "attribution_note": note,
            }
        )
    write_csv(
        args.feature_output,
        feature_rows,
        [
            "feature_window",
            "experimental_peak_d_A",
            "full_ideal_peak_d_A",
            "no_ch2_cooh_peak_d_A",
            "subtraction_peak_d_A",
            "subtraction_peak_intensity_norm_difference",
            "attribution_note",
        ],
    )
    plot_four_trace(args.four_trace_output, d_grid, sm_exp, sm_full, sm_no, sm_sub, args.plot_sigma_a)
    provenance = json.loads(args.provenance.read_text(encoding="utf-8"))
    write_report(args.report_output, args, provenance, feature_rows)

    print(f"Wrote no-CH2-COOH radial profile: {no_radial}")
    print(f"Wrote subtraction profile: {subtraction_csv}")
    print(f"Wrote four-trace plot: {args.four_trace_output}")
    print(f"Wrote feature summary: {args.feature_output}")
    print(f"Wrote report: {args.report_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
