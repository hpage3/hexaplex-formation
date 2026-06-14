#!/usr/bin/env python3
"""Summarize length/twist diffraction sensitivity without inventing twist models."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from hexaplex_formation.pdb_utils import chain_ids, load_pdb_atoms  # noqa: E402
from hexaplex_formation.scattering import q_from_d  # noqa: E402


WINDOWS = {
    "d_4p5_5p0": ("4.5-5.0 A", 4.5, 5.0),
    "d_3p4": ("3.4 A", 3.25, 3.55),
    "d_3p0": ("3.0 A", 2.90, 3.10),
    "d_4p1_8p4": ("4.1-8.4 A", 4.1, 8.4),
}

DEFAULT_LENGTH_MODELS = [
    ("length_4_central_twist_30", "central4_units", 4, "outputs/mini_hexaplex/structures/mini_hexaplex_central4_units.pdb"),
    ("length_7_central_twist_30", "central7_units", 7, "outputs/mini_hexaplex/structures/mini_hexaplex_central7_units.pdb"),
    ("length_8_central_twist_30", "central8_units", 8, "outputs/mini_hexaplex/structures/mini_hexaplex_central8_units.pdb"),
    ("length_12_central_twist_30", "central12_units", 12, "outputs/mini_hexaplex/structures/mini_hexaplex_central12_units.pdb"),
    (
        "full_length_twist_30",
        "full_length_baseline",
        15,
        "outputs/intermediates/ai_candidate_inputs/full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain.pdb",
    ),
]

PLANNED_TWISTS = [24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0]

FEATURE_FIELDNAMES = [
    "model_id",
    "length_units",
    "twist_deg",
    "atom_count",
    "chain_count",
    "q_Ainv_at_max_in_4p5_5p0_A_window",
    "d_A_at_max_in_4p5_5p0_A_window",
    "max_intensity_4p5_5p0_A",
    "mean_intensity_4p5_5p0_A",
    "integrated_intensity_4p5_5p0_A",
    "fwhm_q_Ainv_4p5_5p0_A",
    "width_proxy_q_Ainv_4p5_5p0_A",
    "second_moment_width_q_Ainv_4p5_5p0_A",
    "q_Ainv_at_max_in_3p4_A_window",
    "d_A_at_max_in_3p4_A_window",
    "max_intensity_3p4_A",
    "mean_intensity_3p4_A",
    "integrated_intensity_3p4_A",
    "fwhm_q_Ainv_3p4_A",
    "width_proxy_q_Ainv_3p4_A",
    "q_Ainv_at_max_in_3p0_A_window",
    "d_A_at_max_in_3p0_A_window",
    "max_intensity_3p0_A",
    "mean_intensity_3p0_A",
    "integrated_intensity_3p0_A",
    "fwhm_q_Ainv_3p0_A",
    "width_proxy_q_Ainv_3p0_A",
    "ratio_integrated_4p5_5p0_A_to_3p4_A",
    "notes",
    "warnings",
]

WIDTH_FIELDNAMES = [
    "model_id",
    "length_units",
    "twist_deg",
    "window_name",
    "d_min_A",
    "d_max_A",
    "q_min_Ainv",
    "q_max_Ainv",
    "point_count",
    "q_Ainv_at_max",
    "d_A_at_max",
    "max_intensity",
    "integrated_intensity_sum",
    "integrated_intensity_trapezoid_q",
    "fwhm_q_Ainv",
    "half_max_sample_width_q_Ainv",
    "equivalent_width_q_Ainv",
    "second_moment_width_q_Ainv",
    "width_note",
]

MANIFEST_FIELDNAMES = [
    "model_id",
    "length_units",
    "twist_deg",
    "atom_count",
    "chain_count",
    "model_path",
    "radial_profile_path",
    "native_detector_npy",
    "model_status",
    "source",
    "notes",
    "warnings",
]


@dataclass(frozen=True)
class Model:
    model_id: str
    radial_id: str
    length_units: int
    twist_deg: float
    model_path: Path
    radial_profile_path: Path
    native_detector_npy: Path
    status: str
    source: str
    notes: str = ""
    warnings: str = ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile-dir", type=Path, default=Path("outputs/mini_hexaplex/radial_profiles"))
    parser.add_argument("--diffraction-dir", type=Path, default=Path("outputs/mini_hexaplex/diffraction"))
    parser.add_argument("--plot-dir", type=Path, default=Path("outputs/length_twist_diffraction/plots"))
    parser.add_argument("--builder", type=Path, default=None, help="Optional official PNAB/proto-nucleic builder executable or directory.")
    parser.add_argument("--summary-csv", type=Path, default=Path("outputs/metrics/length_twist_feature_summary.csv"))
    parser.add_argument("--manifest-csv", type=Path, default=Path("outputs/metrics/length_twist_model_manifest.csv"))
    parser.add_argument("--widths-csv", type=Path, default=Path("outputs/metrics/length_twist_peak_widths.csv"))
    parser.add_argument(
        "--contact-summary-csv",
        type=Path,
        default=Path("outputs/metrics/length_twist_contact_network_summary.csv"),
    )
    parser.add_argument("--report", type=Path, default=Path("outputs/reports/length_twist_diffraction_sensitivity_report.md"))
    parser.add_argument("--include-planned-twist-rows", action="store_true", default=True)
    return parser.parse_args()


def safe_float(value: object) -> float | None:
    try:
        parsed = float(str(value))
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def fmt(value: float | None, digits: int = 6) -> str:
    if value is None or not math.isfinite(value):
        return ""
    return f"{value:.{digits}f}"


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def intensity_from_row(row: dict[str, str]) -> float | None:
    for key in ("intensity", "mean_intensity", "intensity_mean"):
        value = safe_float(row.get(key))
        if value is not None:
            return value
    return None


def read_profile(path: Path) -> list[dict[str, float]]:
    points: list[dict[str, float]] = []
    for row in read_csv_rows(path):
        q_value = safe_float(row.get("q_Ainv") or row.get("q_center_inv_angstrom"))
        d_value = safe_float(row.get("d_A") or row.get("d_center_angstrom"))
        intensity = intensity_from_row(row)
        sample_count = safe_float(row.get("sample_count") or row.get("pixel_count") or "1")
        if q_value is None or d_value is None or intensity is None:
            continue
        if sample_count is not None and sample_count <= 0:
            continue
        points.append({"q_Ainv": q_value, "d_A": d_value, "intensity": intensity})
    return sorted(points, key=lambda point: point["q_Ainv"])


def select_window(points: list[dict[str, float]], d_min: float, d_max: float) -> list[dict[str, float]]:
    return [point for point in points if d_min <= point["d_A"] <= d_max]


def trapezoid_integral_q(points: list[dict[str, float]]) -> float | None:
    ordered = sorted(points, key=lambda point: point["q_Ainv"])
    if len(ordered) < 2:
        return None
    total = 0.0
    for left, right in zip(ordered[:-1], ordered[1:]):
        total += (right["q_Ainv"] - left["q_Ainv"]) * (left["intensity"] + right["intensity"]) / 2.0
    return total


def half_max_crossing_q(left: dict[str, float], right: dict[str, float], half_max: float) -> float:
    i0 = left["intensity"]
    i1 = right["intensity"]
    if i1 == i0:
        return (left["q_Ainv"] + right["q_Ainv"]) / 2.0
    fraction = (half_max - i0) / (i1 - i0)
    return left["q_Ainv"] + fraction * (right["q_Ainv"] - left["q_Ainv"])


def estimate_fwhm_q(points: list[dict[str, float]], max_point: dict[str, float]) -> tuple[float | None, float | None, str]:
    ordered = sorted(points, key=lambda point: point["q_Ainv"])
    if len(ordered) < 3 or max_point not in ordered:
        return None, None, "fewer than three profile points; FWHM unavailable"
    max_index = ordered.index(max_point)
    half_max = max_point["intensity"] / 2.0
    above = [point for point in ordered if point["intensity"] >= half_max]
    sample_width = (max(point["q_Ainv"] for point in above) - min(point["q_Ainv"] for point in above)) if above else None

    left_cross = None
    for index in range(max_index - 1, -1, -1):
        if ordered[index]["intensity"] <= half_max <= ordered[index + 1]["intensity"]:
            left_cross = half_max_crossing_q(ordered[index], ordered[index + 1], half_max)
            break
    right_cross = None
    for index in range(max_index, len(ordered) - 1):
        if ordered[index]["intensity"] >= half_max >= ordered[index + 1]["intensity"]:
            right_cross = half_max_crossing_q(ordered[index], ordered[index + 1], half_max)
            break
    if left_cross is None or right_cross is None:
        return None, sample_width, "half-maximum crossing is not bracketed on both sides inside the scored window"
    return abs(right_cross - left_cross), sample_width, "linear-interpolated FWHM within the scored simulated radial window"


def second_moment_width_q(points: list[dict[str, float]]) -> float | None:
    weights = [max(0.0, point["intensity"]) for point in points]
    total = sum(weights)
    if total <= 0:
        return None
    mean_q = sum(point["q_Ainv"] * weight for point, weight in zip(points, weights)) / total
    variance = sum(weight * (point["q_Ainv"] - mean_q) ** 2 for point, weight in zip(points, weights)) / total
    return math.sqrt(max(0.0, variance))


def window_metrics(points: list[dict[str, float]], d_min: float, d_max: float) -> dict[str, float | str | None]:
    selected = select_window(points, d_min, d_max)
    if not selected:
        return {
            "point_count": 0,
            "max_intensity": None,
            "mean_intensity": None,
            "integrated_intensity": None,
            "q_at_max": None,
            "d_at_max": None,
            "trapz_q": None,
            "fwhm_q": None,
            "half_max_sample_width_q": None,
            "equivalent_width_q": None,
            "second_moment_width_q": None,
            "width_note": "empty scored d-spacing window",
        }
    max_point = max(selected, key=lambda point: point["intensity"])
    total = sum(point["intensity"] for point in selected)
    trapz = trapezoid_integral_q(selected)
    fwhm, sample_width, width_note = estimate_fwhm_q(selected, max_point)
    equivalent_width = trapz / max_point["intensity"] if trapz is not None and max_point["intensity"] > 0 else None
    return {
        "point_count": len(selected),
        "max_intensity": max_point["intensity"],
        "mean_intensity": total / len(selected),
        "integrated_intensity": total,
        "q_at_max": max_point["q_Ainv"],
        "d_at_max": max_point["d_A"],
        "trapz_q": trapz,
        "fwhm_q": fwhm,
        "half_max_sample_width_q": sample_width,
        "equivalent_width_q": equivalent_width,
        "second_moment_width_q": second_moment_width_q(selected),
        "width_note": width_note,
    }


def model_atom_chain_counts(path: Path) -> tuple[int | None, int | None, str]:
    if not path.exists():
        return None, None, f"missing model PDB: {path}"
    atoms = load_pdb_atoms(path)
    return len(atoms), len(chain_ids(atoms)), ""


def find_builder(explicit: Path | None) -> tuple[bool, str]:
    if explicit is not None:
        return explicit.exists(), str(explicit)
    candidates = []
    for root in (REPO_ROOT, REPO_ROOT.parent):
        for pattern in ("*PNAB*", "*pnab*", "*Proto*Nucleic*", "*proto*nucleic*", "*builder*"):
            candidates.extend(root.glob(pattern))
    existing = [path for path in candidates if path.exists()]
    return bool(existing), "; ".join(str(path) for path in existing)


def build_models(args: argparse.Namespace, builder_available: bool) -> list[Model]:
    models: list[Model] = []
    for model_id, radial_id, length_units, model_path in DEFAULT_LENGTH_MODELS:
        models.append(
            Model(
                model_id=model_id,
                radial_id=radial_id,
                length_units=length_units,
                twist_deg=30.0,
                model_path=Path(model_path),
                radial_profile_path=args.profile_dir / f"{radial_id}_radial.csv",
                native_detector_npy=args.diffraction_dir / f"{radial_id}.npy",
                status="available" if (args.profile_dir / f"{radial_id}_radial.csv").exists() else "missing_radial_profile",
                source="current 30 degree baseline length sweep",
            )
        )
    if args.include_planned_twist_rows:
        for twist in PLANNED_TWISTS:
            if twist == 30.0:
                continue
            models.append(
                Model(
                    model_id=f"full_length_twist_{int(twist):02d}_planned",
                    radial_id=f"full_length_twist_{int(twist):02d}",
                    length_units=15,
                    twist_deg=twist,
                    model_path=Path(f"outputs/length_twist_diffraction/structures/full_length_twist_{int(twist):02d}.pdb"),
                    radial_profile_path=args.profile_dir / f"full_length_twist_{int(twist):02d}_radial.csv",
                    native_detector_npy=args.diffraction_dir / f"full_length_twist_{int(twist):02d}.npy",
                    status="pending_official_builder" if not builder_available else "pending_generation",
                    source="planned official PNAB/proto-nucleic builder twist sweep",
                    notes="placeholder row; coordinates not generated by this workflow",
                    warnings="" if builder_available else "official builder not found locally",
                )
            )
    return models


def manifest_rows(models: list[Model]) -> list[dict[str, str]]:
    rows = []
    for model in models:
        atom_count, chain_count, warning = model_atom_chain_counts(model.model_path)
        warnings = "; ".join(item for item in [model.warnings, warning] if item)
        rows.append(
            {
                "model_id": model.model_id,
                "length_units": str(model.length_units),
                "twist_deg": fmt(model.twist_deg, 2),
                "atom_count": str(atom_count) if atom_count is not None else "",
                "chain_count": str(chain_count) if chain_count is not None else "",
                "model_path": str(model.model_path),
                "radial_profile_path": str(model.radial_profile_path),
                "native_detector_npy": str(model.native_detector_npy),
                "model_status": model.status,
                "source": model.source,
                "notes": model.notes,
                "warnings": warnings,
            }
        )
    return rows


def feature_and_width_rows(models: list[Model]) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    features: list[dict[str, str]] = []
    widths: list[dict[str, str]] = []
    for model in models:
        atom_count, chain_count, warning = model_atom_chain_counts(model.model_path)
        points = read_profile(model.radial_profile_path) if model.radial_profile_path.exists() else []
        metrics_by_window = {key: window_metrics(points, d_min, d_max) for key, (_label, d_min, d_max) in WINDOWS.items()}
        primary = metrics_by_window["d_4p5_5p0"]
        d3p4 = metrics_by_window["d_3p4"]
        d3p0 = metrics_by_window["d_3p0"]
        primary_integrated = safe_float(primary["integrated_intensity"])
        d3p4_integrated = safe_float(d3p4["integrated_intensity"])
        ratio = primary_integrated / d3p4_integrated if primary_integrated is not None and d3p4_integrated else None
        notes = model.notes
        warnings = "; ".join(item for item in [model.warnings, warning] if item)
        if not points:
            warnings = "; ".join(item for item in [warnings, f"missing radial profile: {model.radial_profile_path}"] if item)
        features.append(
            {
                "model_id": model.model_id,
                "length_units": str(model.length_units),
                "twist_deg": fmt(model.twist_deg, 2),
                "atom_count": str(atom_count) if atom_count is not None else "",
                "chain_count": str(chain_count) if chain_count is not None else "",
                "q_Ainv_at_max_in_4p5_5p0_A_window": fmt(safe_float(primary["q_at_max"])),
                "d_A_at_max_in_4p5_5p0_A_window": fmt(safe_float(primary["d_at_max"])),
                "max_intensity_4p5_5p0_A": fmt(safe_float(primary["max_intensity"])),
                "mean_intensity_4p5_5p0_A": fmt(safe_float(primary["mean_intensity"])),
                "integrated_intensity_4p5_5p0_A": fmt(primary_integrated),
                "fwhm_q_Ainv_4p5_5p0_A": fmt(safe_float(primary["fwhm_q"])),
                "width_proxy_q_Ainv_4p5_5p0_A": fmt(safe_float(primary["equivalent_width_q"])),
                "second_moment_width_q_Ainv_4p5_5p0_A": fmt(safe_float(primary["second_moment_width_q"])),
                "q_Ainv_at_max_in_3p4_A_window": fmt(safe_float(d3p4["q_at_max"])),
                "d_A_at_max_in_3p4_A_window": fmt(safe_float(d3p4["d_at_max"])),
                "max_intensity_3p4_A": fmt(safe_float(d3p4["max_intensity"])),
                "mean_intensity_3p4_A": fmt(safe_float(d3p4["mean_intensity"])),
                "integrated_intensity_3p4_A": fmt(d3p4_integrated),
                "fwhm_q_Ainv_3p4_A": fmt(safe_float(d3p4["fwhm_q"])),
                "width_proxy_q_Ainv_3p4_A": fmt(safe_float(d3p4["equivalent_width_q"])),
                "q_Ainv_at_max_in_3p0_A_window": fmt(safe_float(d3p0["q_at_max"])),
                "d_A_at_max_in_3p0_A_window": fmt(safe_float(d3p0["d_at_max"])),
                "max_intensity_3p0_A": fmt(safe_float(d3p0["max_intensity"])),
                "mean_intensity_3p0_A": fmt(safe_float(d3p0["mean_intensity"])),
                "integrated_intensity_3p0_A": fmt(safe_float(d3p0["integrated_intensity"])),
                "fwhm_q_Ainv_3p0_A": fmt(safe_float(d3p0["fwhm_q"])),
                "width_proxy_q_Ainv_3p0_A": fmt(safe_float(d3p0["equivalent_width_q"])),
                "ratio_integrated_4p5_5p0_A_to_3p4_A": fmt(ratio),
                "notes": notes,
                "warnings": warnings,
            }
        )
        for key, (label, d_min, d_max) in WINDOWS.items():
            metric = metrics_by_window[key]
            widths.append(
                {
                    "model_id": model.model_id,
                    "length_units": str(model.length_units),
                    "twist_deg": fmt(model.twist_deg, 2),
                    "window_name": label,
                    "d_min_A": fmt(d_min, 2),
                    "d_max_A": fmt(d_max, 2),
                    "q_min_Ainv": fmt(q_from_d(d_max)),
                    "q_max_Ainv": fmt(q_from_d(d_min)),
                    "point_count": str(metric["point_count"]),
                    "q_Ainv_at_max": fmt(safe_float(metric["q_at_max"])),
                    "d_A_at_max": fmt(safe_float(metric["d_at_max"])),
                    "max_intensity": fmt(safe_float(metric["max_intensity"])),
                    "integrated_intensity_sum": fmt(safe_float(metric["integrated_intensity"])),
                    "integrated_intensity_trapezoid_q": fmt(safe_float(metric["trapz_q"])),
                    "fwhm_q_Ainv": fmt(safe_float(metric["fwhm_q"])),
                    "half_max_sample_width_q_Ainv": fmt(safe_float(metric["half_max_sample_width_q"])),
                    "equivalent_width_q_Ainv": fmt(safe_float(metric["equivalent_width_q"])),
                    "second_moment_width_q_Ainv": fmt(safe_float(metric["second_moment_width_q"])),
                    "width_note": str(metric["width_note"]),
                }
            )
    return features, widths


def write_contact_subset(path: Path) -> None:
    source = Path("outputs/metrics/seed_contact_network_summary.csv")
    rows = []
    wanted = {"central4_units", "central7_units", "central8_units", "central12_units"}
    for row in read_csv_rows(source):
        if row.get("variant_id") in wanted:
            row = dict(row)
            row["model_id"] = f"length_{row.get('units_per_chain')}_central_twist_30"
            rows.append(row)
    if rows:
        fieldnames = ["model_id"] + [name for name in rows[0].keys() if name != "model_id"]
        write_csv(path, rows, fieldnames)


def write_plots(features: list[dict[str, str]], profile_dir: Path, plot_dir: Path) -> list[Path]:
    available = [row for row in features if row["twist_deg"] == "30.00" and row["integrated_intensity_4p5_5p0_A"]]
    if not available:
        return []
    plot_dir.mkdir(parents=True, exist_ok=True)
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover
        print(f"WARNING: matplotlib unavailable; writing SVG fallback plots: {exc}", file=sys.stderr)
        return write_svg_fallback_plots(available, profile_dir, plot_dir)

    radial_by_model = {
        model_id: read_profile(profile_dir / radial_name)
        for model_id, radial_name in {
            "length_4_central_twist_30": "central4_units_radial.csv",
            "length_7_central_twist_30": "central7_units_radial.csv",
            "length_8_central_twist_30": "central8_units_radial.csv",
            "length_12_central_twist_30": "central12_units_radial.csv",
            "full_length_twist_30": "full_length_baseline_radial.csv",
        }.items()
        if (profile_dir / radial_name).exists()
    }
    outputs: list[Path] = []

    def normalized(values: list[float]) -> list[float]:
        max_value = max(values, default=0.0)
        return [value / max_value for value in values] if max_value > 0 else values

    def overlay(path: Path, title: str, d_range: tuple[float, float] | None = None) -> None:
        fig, ax = plt.subplots(figsize=(9.0, 5.2), dpi=180)
        for model_id, points in radial_by_model.items():
            selected = points if d_range is None else select_window(points, *d_range)
            if not selected:
                continue
            x_values = [point["d_A"] for point in selected]
            y_values = normalized([point["intensity"] for point in selected])
            ax.plot(x_values, y_values, linewidth=1.2, marker="o" if d_range == (4.5, 5.0) else None, markersize=2.5, label=model_id)
        ax.invert_xaxis()
        ax.set_xlabel("d-spacing (A)")
        ax.set_ylabel("normalized radial intensity")
        ax.set_title(title)
        ax.grid(True, alpha=0.25, linewidth=0.6)
        ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(path)
        plt.close(fig)
        outputs.append(path)

    overlay(plot_dir / "length_twist_length_overlay_d_profile.png", "Length sweep at 30 degree twist")
    overlay(plot_dir / "length_twist_length_overlay_d_4p1_8p4_zoom.png", "Length sweep, d = 4.1-8.4 A", (4.1, 8.4))
    overlay(plot_dir / "length_twist_length_overlay_d_4p5_5p0_zoom.png", "Length sweep, d = 4.5-5.0 A", (4.5, 5.0))

    def line_plot(path: Path, y_key: str, ylabel: str, title: str) -> None:
        points = sorted(
            [(safe_float(row["length_units"]), safe_float(row[y_key])) for row in available],
            key=lambda item: item[0] or 0.0,
        )
        fig, ax = plt.subplots(figsize=(7.2, 4.6), dpi=180)
        ax.plot([point[0] for point in points], [point[1] for point in points], marker="o", linewidth=1.4)
        ax.set_xlabel("units per chain")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(True, alpha=0.25, linewidth=0.6)
        fig.tight_layout()
        fig.savefig(path)
        plt.close(fig)
        outputs.append(path)

    line_plot(
        plot_dir / "length_twist_width_proxy_vs_length.png",
        "width_proxy_q_Ainv_4p5_5p0_A",
        "simulated equivalent-width proxy in q (A^-1)",
        "4.5-5.0 A simulated width proxy versus length",
    )
    line_plot(
        plot_dir / "length_twist_4p5_5p0_to_3p4_ratio_vs_length.png",
        "ratio_integrated_4p5_5p0_A_to_3p4_A",
        "integrated intensity ratio",
        "4.5-5.0 A / 3.4 A response ratio versus length",
    )
    return outputs


def svg_polyline(points: list[tuple[float, float]], color: str) -> str:
    if not points:
        return ""
    joined = " ".join(f"{x:.2f},{y:.2f}" for x, y in points)
    return f'<polyline points="{joined}" fill="none" stroke="{color}" stroke-width="1.7" />'


def write_svg(path: Path, title: str, body: list[str], width: int = 900, height: int = 520) -> None:
    path.write_text(
        "\n".join(
            [
                f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
                '<rect width="100%" height="100%" fill="white" />',
                f'<text x="24" y="30" font-family="Arial, sans-serif" font-size="18">{title}</text>',
                *body,
                "</svg>",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def scale_points(
    points: list[tuple[float, float]],
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
    invert_x: bool = False,
) -> list[tuple[float, float]]:
    left, top, width, height = 70.0, 55.0, 760.0, 380.0
    if x_max == x_min:
        x_max = x_min + 1.0
    if y_max == y_min:
        y_max = y_min + 1.0
    scaled = []
    for x_value, y_value in points:
        x_fraction = (x_value - x_min) / (x_max - x_min)
        if invert_x:
            x_fraction = 1.0 - x_fraction
        y_fraction = (y_value - y_min) / (y_max - y_min)
        scaled.append((left + x_fraction * width, top + (1.0 - y_fraction) * height))
    return scaled


def write_svg_fallback_plots(features: list[dict[str, str]], profile_dir: Path, plot_dir: Path) -> list[Path]:
    plot_dir.mkdir(parents=True, exist_ok=True)
    colors = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e", "#17becf"]
    radial_files = {
        "length_4_central_twist_30": "central4_units_radial.csv",
        "length_7_central_twist_30": "central7_units_radial.csv",
        "length_8_central_twist_30": "central8_units_radial.csv",
        "length_12_central_twist_30": "central12_units_radial.csv",
        "full_length_twist_30": "full_length_baseline_radial.csv",
    }
    radial_by_model = {
        model_id: read_profile(profile_dir / filename)
        for model_id, filename in radial_files.items()
        if (profile_dir / filename).exists()
    }
    outputs: list[Path] = []

    def normalized(values: list[float]) -> list[float]:
        max_value = max(values, default=0.0)
        return [value / max_value for value in values] if max_value > 0 else values

    def axes_body(x_label: str, y_label: str) -> list[str]:
        return [
            '<line x1="70" y1="435" x2="830" y2="435" stroke="#333" stroke-width="1" />',
            '<line x1="70" y1="55" x2="70" y2="435" stroke="#333" stroke-width="1" />',
            f'<text x="420" y="492" font-family="Arial, sans-serif" font-size="13">{x_label}</text>',
            f'<text x="15" y="250" transform="rotate(-90 15,250)" font-family="Arial, sans-serif" font-size="13">{y_label}</text>',
        ]

    def overlay_svg(path: Path, title: str, d_range: tuple[float, float] | None = None) -> None:
        series = []
        for index, (model_id, points) in enumerate(radial_by_model.items()):
            selected = points if d_range is None else select_window(points, *d_range)
            if selected:
                y_values = normalized([point["intensity"] for point in selected])
                series.append((model_id, list(zip([point["d_A"] for point in selected], y_values)), colors[index % len(colors)]))
        all_points = [point for _label, points, _color in series for point in points]
        if not all_points:
            return
        x_values = [point[0] for point in all_points]
        body = axes_body("d-spacing (A)", "normalized radial intensity")
        for label, points, color in series:
            body.append(svg_polyline(scale_points(points, min(x_values), max(x_values), 0.0, 1.0, invert_x=True), color))
            legend_y = 62 + 18 * len([item for item in body if "legend" in item])
            body.append(f'<text data-kind="legend" x="640" y="{legend_y}" font-family="Arial, sans-serif" font-size="11" fill="{color}">{label}</text>')
        write_svg(path, title, body)
        outputs.append(path)

    overlay_svg(plot_dir / "length_twist_length_overlay_d_profile.svg", "Length sweep at 30 degree twist")
    overlay_svg(plot_dir / "length_twist_length_overlay_d_4p1_8p4_zoom.svg", "Length sweep, d = 4.1-8.4 A", (4.1, 8.4))
    overlay_svg(plot_dir / "length_twist_length_overlay_d_4p5_5p0_zoom.svg", "Length sweep, d = 4.5-5.0 A", (4.5, 5.0))

    def line_svg(path: Path, y_key: str, y_label: str, title: str) -> None:
        points = sorted(
            [
                (safe_float(row["length_units"]), safe_float(row[y_key]))
                for row in features
                if safe_float(row["length_units"]) is not None and safe_float(row[y_key]) is not None
            ],
            key=lambda item: item[0] or 0.0,
        )
        if not points:
            return
        x_values = [point[0] or 0.0 for point in points]
        y_values = [point[1] or 0.0 for point in points]
        scaled = scale_points([(x or 0.0, y or 0.0) for x, y in points], min(x_values), max(x_values), min(y_values), max(y_values))
        body = axes_body("units per chain", y_label)
        body.append(svg_polyline(scaled, "#1f77b4"))
        for x, y in scaled:
            body.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="3.2" fill="#1f77b4" />')
        write_svg(path, title, body)
        outputs.append(path)

    line_svg(
        plot_dir / "length_twist_width_proxy_vs_length.svg",
        "width_proxy_q_Ainv_4p5_5p0_A",
        "simulated equivalent-width proxy in q (A^-1)",
        "4.5-5.0 A simulated width proxy versus length",
    )
    line_svg(
        plot_dir / "length_twist_4p5_5p0_to_3p4_ratio_vs_length.svg",
        "ratio_integrated_4p5_5p0_A_to_3p4_A",
        "integrated intensity ratio",
        "4.5-5.0 A / 3.4 A response ratio versus length",
    )
    return outputs


def q_window_table() -> list[dict[str, str]]:
    rows = []
    for _key, (label, d_min, d_max) in WINDOWS.items():
        rows.append(
            {
                "feature": label,
                "d_min_A": fmt(d_min, 2),
                "d_max_A": fmt(d_max, 2),
                "q_min_Ainv": fmt(q_from_d(d_max)),
                "q_max_Ainv": fmt(q_from_d(d_min)),
            }
        )
    return rows


def markdown_table(rows: list[dict[str, str]], columns: list[str]) -> str:
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join("---" for _ in columns) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(row.get(column, "") for column in columns) + " |")
    return "\n".join(lines)


def write_builder_stub(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        """#!/usr/bin/env bash
set -euo pipefail

cat >&2 <<'EOF'
The official Proto-Nucleic Acids Building / PNAB program was not found locally.

Provide the official builder and baseline input parameters, then replace the
PNAB_BUILDER placeholder below with the actual executable and option names.

Required inputs:
- current 30 degree six-strand baseline parameter file from the original builder
- sequence/residue specification matching full_hexaplex_anti_parallel_30deg_ideal
- strand count fixed at six
- helical twist values, for example: 24,26,28,30,32,34,36
- output PDB path for each generated twist variant

Proposed command template:
  PNAB_BUILDER --input baseline_builder_parameters.json \\
    --strand-count 6 \\
    --helical-twist-deg ${TWIST_DEG} \\
    --output-pdb outputs/length_twist_diffraction/structures/full_length_twist_${TWIST_DEG}.pdb

This stub intentionally exits without generating coordinates.
EOF
exit 1
""",
        encoding="utf-8",
    )
    path.chmod(0o755)


def write_report(
    path: Path,
    features: list[dict[str, str]],
    manifest: list[dict[str, str]],
    plot_paths: list[Path],
    builder_available: bool,
    builder_location: str,
) -> None:
    length_rows = [row for row in features if row["twist_deg"] == "30.00" and row["integrated_intensity_4p5_5p0_A"]]
    planned = [row for row in manifest if row["model_status"].startswith("pending")]
    lines = [
        "# Length/Twist Diffraction Sensitivity Report",
        "",
        "## Purpose",
        "",
        "This workflow extends the controlled forward-model sensitivity study for stack length and helical twist. It compares simulated radial profiles in q/d space and does not attempt structure determination.",
        "",
        "## Experimental Target Windows",
        "",
        "The q convention is q = 2*pi/d. Native detector/image-plate `.npy` outputs remain detector-space arrays and are not compared directly to experimental d-spacings.",
        "",
        markdown_table(q_window_table(), ["feature", "d_min_A", "d_max_A", "q_min_Ainv", "q_max_Ainv"]),
        "",
        "## Baseline Model",
        "",
        "- Baseline PDB: outputs/intermediates/ai_candidate_inputs/full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain.pdb",
        "- Baseline twist: 30 degrees between adjacent stacked hexads.",
        "- Current length variants are coordinate truncations from the 30 degree baseline; they are not relaxed structures.",
        "",
        "## Length Sweep At 30 Degrees",
        "",
        markdown_table(
            length_rows,
            [
                "model_id",
                "length_units",
                "d_A_at_max_in_4p5_5p0_A_window",
                "integrated_intensity_4p5_5p0_A",
                "width_proxy_q_Ainv_4p5_5p0_A",
                "integrated_intensity_3p4_A",
                "integrated_intensity_3p0_A",
                "ratio_integrated_4p5_5p0_A_to_3p4_A",
            ],
        ),
        "",
        "## Helical Twist Sweep",
        "",
    ]
    if builder_available:
        lines.append(f"- Official/proto builder candidate found: {builder_location}")
        lines.append("- Twist generation is still staged through the official builder integration path so parameters can be audited before running.")
    else:
        lines.append("- Official Proto-Nucleic Acids Building / PNAB program was not found in this repo or the nearby parent directory.")
        lines.append("- Twist variants are represented as pending manifest rows only. No chemically questionable twist models were fabricated.")
    if planned:
        lines.extend(["", markdown_table(planned, ["model_id", "length_units", "twist_deg", "model_status", "notes", "warnings"])])
    lines.extend(
        [
            "",
            "## Feature Definition And Width Metric",
            "",
            "- The primary feature is the simulated 4.5-5.0 A radial-window response.",
            "- Integrated intensity is reported as the sum of radial-bin intensities inside the scored window, matching the prior length-response convention.",
            "- `width_proxy_q_Ainv` is the trapezoid-integrated intensity in q divided by the local maximum inside the window.",
            "- FWHM is reported only when both half-maximum crossings are bracketed inside the scored simulated window.",
            "- These are simulated width proxies, not Scherrer/domain-size estimates.",
            "",
            "## Conservative Interpretation",
            "",
            "- Nick's conceptual expectation is represented here as a testable trend: longer ordered stacks should produce more defined simulated diffraction features.",
            "- The 3.4 A window remains a reference stacking-associated feature in the current profiles.",
            "- The 3.0 A and 4.5-5.0 A windows are tracked as compatibility metrics, not as fitted structural constraints.",
            "- The Emory images are fiber-like/oriented, while these simulations are powder-averaged; this mismatch limits direct inference.",
            "- A lower bound on coherent stack length should be treated as approximate and model-dependent.",
            "",
            "## Plots",
            "",
        ]
    )
    lines.extend(f"- {plot_path}" for plot_path in plot_paths) if plot_paths else lines.append("- Plot generation skipped or no profiles available.")
    lines.extend(
        [
            "",
            "## Outputs",
            "",
            "- outputs/metrics/length_twist_model_manifest.csv",
            "- outputs/metrics/length_twist_feature_summary.csv",
            "- outputs/metrics/length_twist_peak_widths.csv",
            "- outputs/metrics/length_twist_contact_network_summary.csv when contact-network inputs are available",
            "- scripts/run_pnab_twist_builder_template.sh",
            "",
            "## Next Steps",
            "",
            "1. Provide the official builder and exact baseline builder parameter files from the original publication/SI.",
            "2. Generate single-parameter twist variants first, preferably full length and/or 8-unit.",
            "3. Run the same detector simulation and radial averaging on official twist variants.",
            "4. Only then run a small length-by-twist grid such as lengths 4, 7, 8, 12, full and twists 28, 30, 32 degrees.",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_metadata(path: Path, builder_available: bool, builder_location: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            {
                "q_convention": "q_Ainv = 2*pi/d_A",
                "native_npy_convention": "detector/image-plate intensity in mm coordinates; unchanged",
                "builder_available": builder_available,
                "builder_location": builder_location,
                "target_windows": q_window_table(),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )


def main() -> int:
    args = parse_args()
    builder_available, builder_location = find_builder(args.builder)
    models = build_models(args, builder_available)
    manifest = manifest_rows(models)
    features, widths = feature_and_width_rows(models)
    write_csv(args.manifest_csv, manifest, MANIFEST_FIELDNAMES)
    write_csv(args.summary_csv, features, FEATURE_FIELDNAMES)
    write_csv(args.widths_csv, widths, WIDTH_FIELDNAMES)
    write_contact_subset(args.contact_summary_csv)
    plot_paths = write_plots(features, args.profile_dir, args.plot_dir)
    write_builder_stub(Path("scripts/run_pnab_twist_builder_template.sh"))
    write_report(args.report, features, manifest, plot_paths, builder_available, builder_location)
    write_metadata(args.summary_csv.with_suffix(".metadata.json"), builder_available, builder_location)
    print(f"Wrote {args.manifest_csv}")
    print(f"Wrote {args.summary_csv}")
    print(f"Wrote {args.widths_csv}")
    print(f"Wrote {args.report}")
    for path in plot_paths:
        print(f"Wrote {path}")
    if not builder_available:
        print("Official PNAB/proto-nucleic builder not found; twist rows are placeholders only.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
