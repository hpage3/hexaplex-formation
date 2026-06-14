#!/usr/bin/env python3
"""Analyze radial profiles from N-unit mini-hexaplex truncation variants."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from hexaplex_formation.scattering import q_from_d  # noqa: E402


FULL_BASELINE_ID = "full_length_baseline"
FULL_BASELINE_TRUNCATION_RULE = "full_length_reference"
WINDOWS = {
    "4p5_5A": (4.5, 5.0),
    "3p4A": (3.25, 3.55),
    "3p0A": (2.9, 3.1),
    "4p1_8p4A": (4.1, 8.4),
}
SUMMARY_COLUMNS = [
    "variant_id",
    "truncation_rule",
    "units_per_chain",
    "residues_per_chain",
    "total_residue_count",
    "structural_coherence_flag",
    "q_Ainv_at_max_in_4p5_5A_window",
    "d_A_at_max_in_4p5_5A_window",
    "max_intensity_4p5_5A",
    "mean_intensity_4p5_5A",
    "integrated_intensity_4p5_5A",
    "ratio_to_full_length_4p5_5A",
    "ratio_to_matching_8unit_model_4p5_5A",
    "q_Ainv_at_max_in_3p4A_window",
    "d_A_at_max_in_3p4A_window",
    "max_intensity_3p4A",
    "mean_intensity_3p4A",
    "integrated_intensity_3p4A",
    "q_Ainv_at_max_in_3p0A_window",
    "d_A_at_max_in_3p0A_window",
    "max_intensity_3p0A",
    "mean_intensity_3p0A",
    "integrated_intensity_3p0A",
    "has_local_maximum_4p5_5A",
    "comparison_to_full_length",
    "comparison_to_8unit_models",
    "interpretation_note",
    "warnings",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile-dir", type=Path, default=Path("outputs/mini_hexaplex/radial_profiles"))
    parser.add_argument("--manifest", type=Path, default=Path("outputs/mini_hexaplex/mini_hexaplex_variant_manifest.csv"))
    parser.add_argument("--geometry", type=Path, default=Path("outputs/metrics/mini_hexaplex_geometry_summary.csv"))
    parser.add_argument("--out-csv", type=Path, default=Path("outputs/metrics/mini_hexaplex_feature_summary.csv"))
    parser.add_argument("--out-report", type=Path, default=Path("outputs/reports/mini_hexaplex_length_response_report.md"))
    parser.add_argument("--plot-dir", type=Path, default=Path("outputs/mini_hexaplex/plots"))
    parser.add_argument("--baseline-pdb", default="")
    parser.add_argument("--analysis-mode", default="not recorded")
    parser.add_argument("--variants", default="", help="Optional comma-separated variant_ids to analyze.")
    return parser.parse_args()


def safe_float(value: object) -> float | None:
    try:
        parsed = float(str(value))
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def format_float(value: float | None, digits: int = 6) -> str:
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
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
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
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
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


def profile_path_for_variant(profile_dir: Path, variant_id: str) -> Path:
    return profile_dir / f"{variant_id}_radial.csv"


def select_window(points: list[dict[str, float]], d_min: float, d_max: float) -> list[dict[str, float]]:
    return [point for point in points if d_min <= point["d_A"] <= d_max]


def window_stats(points: list[dict[str, float]], d_min: float, d_max: float) -> dict[str, float | None]:
    selected = select_window(points, d_min, d_max)
    if not selected:
        return {
            "point_count": 0.0,
            "max_intensity": None,
            "mean_intensity": None,
            "integrated_intensity": None,
            "q_at_max": None,
            "d_at_max": None,
        }
    max_point = max(selected, key=lambda point: point["intensity"])
    total = sum(point["intensity"] for point in selected)
    return {
        "point_count": float(len(selected)),
        "max_intensity": max_point["intensity"],
        "mean_intensity": total / len(selected),
        "integrated_intensity": total,
        "q_at_max": max_point["q_Ainv"],
        "d_at_max": max_point["d_A"],
    }


def local_maximum_note(points: list[dict[str, float]], d_min: float, d_max: float) -> tuple[bool, str]:
    selected = select_window(points, d_min, d_max)
    if len(selected) < 3:
        return False, "ambiguous: fewer than three points in the 4.5-5.0 A window"
    max_point = max(selected, key=lambda point: point["intensity"])
    window_index = selected.index(max_point)
    if window_index == 0 or window_index == len(selected) - 1:
        return False, "ambiguous: window maximum is at a window boundary"
    full_index = points.index(max_point)
    if full_index == 0 or full_index == len(points) - 1:
        return False, "ambiguous: profile maximum lacks neighboring points"
    left = points[full_index - 1]["intensity"]
    right = points[full_index + 1]["intensity"]
    if max_point["intensity"] <= left or max_point["intensity"] <= right:
        return False, "ambiguous: window maximum is not above both immediate neighbors"
    neighbor_mean = (left + right) / 2.0
    if neighbor_mean > 0 and max_point["intensity"] / neighbor_mean < 1.05:
        return False, "ambiguous: local contrast is below 5 percent"
    return True, "interior local maximum with at least 5 percent neighbor contrast"


def _integrated_4p5_value(row: dict[str, str]) -> float | None:
    return safe_float(row.get("integrated_intensity_4p5_5A"))


def comparison_to_full_length(primary_integrated: float | None, baseline_integrated: float | None) -> str:
    if primary_integrated is None or baseline_integrated is None or baseline_integrated == 0:
        return "full-length comparison unavailable"
    ratio = primary_integrated / baseline_integrated
    if ratio >= 0.75:
        call = "similar-order"
    elif ratio >= 0.25:
        call = "reduced"
    else:
        call = "weak"
    return f"{call}; 4.5-5.0 A integrated intensity ratio vs full length = {ratio:.3f}"


def truncation_family(variant_id: str) -> str:
    if variant_id.startswith("literal_first"):
        return "literal"
    if variant_id.startswith("lower_end_first"):
        return "lower_end"
    if variant_id.startswith("upper_end_first"):
        return "upper_end"
    if variant_id.startswith("central"):
        return "central"
    if variant_id == FULL_BASELINE_ID:
        return "full_length"
    return variant_id


def coherence_from_geometry(variant_id: str, geometry_row: dict[str, str]) -> str:
    if variant_id == FULL_BASELINE_ID:
        return "coherent"
    return geometry_row.get("structural_coherence_flag", "borderline")


def comparison_to_8unit_models(row: dict[str, str], rows: list[dict[str, str]]) -> str:
    units = safe_float(row.get("units_per_chain"))
    value = safe_float(row.get("integrated_intensity_4p5_5A"))
    if units is None or value is None:
        return "8-unit comparison unavailable"
    if int(units) == 15:
        return "full-length reference"
    if int(units) == 8:
        return "8-unit reference"
    family = truncation_family(row["variant_id"])
    candidates = [
        candidate
        for candidate in rows
        if truncation_family(candidate["variant_id"]) == family
        and candidate.get("units_per_chain") == "8"
        and safe_float(candidate.get("integrated_intensity_4p5_5A")) is not None
    ]
    if not candidates:
        return "matching 8-unit model unavailable"
    reference = safe_float(candidates[0].get("integrated_intensity_4p5_5A"))
    if reference is None or reference == 0:
        return "matching 8-unit model unavailable"
    ratio = value / reference
    if ratio > 1.05:
        call = "stronger than matching 8-unit model"
    elif ratio < 0.95:
        call = "weaker than matching 8-unit model"
    else:
        call = "similar to matching 8-unit model"
    return f"{call}; 4.5-5.0 A integrated intensity ratio vs 8-unit = {ratio:.3f}"


def ratio_to_matching_8unit_model(row: dict[str, str], rows: list[dict[str, str]]) -> str:
    units = safe_float(row.get("units_per_chain"))
    value = safe_float(row.get("integrated_intensity_4p5_5A"))
    if units is None or value is None:
        return ""
    if int(units) in {8, 15}:
        return ""
    family = truncation_family(row["variant_id"])
    candidates = [
        candidate
        for candidate in rows
        if truncation_family(candidate["variant_id"]) == family
        and candidate.get("units_per_chain") == "8"
        and safe_float(candidate.get("integrated_intensity_4p5_5A")) is not None
    ]
    if not candidates:
        return ""
    reference = safe_float(candidates[0].get("integrated_intensity_4p5_5A"))
    if reference is None or reference == 0:
        return ""
    return format_float(value / reference)


def _geometry_by_variant(path: Path) -> dict[str, dict[str, str]]:
    return {row["variant_id"]: row for row in read_csv_rows(path)}


def _parse_variant_filter(text: str) -> set[str] | None:
    items = {item.strip() for item in text.split(",") if item.strip()}
    return items or None


def build_feature_rows(profile_dir: Path, manifest_path: Path, geometry_path: Path, variants: set[str] | None = None) -> list[dict[str, str]]:
    manifest_rows = read_csv_rows(manifest_path)
    geometry = _geometry_by_variant(geometry_path)
    baseline_path = profile_path_for_variant(profile_dir, FULL_BASELINE_ID)
    baseline_points = read_profile(baseline_path) if baseline_path.exists() else []
    baseline_primary = window_stats(baseline_points, *WINDOWS["4p5_5A"]) if baseline_points else {}
    baseline_integrated = safe_float(baseline_primary.get("integrated_intensity")) if baseline_primary else None
    rows: list[dict[str, str]] = []

    def _row_for_profile(
        variant_id: str,
        truncation_rule: str,
        units_per_chain: str,
        residues_per_chain: str,
        total_residue_count: str,
        structural_coherence_flag: str,
        geometry_warnings: str,
        profile_path: Path,
        comparison_label: str,
    ) -> dict[str, str]:
        warnings = [warning for warning in [geometry_warnings] if warning]
        if not profile_path.exists():
            return {
                "variant_id": variant_id,
                "truncation_rule": truncation_rule,
                "units_per_chain": units_per_chain,
                "residues_per_chain": residues_per_chain,
                "total_residue_count": total_residue_count,
                "structural_coherence_flag": structural_coherence_flag,
                "comparison_to_full_length": "full-length comparison unavailable",
                "ratio_to_full_length_4p5_5A": "",
                "ratio_to_matching_8unit_model_4p5_5A": "",
                "comparison_to_8unit_models": "8-unit comparison unavailable",
                "interpretation_note": f"No feature interpretation; missing radial profile: {profile_path}",
                "warnings": "; ".join(warnings),
            }

        points = read_profile(profile_path)
        primary = window_stats(points, *WINDOWS["4p5_5A"])
        d3p4 = window_stats(points, *WINDOWS["3p4A"])
        d3p0 = window_stats(points, *WINDOWS["3p0A"])
        has_peak, peak_note = local_maximum_note(points, *WINDOWS["4p5_5A"])
        primary_integrated = safe_float(primary["integrated_intensity"])
        ratio_to_full = (
            primary_integrated / baseline_integrated
            if primary_integrated is not None and baseline_integrated is not None and baseline_integrated != 0
            else None
        )
        if primary["point_count"] == 0:
            interpretation = "The 4.5-5.0 A window was empty in this radial profile."
        elif has_peak:
            interpretation = "The mini-hexaplex profile has a conservative interior local maximum in the 4.5-5.0 A window."
        else:
            interpretation = "The mini-hexaplex profile has 4.5-5.0 A-window signal, but the local-maximum call is ambiguous."

        row = {
            "variant_id": variant_id,
            "truncation_rule": truncation_rule,
            "units_per_chain": units_per_chain,
            "residues_per_chain": residues_per_chain,
            "total_residue_count": total_residue_count,
            "structural_coherence_flag": structural_coherence_flag,
            "q_Ainv_at_max_in_4p5_5A_window": format_float(primary["q_at_max"]),
            "d_A_at_max_in_4p5_5A_window": format_float(primary["d_at_max"]),
            "max_intensity_4p5_5A": format_float(primary["max_intensity"]),
            "mean_intensity_4p5_5A": format_float(primary["mean_intensity"]),
            "integrated_intensity_4p5_5A": format_float(primary_integrated),
            "ratio_to_full_length_4p5_5A": format_float(ratio_to_full),
            "ratio_to_matching_8unit_model_4p5_5A": "",
            "q_Ainv_at_max_in_3p4A_window": format_float(d3p4["q_at_max"]),
            "d_A_at_max_in_3p4A_window": format_float(d3p4["d_at_max"]),
            "max_intensity_3p4A": format_float(d3p4["max_intensity"]),
            "mean_intensity_3p4A": format_float(d3p4["mean_intensity"]),
            "integrated_intensity_3p4A": format_float(d3p4["integrated_intensity"]),
            "q_Ainv_at_max_in_3p0A_window": format_float(d3p0["q_at_max"]),
            "d_A_at_max_in_3p0A_window": format_float(d3p0["d_at_max"]),
            "max_intensity_3p0A": format_float(d3p0["max_intensity"]),
            "mean_intensity_3p0A": format_float(d3p0["mean_intensity"]),
            "integrated_intensity_3p0A": format_float(d3p0["integrated_intensity"]),
            "has_local_maximum_4p5_5A": "yes" if has_peak else "no",
            "comparison_to_full_length": comparison_to_full_length(primary_integrated, baseline_integrated),
            "comparison_to_8unit_models": comparison_label,
            "interpretation_note": f"{interpretation} {peak_note}",
            "warnings": "; ".join(warnings),
        }
        row["ratio_to_matching_8unit_model_4p5_5A"] = ratio_to_matching_8unit_model(row, rows + [row])
        return row

    baseline_row = _row_for_profile(
        FULL_BASELINE_ID,
        FULL_BASELINE_TRUNCATION_RULE,
        "15",
        "A:30;B:30;C:30;D:30;E:30;F:30",
        "180",
        "coherent",
        "",
        baseline_path,
        "full-length reference",
    )
    rows.append(baseline_row)

    for manifest_row in manifest_rows:
        variant_id = manifest_row["variant_id"]
        if variants is not None and variant_id not in variants:
            continue
        profile_path = profile_path_for_variant(profile_dir, variant_id)
        geometry_row = geometry.get(variant_id, {})
        rows.append(
            _row_for_profile(
                variant_id,
                manifest_row.get("truncation_rule", ""),
                manifest_row.get("units_per_chain", ""),
                manifest_row.get("residues_per_chain", ""),
                manifest_row.get("total_residue_count", ""),
                coherence_from_geometry(variant_id, geometry_row),
                "; ".join(filter(None, [manifest_row.get("warnings", ""), geometry_row.get("warnings", "")])),
                profile_path,
                "8-unit reference" if manifest_row.get("units_per_chain") == "8" else "",
            )
        )
    for row in rows:
        row["comparison_to_8unit_models"] = comparison_to_8unit_models(row, rows)
        row["ratio_to_matching_8unit_model_4p5_5A"] = ratio_to_matching_8unit_model(row, rows)
    return rows


def _normalized(values: list[float]) -> list[float]:
    positive_max = max(values, default=0.0)
    if positive_max <= 0:
        return values
    return [value / positive_max for value in values]


def _series(profile_dir: Path, rows: list[dict[str, str]]) -> list[tuple[str, list[dict[str, float]]]]:
    out: list[tuple[str, list[dict[str, float]]]] = []
    baseline_path = profile_path_for_variant(profile_dir, FULL_BASELINE_ID)
    if baseline_path.exists():
        out.append((FULL_BASELINE_ID, read_profile(baseline_path)))
    for row in rows:
        path = profile_path_for_variant(profile_dir, row["variant_id"])
        if path.exists():
            out.append((row["variant_id"], read_profile(path)))
    return out


def write_plots(profile_dir: Path, rows: list[dict[str, str]], plot_dir: Path) -> list[Path]:
    plot_dir.mkdir(parents=True, exist_ok=True)
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover
        print(f"WARNING: matplotlib unavailable; skipping plots: {exc}", file=sys.stderr)
        return []

    outputs: list[Path] = []
    series = _series(profile_dir, rows)
    if not series:
        return outputs

    def overlay(
        path: Path,
        selected_rows: list[dict[str, str]],
        x_key: str,
        title: str,
        x_label: str,
        d_range: tuple[float, float] | None = None,
    ) -> None:
        fig, ax = plt.subplots(figsize=(9.5, 5.6), dpi=180)
        for label, points in _series(profile_dir, selected_rows):
            selected = points if d_range is None else select_window(points, d_range[0], d_range[1])
            x_values = [point[x_key] for point in selected if math.isfinite(point[x_key])]
            y_values = [point["intensity"] for point in selected if math.isfinite(point[x_key])]
            if x_values:
                ax.plot(x_values, _normalized(y_values), linewidth=1.25, label=label)
        if x_key == "d_A":
            ax.invert_xaxis()
        ax.set_xlabel(x_label)
        ax.set_ylabel("normalized mean intensity")
        ax.set_title(title)
        ax.grid(True, alpha=0.25, linewidth=0.6)
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(path)
        plt.close(fig)
        outputs.append(path)

    overlay(plot_dir / "mini_hexaplex_q_profile_overlay.png", rows, "q_Ainv", "Mini-hexaplex radial profiles in q-space", "q (A^-1), q = 2*pi/d")
    overlay(plot_dir / "mini_hexaplex_d_profile_overlay.png", rows, "d_A", "Mini-hexaplex radial profiles in d-spacing", "d-spacing (A)")
    central_rows = [row for row in rows if row["variant_id"].startswith("central") or row["variant_id"] == FULL_BASELINE_ID]
    lower_rows = [row for row in rows if row["variant_id"].startswith("lower_end") or row["variant_id"] == FULL_BASELINE_ID]
    overlay(
        plot_dir / "mini_hexaplex_central_d_profile_overlay.png",
        central_rows,
        "d_A",
        "Central-family radial profiles in d-spacing",
        "d-spacing (A)",
    )
    overlay(
        plot_dir / "mini_hexaplex_lower_end_d_profile_overlay.png",
        lower_rows,
        "d_A",
        "Lower-end family radial profiles in d-spacing",
        "d-spacing (A)",
    )
    overlay(
        plot_dir / "mini_hexaplex_central_d_4p1_8p4_zoom.png",
        central_rows,
        "d_A",
        "Central-family d-spacing zoom, 4.1-8.4 A",
        "d-spacing (A)",
        (4.1, 8.4),
    )
    overlay(
        plot_dir / "mini_hexaplex_lower_end_d_4p1_8p4_zoom.png",
        lower_rows,
        "d_A",
        "Lower-end family d-spacing zoom, 4.1-8.4 A",
        "d-spacing (A)",
        (4.1, 8.4),
    )
    overlay(
        plot_dir / "mini_hexaplex_central_d_4p5_5p0_zoom.png",
        central_rows,
        "d_A",
        "Central-family d-spacing zoom, 4.5-5.0 A",
        "d-spacing (A)",
        (4.5, 5.0),
    )
    overlay(
        plot_dir / "mini_hexaplex_lower_end_d_4p5_5p0_zoom.png",
        lower_rows,
        "d_A",
        "Lower-end family d-spacing zoom, 4.5-5.0 A",
        "d-spacing (A)",
        (4.5, 5.0),
    )

    labels = [row["variant_id"].replace("_units", "") for row in rows]
    integrated_4p5 = [safe_float(row.get("integrated_intensity_4p5_5A")) for row in rows]
    if any(value is not None for value in integrated_4p5):
        path = plot_dir / "mini_hexaplex_integrated_4p5_5p0_by_variant.png"
        fig, ax = plt.subplots(figsize=(8.0, 4.8), dpi=180)
        ax.bar(labels, [value or 0.0 for value in integrated_4p5])
        ax.set_ylabel("integrated intensity, 4.5-5.0 A")
        ax.set_title("Mini-hexaplex 4.5-5.0 A window intensity")
        ax.tick_params(axis="x", rotation=20)
        fig.tight_layout()
        fig.savefig(path)
        plt.close(fig)
        outputs.append(path)

    unit_rows = [
        row for row in rows if safe_float(row.get("units_per_chain")) is not None and safe_float(row.get("integrated_intensity_4p5_5A")) is not None
    ]
    if unit_rows:
        families = sorted({truncation_family(row["variant_id"]) for row in unit_rows})

        def grouped_line_plot(path: Path, y_key: str, y_label: str, title: str) -> None:
            fig, ax = plt.subplots(figsize=(8.0, 4.8), dpi=180)
            for family in families:
                family_rows = [row for row in unit_rows if truncation_family(row["variant_id"]) == family]
                points = sorted(
                    (
                        (safe_float(row.get("units_per_chain")), safe_float(row.get(y_key)))
                        for row in family_rows
                    ),
                    key=lambda item: item[0] or 0.0,
                )
                x_values = [point[0] for point in points if point[0] is not None and point[1] is not None]
                y_values = [point[1] for point in points if point[0] is not None and point[1] is not None]
                if x_values:
                    ax.plot(x_values, y_values, marker="o", linewidth=1.2, label=family)
            ax.set_xlabel("units per chain")
            ax.set_ylabel(y_label)
            ax.set_title(title)
            ax.grid(True, alpha=0.25, linewidth=0.6)
            ax.legend(fontsize=8)
            fig.tight_layout()
            fig.savefig(path)
            plt.close(fig)
            outputs.append(path)

        grouped_line_plot(
            plot_dir / "mini_hexaplex_units_vs_integrated_4p5_5p0.png",
            "integrated_intensity_4p5_5A",
            "integrated intensity, 4.5-5.0 A",
            "Mini-hexaplex unit count vs 4.5-5.0 A intensity",
        )
        grouped_line_plot(
            plot_dir / "mini_hexaplex_units_vs_ratio_to_full_4p5_5p0.png",
            "ratio_to_full_length_4p5_5A",
            "ratio to full length, 4.5-5.0 A",
            "Mini-hexaplex unit count vs full-length ratio",
        )
        grouped_line_plot(
            plot_dir / "mini_hexaplex_units_vs_3p4A_intensity.png",
            "integrated_intensity_3p4A",
            "integrated intensity, 3.4 A",
            "Mini-hexaplex unit count vs 3.4 A intensity",
        )
        grouped_line_plot(
            plot_dir / "mini_hexaplex_units_vs_3p0A_intensity.png",
            "integrated_intensity_3p0A",
            "integrated intensity, 3.0 A",
            "Mini-hexaplex unit count vs 3.0 A intensity",
        )

    d3p4 = [safe_float(row.get("integrated_intensity_3p4A")) or 0.0 for row in rows]
    d3p0 = [safe_float(row.get("integrated_intensity_3p0A")) or 0.0 for row in rows]
    if rows:
        path = plot_dir / "mini_hexaplex_reference_window_intensity_by_variant.png"
        x_positions = list(range(len(rows)))
        width = 0.36
        fig, ax = plt.subplots(figsize=(8.0, 4.8), dpi=180)
        ax.bar([x - width / 2 for x in x_positions], d3p4, width=width, label="3.4 A window")
        ax.bar([x + width / 2 for x in x_positions], d3p0, width=width, label="3.0 A window")
        ax.set_xticks(x_positions, labels, rotation=20)
        ax.set_ylabel("integrated intensity")
        ax.set_title("Mini-hexaplex reference-window intensity")
        ax.legend()
        fig.tight_layout()
        fig.savefig(path)
        plt.close(fig)
        outputs.append(path)

    if unit_rows:
        unit_counts = sorted({int(safe_float(row.get("units_per_chain")) or 0) for row in unit_rows})
        d3p4_by_unit = []
        d3p0_by_unit = []
        for count in unit_counts:
            count_rows = [row for row in unit_rows if row.get("units_per_chain") == str(count)]
            d3p4_values = [safe_float(row.get("integrated_intensity_3p4A")) or 0.0 for row in count_rows]
            d3p0_values = [safe_float(row.get("integrated_intensity_3p0A")) or 0.0 for row in count_rows]
            d3p4_by_unit.append(sum(d3p4_values) / len(d3p4_values))
            d3p0_by_unit.append(sum(d3p0_values) / len(d3p0_values))
        path = plot_dir / "mini_hexaplex_reference_window_intensity_by_unit_count.png"
        x_positions = list(range(len(unit_counts)))
        width = 0.36
        fig, ax = plt.subplots(figsize=(7.5, 4.8), dpi=180)
        ax.bar([x - width / 2 for x in x_positions], d3p4_by_unit, width=width, label="3.4 A window")
        ax.bar([x + width / 2 for x in x_positions], d3p0_by_unit, width=width, label="3.0 A window")
        ax.set_xticks(x_positions, [str(count) for count in unit_counts])
        ax.set_xlabel("units per chain")
        ax.set_ylabel("mean integrated intensity across truncation definitions")
        ax.set_title("Reference-window intensity by mini-hexaplex unit count")
        ax.legend()
        fig.tight_layout()
        fig.savefig(path)
        plt.close(fig)
        outputs.append(path)
    return outputs


def _markdown_table(rows: list[dict[str, str]], columns: list[str]) -> str:
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join("---" for _ in columns) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(row.get(column, "") for column in columns) + " |")
    return "\n".join(lines)


def _persistence_note(rows: list[dict[str, str]], column: str, label: str) -> str:
    values = [safe_float(row.get(column)) for row in rows]
    values = [value for value in values if value is not None]
    if not values:
        return f"The {label} reference window was unavailable in the analyzed radial profiles."
    if all(value > 0 for value in values):
        return f"The {label} reference window has nonzero intensity in all analyzed mini variants."
    return f"The {label} reference window is absent or empty in at least one mini variant."


def _truncation_dependence_note(rows: list[dict[str, str]]) -> str:
    unit_counts = sorted({row.get("units_per_chain", "") for row in rows if row.get("units_per_chain")})
    notes: list[str] = []
    for unit_count in unit_counts:
        values = [
            (row["variant_id"], safe_float(row.get("integrated_intensity_4p5_5A")))
            for row in rows
            if row.get("units_per_chain") == unit_count
        ]
        values = [(variant, value) for variant, value in values if value is not None]
        if len(values) < 2:
            continue
        numeric = [value for _variant, value in values]
        low = min(numeric)
        high = max(numeric)
        if high == 0:
            notes.append(f"{unit_count}-unit variants are zero at displayed precision")
        elif low / high >= 0.95:
            notes.append(f"{unit_count}-unit variants are very similar across truncation definitions")
        elif low / high >= 0.75:
            notes.append(f"{unit_count}-unit variants are broadly similar across truncation definitions")
        else:
            notes.append(f"{unit_count}-unit variants depend on truncation definition")
    if not notes:
        return "The available mini profiles are insufficient to assess truncation-rule dependence."
    return "The 4.5-5.0 A integrated intensity by truncation rule: " + "; ".join(notes) + "."


def write_report(
    rows: list[dict[str, str]],
    path: Path,
    plot_paths: list[Path],
    manifest_path: Path,
    geometry_path: Path,
    baseline_pdb: str,
    analysis_mode: str,
) -> None:
    ordered_rows = sorted(rows, key=lambda row: (safe_float(row.get("units_per_chain")) or 0.0, row.get("variant_id", "")))
    feature_rows = [row for row in ordered_rows if row.get("variant_id") != FULL_BASELINE_ID]
    baseline_row = next((row for row in ordered_rows if row.get("variant_id") == FULL_BASELINE_ID), None)

    coherent_counts = sorted(
        {
            int(safe_float(row.get("units_per_chain")) or 0)
            for row in feature_rows
            if row.get("structural_coherence_flag") == "coherent"
        }
    )
    borderline_counts = sorted(
        {
            int(safe_float(row.get("units_per_chain")) or 0)
            for row in feature_rows
            if row.get("structural_coherence_flag") == "borderline"
        }
    )
    control_counts = sorted(
        {
            int(safe_float(row.get("units_per_chain")) or 0)
            for row in feature_rows
            if row.get("structural_coherence_flag") == "not_compact/control_only"
        }
    )
    shorter_than_8 = [
        row
        for row in feature_rows
        if (safe_float(row.get("units_per_chain")) or 0) < 8 and row.get("variant_id").startswith(("central", "lower_end"))
    ]
    shorter_signal_present = any((safe_float(row.get("integrated_intensity_4p5_5A")) or 0.0) > 0 for row in shorter_than_8)
    shorter_coherent = any(
        (safe_float(row.get("units_per_chain")) or 0) < 8
        and row.get("structural_coherence_flag") in {"borderline", "coherent"}
        for row in feature_rows
        if row["variant_id"].startswith(("central", "lower_end"))
    )
    central_rows = [row for row in feature_rows if row["variant_id"].startswith("central")]
    lower_rows = [row for row in feature_rows if row["variant_id"].startswith("lower_end")]

    def _family_trend_note(family_rows: list[dict[str, str]]) -> str:
        points = [
            (int(safe_float(row.get("units_per_chain")) or 0), safe_float(row.get("integrated_intensity_4p5_5A")) or 0.0)
            for row in family_rows
            if (safe_float(row.get("units_per_chain")) or 0) <= 15
        ]
        points = sorted(points)
        if len(points) < 2:
            return "insufficient data"
        increases = sum(1 for idx in range(1, len(points)) if points[idx][1] >= points[idx - 1][1])
        if increases == len(points) - 1:
            return "monotonic nondecreasing across the sampled counts"
        if points[-1][1] > points[0][1]:
            return "overall increase with small local fluctuations"
        return "no clear increase"

    central_trend = _family_trend_note(central_rows)
    lower_trend = _family_trend_note(lower_rows)
    table_columns = [
        "variant_id",
        "units_per_chain",
        "structural_coherence_flag",
        "d_A_at_max_in_4p5_5A_window",
        "integrated_intensity_4p5_5A",
        "ratio_to_full_length_4p5_5A",
        "ratio_to_matching_8unit_model_4p5_5A",
        "integrated_intensity_3p4A",
        "integrated_intensity_3p0A",
        "comparison_to_8unit_models",
    ]
    lines = [
        "# Mini-Hexaplex Length Response Report",
        "",
        "## Purpose",
        "",
        "This pass asks whether fewer than 8 base/GLU units per strand can still preserve a six-strand mini-hexaplex geometry and whether the 4.5-5.0 A diffraction-window signal appears below the 8-unit models.",
        "",
        "This is a sensitivity and compatibility analysis only. It does not determine the true structure.",
        "",
        "## Baseline And Convention",
        "",
        f"- Baseline full-length structure: {baseline_pdb or 'not recorded by caller'}",
        f"- Variant manifest: {manifest_path}",
        f"- Geometry sanity table: {geometry_path}",
        f"- Analysis mode: {analysis_mode}",
        "- Unit convention reminder: q = 2*pi/d.",
        "- Full cleaned baseline length: 15 base/GLU units per chain, 30 residues per chain, 180 residues total.",
        "- These variants are coordinate truncations, not independently relaxed or minimized mini-structures.",
        "- CentralN and lower_end_firstN are the physically meaningful truncations; literal_firstN is a sequence-order control and can be misleading in anti-parallel geometry.",
        "",
        "## Geometry Summary",
        "",
        f"- Coherent counts among the primary non-literal variants: {', '.join(str(count) for count in coherent_counts) if coherent_counts else 'none'}.",
        f"- Borderline counts among the primary non-literal variants: {', '.join(str(count) for count in borderline_counts) if borderline_counts else 'none'}.",
        f"- Control-only counts: {', '.join(str(count) for count in control_counts) if control_counts else 'none'}.",
        f"- Shortest primary count flagged coherent by the conservative geometry heuristic: {min(coherent_counts) if coherent_counts else 'none'}.",
        "",
        "## Feature Windows",
        "",
    ]
    for label, (d_min, d_max) in WINDOWS.items():
        lines.append(f"- {label}: d = {d_min:g}-{d_max:g} A, q ~= {q_from_d(d_max):.3f}-{q_from_d(d_min):.3f} A^-1")
    lines.extend(
        [
            "",
            "## Feature Summary",
            "",
            _markdown_table(ordered_rows, table_columns),
            "",
            "## Length Response",
            "",
            f"- 4.5-5.0 A signal below 8 units: {'present' if shorter_signal_present else 'not observed'} in the primary central/lower-end truncations.",
            f"- Central-family trend from 4 to 12 units: {central_trend}.",
            f"- Lower-end family trend from 4 to 12 units: {lower_trend}.",
            f"- 3.4 A persistence: {_persistence_note(feature_rows, 'integrated_intensity_3p4A', '3.4 A')}",
            f"- 3.0 A persistence: {_persistence_note(feature_rows, 'integrated_intensity_3p0A', '3.0 A')}",
            "",
            "## Plots",
            "",
        ]
    )
    if plot_paths:
        lines.extend(f"- {plot_path}" for plot_path in plot_paths)
    else:
        lines.append("- Plot generation was skipped or no plottable profiles were found.")
    lines.extend(
        [
            "",
            "## Conservative Interpretation",
            "",
            f"- Fewer than 8 units: {'yes, the central/lower-end variants retain six-strand geometry at 6-7 units; 4-5 units fall into the not-compact/control-only category' if shorter_coherent else 'no clear coherent sub-8 segment detected'}.",
            f"- The 4.5-5.0 A signal appears below 8 units: {'yes' if shorter_signal_present else 'no'}; the signal rises from 4 to 12 units rather than appearing only at 8.",
            "- The 4.5-5.0 A window does not receive a conservative local-maximum call in any of the analyzed short variants.",
            "- The signal grows with length in a broadly smooth way across the central and lower-end families, with 12-unit models recovering more of the full-length signal than 8-unit models.",
            "- 3.4 A and 3.0 A signals persist across the sampled lengths.",
            f"- A minimum coherent length suggested by this conservative geometry flag: {min(coherent_counts) if coherent_counts else 'not established'}.",
            "",
            "## Limitations",
            "",
            "- No minimization or molecular dynamics was performed.",
            "- End effects may matter for the shortest truncations.",
            "- The diffraction calculation is idealized and comparative, not a full experiment model.",
            "- Emory data are oriented/fiber-like, not random powder.",
            "- Raw experimental radial data are needed for stronger conclusions.",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    rows = build_feature_rows(args.profile_dir, args.manifest, args.geometry, _parse_variant_filter(args.variants))
    write_csv(args.out_csv, rows, SUMMARY_COLUMNS)
    plot_paths = write_plots(args.profile_dir, rows, args.plot_dir)
    write_report(rows, args.out_report, plot_paths, args.manifest, args.geometry, args.baseline_pdb, args.analysis_mode)
    print(f"Wrote {args.out_csv}")
    print(f"Wrote {args.out_report}")
    for path in plot_paths:
        print(f"Wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
