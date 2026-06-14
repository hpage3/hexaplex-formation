#!/usr/bin/env python3
"""Analyze radial diffraction profiles from base-length variants."""

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


MAIN_WINDOW = (4.5, 5.0)
REFERENCE_WINDOWS = {
    "d_4p5_5p0": (4.5, 5.0),
    "d_4p1": (3.95, 4.25),
    "d_5p5": (5.35, 5.75),
    "d_7p0": (6.75, 7.25),
    "d_8p4": (8.1, 8.7),
    "d_3p4": (3.25, 3.55),
    "d_3p0": (2.9, 3.1),
}
SUMMARY_COLUMNS = [
    "variant_id",
    "scale_factor",
    "geometry_warning",
    "q_Ainv_at_max_in_4p5_5A_window",
    "d_A_at_max_in_4p5_5A_window",
    "max_intensity_4p5_5A",
    "mean_intensity_4p5_5A",
    "integrated_intensity_4p5_5A",
    "max_intensity_3p4A",
    "mean_intensity_3p4A",
    "integrated_intensity_3p4A",
    "max_intensity_3p0A",
    "mean_intensity_3p0A",
    "integrated_intensity_3p0A",
    "has_local_maximum_4p5_5A",
    "local_maximum_note",
    "interpretation_note",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--profile-dir",
        type=Path,
        default=Path("outputs/base_length_sweep/radial_profiles"),
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("outputs/base_length_sweep/structures/base_length_variant_manifest.csv"),
    )
    parser.add_argument(
        "--geometry",
        type=Path,
        default=Path("outputs/metrics/base_length_variant_geometry.csv"),
    )
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=Path("outputs/metrics/base_length_sweep_feature_summary.csv"),
    )
    parser.add_argument(
        "--out-report",
        type=Path,
        default=Path("outputs/reports/base_length_sweep_report.md"),
    )
    parser.add_argument(
        "--plot-dir",
        type=Path,
        default=Path("outputs/base_length_sweep/plots"),
    )
    parser.add_argument("--scales", default="", help="Optional comma-separated scale factors to analyze.")
    parser.add_argument("--baseline-pdb", default="")
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


def scale_to_token(scale: float | str) -> str:
    return f"{float(scale):.2f}".replace(".", "p")


def variant_id_for_scale(scale: float | str) -> str:
    return f"hexaplex_base_length_scale_{scale_to_token(scale)}"


def parse_scale_list(text: str) -> list[float]:
    values: list[float] = []
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        values.append(float(item))
    return values


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


def profile_path_for_variant(profile_dir: Path, variant_id: str) -> Path:
    return profile_dir / f"{variant_id}_radial.csv"


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


def select_window(points: list[dict[str, float]], d_min: float, d_max: float) -> list[dict[str, float]]:
    if d_max < d_min:
        raise ValueError("d_max must be greater than or equal to d_min")
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

    ordered = points
    full_index = ordered.index(max_point)
    if full_index == 0 or full_index == len(ordered) - 1:
        return False, "ambiguous: profile maximum lacks neighboring points"
    left = ordered[full_index - 1]["intensity"]
    right = ordered[full_index + 1]["intensity"]
    if max_point["intensity"] <= left or max_point["intensity"] <= right:
        return False, "ambiguous: window maximum is not above both immediate neighbors"

    neighbor_mean = (left + right) / 2.0
    if neighbor_mean > 0 and max_point["intensity"] / neighbor_mean < 1.05:
        return False, "ambiguous: local contrast is below 5 percent"
    return True, "interior local maximum with at least 5 percent neighbor contrast"


def geometry_by_variant(path: Path) -> dict[str, dict[str, str]]:
    return {row["variant_id"]: row for row in read_csv_rows(path)}


def manifest_rows_for_scales(path: Path, scales: list[float]) -> list[dict[str, str]]:
    rows = read_csv_rows(path)
    if not scales:
        return rows
    wanted = {variant_id_for_scale(scale) for scale in scales}
    return [row for row in rows if row.get("variant_id") in wanted]


def build_feature_rows(
    profile_dir: Path,
    manifest_path: Path,
    geometry_path: Path,
    scales: list[float],
) -> list[dict[str, str]]:
    geometry = geometry_by_variant(geometry_path)
    rows: list[dict[str, str]] = []
    for manifest_row in manifest_rows_for_scales(manifest_path, scales):
        variant_id = manifest_row["variant_id"]
        profile_path = profile_path_for_variant(profile_dir, variant_id)
        scale_factor = manifest_row.get("scale_factor", "")
        geometry_row = geometry.get(variant_id, {})
        geometry_warning = geometry_row.get("warnings", "")
        if not profile_path.exists():
            rows.append(
                {
                    "variant_id": variant_id,
                    "scale_factor": scale_factor,
                    "geometry_warning": geometry_warning,
                    "local_maximum_note": f"missing radial profile: {profile_path}",
                    "interpretation_note": "No feature interpretation; radial profile was unavailable.",
                }
            )
            continue

        points = read_profile(profile_path)
        primary = window_stats(points, *MAIN_WINDOW)
        d3p4 = window_stats(points, *REFERENCE_WINDOWS["d_3p4"])
        d3p0 = window_stats(points, *REFERENCE_WINDOWS["d_3p0"])
        has_peak, peak_note = local_maximum_note(points, *MAIN_WINDOW)
        if primary["point_count"] == 0:
            interpretation = "The 4.5-5.0 A window was empty in this radial profile."
        elif has_peak:
            interpretation = "A conservative interior local maximum was detected in the 4.5-5.0 A window."
        else:
            interpretation = "The 4.5-5.0 A window has signal, but the local-maximum call is ambiguous."

        rows.append(
            {
                "variant_id": variant_id,
                "scale_factor": scale_factor,
                "geometry_warning": geometry_warning,
                "q_Ainv_at_max_in_4p5_5A_window": format_float(primary["q_at_max"]),
                "d_A_at_max_in_4p5_5A_window": format_float(primary["d_at_max"]),
                "max_intensity_4p5_5A": format_float(primary["max_intensity"]),
                "mean_intensity_4p5_5A": format_float(primary["mean_intensity"]),
                "integrated_intensity_4p5_5A": format_float(primary["integrated_intensity"]),
                "max_intensity_3p4A": format_float(d3p4["max_intensity"]),
                "mean_intensity_3p4A": format_float(d3p4["mean_intensity"]),
                "integrated_intensity_3p4A": format_float(d3p4["integrated_intensity"]),
                "max_intensity_3p0A": format_float(d3p0["max_intensity"]),
                "mean_intensity_3p0A": format_float(d3p0["mean_intensity"]),
                "integrated_intensity_3p0A": format_float(d3p0["integrated_intensity"]),
                "has_local_maximum_4p5_5A": "yes" if has_peak else "no",
                "local_maximum_note": peak_note,
                "interpretation_note": interpretation,
            }
        )
    return rows


def _profile_series(profile_dir: Path, rows: list[dict[str, str]]) -> list[tuple[str, list[dict[str, float]]]]:
    series = []
    for row in rows:
        path = profile_path_for_variant(profile_dir, row["variant_id"])
        if path.exists():
            series.append((row["scale_factor"], read_profile(path)))
    return series


def _normalized(values: list[float]) -> list[float]:
    positive_max = max(values, default=0.0)
    if positive_max <= 0:
        return values
    return [value / positive_max for value in values]


def write_plots(profile_dir: Path, rows: list[dict[str, str]], plot_dir: Path) -> list[Path]:
    plot_dir.mkdir(parents=True, exist_ok=True)
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover - depends on optional runtime dependency
        print(f"WARNING: matplotlib unavailable; skipping plots: {exc}", file=sys.stderr)
        return []

    outputs: list[Path] = []
    series = _profile_series(profile_dir, rows)
    if not series:
        return outputs

    def overlay(path: Path, x_key: str, title: str, x_label: str, d_range: tuple[float, float] | None = None) -> None:
        fig, ax = plt.subplots(figsize=(9.5, 5.6), dpi=180)
        for scale, points in series:
            selected = points
            if d_range is not None:
                selected = select_window(points, d_range[0], d_range[1])
            x_values = [point[x_key] for point in selected if math.isfinite(point[x_key])]
            y_values = [point["intensity"] for point in selected if math.isfinite(point[x_key])]
            if not x_values:
                continue
            ax.plot(x_values, _normalized(y_values), linewidth=1.2, label=f"scale {scale}")
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

    overlay(
        plot_dir / "base_length_sweep_q_profile_overlay.png",
        "q_Ainv",
        "Base-length sweep radial profiles in q-space",
        "q (A^-1), convention q = 2*pi/d",
    )
    overlay(
        plot_dir / "base_length_sweep_d_profile_overlay.png",
        "d_A",
        "Base-length sweep radial profiles in d-spacing",
        "d-spacing (A)",
    )
    overlay(
        plot_dir / "base_length_sweep_d_4p1_8p4_zoom.png",
        "d_A",
        "Base-length sweep d-spacing zoom, 4.1-8.4 A",
        "d-spacing (A)",
        d_range=(4.1, 8.4),
    )
    overlay(
        plot_dir / "base_length_sweep_d_4p5_5p0_zoom.png",
        "d_A",
        "Base-length sweep d-spacing zoom, 4.5-5.0 A",
        "d-spacing (A)",
        d_range=(4.5, 5.0),
    )

    scale_values = []
    integrated_values = []
    peak_scales = []
    peak_d_values = []
    for row in rows:
        scale = safe_float(row.get("scale_factor"))
        integrated = safe_float(row.get("integrated_intensity_4p5_5A"))
        peak_d = safe_float(row.get("d_A_at_max_in_4p5_5A_window"))
        if scale is not None and integrated is not None:
            scale_values.append(scale)
            integrated_values.append(integrated)
        if scale is not None and peak_d is not None and row.get("has_local_maximum_4p5_5A") == "yes":
            peak_scales.append(scale)
            peak_d_values.append(peak_d)

    if scale_values:
        path = plot_dir / "base_length_scale_vs_integrated_4p5_5p0.png"
        fig, ax = plt.subplots(figsize=(7.5, 4.8), dpi=180)
        ax.plot(scale_values, integrated_values, marker="o", linewidth=1.2)
        ax.set_xlabel("base-length scale factor")
        ax.set_ylabel("integrated intensity, 4.5-5.0 A window")
        ax.set_title("Base-length scale vs 4.5-5.0 A window intensity")
        ax.grid(True, alpha=0.25, linewidth=0.6)
        fig.tight_layout()
        fig.savefig(path)
        plt.close(fig)
        outputs.append(path)

    if peak_scales:
        path = plot_dir / "base_length_scale_vs_peak_d_4p5_5p0.png"
        fig, ax = plt.subplots(figsize=(7.5, 4.8), dpi=180)
        ax.plot(peak_scales, peak_d_values, marker="o", linewidth=1.2)
        ax.set_xlabel("base-length scale factor")
        ax.set_ylabel("peak d-spacing in 4.5-5.0 A window")
        ax.set_title("Base-length scale vs conservative 4.5-5.0 A peak position")
        ax.grid(True, alpha=0.25, linewidth=0.6)
        fig.tight_layout()
        fig.savefig(path)
        plt.close(fig)
        outputs.append(path)
    return outputs


def _markdown_table(rows: list[dict[str, str]], columns: list[str]) -> str:
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(row.get(column, "") for column in columns) + " |")
    return "\n".join(lines)


def _trend_note(rows: list[dict[str, str]]) -> str:
    pairs = [
        (safe_float(row.get("scale_factor")), safe_float(row.get("integrated_intensity_4p5_5A")))
        for row in rows
    ]
    pairs = [(scale, intensity) for scale, intensity in pairs if scale is not None and intensity is not None]
    if len(pairs) < 2:
        return "The available profiles are insufficient to assess a scale trend."
    pairs.sort()
    first = pairs[0][1]
    last = pairs[-1][1]
    if last > first:
        return "The integrated 4.5-5.0 A-window intensity increases across the analyzed scale range."
    if last < first:
        return "The integrated 4.5-5.0 A-window intensity decreases across the analyzed scale range."
    return "The integrated 4.5-5.0 A-window intensity is unchanged at the displayed precision."


def _reference_persistence_note(rows: list[dict[str, str]], column: str, label: str) -> str:
    values = [safe_float(row.get(column)) for row in rows]
    values = [value for value in values if value is not None]
    if not values:
        return f"The {label} reference window was unavailable in the analyzed radial profiles."
    if all(value > 0 for value in values):
        return f"The {label} reference window has nonzero intensity in all analyzed variants."
    return f"The {label} reference window is absent or empty in at least one analyzed variant."


def write_report(
    rows: list[dict[str, str]],
    path: Path,
    plot_paths: list[Path],
    manifest_path: Path,
    geometry_path: Path,
    baseline_pdb: str,
) -> None:
    manifest = read_csv_rows(manifest_path)
    geometry = geometry_by_variant(geometry_path)
    analyzed_ids = {row["variant_id"] for row in rows}
    analyzed_manifest = [row for row in manifest if row.get("variant_id") in analyzed_ids]
    transformed_counts = sorted({row.get("transformed_atom_count", "") for row in analyzed_manifest if row.get("transformed_atom_count")})
    fixed_counts = sorted({row.get("fixed_atom_count", "") for row in analyzed_manifest if row.get("fixed_atom_count")})
    scale_values = [row.get("scale_factor", "") for row in rows]
    warnings = [geometry[row["variant_id"]].get("warnings", "") for row in rows if geometry.get(row["variant_id"], {}).get("warnings")]
    has_peak = any(row.get("has_local_maximum_4p5_5A") == "yes" for row in rows)
    lines = [
        "# Base-Length Diffraction Sweep Report",
        "",
        "## Purpose",
        "",
        "This controlled sweep asks whether changing local CYP/MEP base/hexad-arm length in a six-strand hexaplex model generates or strengthens a simulated diffraction feature near d ~= 4.5-5.0 A.",
        "",
        "This is a computational sensitivity study only. It does not determine the structure.",
        "",
        "## Scientific cautions",
        "",
        "- The d ~= 4.5 A feature is reciprocal-space-like, not a literal atom-distance assignment.",
        "- The unit convention is q = 2*pi/d.",
        "- The simulation is powder-averaged and simplified relative to oriented/fiber-like experimental arcs.",
        "- The variants are idealized and depend on the atom-selection and local-anchor transformation rule.",
        "- CYP/MEP candidate arm atoms are axis-facing in the fitted-axis inspection; the transform scales local anchor-to-atom vectors, not global outward radial vectors.",
        "- Raw experimental radial data are not included; comparisons are limited to q or d-spacing regions.",
        "",
        "## Inputs and variants",
        "",
        f"- Baseline structure: {baseline_pdb or 'not recorded by caller'}",
        f"- Variant manifest: {manifest_path}",
        f"- Geometry sanity table: {geometry_path}",
        f"- Scale factors analyzed: {', '.join(scale_values)}",
        f"- Transformed atom count(s): {', '.join(transformed_counts) if transformed_counts else 'not available'}",
        f"- Fixed atom count(s): {', '.join(fixed_counts) if fixed_counts else 'not available'}",
        "- Scale 1.20 is excluded from default interpretation because geometry checks found suspicious heavy-atom overlaps below 1.00 A. It should be treated only as a stress-test variant if explicitly included.",
        "",
        "Operationally, base/hexad-arm length is the local distance from each CYP/MEP residue anchor to selected non-backbone candidate arm atoms. GLU atoms and CYP/MEP backbone-like atoms remain fixed.",
        "",
        "## Feature windows",
        "",
    ]
    for label, (d_min, d_max) in REFERENCE_WINDOWS.items():
        lines.append(f"- {label}: d = {d_min:g}-{d_max:g} A, q ~= {q_from_d(d_max):.3f}-{q_from_d(d_min):.3f} A^-1")
    lines.extend(
        [
            "",
            "## Summary table",
            "",
            _markdown_table(
                rows,
                [
                    "variant_id",
                    "scale_factor",
                    "geometry_warning",
                    "d_A_at_max_in_4p5_5A_window",
                    "integrated_intensity_4p5_5A",
                    "has_local_maximum_4p5_5A",
                    "integrated_intensity_3p4A",
                    "integrated_intensity_3p0A",
                ],
            ),
            "",
            "## Plots",
            "",
        ]
    )
    if plot_paths:
        lines.extend(f"- {path}" for path in plot_paths)
    else:
        lines.append("- Plot generation was skipped or no plottable profiles were found.")
    lines.extend(
        [
            "",
            "## Conservative interpretation",
            "",
            f"- 4.5-5.0 A local maximum call: {'at least one analyzed variant has a conservative local maximum' if has_peak else 'no analyzed variant has a conservative local maximum; window signal should be treated as ambiguous'}.",
            f"- {_trend_note(rows)}",
            f"- {_reference_persistence_note(rows, 'integrated_intensity_3p4A', '3.4 A')}",
            f"- {_reference_persistence_note(rows, 'integrated_intensity_3p0A', '3.0 A')}",
            "",
            "In this simplified comparative dataset, a change in the 4.5-5.0 A-window score with base-length scale is consistent with sensitivity to local CYP/MEP arm geometry. It does not establish causality or structural identity.",
            "",
            "## Geometry warnings",
            "",
        ]
    )
    if warnings:
        lines.extend(f"- {warning}" for warning in warnings)
    else:
        lines.append("- No geometry warnings were reported for the analyzed variants.")
    lines.extend(
        [
            "",
            "## Limitations",
            "",
            "- Powder-averaged simulation is not equivalent to Emory oriented/fiber-like arcs.",
            "- Structural variants are idealized local-anchor perturbations.",
            "- Results depend on the selected CYP/MEP transformable atoms.",
            "- Raw experimental radial profiles are unavailable, so this is not a direct fit.",
            "- Scale 1.20 has a known overlap warning and is not part of the default main interpretation.",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    scales = parse_scale_list(args.scales)
    rows = build_feature_rows(args.profile_dir, args.manifest, args.geometry, scales)
    write_csv(args.out_csv, rows, SUMMARY_COLUMNS)
    plot_paths = write_plots(args.profile_dir, rows, args.plot_dir)
    write_report(rows, args.out_report, plot_paths, args.manifest, args.geometry, args.baseline_pdb)
    print(f"Wrote {args.out_csv}")
    print(f"Wrote {args.out_report}")
    for path in plot_paths:
        print(f"Wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
