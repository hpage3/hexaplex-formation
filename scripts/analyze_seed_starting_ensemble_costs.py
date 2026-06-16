#!/usr/bin/env python3
"""Audit loose-start costs against formed seed endpoint order parameters."""

from __future__ import annotations

import argparse
import csv
import math
from collections import OrderedDict
from pathlib import Path


FEATURE_GROUPS: OrderedDict[str, list[str]] = OrderedDict(
    [
        (
            "geometric",
            [
                "RMSD_to_formed_seed_A",
                "compactness_score",
                "radius_of_gyration_A",
                "radial_extent_A",
                "axial_extent_A",
                "helical_axis_alignment_score",
            ],
        ),
        (
            "contact_recovery",
            [
                "contact_fraction_vs_target",
                "CYP_MEP_contact_fraction_vs_target",
                "backbone_contact_fraction_vs_target",
            ],
        ),
        (
            "register_phase",
            [
                "axial_register_score",
                "angular_phase_order_score",
                "refined_angular_phase_score",
            ],
        ),
    ]
)

COMPONENT_COLUMNS = [
    "unit_count",
    "start_class",
    "feature_group",
    "feature_name",
    "loose_sample_count",
    "formed_reference_sample_count",
    "loose_mean",
    "formed_reference_mean",
    "formed_reference_std",
    "effective_reference_std",
    "std_floor_used",
    "signed_z_deviation",
    "abs_z_deviation",
]

SUMMARY_COLUMNS = [
    "unit_count",
    "start_class",
    "n_features_total",
    "n_geometric_features",
    "n_contact_features",
    "n_register_phase_features",
    "geometric_cost",
    "contact_recovery_cost",
    "register_phase_cost",
    "overall_cost",
    "dominant_cost_group",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-csv", type=Path, default=Path("outputs/metrics/seed_formation_order_parameters.csv"))
    parser.add_argument("--components-csv", type=Path, default=Path("outputs/metrics/seed_starting_ensemble_cost_components.csv"))
    parser.add_argument("--summary-csv", type=Path, default=Path("outputs/metrics/seed_starting_ensemble_cost_summary.csv"))
    parser.add_argument("--report", type=Path, default=Path("outputs/reports/seed_starting_ensemble_cost_report.md"))
    parser.add_argument("--plot-dir", type=Path, default=Path("outputs/plots/seed_starting_ensemble_costs"))
    parser.add_argument("--reference-class", default="formed_perturbed")
    parser.add_argument("--unit-counts", default="")
    parser.add_argument("--std-floor", type=float, default=1e-3)
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


def mean(values: list[float]) -> float | None:
    return sum(values) / len(values) if values else None


def sample_std(values: list[float]) -> float | None:
    if not values:
        return None
    if len(values) == 1:
        return 0.0
    avg = sum(values) / len(values)
    return math.sqrt(sum((value - avg) ** 2 for value in values) / (len(values) - 1))


def rms(values: list[float]) -> float | None:
    if not values:
        return None
    return math.sqrt(sum(value * value for value in values) / len(values))


def read_csv_rows(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    if not path.exists():
        raise FileNotFoundError(f"Input CSV not found: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        return rows, list(reader.fieldnames or [])


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_unit_filter(value: str) -> set[int] | None:
    if not value.strip():
        return None
    return {int(token.strip()) for token in value.split(",") if token.strip()}


def row_unit(row: dict[str, str]) -> int | None:
    try:
        return int(row.get("unit_count", ""))
    except ValueError:
        return None


def finite_values(rows: list[dict[str, str]], column: str) -> list[float]:
    values = [safe_float(row.get(column, "")) for row in rows]
    return [value for value in values if value is not None]


def all_feature_columns(fieldnames: list[str]) -> OrderedDict[str, list[str]]:
    available = OrderedDict()
    for group, candidates in FEATURE_GROUPS.items():
        available[group] = [column for column in candidates if column in fieldnames]
    return available


def effective_std(mean_value: float, std_value: float, std_floor: float) -> tuple[float, bool]:
    floor = max(std_floor, abs(mean_value) * std_floor)
    if std_value < floor:
        return floor, True
    return std_value, False


def build_component_rows(
    rows: list[dict[str, str]],
    fieldnames: list[str],
    reference_class: str,
    unit_filter: set[int] | None,
    std_floor: float,
) -> tuple[list[dict[str, str]], dict[str, list[str]]]:
    available = all_feature_columns(fieldnames)
    units = sorted({unit for row in rows if (unit := row_unit(row)) is not None and (unit_filter is None or unit in unit_filter)})
    start_classes = sorted({row.get("ensemble_type", "") for row in rows if row.get("ensemble_type", "") and row.get("ensemble_type") != reference_class})
    skipped: dict[str, list[str]] = {group: [] for group in FEATURE_GROUPS}
    component_rows: list[dict[str, str]] = []

    for group, candidates in FEATURE_GROUPS.items():
        missing = [feature for feature in candidates if feature not in fieldnames]
        if missing:
            skipped[group].extend(f"{feature}: missing column" for feature in missing)

    for unit in units:
        reference_rows = [row for row in rows if row_unit(row) == unit and row.get("ensemble_type") == reference_class]
        for start_class in start_classes:
            loose_rows = [row for row in rows if row_unit(row) == unit and row.get("ensemble_type") == start_class]
            if not reference_rows or not loose_rows:
                continue
            for group, features in available.items():
                for feature in features:
                    formed_values = finite_values(reference_rows, feature)
                    loose_values = finite_values(loose_rows, feature)
                    if not formed_values:
                        skipped[group].append(f"{feature}: no formed reference values for unit {unit}")
                        continue
                    if not loose_values:
                        skipped[group].append(f"{feature}: no {start_class} values for unit {unit}")
                        continue
                    formed_mean = mean(formed_values)
                    formed_std = sample_std(formed_values)
                    loose_mean = mean(loose_values)
                    if formed_mean is None or formed_std is None or loose_mean is None:
                        continue
                    std_eff, floor_used = effective_std(formed_mean, formed_std, std_floor)
                    signed_z = (loose_mean - formed_mean) / std_eff
                    component_rows.append(
                        {
                            "unit_count": str(unit),
                            "start_class": start_class,
                            "feature_group": group,
                            "feature_name": feature,
                            "loose_sample_count": str(len(loose_values)),
                            "formed_reference_sample_count": str(len(formed_values)),
                            "loose_mean": format_float(loose_mean),
                            "formed_reference_mean": format_float(formed_mean),
                            "formed_reference_std": format_float(formed_std),
                            "effective_reference_std": format_float(std_eff),
                            "std_floor_used": "true" if floor_used else "false",
                            "signed_z_deviation": format_float(signed_z),
                            "abs_z_deviation": format_float(abs(signed_z)),
                        }
                    )
    for group in skipped:
        skipped[group] = sorted(set(skipped[group]))
    return component_rows, skipped


def build_summary_rows(component_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    keys = sorted({(int(row["unit_count"]), row["start_class"]) for row in component_rows})
    summary_rows: list[dict[str, str]] = []
    for unit, start_class in keys:
        selected = [row for row in component_rows if int(row["unit_count"]) == unit and row["start_class"] == start_class]
        costs: dict[str, float | None] = {}
        counts: dict[str, int] = {}
        for group in FEATURE_GROUPS:
            values = [safe_float(row["abs_z_deviation"]) for row in selected if row["feature_group"] == group]
            finite = [value for value in values if value is not None]
            costs[group] = rms(finite)
            counts[group] = len(finite)
        all_values = [safe_float(row["abs_z_deviation"]) for row in selected]
        all_finite = [value for value in all_values if value is not None]
        available_costs = {group: value for group, value in costs.items() if value is not None}
        dominant = max(available_costs, key=lambda group: available_costs[group]) if available_costs else ""
        summary_rows.append(
            {
                "unit_count": str(unit),
                "start_class": start_class,
                "n_features_total": str(len(all_finite)),
                "n_geometric_features": str(counts["geometric"]),
                "n_contact_features": str(counts["contact_recovery"]),
                "n_register_phase_features": str(counts["register_phase"]),
                "geometric_cost": format_float(costs["geometric"]),
                "contact_recovery_cost": format_float(costs["contact_recovery"]),
                "register_phase_cost": format_float(costs["register_phase"]),
                "overall_cost": format_float(rms(all_finite)),
                "dominant_cost_group": dominant,
            }
        )
    return summary_rows


def sort_summary_for_ranking(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return sorted(rows, key=lambda row: (safe_float(row.get("overall_cost")) is None, safe_float(row.get("overall_cost")) or float("inf")))


def svg_cost_plot(summary_rows: list[dict[str, str]], path: Path) -> None:
    if not summary_rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = sorted(summary_rows, key=lambda row: (int(row["unit_count"]), row["start_class"]))
    values = [safe_float(row["overall_cost"]) or 0.0 for row in rows]
    max_value = max(values + [1.0])
    width = max(760, 115 * len(rows) + 110)
    height = 430
    left = 70
    bottom = 330
    plot_h = 250
    bar_w = 44
    gap = 34
    colors = {
        "loose_initial": "#4c78a8",
        "angular_randomized_loose_initial": "#f58518",
        "radially_separated": "#54a24b",
        "axially_misregistered": "#b279a2",
    }
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width / 2:.0f}" y="28" text-anchor="middle" font-family="Arial" font-size="17">Seed starting-ensemble cost audit</text>',
        f'<line x1="{left}" y1="{bottom - plot_h}" x2="{left}" y2="{bottom}" stroke="#333"/>',
        f'<line x1="{left}" y1="{bottom}" x2="{width - 24}" y2="{bottom}" stroke="#333"/>',
        f'<text x="18" y="{bottom - plot_h / 2:.0f}" text-anchor="middle" font-family="Arial" font-size="12" transform="rotate(-90 18,{bottom - plot_h / 2:.0f})">overall cost</text>',
    ]
    for index, row in enumerate(rows):
        value = safe_float(row["overall_cost"]) or 0.0
        h = 0.0 if max_value <= 0 else value / max_value * (plot_h - 18)
        x = left + 26 + index * (bar_w + gap)
        y = bottom - h
        label = f"{row['unit_count']} {row['start_class'].replace('_', ' ')}"
        color = colors.get(row["start_class"], "#54a24b")
        parts.append(f'<rect x="{x}" y="{y:.2f}" width="{bar_w}" height="{h:.2f}" fill="{color}"/>')
        parts.append(f'<text x="{x + bar_w / 2:.0f}" y="{y - 6:.2f}" text-anchor="middle" font-family="Arial" font-size="10">{value:.2f}</text>')
        parts.append(f'<text x="{x + bar_w / 2:.0f}" y="{bottom + 14}" text-anchor="end" font-family="Arial" font-size="9" transform="rotate(-40 {x + bar_w / 2:.0f},{bottom + 14})">{label}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def markdown_table(rows: list[dict[str, str]], columns: list[str]) -> list[str]:
    lines = ["| " + " | ".join(columns) + " |", "|" + "|".join("---" for _ in columns) + "|"]
    for row in rows:
        lines.append("| " + " | ".join(row.get(column, "") for column in columns) + " |")
    return lines


def interpret_radial_class(summary_rows: list[dict[str, str]]) -> list[str]:
    lines: list[str] = []
    units = sorted({row["unit_count"] for row in summary_rows})
    by_key = {(row["unit_count"], row["start_class"]): row for row in summary_rows}
    for unit in units:
        radial = by_key.get((unit, "radially_separated"))
        loose = by_key.get((unit, "loose_initial"))
        if radial is None:
            continue
        radial_cost = safe_float(radial.get("overall_cost"))
        loose_cost = safe_float(loose.get("overall_cost")) if loose else None
        if radial_cost is not None and loose_cost is not None:
            direction = "higher" if radial_cost > loose_cost else "lower"
            comparison = f"{direction} than `loose_initial` ({format_float(loose_cost, 3)})"
        else:
            comparison = "not directly comparable to `loose_initial` in this run"
        lines.append(
            f"- Unit {unit} `radially_separated`: overall cost {radial.get('overall_cost', '')}, "
            f"{comparison}; dominant group `{radial.get('dominant_cost_group', '')}`."
        )
    return lines


def compare_start_class(
    summary_rows: list[dict[str, str]],
    target_class: str,
    reference_classes: list[str],
    unit_filter: set[str] | None = None,
) -> list[str]:
    lines: list[str] = []
    units = sorted({row["unit_count"] for row in summary_rows if unit_filter is None or row["unit_count"] in unit_filter})
    by_key = {(row["unit_count"], row["start_class"]): row for row in summary_rows}
    for unit in units:
        target = by_key.get((unit, target_class))
        if target is None:
            continue
        target_cost = safe_float(target.get("overall_cost"))
        comparisons: list[str] = []
        for reference_class in reference_classes:
            reference = by_key.get((unit, reference_class))
            reference_cost = safe_float(reference.get("overall_cost")) if reference else None
            if target_cost is None or reference_cost is None:
                comparisons.append(f"not directly comparable to `{reference_class}`")
                continue
            direction = "higher" if target_cost > reference_cost else "lower"
            comparisons.append(f"{direction} than `{reference_class}` ({format_float(reference_cost, 3)})")
        lines.append(
            f"- Unit {unit} `{target_class}`: overall cost {target.get('overall_cost', '')}, "
            + "; ".join(comparisons)
            + f"; dominant group `{target.get('dominant_cost_group', '')}`."
        )
    return lines


def write_report(
    path: Path,
    args: argparse.Namespace,
    component_rows: list[dict[str, str]],
    summary_rows: list[dict[str, str]],
    skipped: dict[str, list[str]],
    plot_paths: list[Path],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    used_by_group: dict[str, list[str]] = {}
    for group in FEATURE_GROUPS:
        used_by_group[group] = sorted({row["feature_name"] for row in component_rows if row["feature_group"] == group})
    focus = [row for row in summary_rows if row["unit_count"] in {"6", "7"}]
    ranking = sort_summary_for_ranking(summary_rows)
    lines = [
        "# Seed Starting-Ensemble Cost Audit",
        "",
        "## Purpose",
        "",
        "This audit compares synthetic loose-start endpoint classes against formed-perturbed endpoints in existing order-parameter space. It is intended to identify which cost groups dominate the distance from each loose-start class to the formed endpoint distribution.",
        "",
        "This analysis is hypothesis-generating. It is not an atomistic Schrodinger bridge, not molecular dynamics, and not evidence of a physical nucleation pathway.",
        "",
        "The `radially_separated` class is a dry-down/concentration-compaction proxy: it keeps rough angular phase and small axial offsets while moving chains farther apart radially. It does not model solvent evaporation, crystallization, molecular dynamics, or a full atomistic Schrodinger bridge.",
        "",
        "The `axially_misregistered` class is a synthetic register-perturbation proxy. It starts from the loose baseline and applies symmetric chain-specific offsets along the fitted stack axis, asking whether axial register disruption gives a distinct cost signature from radial separation.",
        "",
        "## Inputs Used",
        "",
        f"- Order-parameter CSV: `{args.input_csv}`",
        f"- Reference class: `{args.reference_class}`",
        "",
        "## Cost Definition",
        "",
        "For each unit count, the formed-perturbed rows define a reference distribution. For each available feature, the audit computes the formed-reference mean and standard deviation, compares the loose-start class mean against that formed mean, and records signed and absolute standardized deviations.",
        "",
        "Group costs are root-mean-square absolute standardized deviations across available features. Lower cost means closer to the formed-perturbed endpoint distribution in the selected exploratory order parameters.",
        "",
        "A small standard-deviation floor is used when the formed-reference standard deviation is zero or tiny, so constant synthetic reference features do not crash the audit.",
        "",
        "## Feature Groups Used",
        "",
    ]
    for group, features in used_by_group.items():
        lines.append(f"- `{group}`: " + (", ".join(f"`{feature}`" for feature in features) if features else "no available features"))
    if any(skipped.values()):
        lines.extend(["", "Skipped or unavailable feature notes:"])
        for group, notes in skipped.items():
            if notes:
                shown = "; ".join(notes[:8])
                suffix = " ..." if len(notes) > 8 else ""
                lines.append(f"- `{group}`: {shown}{suffix}")
    lines.extend(["", "## 6-vs-7 Focused Summary", ""])
    lines.extend(markdown_table(focus, ["unit_count", "start_class", "geometric_cost", "contact_recovery_cost", "register_phase_cost", "overall_cost", "dominant_cost_group"]))
    lines.extend(["", "## Ranking By Overall Cost", ""])
    lines.extend(markdown_table(ranking, ["unit_count", "start_class", "overall_cost", "dominant_cost_group"]))
    lines.extend(["", "## Cost Decomposition", ""])
    for row in ranking:
        lines.append(
            f"- Unit {row['unit_count']} `{row['start_class']}`: overall cost {row['overall_cost']}; "
            f"dominant group `{row['dominant_cost_group']}`."
        )
    radial_lines = interpret_radial_class(summary_rows)
    if radial_lines:
        lines.extend(["", "## Radial-Separation Readout", ""])
        lines.extend(radial_lines)
        lines.append("")
        lines.append("This readout asks whether a radial compaction proxy changes the cost relative to the older loose baseline. It should be interpreted as a geometric endpoint comparison, not as a simulation of the lab dry-down process.")
    axial_lines = compare_start_class(
        summary_rows,
        "axially_misregistered",
        ["loose_initial", "radially_separated", "angular_randomized_loose_initial"],
        {"6", "7"},
    )
    if axial_lines:
        lines.extend(["", "## Axial-Misregistration Readout", ""])
        lines.extend(axial_lines)
        lines.append("")
        lines.append("This readout asks whether stack-axis offsets create a distinct register or phase cost signature. It should be interpreted as a synthetic endpoint comparison, not as a physical trajectory.")
    lines.extend(
        [
            "",
            "## Conservative Interpretation",
            "",
            "The cost values are order-parameter distance proxies. They indicate how far synthetic loose-start classes are from formed-perturbed endpoints in selected summary coordinates. They do not establish a physical mechanism or a validated assembly route.",
            "",
            "A lower cost class is closer to the formed-perturbed endpoint distribution in this feature set. A dominant cost group suggests which type of coordinate difference is largest: geometry, contact recovery, or register/phase.",
            "",
            "For the current 6-vs-7 focus, this audit should be read alongside the known endpoint construction: mini-hexaplex endpoints are coordinate-derived fragments cut from already-formed models, so formed geometry is present by construction.",
            "",
            "## Recommended Next Step",
            "",
            "Keep the start-class set limited to `loose_initial`, `angular_randomized_loose_initial`, `radially_separated`, and `axially_misregistered` until the axial-register readout is reviewed. Additional start classes should be added only one at a time with the same cost-decomposition checks.",
        ]
    )
    if plot_paths:
        lines.extend(["", "## Plots", ""])
        for plot_path in plot_paths:
            lines.append(f"- `{plot_path}`")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(args: argparse.Namespace) -> dict[str, object]:
    if args.std_floor <= 0:
        raise ValueError("--std-floor must be greater than zero")
    rows, fieldnames = read_csv_rows(args.input_csv)
    unit_filter = parse_unit_filter(args.unit_counts)
    component_rows, skipped = build_component_rows(rows, fieldnames, args.reference_class, unit_filter, args.std_floor)
    summary_rows = build_summary_rows(component_rows)
    write_csv(args.components_csv, component_rows, COMPONENT_COLUMNS)
    write_csv(args.summary_csv, summary_rows, SUMMARY_COLUMNS)
    plot_path = args.plot_dir / "seed_starting_ensemble_overall_costs.svg"
    svg_cost_plot(summary_rows, plot_path)
    plot_paths = [plot_path] if plot_path.exists() else []
    write_report(args.report, args, component_rows, summary_rows, skipped, plot_paths)
    return {
        "component_rows": len(component_rows),
        "summary_rows": len(summary_rows),
        "plot_paths": plot_paths,
    }


def main() -> None:
    result = run(parse_args())
    print(f"Wrote {result['component_rows']} cost component rows")
    print(f"Wrote {result['summary_rows']} cost summary rows")
    print(f"Wrote {len(result['plot_paths'])} plots")


if __name__ == "__main__":
    main()
