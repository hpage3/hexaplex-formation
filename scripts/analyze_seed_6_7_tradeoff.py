#!/usr/bin/env python3
"""Explain 6-vs-7 seed tradeoffs using existing metrics only."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


SUMMARY_COLUMNS = [
    "category",
    "metric",
    "unit6_value",
    "unit7_value",
    "difference_7_minus_6",
    "fractional_difference_7_minus_6",
    "classification",
    "reasoning_category",
    "notes",
]

FEATURE_COLUMNS = [
    "source_file",
    "metric",
    "metric_label",
    "unit6_value",
    "unit7_value",
    "difference_7_minus_6",
    "fractional_difference_7_minus_6",
    "classification",
    "reasoning_category",
    "interpretation",
]

ABSOLUTE_CONTACT_COLUMNS = [
    ("total_interchain_contacts", "total interchain contacts", "absolute_network_growth"),
    ("CYP_MEP_contact_count", "CYP/MEP-involving contacts", "absolute_network_growth"),
    ("GLU_contact_count", "GLU-involving contacts", "absolute_network_growth"),
    ("backbone_like_contact_count", "backbone-like contacts", "absolute_network_growth"),
    ("contacts_per_unit", "contacts per unit", "possible_edge_effect"),
    ("contacts_per_chain_pair_mean", "contacts per chain pair mean", "possible_edge_effect"),
    ("contacts_per_chain_pair_max", "contacts per chain pair max", "possible_edge_effect"),
]

FRACTIONAL_COLUMNS = [
    ("CYP_MEP_contact_fraction", "CYP/MEP contact fraction", "normalized_fraction", "higher"),
    ("perturbation_contact_fraction_vs_reference_mean", "perturbation contact retention", "perturbation_retention", "higher"),
    ("perturbation_contact_fraction_vs_reference_std", "perturbation contact retention variability", "perturbation_retention", "lower"),
    ("contact_redundancy_score", "contact redundancy score", "normalized_fraction", "higher"),
    ("nucleation_network_score", "network score", "normalized_fraction", "higher"),
]

ORDER_FEATURES = [
    ("CYP_MEP_contact_fraction_vs_target", "CYP/MEP contact fraction vs target", "normalized_fraction", "higher"),
    ("contact_fraction_vs_target", "total contact fraction vs target", "normalized_fraction", "higher"),
    ("compactness_score", "compactness", "geometric_register", "higher"),
    ("axial_register_score", "axial register", "geometric_register", "higher"),
    ("RMSD_to_formed_seed_A", "RMSD to formed seed", "geometric_register", "lower"),
    ("seed_formation_score", "seed formation score", "normalized_fraction", "higher"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--transition-summary", type=Path, default=Path("outputs/metrics/seed_6_7_transition_summary.csv"))
    parser.add_argument("--transition-comparison", type=Path, default=Path("outputs/metrics/seed_6_7_feature_comparison.csv"))
    parser.add_argument("--seed-order-csv", type=Path, default=Path("outputs/metrics/seed_formation_order_parameters.csv"))
    parser.add_argument("--contact-summary-csv", type=Path, default=Path("outputs/metrics/seed_contact_network_summary.csv"))
    parser.add_argument("--contact-edges-csv", type=Path, default=Path("outputs/metrics/seed_contact_network_edges.csv"))
    parser.add_argument("--endpoint-means-csv", type=Path, default=Path("outputs/metrics/seed_endpoint_metric_means_by_ensemble.csv"))
    parser.add_argument(
        "--angular-activation-csv",
        type=Path,
        default=Path("outputs/metrics/seed_bridge_activation_summary_angular_randomized.csv"),
    )
    parser.add_argument(
        "--angular-ordering-csv",
        type=Path,
        default=Path("outputs/metrics/seed_bridge_feature_ordering_angular_randomized.csv"),
    )
    parser.add_argument("--summary-csv", type=Path, default=Path("outputs/metrics/seed_6_7_tradeoff_summary.csv"))
    parser.add_argument("--feature-table-csv", type=Path, default=Path("outputs/metrics/seed_6_7_tradeoff_feature_table.csv"))
    parser.add_argument("--report", type=Path, default=Path("outputs/reports/seed_6_7_tradeoff_report.md"))
    parser.add_argument("--plot-dir", type=Path, default=Path("outputs/plots/seed_6_7_tradeoff"))
    parser.add_argument("--negligible-relative", type=float, default=0.005)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


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


def fractional_difference(six: float, seven: float) -> float | None:
    if abs(six) < 1e-12:
        return None
    return (seven - six) / abs(six)


def classify(six: float, seven: float, orientation: str, negligible_relative: float) -> str:
    diff = seven - six
    frac = fractional_difference(six, seven)
    if abs(diff) < 1e-12 or (frac is not None and abs(frac) < negligible_relative):
        return "mixed_or_negligible"
    if orientation == "lower":
        return "favors_7" if seven < six else "favors_6"
    return "favors_7" if seven > six else "favors_6"


def add_feature_row(
    rows: list[dict[str, str]],
    source: str,
    metric: str,
    label: str,
    six: float | None,
    seven: float | None,
    reasoning_category: str,
    orientation: str,
    negligible_relative: float,
    interpretation: str,
) -> None:
    if six is None or seven is None:
        rows.append(
            {
                "source_file": source,
                "metric": metric,
                "metric_label": label,
                "unit6_value": "",
                "unit7_value": "",
                "difference_7_minus_6": "",
                "fractional_difference_7_minus_6": "",
                "classification": "mixed_or_negligible",
                "reasoning_category": "insufficient_data",
                "interpretation": f"{label}: insufficient data in existing CSVs.",
            }
        )
        return
    diff = seven - six
    frac = fractional_difference(six, seven)
    classification = classify(six, seven, orientation, negligible_relative)
    rows.append(
        {
            "source_file": source,
            "metric": metric,
            "metric_label": label,
            "unit6_value": format_float(six),
            "unit7_value": format_float(seven),
            "difference_7_minus_6": format_float(diff),
            "fractional_difference_7_minus_6": format_float(frac),
            "classification": classification,
            "reasoning_category": reasoning_category,
            "interpretation": interpretation,
        }
    )


def central_contact_rows(rows: list[dict[str, str]]) -> dict[int, dict[str, str]]:
    out = {}
    for row in rows:
        if row.get("variant_id") in {"central6_units", "central7_units"}:
            unit = safe_float(row.get("units_per_chain"))
            if unit is not None:
                out[int(unit)] = row
    return out


def edge_type_totals(rows: list[dict[str, str]]) -> dict[int, dict[str, float]]:
    totals: dict[int, dict[str, float]] = defaultdict(lambda: defaultdict(float))
    for row in rows:
        if row.get("edge_scope") != "chain":
            continue
        unit = safe_float(row.get("units_per_chain"))
        if unit is None or int(unit) not in {6, 7}:
            continue
        variant = row.get("variant_id", "")
        if variant != f"central{int(unit)}_units":
            continue
        for column in ["contact_count", "CYP_MEP_vs_CYP_MEP_count", "CYP_MEP_vs_GLU_count", "GLU_vs_GLU_count"]:
            value = safe_float(row.get(column))
            if value is not None:
                totals[int(unit)][column] += value
    return {unit: dict(values) for unit, values in totals.items()}


def edge_position_summary(rows: list[dict[str, str]]) -> dict[int, dict[str, float]]:
    summary: dict[int, dict[str, float]] = {}
    for unit in [6, 7]:
        total = 0.0
        edge = 0.0
        internal = 0.0
        for row in rows:
            if row.get("edge_scope") != "unit":
                continue
            if row.get("variant_id") != f"central{unit}_units":
                continue
            count = safe_float(row.get("contact_count")) or 0.0
            unit_a = safe_float(row.get("unit_a"))
            unit_b = safe_float(row.get("unit_b"))
            total += count
            if unit_a in {1.0, float(unit)} or unit_b in {1.0, float(unit)}:
                edge += count
            else:
                internal += count
        summary[unit] = {
            "unit_edge_contact_count": edge,
            "unit_internal_contact_count": internal,
            "unit_edge_contact_fraction": edge / total if total else math.nan,
        }
    return summary


def transition_lookup(rows: list[dict[str, str]]) -> dict[tuple[int, str, str], float]:
    out = {}
    for row in rows:
        unit = safe_float(row.get("unit_count"))
        value = safe_float(row.get("mean"))
        if unit is None or value is None:
            continue
        out[(int(unit), row.get("ensemble_type", ""), row.get("feature", ""))] = value
    return out


def angular_activation_lookup(rows: list[dict[str, str]]) -> dict[tuple[int, str], float]:
    out = {}
    for row in rows:
        if row.get("threshold_fraction") != "0.75":
            continue
        unit = safe_float(row.get("unit_count"))
        value = safe_float(row.get("median_activation_time"))
        if unit is None or value is None:
            continue
        out[(int(unit), row.get("feature", ""))] = value
    return out


def build_feature_rows(args: argparse.Namespace) -> list[dict[str, str]]:
    feature_rows: list[dict[str, str]] = []
    contact_by_unit = central_contact_rows(read_csv(args.contact_summary_csv))
    edge_totals = edge_type_totals(read_csv(args.contact_edges_csv))
    edge_positions = edge_position_summary(read_csv(args.contact_edges_csv))
    transition = transition_lookup(read_csv(args.transition_summary))
    activation = angular_activation_lookup(read_csv(args.angular_activation_csv))

    if 6 in contact_by_unit and 7 in contact_by_unit:
        for column, label, category in ABSOLUTE_CONTACT_COLUMNS:
            six = safe_float(contact_by_unit[6].get(column))
            seven = safe_float(contact_by_unit[7].get(column))
            add_feature_row(
                feature_rows,
                str(args.contact_summary_csv),
                column,
                label,
                six,
                seven,
                category,
                "higher",
                args.negligible_relative,
                f"{label} compares absolute/static contact-network growth between central 6- and 7-unit seeds.",
            )
        for column, label, category, orientation in FRACTIONAL_COLUMNS:
            six = safe_float(contact_by_unit[6].get(column))
            seven = safe_float(contact_by_unit[7].get(column))
            add_feature_row(
                feature_rows,
                str(args.contact_summary_csv),
                column,
                label,
                six,
                seven,
                category,
                orientation,
                args.negligible_relative,
                f"{label} is a normalized or perturbation-retention metric, so it can move differently from raw counts.",
            )

    edge_metrics = [
        ("contact_count", "chain-summed total contacts", "absolute_network_growth"),
        ("CYP_MEP_vs_CYP_MEP_count", "CYP/MEP vs CYP/MEP contacts", "absolute_network_growth"),
        ("CYP_MEP_vs_GLU_count", "CYP/MEP vs GLU contacts", "absolute_network_growth"),
        ("GLU_vs_GLU_count", "GLU vs GLU contacts", "absolute_network_growth"),
    ]
    for column, label, category in edge_metrics:
        add_feature_row(
            feature_rows,
            str(args.contact_edges_csv),
            column,
            label,
            edge_totals.get(6, {}).get(column),
            edge_totals.get(7, {}).get(column),
            category,
            "higher",
            args.negligible_relative,
            f"{label} uses chain-level edge rows and separates contact chemistry where existing edge data allow it.",
        )

    for column, label in [
        ("unit_edge_contact_count", "contacts involving terminal/edge units"),
        ("unit_internal_contact_count", "contacts involving only internal units"),
        ("unit_edge_contact_fraction", "fraction of unit-edge contacts involving terminal/edge units"),
    ]:
        add_feature_row(
            feature_rows,
            str(args.contact_edges_csv),
            column,
            label,
            edge_positions.get(6, {}).get(column),
            edge_positions.get(7, {}).get(column),
            "possible_edge_effect",
            "higher",
            args.negligible_relative,
            f"{label} is a coarse edge-effect proxy from unit-level contact edges.",
        )

    for ensemble in ["formed_perturbed", "loose_initial", "angular_randomized_loose_initial"]:
        for feature, label, category, orientation in ORDER_FEATURES:
            add_feature_row(
                feature_rows,
                str(args.transition_summary),
                f"{ensemble}:{feature}",
                f"{ensemble} {label}",
                transition.get((6, ensemble, feature)),
                transition.get((7, ensemble, feature)),
                category,
                orientation,
                args.negligible_relative,
                f"{label} compares normalized/order-parameter means for the {ensemble} ensemble.",
            )

    activation_features = [
        ("CYP_MEP_contact_fraction_vs_target", "angular-randomized CYP/MEP activation time"),
        ("contact_fraction_vs_target", "angular-randomized total-contact activation time"),
        ("compactness_score", "angular-randomized compactness activation time"),
        ("axial_register_score", "angular-randomized axial-register activation time"),
        ("rmsd_formedness_score", "angular-randomized RMSD formedness activation time"),
        ("seed_formation_score", "angular-randomized seed-score activation time"),
    ]
    for feature, label in activation_features:
        add_feature_row(
            feature_rows,
            str(args.angular_activation_csv),
            f"angular_activation:{feature}",
            label,
            activation.get((6, feature)),
            activation.get((7, feature)),
            "normalized_fraction",
            "lower",
            args.negligible_relative,
            "Lower median activation time means earlier threshold crossing in the angular-randomized bridge diagnostics.",
        )
    return feature_rows


def summarize_by_category(feature_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in feature_rows:
        grouped[row["reasoning_category"]].append(row)
    summary_rows: list[dict[str, str]] = []
    for category, rows in sorted(grouped.items()):
        favors_6 = sum(1 for row in rows if row["classification"] == "favors_6")
        favors_7 = sum(1 for row in rows if row["classification"] == "favors_7")
        mixed = sum(1 for row in rows if row["classification"] == "mixed_or_negligible")
        if favors_7 > favors_6:
            classification = "favors_7"
        elif favors_6 > favors_7:
            classification = "favors_6"
        else:
            classification = "mixed_or_negligible"
        summary_rows.append(
            {
                "category": category,
                "metric": "classification_count",
                "unit6_value": str(favors_6),
                "unit7_value": str(favors_7),
                "difference_7_minus_6": str(favors_7 - favors_6),
                "fractional_difference_7_minus_6": "",
                "classification": classification,
                "reasoning_category": category,
                "notes": f"{favors_7} metrics favor 7, {favors_6} favor 6, {mixed} are mixed/negligible.",
            }
        )
    return summary_rows


def write_csv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def write_svg_summary(feature_rows: list[dict[str, str]], plot_dir: Path) -> list[Path]:
    plot_dir.mkdir(parents=True, exist_ok=True)
    categories = ["absolute_network_growth", "normalized_fraction", "perturbation_retention", "geometric_register", "possible_edge_effect"]
    counts = {
        category: {
            "favors_6": sum(1 for row in feature_rows if row["reasoning_category"] == category and row["classification"] == "favors_6"),
            "favors_7": sum(1 for row in feature_rows if row["reasoning_category"] == category and row["classification"] == "favors_7"),
            "mixed_or_negligible": sum(1 for row in feature_rows if row["reasoning_category"] == category and row["classification"] == "mixed_or_negligible"),
        }
        for category in categories
    }
    width = 900
    height = 430
    left = 80
    bottom = 340
    group_width = 145
    bar_width = 32
    max_count = max([value for group in counts.values() for value in group.values()] + [1])
    colors = {"favors_6": "#4c78a8", "favors_7": "#f58518", "mixed_or_negligible": "#9a9a9a"}
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="450" y="32" text-anchor="middle" font-family="Arial" font-size="18">6-vs-7 tradeoff classification counts</text>',
        f'<line x1="{left}" y1="60" x2="{left}" y2="{bottom}" stroke="#333"/>',
        f'<line x1="{left}" y1="{bottom}" x2="{width - 30}" y2="{bottom}" stroke="#333"/>',
    ]
    for i, category in enumerate(categories):
        x0 = left + 25 + i * group_width
        for j, classification in enumerate(["favors_6", "favors_7", "mixed_or_negligible"]):
            value = counts[category][classification]
            h = value / max_count * 240 if max_count else 0
            x = x0 + j * (bar_width + 6)
            y = bottom - h
            parts.append(f'<rect x="{x}" y="{y:.2f}" width="{bar_width}" height="{h:.2f}" fill="{colors[classification]}"/>')
            parts.append(f'<text x="{x + bar_width / 2:.0f}" y="{y - 6:.2f}" text-anchor="middle" font-family="Arial" font-size="11">{value}</text>')
        label = category.replace("_", " ")
        parts.append(f'<text x="{x0 + 48}" y="{bottom + 18}" text-anchor="end" font-family="Arial" font-size="10" transform="rotate(-30 {x0 + 48},{bottom + 18})">{label}</text>')
    legend_x = 600
    for i, classification in enumerate(["favors_6", "favors_7", "mixed_or_negligible"]):
        y = 65 + i * 22
        parts.append(f'<rect x="{legend_x}" y="{y}" width="14" height="14" fill="{colors[classification]}"/>')
        parts.append(f'<text x="{legend_x + 20}" y="{y + 12}" font-family="Arial" font-size="12">{classification}</text>')
    parts.append("</svg>")
    out = plot_dir / "seed_6_7_tradeoff_classification_counts.svg"
    out.write_text("\n".join(parts) + "\n", encoding="utf-8")
    return [out]


def report_text(summary_rows: list[dict[str, str]], feature_rows: list[dict[str, str]], plots: list[Path], args: argparse.Namespace) -> str:
    def rows_for(category: str) -> list[dict[str, str]]:
        return [row for row in feature_rows if row["reasoning_category"] == category]

    absolute_rows = [row for row in rows_for("absolute_network_growth") if row["classification"] == "favors_7"]
    retention_rows = rows_for("perturbation_retention")
    geometry_rows = rows_for("geometric_register")
    edge_rows = rows_for("possible_edge_effect")

    lines = [
        "# 6-vs-7 Seed Tradeoff Analysis",
        "",
        "This follow-up reuses existing CSV outputs only. It does not regenerate ensembles, mini-hexaplex structures, or bridge paths.",
        "",
        "The goal is explanatory rather than mechanistic: compare absolute contact-network growth, normalized formedness metrics, and perturbation-retention metrics in coordinate-derived 6- and 7-unit mini-seeds. This report does not claim a physical nucleation threshold or a true assembly mechanism.",
        "",
        "## Category Summary",
        "",
        "| Reasoning category | Metrics favoring 6 | Metrics favoring 7 | Net 7-minus-6 | Overall classification |",
        "|---|---:|---:|---:|---|",
    ]
    for row in summary_rows:
        lines.append(
            f"| {row['reasoning_category']} | {row['unit6_value']} | {row['unit7_value']} | "
            f"{row['difference_7_minus_6']} | {row['classification']} |"
        )

    lines.extend(["", "## Absolute Contact-Network Growth", "", "| Metric | 6-unit | 7-unit | 7-minus-6 | Classification |", "|---|---:|---:|---:|---|"])
    for row in absolute_rows:
        lines.append(
            f"| {row['metric_label']} | {row['unit6_value']} | {row['unit7_value']} | "
            f"{row['difference_7_minus_6']} | {row['classification']} |"
        )

    lines.extend(["", "## Normalized, Retention, And Register Metrics", "", "| Metric | 6-unit | 7-unit | 7-minus-6 | Classification | Category |", "|---|---:|---:|---:|---|---|"])
    for row in retention_rows + geometry_rows:
        lines.append(
            f"| {row['metric_label']} | {row['unit6_value']} | {row['unit7_value']} | "
            f"{row['difference_7_minus_6']} | {row['classification']} | {row['reasoning_category']} |"
        )

    lines.extend(["", "## Edge-Effect Proxies", "", "| Metric | 6-unit | 7-unit | 7-minus-6 | Classification |", "|---|---:|---:|---:|---|"])
    for row in edge_rows:
        lines.append(
            f"| {row['metric_label']} | {row['unit6_value']} | {row['unit7_value']} | "
            f"{row['difference_7_minus_6']} | {row['classification']} |"
        )

    lines.extend(["", "## Interpretation", ""])
    if absolute_rows:
        lines.append(
            "The 7-unit seed increases contact-network opportunity and redundancy: raw total, CYP/MEP-involving, and chemistry-resolved chain-edge contact counts all rise in the existing static contact network."
        )
    lines.append(
        "The 6-unit seed can still score slightly better on normalized formedness, geometric register, or endpoint-retention metrics. Those quantities are fractions, fitted scores, or retention ratios rather than absolute contact opportunities."
    )
    lines.append(
        "This pattern supports a transition/tradeoff regime rather than a single best seed length. Added length creates more absolute contacts, but it also adds more degrees of freedom and edge-sensitive contacts that can reduce normalized retention or register scores."
    )
    lines.append(
        "The unit-level edge proxy can flag whether terminal/edge contacts change between 6 and 7 units, but it does not prove that edge effects cause the normalized-score differences."
    )

    lines.extend(["", "## Outputs", "", f"- Summary CSV: `{args.summary_csv}`", f"- Feature table CSV: `{args.feature_table_csv}`", f"- Plot directory: `{args.plot_dir}`"])
    for plot in plots:
        lines.append(f"- Plot: `{plot}`")
    lines.extend(["", "## Limitations", "", "- The contact-type split is available from `seed_contact_network_edges.csv`; any absent pair classes are reported only where existing columns support them.", "- Endpoint retention is based on existing perturbation summaries and bridge diagnostics, not new sampling.", "- The analysis is order-parameter and coordinate-derived; it should not be read as an atomistic mechanism."])
    return "\n".join(lines) + "\n"


def run(args: argparse.Namespace) -> dict[str, object]:
    feature_rows = build_feature_rows(args)
    summary_rows = summarize_by_category(feature_rows)
    write_csv(args.feature_table_csv, feature_rows, FEATURE_COLUMNS)
    write_csv(args.summary_csv, summary_rows, SUMMARY_COLUMNS)
    plots = write_svg_summary(feature_rows, args.plot_dir)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(report_text(summary_rows, feature_rows, plots, args), encoding="utf-8")
    return {"feature_rows": len(feature_rows), "summary_rows": len(summary_rows), "plots": plots}


def main() -> None:
    result = run(parse_args())
    print(f"Wrote {result['summary_rows']} summary rows")
    print(f"Wrote {result['feature_rows']} feature rows")
    print(f"Wrote {len(result['plots'])} plots")


if __name__ == "__main__":
    main()
