#!/usr/bin/env python3
"""Focused 6-vs-7 unit seed-transition analysis from existing outputs."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


PRIMARY_FEATURES = [
    ("CYP_MEP_contact_fraction_vs_target", "CYP/MEP contact formation", "higher"),
    ("contact_fraction_vs_target", "total contact formation", "higher"),
    ("compactness_score", "compactness", "higher"),
    ("axial_register_score", "axial register", "higher"),
    ("RMSD_to_formed_seed_A", "RMSD formedness", "lower"),
    ("seed_formation_score", "seed formedness score", "higher"),
]

ENSEMBLES = ["loose_initial", "angular_randomized_loose_initial", "formed_perturbed"]

SUMMARY_COLUMNS = [
    "unit_count",
    "ensemble_type",
    "feature",
    "feature_label",
    "sample_count",
    "mean",
    "std",
    "orientation",
]

COMPARISON_COLUMNS = [
    "metric",
    "feature",
    "ensemble_type",
    "unit6_value",
    "unit7_value",
    "difference_7_minus_6",
    "standardized_effect_size",
    "direction_of_change",
    "cautious_interpretation",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--seed-order-csv",
        type=Path,
        default=Path("outputs/metrics/seed_formation_order_parameters.csv"),
    )
    parser.add_argument(
        "--contact-summary-csv",
        type=Path,
        default=Path("outputs/metrics/seed_contact_network_summary.csv"),
    )
    parser.add_argument(
        "--bridge-activation-csv",
        type=Path,
        default=Path("outputs/metrics/seed_bridge_activation_summary.csv"),
    )
    parser.add_argument(
        "--bridge-ordering-csv",
        type=Path,
        default=Path("outputs/metrics/seed_bridge_feature_ordering.csv"),
    )
    parser.add_argument(
        "--angular-bridge-activation-csv",
        type=Path,
        default=Path("outputs/metrics/seed_bridge_activation_summary_angular_randomized.csv"),
    )
    parser.add_argument(
        "--angular-bridge-ordering-csv",
        type=Path,
        default=Path("outputs/metrics/seed_bridge_feature_ordering_angular_randomized.csv"),
    )
    parser.add_argument(
        "--mini-geometry-csv",
        type=Path,
        default=Path("outputs/metrics/mini_hexaplex_geometry_summary.csv"),
    )
    parser.add_argument(
        "--mini-helicity-csv",
        type=Path,
        default=Path("outputs/metrics/mini_hexaplex_helicity_summary.csv"),
    )
    parser.add_argument(
        "--summary-csv",
        type=Path,
        default=Path("outputs/metrics/seed_6_7_transition_summary.csv"),
    )
    parser.add_argument(
        "--comparison-csv",
        type=Path,
        default=Path("outputs/metrics/seed_6_7_feature_comparison.csv"),
    )
    parser.add_argument(
        "--plot-dir",
        type=Path,
        default=Path("outputs/plots/seed_6_7_transition"),
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=Path("outputs/reports/seed_6_7_transition_report.md"),
    )
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def safe_float(value: object) -> float | None:
    try:
        result = float(str(value))
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def format_float(value: float | None, digits: int = 6) -> str:
    if value is None or not math.isfinite(value):
        return ""
    return f"{value:.{digits}f}"


def mean(values: list[float]) -> float | None:
    if not values:
        return None
    return sum(values) / len(values)


def sample_std(values: list[float]) -> float | None:
    if len(values) < 2:
        return 0.0 if values else None
    avg = mean(values)
    assert avg is not None
    return math.sqrt(sum((value - avg) ** 2 for value in values) / (len(values) - 1))


def pooled_std(std_a: float | None, n_a: int, std_b: float | None, n_b: int) -> float | None:
    if std_a is None or std_b is None or n_a < 2 or n_b < 2:
        return None
    numerator = ((n_a - 1) * std_a**2) + ((n_b - 1) * std_b**2)
    denominator = n_a + n_b - 2
    if denominator <= 0:
        return None
    pooled = math.sqrt(numerator / denominator)
    return pooled if pooled > 0 else None


def summarize_seed_order(rows: list[dict[str, str]]) -> tuple[list[dict[str, str]], dict[tuple[int, str, str], dict[str, float]]]:
    grouped: dict[tuple[int, str, str], list[float]] = defaultdict(list)
    for row in rows:
        try:
            unit_count = int(row.get("unit_count", ""))
        except ValueError:
            continue
        if unit_count not in {6, 7}:
            continue
        ensemble = row.get("ensemble_type", "")
        if ensemble not in ENSEMBLES:
            continue
        for feature, _, _ in PRIMARY_FEATURES:
            value = safe_float(row.get(feature))
            if value is not None:
                grouped[(unit_count, ensemble, feature)].append(value)

    summary_rows: list[dict[str, str]] = []
    stats: dict[tuple[int, str, str], dict[str, float]] = {}
    for (unit_count, ensemble, feature), values in sorted(grouped.items()):
        avg = mean(values)
        std = sample_std(values)
        label, orientation = feature_metadata(feature)
        stats[(unit_count, ensemble, feature)] = {
            "mean": avg if avg is not None else math.nan,
            "std": std if std is not None else math.nan,
            "n": float(len(values)),
        }
        summary_rows.append(
            {
                "unit_count": str(unit_count),
                "ensemble_type": ensemble,
                "feature": feature,
                "feature_label": label,
                "sample_count": str(len(values)),
                "mean": format_float(avg),
                "std": format_float(std),
                "orientation": orientation,
            }
        )
    return summary_rows, stats


def feature_metadata(feature: str) -> tuple[str, str]:
    for name, label, orientation in PRIMARY_FEATURES:
        if name == feature:
            return label, orientation
    return feature, "higher"


def direction(feature: str, difference: float | None) -> str:
    if difference is None or abs(difference) < 1e-12:
        return "similar"
    _, orientation = feature_metadata(feature)
    if orientation == "lower":
        return "7-unit stronger" if difference < 0 else "6-unit stronger"
    return "7-unit stronger" if difference > 0 else "6-unit stronger"


def interpretation_for(feature: str, ensemble: str, change: str) -> str:
    label, orientation = feature_metadata(feature)
    if change == "similar":
        return f"{label} is similar between 6 and 7 units in {ensemble}."
    if change == "7-unit stronger":
        qualifier = "lower" if orientation == "lower" else "higher"
        return f"7 units show {qualifier} {label} in {ensemble}, consistent with stronger organization in this existing dataset."
    qualifier = "lower" if orientation == "lower" else "higher"
    return f"6 units show {qualifier} {label} in {ensemble}; this metric does not support a stronger 7-unit signal by itself."


def feature_comparisons(stats: dict[tuple[int, str, str], dict[str, float]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for ensemble in ENSEMBLES:
        for feature, _, _ in PRIMARY_FEATURES:
            six = stats.get((6, ensemble, feature))
            seven = stats.get((7, ensemble, feature))
            if not six or not seven:
                continue
            six_mean = six["mean"]
            seven_mean = seven["mean"]
            difference = seven_mean - six_mean
            psd = pooled_std(six["std"], int(six["n"]), seven["std"], int(seven["n"]))
            effect = difference / psd if psd else None
            change = direction(feature, difference)
            rows.append(
                {
                    "metric": "7-minus-6 ensemble mean",
                    "feature": feature,
                    "ensemble_type": ensemble,
                    "unit6_value": format_float(six_mean),
                    "unit7_value": format_float(seven_mean),
                    "difference_7_minus_6": format_float(difference),
                    "standardized_effect_size": format_float(effect),
                    "direction_of_change": change,
                    "cautious_interpretation": interpretation_for(feature, ensemble, change),
                }
            )

    for feature, _, _ in PRIMARY_FEATURES:
        for unit_count in [6, 7]:
            loose = stats.get((unit_count, "loose_initial", feature))
            formed = stats.get((unit_count, "formed_perturbed", feature))
            if not loose or not formed:
                continue
            separation = formed["mean"] - loose["mean"]
            if feature == "RMSD_to_formed_seed_A":
                separation = loose["mean"] - formed["mean"]
            rows.append(
                {
                    "metric": "formed-minus-loose separation",
                    "feature": feature,
                    "ensemble_type": f"{unit_count}-unit formed_perturbed vs loose_initial",
                    "unit6_value": format_float(separation if unit_count == 6 else None),
                    "unit7_value": format_float(separation if unit_count == 7 else None),
                    "difference_7_minus_6": "",
                    "standardized_effect_size": "",
                    "direction_of_change": "larger separation is more formed-like contrast",
                    "cautious_interpretation": (
                        f"{unit_count} units show a formed-vs-loose contrast of {format_float(separation)} "
                        f"for {feature_metadata(feature)[0]}."
                    ),
                }
            )
    return rows


def contact_network_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    by_unit = {}
    for row in rows:
        unit = row.get("units_per_chain")
        if unit in {"6", "7"} and row.get("variant_id") == f"central{unit}_units":
            by_unit[int(unit)] = row
    metrics = [
        ("CYP_MEP_contact_count", "CYP/MEP-involving contacts"),
        ("CYP_MEP_contact_fraction", "CYP/MEP contact fraction"),
        ("GLU_contact_count", "GLU-involving contacts"),
        ("backbone_like_contact_count", "backbone-like contacts"),
        ("total_interchain_contacts", "total interchain contacts"),
        ("contacts_per_unit", "contacts per unit"),
        ("perturbation_contact_fraction_vs_reference_mean", "contact retention under perturbation"),
        ("perturbation_chain_graph_connected_probability", "interchain graph retention probability"),
    ]
    comparison_rows: list[dict[str, str]] = []
    if 6 not in by_unit or 7 not in by_unit:
        return comparison_rows
    for column, label in metrics:
        six = safe_float(by_unit[6].get(column))
        seven = safe_float(by_unit[7].get(column))
        if six is None or seven is None:
            continue
        diff = seven - six
        change = "7-unit stronger" if diff > 0 else "6-unit stronger" if diff < 0 else "similar"
        comparison_rows.append(
            {
                "metric": label,
                "feature": column,
                "ensemble_type": "central seed contact network",
                "unit6_value": format_float(six),
                "unit7_value": format_float(seven),
                "difference_7_minus_6": format_float(diff),
                "standardized_effect_size": "",
                "direction_of_change": change,
                "cautious_interpretation": (
                    f"Static contact-network summary gives {format_float(six)} for 6 units and "
                    f"{format_float(seven)} for 7 units."
                ),
            }
        )
    return comparison_rows


def activation_rows(rows: list[dict[str, str]], label: str) -> list[dict[str, str]]:
    selected = [
        "CYP_MEP_contact_fraction_vs_target",
        "contact_fraction_vs_target",
        "compactness_score",
        "axial_register_score",
        "rmsd_formedness_score",
        "seed_formation_score",
    ]
    primary_threshold = "0.75"
    output: list[dict[str, str]] = []
    by_key = {}
    for row in rows:
        if row.get("unit_count") in {"6", "7"} and row.get("feature") in selected:
            if row.get("threshold_fraction") == primary_threshold:
                by_key[(int(row["unit_count"]), row["feature"])] = row
    for feature in selected:
        six = by_key.get((6, feature))
        seven = by_key.get((7, feature))
        if not six or not seven:
            continue
        six_time = safe_float(six.get("median_activation_time"))
        seven_time = safe_float(seven.get("median_activation_time"))
        if six_time is None or seven_time is None:
            continue
        diff = seven_time - six_time
        change = "7-unit earlier" if diff < 0 else "6-unit earlier" if diff > 0 else "similar"
        output.append(
            {
                "metric": f"{label} median activation time",
                "feature": feature,
                "ensemble_type": "bridge ordering",
                "unit6_value": format_float(six_time),
                "unit7_value": format_float(seven_time),
                "difference_7_minus_6": format_float(diff),
                "standardized_effect_size": "",
                "direction_of_change": change,
                "cautious_interpretation": (
                    f"At the 0.75 formed-fraction threshold, {feature} activates at median t="
                    f"{format_float(six_time)} for 6 units and {format_float(seven_time)} for 7 units."
                ),
            }
        )
    return output


def mini_geometry_rows(geometry: list[dict[str, str]], helicity: list[dict[str, str]]) -> list[dict[str, str]]:
    comparisons: list[dict[str, str]] = []
    sources = [
        (geometry, "mini-hexaplex geometry", ["axial_extent_A", "radial_extent_mean_A", "suspicious_overlap_count"]),
        (helicity, "mini-hexaplex helicity", ["mean_helical_r2", "mean_twist_per_unit_deg", "helical_coherence_score"]),
    ]
    for rows, metric_label, columns in sources:
        central = {}
        for row in rows:
            variant = row.get("variant_id", "")
            if variant in {"central6_units", "central7_units"}:
                central[int(row["units_per_chain"])] = row
        if 6 not in central or 7 not in central:
            continue
        for column in columns:
            six = safe_float(central[6].get(column))
            seven = safe_float(central[7].get(column))
            if six is None or seven is None:
                continue
            diff = seven - six
            comparisons.append(
                {
                    "metric": metric_label,
                    "feature": column,
                    "ensemble_type": "central mini-hexaplex reference",
                    "unit6_value": format_float(six),
                    "unit7_value": format_float(seven),
                    "difference_7_minus_6": format_float(diff),
                    "standardized_effect_size": "",
                    "direction_of_change": "7-minus-6 reported",
                    "cautious_interpretation": f"Reference central mini-hexaplex {column} changes by {format_float(diff)} from 6 to 7 units.",
                }
            )
    return comparisons


def write_csv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def plot_feature_bars(summary_rows: list[dict[str, str]], plot_dir: Path) -> list[Path]:
    by_feature: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in summary_rows:
        by_feature[row["feature"]].append(row)
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return write_svg_feature_bars(by_feature, plot_dir)
    plot_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for feature, rows in by_feature.items():
        labels = []
        values = []
        errors = []
        for ensemble in ENSEMBLES:
            for unit in ["6", "7"]:
                match = next((row for row in rows if row["ensemble_type"] == ensemble and row["unit_count"] == unit), None)
                if not match:
                    continue
                labels.append(f"{unit} {ensemble.replace('_', ' ')}")
                values.append(float(match["mean"]))
                errors.append(float(match["std"]))
        if not values:
            continue
        fig, ax = plt.subplots(figsize=(9, 4.8))
        ax.bar(range(len(values)), values, yerr=errors, color=["#4c78a8", "#f58518"] * 3, capsize=3)
        ax.set_xticks(range(len(values)))
        ax.set_xticklabels(labels, rotation=30, ha="right")
        ax.set_ylabel(feature)
        ax.set_title(feature_metadata(feature)[0])
        fig.tight_layout()
        out_path = plot_dir / f"{feature}_6_vs_7.png"
        fig.savefig(out_path, dpi=180)
        plt.close(fig)
        written.append(out_path)
    return written


def write_svg_feature_bars(by_feature: dict[str, list[dict[str, str]]], plot_dir: Path) -> list[Path]:
    plot_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for feature, rows in by_feature.items():
        entries = []
        for ensemble in ENSEMBLES:
            for unit in ["6", "7"]:
                match = next((row for row in rows if row["ensemble_type"] == ensemble and row["unit_count"] == unit), None)
                if match:
                    entries.append((f"{unit} {ensemble}", float(match["mean"])))
        if not entries:
            continue
        max_value = max(value for _, value in entries)
        min_value = min(value for _, value in entries)
        span = max(max_value - min_value, 1e-9)
        width = 860
        height = 420
        left = 80
        bottom = 340
        bar_width = 90
        gap = 35
        colors = ["#4c78a8", "#f58518"]
        parts = [
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
            '<rect width="100%" height="100%" fill="white"/>',
            f'<text x="{width / 2:.0f}" y="32" text-anchor="middle" font-family="Arial" font-size="18">{feature_metadata(feature)[0]}</text>',
            f'<line x1="{left}" y1="60" x2="{left}" y2="{bottom}" stroke="#333"/>',
            f'<line x1="{left}" y1="{bottom}" x2="{width - 30}" y2="{bottom}" stroke="#333"/>',
        ]
        for index, (label, value) in enumerate(entries):
            normalized = (value - min_value) / span if span else 0.0
            bar_height = 35 + normalized * 230
            x = left + 25 + index * (bar_width + gap)
            y = bottom - bar_height
            parts.append(f'<rect x="{x}" y="{y:.2f}" width="{bar_width}" height="{bar_height:.2f}" fill="{colors[index % 2]}"/>')
            parts.append(f'<text x="{x + bar_width / 2:.0f}" y="{y - 8:.2f}" text-anchor="middle" font-family="Arial" font-size="11">{value:.3f}</text>')
            parts.append(
                f'<text x="{x + bar_width / 2:.0f}" y="{bottom + 16}" text-anchor="end" font-family="Arial" font-size="10" transform="rotate(-35 {x + bar_width / 2:.0f},{bottom + 16})">{label}</text>'
            )
        parts.append(f'<text x="16" y="205" text-anchor="middle" font-family="Arial" font-size="12" transform="rotate(-90 16,205)">{feature}</text>')
        parts.append("</svg>")
        out_path = plot_dir / f"{feature}_6_vs_7.svg"
        out_path.write_text("\n".join(parts) + "\n", encoding="utf-8")
        written.append(out_path)
    return written


def report_text(
    summary_rows: list[dict[str, str]],
    comparison_rows: list[dict[str, str]],
    ordering_rows: list[dict[str, str]],
    plot_paths: list[Path],
    args: argparse.Namespace,
) -> str:
    top_contact = [row for row in comparison_rows if row["metric"] in {"CYP/MEP-involving contacts", "contact retention under perturbation", "total interchain contacts"}]
    primary = [
        row
        for row in comparison_rows
        if row["metric"] == "7-minus-6 ensemble mean" and row["ensemble_type"] == "formed_perturbed"
    ]

    lines = [
        "# Focused 6-vs-7 Seed Transition Analysis",
        "",
        "This report reuses existing mini-hexaplex, seed-formation, contact-network, and bridge-ordering outputs. No ensembles or structures are regenerated.",
        "",
        "The framing is cautious: this is an SB-inspired/order-parameter comparison, not a physical proof of nucleation and not an atomistic Schrodinger bridge. The defensible question is whether the existing 6- and 7-unit data show stronger, earlier, or more perturbation-resistant structural organization in the 7-unit seed than the 6-unit seed.",
        "",
        "## Inputs",
        "",
        f"- Seed order parameters: `{args.seed_order_csv}`",
        f"- Contact-network summary: `{args.contact_summary_csv}`",
        f"- Bridge activation summary: `{args.bridge_activation_csv}`",
        f"- Angular-randomized bridge activation summary: `{args.angular_bridge_activation_csv}`",
        f"- Mini-hexaplex geometry/helicity summaries: `{args.mini_geometry_csv}`, `{args.mini_helicity_csv}`",
        "",
        "## Primary Formed-Perturbed 7-vs-6 Comparison",
        "",
        "| Feature | 6-unit mean | 7-unit mean | 7-minus-6 | Direction |",
        "|---|---:|---:|---:|---|",
    ]
    for row in primary:
        lines.append(
            f"| {row['feature']} | {row['unit6_value']} | {row['unit7_value']} | "
            f"{row['difference_7_minus_6']} | {row['direction_of_change']} |"
        )

    lines.extend(["", "## Contact-Network And Retention Summary", "", "| Metric | 6-unit | 7-unit | 7-minus-6 | Direction |", "|---|---:|---:|---:|---|"])
    for row in top_contact:
        lines.append(
            f"| {row['metric']} | {row['unit6_value']} | {row['unit7_value']} | "
            f"{row['difference_7_minus_6']} | {row['direction_of_change']} |"
        )

    lines.extend(["", "## Bridge Ordering", "", "| Source | Feature | 6-unit median t | 7-unit median t | Direction |", "|---|---|---:|---:|---|"])
    for row in ordering_rows:
        lines.append(
            f"| {row['metric']} | {row['feature']} | {row['unit6_value']} | "
            f"{row['unit7_value']} | {row['direction_of_change']} |"
        )

    lines.extend(["", "## Interpretation", ""])
    signals_7 = sum(1 for row in primary if row["direction_of_change"] == "7-unit stronger")
    signals_6 = sum(1 for row in primary if row["direction_of_change"] == "6-unit stronger")
    if signals_7 > signals_6:
        lines.append(
            "Across the formed-perturbed ensemble means, more primary features favor the 7-unit seed than the 6-unit seed. "
            "Together with the contact-network increase in total and CYP/MEP-involving contacts, this supports the cautious interpretation that the 7-unit seed shows stronger structural organization in the existing coordinate-derived data."
        )
    else:
        lines.append(
            "The primary formed-perturbed ensemble means are mixed. The existing data do not support a simple universal 7-unit strengthening across every feature."
        )
    lines.append(
        "The contact-network summary separates CYP/MEP-involving, GLU-involving, backbone-like, interchain, per-unit, and perturbation-retention measures where those columns are available. It does not provide a full CYP/MEP-CYP/MEP versus CYP/MEP-GLU pair-type decomposition, so that finer split should not be overclaimed here."
    )
    lines.append(
        "Activation times are order-parameter diagnostics along previously generated paired paths. Earlier activation should be read as earlier threshold crossing in this constructed feature space, not as a mechanistic assembly step."
    )

    lines.extend(["", "## Output Files", "", f"- Summary CSV: `{args.summary_csv}`", f"- Feature comparison CSV: `{args.comparison_csv}`", f"- Plot directory: `{args.plot_dir}`"])
    for path in plot_paths:
        lines.append(f"- Plot: `{path}`")
    lines.extend(["", "## Limitations", "", "- Existing ensembles have finite sample counts and are coordinate-derived perturbations.", "- The analysis compares order parameters and static contact summaries; it does not establish a physical nucleation threshold.", "- Contact classes are limited by the columns available in the prior contact-network CSV.", "- Bridge ordering is SB-inspired feature interpolation/ordering, not atomistic Schrodinger-bridge dynamics."])
    return "\n".join(lines) + "\n"


def run(args: argparse.Namespace) -> dict[str, object]:
    seed_rows = read_csv(args.seed_order_csv)
    contact_rows = read_csv(args.contact_summary_csv)
    bridge_activation = read_csv(args.bridge_activation_csv)
    angular_activation = read_csv(args.angular_bridge_activation_csv)
    geometry_rows = read_csv(args.mini_geometry_csv)
    helicity_rows = read_csv(args.mini_helicity_csv)

    summary_rows, stats = summarize_seed_order(seed_rows)
    comparison_rows = feature_comparisons(stats)
    comparison_rows.extend(contact_network_rows(contact_rows))
    comparison_rows.extend(activation_rows(bridge_activation, "loose-initial bridge"))
    angular_rows = activation_rows(angular_activation, "angular-randomized bridge")
    comparison_rows.extend(angular_rows)
    comparison_rows.extend(mini_geometry_rows(geometry_rows, helicity_rows))

    write_csv(args.summary_csv, summary_rows, SUMMARY_COLUMNS)
    write_csv(args.comparison_csv, comparison_rows, COMPARISON_COLUMNS)
    plot_paths = plot_feature_bars(summary_rows, args.plot_dir)

    ordering_rows = [
        row
        for row in comparison_rows
        if row["metric"] in {"loose-initial bridge median activation time", "angular-randomized bridge median activation time"}
    ]
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(report_text(summary_rows, comparison_rows, ordering_rows, plot_paths, args), encoding="utf-8")
    return {
        "summary_rows": len(summary_rows),
        "comparison_rows": len(comparison_rows),
        "plot_paths": plot_paths,
    }


def main() -> None:
    result = run(parse_args())
    print(f"Wrote {result['summary_rows']} summary rows")
    print(f"Wrote {result['comparison_rows']} comparison rows")
    print(f"Wrote {len(result['plot_paths'])} plots")


if __name__ == "__main__":
    main()
