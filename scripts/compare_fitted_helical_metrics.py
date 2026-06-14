#!/usr/bin/env python3
"""Compare fitted-axis helical order summaries."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


REPORT_COLUMNS = [
    "model",
    "residue_count",
    "mean_radius_fitted",
    "axial_span",
    "angular_coverage_rad",
    "approximate_turns",
    "approximate_pitch_per_turn",
    "z_axis_angle_degrees",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics-dir", type=Path, default=Path("outputs/metrics/fitted_helical_order"))
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=Path("outputs/metrics/fitted_helical_order/fitted_helical_comparison.csv"),
    )
    parser.add_argument("--out-md", type=Path, default=Path("outputs/reports/fitted_helical_order_report.md"))
    return parser.parse_args()


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def model_name_from_summary(path: Path) -> str:
    suffix = "_fitted_helical_order_summary"
    stem = path.stem
    return stem[: -len(suffix)] if stem.endswith(suffix) else stem


def aggregate_summaries(metrics_dir: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for path in sorted(metrics_dir.glob("*_fitted_helical_order_summary.csv")):
        summary_rows = _rows(path)
        if not summary_rows:
            continue
        row = dict(summary_rows[0])
        row["model"] = model_name_from_summary(path)
        rows.append(row)
    return rows


def write_csv(rows: list[dict[str, str]], path: Path) -> None:
    if not rows:
        fieldnames = ["model"]
    else:
        fieldnames = ["model", *[column for column in rows[0] if column != "model"]]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _as_float(row: dict[str, str] | None, column: str) -> float:
    if row is None:
        return 0.0
    return float(row.get(column, "") or 0.0)


def _find_model(rows: list[dict[str, str]], model: str) -> dict[str, str] | None:
    return next((row for row in rows if row["model"] == model), None)


def report_notes(rows: list[dict[str, str]]) -> list[str]:
    notes: list[str] = []
    scaffold = _find_model(rows, "hexaplex_scaffold_only_complement_heavy_deduped")
    block1 = _find_model(rows, "scaffold_blocks_1_heavy_deduped")

    if scaffold is not None:
        notes.append(
            "For the scaffold-only complement, the fitted principal axis is "
            f"{_as_float(scaffold, 'z_axis_angle_degrees'):.3f} degrees from the global z-axis."
        )
    else:
        notes.append("The scaffold-only complement fitted-axis summary was not present.")

    if block1 is not None:
        coverage = _as_float(block1, "angular_coverage_rad")
        if coverage > 5.5:
            notes.append(
                "Scaffold block 1 still shows near-full angular coverage under the fitted-axis metric "
                f"({coverage:.3f} rad), supporting the current individual folded/twisted path interpretation."
            )
        else:
            notes.append(
                "Scaffold block 1 does not show near-full angular coverage under the fitted-axis metric "
                f"({coverage:.3f} rad), which weakens the individual folded/twisted path interpretation."
            )
    else:
        notes.append("Scaffold block 1 fitted-axis summary was not present.")

    notes.append(
        "These fitted-axis results are geometric descriptors, not proof of a physical helix or temporal assembly order."
    )
    notes.append("The candidate block mapping still requires validation against PyMOL colored strand paths.")
    return notes


def _markdown_table(rows: list[dict[str, str]], columns: list[str]) -> str:
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(row.get(column, "") for column in columns) + " |")
    return "\n".join(lines)


def write_markdown_report(rows: list[dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Fitted-Axis Helical Order Report",
        "",
        "## Scientific cautions",
        "",
        "The fitted-axis helical metrics are geometric descriptors computed from residue centroids. They do not prove a physical helix, formation pathway, or molecular dynamics. The d ~= 4.5 A diffraction feature remains a reciprocal-space scaffold signature in the working hypothesis, not a literal atom-contact assignment.",
        "",
        "The current block map is a candidate contiguous-residue map. Fitted-axis results help test whether candidate blocks behave like folded or twisted paths, but the mapping still requires validation against PyMOL colored strand paths.",
        "",
        "Unwrapped angle span, approximate turns, and pitch-like values depend on residue order, so they should be interpreted as order-aware descriptors rather than invariant helix parameters.",
        "",
        "## Comparison table",
        "",
        _markdown_table(rows, REPORT_COLUMNS),
        "",
        "## Automatically generated notes",
        "",
    ]
    lines.extend(f"- {note}" for note in report_notes(rows))
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    rows = aggregate_summaries(args.metrics_dir)
    write_csv(rows, args.out_csv)
    write_markdown_report(rows, args.out_md)
    print(f"Wrote {args.out_csv} with {len(rows)} model row(s)")
    print(f"Wrote {args.out_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
