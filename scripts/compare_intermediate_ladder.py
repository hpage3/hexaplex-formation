#!/usr/bin/env python3
"""Aggregate and report metrics for candidate intermediate ladder models."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


OUTPUT_COLUMNS = [
    "model_name",
    "atom_mode",
    "included_blocks",
    "includes_hexads",
    "atom_count",
    "residue_count",
    "contact_count_4p5A",
    "motif_GLU_GLU",
    "motif_GLU_any",
    "hbond_candidate_count",
    "mean_radius_xy",
    "z_span",
    "angular_coverage_rad",
]

REPORT_COLUMNS = OUTPUT_COLUMNS


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=Path("outputs/metrics/intermediate_ladder_summary.csv"))
    parser.add_argument("--metrics-dir", type=Path, default=Path("outputs/metrics/intermediate_ladder"))
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=Path("outputs/metrics/intermediate_ladder/intermediate_ladder_comparison.csv"),
    )
    parser.add_argument("--out-md", type=Path, default=Path("outputs/reports/intermediate_ladder_report.md"))
    return parser.parse_args()


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def metric_base_from_name(path_or_name: str | Path) -> str:
    stem = Path(path_or_name).stem
    for suffix in [
        "_contacts_4p5A",
        "_glu_motifs",
        "_helical_order_summary",
        "_hbond_candidates",
    ]:
        if stem.endswith(suffix):
            return stem[: -len(suffix)]
    return stem


def _empty_metric_row(summary_row: dict[str, str]) -> dict[str, str]:
    return {
        "model_name": summary_row["model_name"],
        "atom_mode": summary_row["atom_mode"],
        "included_blocks": summary_row["included_blocks"],
        "includes_hexads": summary_row["includes_hexads"],
        "atom_count": summary_row["atom_count"],
        "residue_count": summary_row["residue_count"],
        "contact_count_4p5A": "0",
        "motif_GLU_GLU": "0",
        "motif_GLU_any": "0",
        "hbond_candidate_count": "0",
        "mean_radius_xy": "",
        "z_span": "",
        "angular_coverage_rad": "",
    }


def _read_motif_counts(path: Path) -> dict[str, str]:
    counts = {"GLU-GLU": "0", "GLU-any": "0"}
    for row in _rows(path):
        motif_type = row.get("motif_type", "")
        if motif_type in counts:
            counts[motif_type] = row.get("count", "0") or "0"
    return counts


def aggregate_ladder_metrics(summary_path: Path, metrics_dir: Path) -> list[dict[str, str]]:
    summary_rows = _rows(summary_path)
    rows_by_model = {row["model_name"]: _empty_metric_row(row) for row in summary_rows}

    for path in sorted(metrics_dir.glob("*_contacts_4p5A.csv")):
        model_name = metric_base_from_name(path)
        if model_name in rows_by_model:
            rows_by_model[model_name]["contact_count_4p5A"] = str(len(_rows(path)))

    for path in sorted(metrics_dir.glob("*_glu_motifs.csv")):
        model_name = metric_base_from_name(path)
        if model_name in rows_by_model:
            counts = _read_motif_counts(path)
            rows_by_model[model_name]["motif_GLU_GLU"] = counts["GLU-GLU"]
            rows_by_model[model_name]["motif_GLU_any"] = counts["GLU-any"]

    for path in sorted(metrics_dir.glob("*_helical_order_summary.csv")):
        model_name = metric_base_from_name(path)
        if model_name not in rows_by_model:
            continue
        rows = _rows(path)
        if not rows:
            continue
        row = rows[0]
        rows_by_model[model_name]["mean_radius_xy"] = row.get("mean_radius_xy", "")
        rows_by_model[model_name]["z_span"] = row.get("z_span", "")
        rows_by_model[model_name]["angular_coverage_rad"] = row.get("angular_coverage_rad", "")

    for path in sorted(metrics_dir.glob("*_hbond_candidates.csv")):
        model_name = metric_base_from_name(path)
        if model_name in rows_by_model:
            rows_by_model[model_name]["hbond_candidate_count"] = str(len(_rows(path)))

    return [rows_by_model[row["model_name"]] for row in summary_rows]


def write_csv(rows: list[dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _as_int(row: dict[str, str], column: str) -> int:
    return int(row.get(column, "") or 0)


def _as_float(row: dict[str, str], column: str) -> float:
    return float(row.get(column, "") or 0.0)


def block_count(row: dict[str, str]) -> int:
    included = row.get("included_blocks", "")
    if not included:
        return 0
    return len([part for part in included.split(",") if part])


def scaffold_only_ladder_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return [
        row
        for row in rows
        if row["model_name"].startswith("scaffold_blocks_")
        and row["atom_mode"] == "heavy_deduped"
        and row["includes_hexads"] == "no"
    ]


def hexads_plus_ladder_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return [
        row
        for row in rows
        if row["model_name"].startswith("hexads_plus_scaffold_blocks_")
        and row["atom_mode"] == "heavy_deduped"
        and row["includes_hexads"] == "yes"
    ]


def _first_row_exceeding(rows: list[dict[str, str]], column: str, threshold: float) -> dict[str, str] | None:
    for row in sorted(rows, key=block_count):
        if _as_float(row, column) > threshold:
            return row
    return None


def report_notes(rows: list[dict[str, str]]) -> list[str]:
    scaffold_rows = scaffold_only_ladder_rows(rows)
    hexad_rows = hexads_plus_ladder_rows(rows)
    notes: list[str] = []

    for threshold in [3.0, 5.5]:
        row = _first_row_exceeding(scaffold_rows, "angular_coverage_rad", threshold)
        if row is None:
            notes.append(f"No scaffold-only candidate intermediate exceeds {threshold:.1f} rad angular coverage.")
        else:
            notes.append(
                f"The first scaffold-only candidate intermediate exceeding {threshold:.1f} rad angular coverage "
                f"uses {block_count(row)} block(s): `{row['model_name']}` "
                f"({_as_float(row, 'angular_coverage_rad'):.6f} rad). This is consistent with emerging coarse "
                "helical order but does not prove temporal assembly order."
            )

    glu_row = next((row for row in sorted(scaffold_rows, key=block_count) if _as_int(row, "motif_GLU_GLU") > 0), None)
    if glu_row is None:
        notes.append("No scaffold-only candidate intermediate has GLU-GLU motifs in this ladder.")
    else:
        notes.append(
            f"The first scaffold-only candidate intermediate with GLU-GLU motifs uses {block_count(glu_row)} "
            f"block(s): `{glu_row['model_name']}` ({_as_int(glu_row, 'motif_GLU_GLU')} motifs). This requires "
            "validation against PyMOL strand mapping and/or simulation."
        )

    paired = []
    hexad_by_blocks = {row["included_blocks"]: row for row in hexad_rows}
    for scaffold_row in scaffold_rows:
        hexad_row = hexad_by_blocks.get(scaffold_row["included_blocks"])
        if hexad_row is not None:
            paired.append((scaffold_row, hexad_row))
    changed = [
        (scaffold, hexad)
        for scaffold, hexad in paired
        if scaffold["contact_count_4p5A"] != hexad["contact_count_4p5A"]
        or scaffold["angular_coverage_rad"] != hexad["angular_coverage_rad"]
        or scaffold["z_span"] != hexad["z_span"]
    ]
    if changed:
        first_scaffold, first_hexad = changed[0]
        notes.append(
            "Adding hexads changes contact counts or helical summaries for at least one matched block count; "
            f"for {block_count(first_scaffold)} block(s), contacts change from "
            f"{first_scaffold['contact_count_4p5A']} to {first_hexad['contact_count_4p5A']}. This is a "
            "candidate intermediate comparison and requires validation against PyMOL strand mapping and/or simulation."
        )
    else:
        notes.append("Adding hexads does not change matched contact counts or helical summaries in this table.")

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
        "# Intermediate Ladder Report",
        "",
        "## Scientific cautions",
        "",
        "The d ~= 4.5 A diffraction feature is a reciprocal-space scaffold signature in the current working hypothesis; it is not a literal 4.5 A atom-contact assignment. Contact maps, GLU motifs, and helical summaries are real-space structural summaries. Hydrogen-bond outputs are rough geometric candidates, not definitive hydrogen bonds.",
        "",
        "The block mapping used here is a candidate contiguous-residue map. These structures are candidate intermediate models for hypothesis testing and do not prove temporal assembly order. Results require validation against PyMOL strand mapping and/or simulation.",
        "",
        "## Ladder model table",
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
    rows = aggregate_ladder_metrics(args.summary, args.metrics_dir)
    write_csv(rows, args.out_csv)
    write_markdown_report(rows, args.out_md)
    print(f"Wrote {args.out_csv} with {len(rows)} model row(s)")
    print(f"Wrote {args.out_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
