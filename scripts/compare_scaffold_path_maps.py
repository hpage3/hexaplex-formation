#!/usr/bin/env python3
"""Compare two scaffold path maps by residue-to-strand assignment."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, OrderedDict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--map-a", type=Path, default=Path("inputs/metadata/scaffold_path_map_candidate.csv"))
    parser.add_argument("--map-b", type=Path, default=Path("inputs/metadata/scaffold_path_map_manual.csv"))
    parser.add_argument("--out-md", type=Path, default=Path("outputs/reports/scaffold_path_map_comparison.md"))
    return parser.parse_args()


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def residue_identity(row: dict[str, str]) -> tuple[str, str, str, str]:
    return (
        row.get("chain_id", "").strip(),
        row.get("residue_name", "").strip(),
        row.get("residue_number", "").strip(),
        row.get("insertion_code", "").strip(),
    )


def identity_label(identity: tuple[str, str, str, str]) -> str:
    chain_id, residue_name, residue_number, insertion_code = identity
    label = f"{residue_name}{residue_number}{insertion_code}"
    return f"{chain_id}:{label}" if chain_id else label


def assignment_index(rows: list[dict[str, str]]) -> dict[tuple[str, str, str, str], dict[str, str]]:
    return {residue_identity(row): row for row in rows}


def strand_sizes(rows: list[dict[str, str]]) -> Counter[str]:
    return Counter(row.get("strand_id", "").strip() for row in rows)


def strand_ranges(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: "OrderedDict[str, list[dict[str, str]]]" = OrderedDict()
    for row in rows:
        grouped.setdefault(row.get("strand_id", "").strip(), []).append(row)
    ranges: list[dict[str, str]] = []
    for strand_id, strand_rows in grouped.items():
        if not strand_id:
            continue
        ranges.append(
            {
                "strand_id": strand_id,
                "strand_label": strand_rows[0].get("strand_label", ""),
                "residue_count": str(len(strand_rows)),
                "first_residue_label": strand_rows[0].get("residue_label", ""),
                "last_residue_label": strand_rows[-1].get("residue_label", ""),
            }
        )
    return ranges


def markdown_table(rows: list[dict[str, str]], columns: list[str]) -> str:
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(row.get(column, "") for column in columns) + " |")
    return "\n".join(lines)


def comparison_report(map_a_path: Path, map_b_path: Path) -> str:
    if not map_b_path.exists():
        return "\n".join(
            [
                "# Scaffold Path Map Comparison",
                "",
                f"- Candidate map: {map_a_path}",
                f"- Manual/PyMOL map: {map_b_path}",
                "- Status: manual/PyMOL map is not available yet.",
                "",
                "Only the candidate map exists. Create `inputs/metadata/scaffold_path_map_manual.csv` from the manual template after PyMOL colored path inspection, then rerun this comparison.",
                "",
                "Different strand labels can represent equivalent maps if labels are permuted; compare residue membership, not only label names.",
                "",
            ]
        )

    rows_a = read_csv_rows(map_a_path)
    rows_b = read_csv_rows(map_b_path)
    by_a = assignment_index(rows_a)
    by_b = assignment_index(rows_b)
    all_identities = sorted(set(by_a) | set(by_b), key=lambda item: (item[0], int(item[2] or 0), item[3], item[1]))

    same = 0
    different_rows: list[dict[str, str]] = []
    for identity in all_identities:
        row_a = by_a.get(identity, {})
        row_b = by_b.get(identity, {})
        strand_a = row_a.get("strand_id", "")
        strand_b = row_b.get("strand_id", "")
        if row_a and row_b and strand_a == strand_b:
            same += 1
        else:
            different_rows.append(
                {
                    "residue": identity_label(identity),
                    "map_a_strand": strand_a,
                    "map_a_label": row_a.get("strand_label", ""),
                    "map_b_strand": strand_b,
                    "map_b_label": row_b.get("strand_label", ""),
                }
            )

    size_rows = []
    sizes_a = strand_sizes(rows_a)
    sizes_b = strand_sizes(rows_b)
    for strand_id in sorted(set(sizes_a) | set(sizes_b)):
        if not strand_id:
            continue
        size_rows.append(
            {
                "strand_id": strand_id,
                "map_a_count": str(sizes_a.get(strand_id, 0)),
                "map_b_count": str(sizes_b.get(strand_id, 0)),
            }
        )

    lines = [
        "# Scaffold Path Map Comparison",
        "",
        f"- Map A: {map_a_path}",
        f"- Map B: {map_b_path}",
        f"- Same residue-to-strand assignments: {same}",
        f"- Different or missing assignments: {len(different_rows)}",
        "",
        "Different strand labels can represent equivalent maps if labels are permuted; compare residue membership and colored path continuity, not only label names.",
        "",
        "## Strand sizes",
        "",
        markdown_table(size_rows, ["strand_id", "map_a_count", "map_b_count"]),
        "",
        "## Map A strand ranges",
        "",
        markdown_table(strand_ranges(rows_a), ["strand_id", "strand_label", "residue_count", "first_residue_label", "last_residue_label"]),
        "",
        "## Map B strand ranges",
        "",
        markdown_table(strand_ranges(rows_b), ["strand_id", "strand_label", "residue_count", "first_residue_label", "last_residue_label"]),
        "",
        "## Differing residue assignments",
        "",
        markdown_table(different_rows[:100], ["residue", "map_a_strand", "map_a_label", "map_b_strand", "map_b_label"]),
        "",
    ]
    if len(different_rows) > 100:
        lines.append(f"Only the first 100 of {len(different_rows)} differing assignments are shown.")
        lines.append("")
    return "\n".join(lines)


def main() -> int:
    args = parse_args()
    args.out_md.parent.mkdir(parents=True, exist_ok=True)
    args.out_md.write_text(comparison_report(args.map_a, args.map_b), encoding="utf-8")
    print(f"Wrote {args.out_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
