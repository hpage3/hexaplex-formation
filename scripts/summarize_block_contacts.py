#!/usr/bin/env python3
"""Summarize block contact decomposition CSVs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


CATEGORY_FIELDNAMES = [
    "contact_category",
    "contact_count",
    "GLU_involved_count",
    "GLU_GLU_count",
    "backbone_O_to_GLU_sidechain_O_count",
    "GLU_sidechain_O_to_GLU_sidechain_O_count",
    "backbone_N_to_GLU_sidechain_O_count",
]

BLOCK_PAIR_FIELDNAMES = [
    "block_pair",
    "contact_count",
    "GLU_involved_count",
    "GLU_GLU_count",
    "min_distance_A",
]

BOOLEAN_COUNT_COLUMNS = {
    "is_GLU_involved": "GLU_involved_count",
    "is_GLU_GLU": "GLU_GLU_count",
    "is_backbone_O_to_GLU_sidechain_O": "backbone_O_to_GLU_sidechain_O_count",
    "is_GLU_sidechain_O_to_GLU_sidechain_O": "GLU_sidechain_O_to_GLU_sidechain_O_count",
    "is_backbone_N_to_GLU_sidechain_O": "backbone_N_to_GLU_sidechain_O_count",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--decomposition-csv", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    return parser.parse_args()


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _yes(row: dict[str, str], column: str) -> bool:
    return row.get(column, "").lower() == "yes"


def _empty_category_row(category: str) -> dict[str, str]:
    row = {field: "0" for field in CATEGORY_FIELDNAMES}
    row["contact_category"] = category
    return row


def summarize_by_category(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    summaries: dict[str, dict[str, str]] = {}
    for row in rows:
        category = row["contact_category"]
        summary = summaries.setdefault(category, _empty_category_row(category))
        summary["contact_count"] = str(int(summary["contact_count"]) + 1)
        for source_column, output_column in BOOLEAN_COUNT_COLUMNS.items():
            if _yes(row, source_column):
                summary[output_column] = str(int(summary[output_column]) + 1)
    return [summaries[key] for key in sorted(summaries)]


def block_pair_for_row(row: dict[str, str]) -> str:
    component_i = row["component_i"]
    component_j = row["component_j"]
    if component_i == "scaffold" and component_j == "scaffold":
        blocks = sorted([row["block_i"], row["block_j"]], key=int)
        return f"{blocks[0]}--{blocks[1]}"
    if component_i == "scaffold" or component_j == "scaffold":
        return "scaffold--hexad_or_other"
    return "hexad_or_other--hexad_or_other"


def summarize_by_block_pair(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    summaries: dict[str, dict[str, str]] = {}
    for row in rows:
        pair = block_pair_for_row(row)
        summary = summaries.setdefault(
            pair,
            {
                "block_pair": pair,
                "contact_count": "0",
                "GLU_involved_count": "0",
                "GLU_GLU_count": "0",
                "min_distance_A": "",
            },
        )
        summary["contact_count"] = str(int(summary["contact_count"]) + 1)
        if _yes(row, "is_GLU_involved"):
            summary["GLU_involved_count"] = str(int(summary["GLU_involved_count"]) + 1)
        if _yes(row, "is_GLU_GLU"):
            summary["GLU_GLU_count"] = str(int(summary["GLU_GLU_count"]) + 1)
        distance = float(row["min_distance_A"])
        if not summary["min_distance_A"] or distance < float(summary["min_distance_A"]):
            summary["min_distance_A"] = f"{distance:.6f}"
    return [summaries[key] for key in sorted(summaries)]


def _block_pair_summary_path(out_path: Path) -> Path:
    if out_path.name.endswith("_block_contact_summary.csv"):
        return out_path.with_name(out_path.name.replace("_block_contact_summary.csv", "_block_pair_summary.csv"))
    return out_path.with_name(f"{out_path.stem}_block_pair_summary.csv")


def write_csv(rows: list[dict[str, str]], path: Path, fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    rows = _rows(args.decomposition_csv)
    write_csv(summarize_by_category(rows), args.out, CATEGORY_FIELDNAMES)
    pair_out = _block_pair_summary_path(args.out)
    write_csv(summarize_by_block_pair(rows), pair_out, BLOCK_PAIR_FIELDNAMES)
    print(f"Wrote {args.out}")
    print(f"Wrote {pair_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
