#!/usr/bin/env python3
"""Convert the legacy contiguous block map into the canonical scaffold path map format."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


FIELDNAMES = [
    "map_name",
    "strand_id",
    "strand_label",
    "residue_index_in_pdb_order",
    "chain_id",
    "residue_name",
    "residue_number",
    "insertion_code",
    "residue_label",
    "source",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("inputs/metadata/strand_map_candidate.csv"))
    parser.add_argument("--out", type=Path, default=Path("inputs/metadata/scaffold_path_map_candidate.csv"))
    parser.add_argument("--map-name", default="contiguous_residue_blocks")
    parser.add_argument("--source", default="generated_from_strand_map_candidate")
    return parser.parse_args()


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def convert_rows(rows: list[dict[str, str]], map_name: str, source: str) -> list[dict[str, str]]:
    converted: list[dict[str, str]] = []
    for row in rows:
        block_id = row.get("block_id", "").strip()
        converted.append(
            {
                "map_name": map_name,
                "strand_id": block_id,
                "strand_label": f"block_{block_id}" if block_id else "",
                "residue_index_in_pdb_order": row.get("residue_index_in_pdb_order", ""),
                "chain_id": row.get("chain_id", ""),
                "residue_name": row.get("residue_name", ""),
                "residue_number": row.get("residue_number", ""),
                "insertion_code": row.get("insertion_code", ""),
                "residue_label": row.get("residue_label", ""),
                "source": source,
            }
        )
    return converted


def write_rows(rows: list[dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    rows = convert_rows(read_csv_rows(args.input), args.map_name, args.source)
    write_rows(rows, args.out)
    print(f"Wrote {args.out} with {len(rows)} residue row(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
