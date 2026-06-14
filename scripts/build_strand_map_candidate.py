#!/usr/bin/env python3
"""Build a candidate scaffold block map from contiguous PDB residue order."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from hexaplex_formation.geometry import group_atoms_by_residue  # noqa: E402
from hexaplex_formation.pdb_utils import load_pdb_atoms  # noqa: E402


FIELDNAMES = [
    "block_id",
    "residue_index_in_pdb_order",
    "chain_id",
    "residue_name",
    "residue_number",
    "insertion_code",
    "residue_label",
]

CAUTION = "Candidate scaffold block map; not yet validated against PyMOL colored strand paths."


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pdb",
        type=Path,
        default=Path("outputs/intermediates/normalized_structures/hexaplex_scaffold_only_complement_heavy_deduped.pdb"),
    )
    parser.add_argument("--out", type=Path, default=Path("inputs/metadata/strand_map_candidate.csv"))
    parser.add_argument("--blocks", type=int, default=6)
    args = parser.parse_args()
    if args.blocks <= 0:
        parser.error("--blocks must be greater than zero")
    return args


def residue_label(chain_id: str, residue_name: str, residue_number: int | None, insertion_code: str) -> str:
    number = "" if residue_number is None else str(residue_number)
    label = f"{residue_name}{number}{insertion_code}"
    return f"{chain_id}:{label}" if chain_id else label


def assign_block(index_zero_based: int, residue_count: int, blocks: int) -> int:
    if residue_count <= 0:
        raise ValueError("residue_count must be greater than zero")
    return (index_zero_based * blocks) // residue_count + 1


def build_rows_from_residue_keys(
    residue_keys: list[tuple[str, str, int | None, str]],
    blocks: int,
) -> list[dict[str, str]]:
    residue_count = len(residue_keys)
    if residue_count == 0:
        return []

    rows: list[dict[str, str]] = []
    for index, (chain_id, residue_name, residue_number, insertion_code) in enumerate(residue_keys, start=1):
        rows.append(
            {
                "block_id": str(assign_block(index - 1, residue_count, blocks)),
                "residue_index_in_pdb_order": str(index),
                "chain_id": chain_id,
                "residue_name": residue_name,
                "residue_number": "" if residue_number is None else str(residue_number),
                "insertion_code": insertion_code,
                "residue_label": residue_label(chain_id, residue_name, residue_number, insertion_code),
            }
        )
    return rows


def residues_per_block(rows: list[dict[str, str]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        block_id = row["block_id"]
        counts[block_id] = counts.get(block_id, 0) + 1
    return counts


def write_map(rows: list[dict[str, str]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_metadata(source_pdb: Path, rows: list[dict[str, str]], blocks: int, out_path: Path) -> None:
    metadata = {
        "source_pdb": str(source_pdb),
        "residue_count": len(rows),
        "blocks": blocks,
        "residues_per_block": residues_per_block(rows),
        "mapping_method": "contiguous_residue_order_candidate",
        "caution": CAUTION,
    }
    metadata_path = out_path.with_name(f"{out_path.stem}.metadata.json")
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    atoms = load_pdb_atoms(args.pdb)
    residue_keys = list(group_atoms_by_residue(atoms).keys())
    rows = build_rows_from_residue_keys(residue_keys, args.blocks)
    write_map(rows, args.out)
    write_metadata(args.pdb, rows, args.blocks, args.out)
    print(f"Wrote {args.out} with {len(rows)} residue row(s)")
    print(f"Wrote {args.out.with_name(f'{args.out.stem}.metadata.json')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
