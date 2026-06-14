#!/usr/bin/env python3
"""Generate a blank scaffold path map template for manual/PyMOL strand assignment."""

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

from convert_strand_map_to_scaffold_path_map import FIELDNAMES  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pdb",
        type=Path,
        default=Path("outputs/intermediates/normalized_structures/hexaplex_scaffold_only_complement_heavy_deduped.pdb"),
    )
    parser.add_argument("--out", type=Path, default=Path("inputs/metadata/scaffold_path_map_manual_template.csv"))
    return parser.parse_args()


def residue_label(chain_id: str, residue_name: str, residue_number: int | None, insertion_code: str) -> str:
    number = "" if residue_number is None else str(residue_number)
    label = f"{residue_name}{number}{insertion_code}"
    return f"{chain_id}:{label}" if chain_id else label


def template_rows_from_pdb(pdb_path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    residue_keys = list(group_atoms_by_residue(load_pdb_atoms(pdb_path)).keys())
    for index, (chain_id, residue_name, residue_number, insertion_code) in enumerate(residue_keys, start=1):
        rows.append(
            {
                "map_name": "manual_pymol_scaffold_path_map",
                "strand_id": "",
                "strand_label": "",
                "residue_index_in_pdb_order": str(index),
                "chain_id": chain_id,
                "residue_name": residue_name,
                "residue_number": "" if residue_number is None else str(residue_number),
                "insertion_code": insertion_code,
                "residue_label": residue_label(chain_id, residue_name, residue_number, insertion_code),
                "source": "manual_template_from_pdb_order",
            }
        )
    return rows


def write_rows(rows: list[dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_metadata(pdb_path: Path, rows: list[dict[str, str]], out_path: Path) -> None:
    metadata = {
        "source_pdb": str(pdb_path),
        "template_csv": str(out_path),
        "residue_count": len(rows),
        "instructions": [
            "Inspect the scaffold path colors in PyMOL.",
            "Fill strand_id and strand_label for each residue.",
            "Save the edited copy as inputs/metadata/scaffold_path_map_manual.csv.",
            "Rerun scripts/run_scaffold_path_map_workflow.sh to validate and compare maps.",
        ],
        "caution": "Manual/PyMOL path assignment validates map consistency, not biological truth or temporal assembly order.",
    }
    metadata_path = out_path.with_name(f"{out_path.stem}.metadata.json")
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    rows = template_rows_from_pdb(args.pdb)
    write_rows(rows, args.out)
    write_metadata(args.pdb, rows, args.out)
    print(f"Wrote {args.out} with {len(rows)} residue row(s)")
    print(f"Wrote {args.out.with_name(f'{args.out.stem}.metadata.json')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
