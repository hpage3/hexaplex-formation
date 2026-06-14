#!/usr/bin/env python3
"""Inventory PDB structures and write conservative structure-level metrics."""

from __future__ import annotations

import csv
import argparse
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from hexaplex_formation.pdb_utils import (  # noqa: E402
    atom_count,
    bounding_box,
    centroid,
    chain_ids,
    load_pdb_atoms,
    residue_count,
    residue_names,
)


FIELDNAMES = [
    "structure_file",
    "atom_count",
    "residue_count",
    "chain_count",
    "chain_ids",
    "residue_names",
    "min_x",
    "max_x",
    "min_y",
    "max_y",
    "min_z",
    "max_z",
    "centroid_x",
    "centroid_y",
    "centroid_z",
]


def _format_float(value: float | None) -> str:
    return "" if value is None else f"{value:.6f}"


def _resolve_repo_path(path: Path) -> Path:
    return path if path.is_absolute() else REPO_ROOT / path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("inputs/structures"),
        help="Directory containing .pdb files.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("outputs/metrics/structure_inventory.csv"),
        help="Output inventory CSV path.",
    )
    return parser.parse_args()


def inventory_structure(path: Path) -> dict[str, str | int]:
    atoms = load_pdb_atoms(path)
    chains = chain_ids(atoms)
    bbox = bounding_box(atoms)
    center = centroid(atoms)

    min_x = max_x = min_y = max_y = min_z = max_z = None
    if bbox is not None:
        min_x, max_x, min_y, max_y, min_z, max_z = bbox

    centroid_x = centroid_y = centroid_z = None
    if center is not None:
        centroid_x, centroid_y, centroid_z = center

    return {
        "structure_file": str(path.relative_to(REPO_ROOT)),
        "atom_count": atom_count(atoms),
        "residue_count": residue_count(atoms),
        "chain_count": len(chains),
        "chain_ids": ";".join(chains),
        "residue_names": ";".join(residue_names(atoms)),
        "min_x": _format_float(min_x),
        "max_x": _format_float(max_x),
        "min_y": _format_float(min_y),
        "max_y": _format_float(max_y),
        "min_z": _format_float(min_z),
        "max_z": _format_float(max_z),
        "centroid_x": _format_float(centroid_x),
        "centroid_y": _format_float(centroid_y),
        "centroid_z": _format_float(centroid_z),
    }


def main() -> int:
    args = parse_args()
    structures_dir = _resolve_repo_path(args.input_dir)
    output_path = _resolve_repo_path(args.out)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    pdb_paths = sorted(structures_dir.glob("*.pdb"))
    rows = [inventory_structure(path) for path in pdb_paths]

    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Scanned {len(pdb_paths)} PDB file(s) in {structures_dir.relative_to(REPO_ROOT)}")
    print(f"Wrote {output_path.relative_to(REPO_ROOT)}")
    for row in rows:
        print(
            f"- {row['structure_file']}: {row['atom_count']} atoms, "
            f"{row['residue_count']} residues, chains={row['chain_ids'] or '(blank)'}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
