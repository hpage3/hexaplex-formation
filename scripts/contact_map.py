#!/usr/bin/env python3
"""Write conservative real-space residue contact candidates from a PDB."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from hexaplex_formation.geometry import group_atoms_by_residue, min_inter_residue_distance  # noqa: E402
from hexaplex_formation.pdb_utils import heavy_atoms, load_pdb_atoms  # noqa: E402


FIELDNAMES = [
    "residue_i",
    "residue_j",
    "chain_i",
    "chain_j",
    "residue_name_i",
    "residue_name_j",
    "residue_number_i",
    "residue_number_j",
    "min_distance_A",
    "contact_cutoff_A",
    "atom_i",
    "atom_j",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--cutoff", type=float, default=4.5)
    parser.add_argument("--heavy-only", action="store_true", help="Exclude hydrogens. This is the default.")
    parser.add_argument("--all-atoms", action="store_true", help="Include hydrogens.")
    parser.add_argument("--residue-filter", default=None, help="Comma-separated residue names to include.")
    parser.add_argument("--min-sequence-separation", type=int, default=0)
    args = parser.parse_args()
    if args.heavy_only and args.all_atoms:
        parser.error("--heavy-only and --all-atoms cannot be supplied together")
    if args.cutoff <= 0:
        parser.error("--cutoff must be greater than zero")
    if args.min_sequence_separation < 0:
        parser.error("--min-sequence-separation must be non-negative")
    return args


def residue_label(key: tuple[str, str, int | None, str]) -> str:
    chain_id, residue_name, residue_number, insertion_code = key
    number = "" if residue_number is None else str(residue_number)
    residue = f"{residue_name}{number}{insertion_code}"
    return f"{chain_id}:{residue}" if chain_id else residue


def _residue_filter(value: str | None) -> set[str] | None:
    if value is None:
        return None
    names = {part.strip().upper() for part in value.split(",") if part.strip()}
    return names or None


def contact_rows(
    atoms,
    cutoff: float,
    use_heavy_only: bool = True,
    residue_filter: set[str] | None = None,
    min_sequence_separation: int = 0,
) -> list[dict[str, str]]:
    if use_heavy_only:
        atoms = heavy_atoms(atoms)
    if residue_filter is not None:
        atoms = [atom for atom in atoms if atom.residue_name.upper() in residue_filter]

    residues = list(group_atoms_by_residue(atoms).items())
    rows: list[dict[str, str]] = []
    for index, (key_i, atoms_i) in enumerate(residues):
        for key_j, atoms_j in residues[index + 1 :]:
            chain_i, residue_name_i, residue_number_i, _icode_i = key_i
            chain_j, residue_name_j, residue_number_j, _icode_j = key_j
            if (
                min_sequence_separation
                and chain_i == chain_j
                and residue_number_i is not None
                and residue_number_j is not None
                and abs(residue_number_i - residue_number_j) < min_sequence_separation
            ):
                continue

            minimum = min_inter_residue_distance(atoms_i, atoms_j, heavy_only=False)
            if minimum is None or minimum[0] > cutoff:
                continue
            min_distance, atom_i, atom_j = minimum
            rows.append(
                {
                    "residue_i": residue_label(key_i),
                    "residue_j": residue_label(key_j),
                    "chain_i": chain_i,
                    "chain_j": chain_j,
                    "residue_name_i": residue_name_i,
                    "residue_name_j": residue_name_j,
                    "residue_number_i": "" if residue_number_i is None else str(residue_number_i),
                    "residue_number_j": "" if residue_number_j is None else str(residue_number_j),
                    "min_distance_A": f"{min_distance:.6f}",
                    "contact_cutoff_A": f"{cutoff:.6f}",
                    "atom_i": atom_i.atom_name,
                    "atom_j": atom_j.atom_name,
                }
            )
    return rows


def main() -> int:
    args = parse_args()
    atoms = load_pdb_atoms(args.pdb)
    rows = contact_rows(
        atoms,
        cutoff=args.cutoff,
        use_heavy_only=not args.all_atoms,
        residue_filter=_residue_filter(args.residue_filter),
        min_sequence_separation=args.min_sequence_separation,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {args.out} with {len(rows)} contact candidate(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
