#!/usr/bin/env python3
"""Compute a simple all-pairs real-space atom-distance histogram."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from hexaplex_formation.pdb_utils import PDBAtom, is_hydrogen, load_pdb_atoms  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb", required=True, type=Path, help="Input PDB file.")
    parser.add_argument("--out", required=True, type=Path, help="Output CSV path.")
    parser.add_argument("--bin-width", type=float, default=0.1, help="Histogram bin width in angstroms.")
    parser.add_argument("--max-distance", type=float, default=30.0, help="Maximum distance in angstroms.")
    parser.add_argument(
        "--residue-filter",
        default=None,
        help="Comma-separated residue names to include, for example GLU,ALA.",
    )
    parser.add_argument(
        "--heavy-only",
        action="store_true",
        help="Exclude hydrogens by element or atom name. This is the default.",
    )
    parser.add_argument(
        "--all-atoms",
        action="store_true",
        help="Include hydrogens in the pair-distance histogram.",
    )
    args = parser.parse_args()
    if args.heavy_only and args.all_atoms:
        parser.error("--heavy-only and --all-atoms cannot be supplied together")
    return args


def _residue_filter(value: str | None) -> set[str] | None:
    if value is None:
        return None
    names = {part.strip().upper() for part in value.split(",") if part.strip()}
    return names or None


def filter_atoms(
    atoms: list[PDBAtom],
    residue_filter: set[str] | None,
    heavy_only: bool,
) -> list[PDBAtom]:
    filtered: list[PDBAtom] = []
    for atom in atoms:
        if residue_filter is not None and atom.residue_name.upper() not in residue_filter:
            continue
        if heavy_only and is_hydrogen(atom):
            continue
        filtered.append(atom)
    return filtered


def pair_distance_histogram(
    atoms: list[PDBAtom],
    bin_width: float,
    max_distance: float,
) -> list[int]:
    if bin_width <= 0:
        raise ValueError("--bin-width must be greater than zero")
    if max_distance <= 0:
        raise ValueError("--max-distance must be greater than zero")

    bin_count = math.ceil(max_distance / bin_width)
    counts = [0] * bin_count

    for i, atom_i in enumerate(atoms):
        for atom_j in atoms[i + 1 :]:
            dx = atom_i.x - atom_j.x
            dy = atom_i.y - atom_j.y
            dz = atom_i.z - atom_j.z
            distance = math.sqrt(dx * dx + dy * dy + dz * dz)
            if distance > max_distance:
                continue
            index = min(int(distance / bin_width), bin_count - 1)
            counts[index] += 1

    return counts


def write_histogram(path: Path, counts: list[int], bin_width: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["bin_start_A", "bin_end_A", "pair_count"])
        writer.writeheader()
        for index, count in enumerate(counts):
            start = index * bin_width
            end = start + bin_width
            writer.writerow(
                {
                    "bin_start_A": f"{start:.6f}",
                    "bin_end_A": f"{end:.6f}",
                    "pair_count": count,
                }
            )


def main() -> int:
    args = parse_args()
    atoms = load_pdb_atoms(args.pdb)
    use_heavy_only = not args.all_atoms
    filtered_atoms = filter_atoms(atoms, _residue_filter(args.residue_filter), use_heavy_only)
    counts = pair_distance_histogram(filtered_atoms, args.bin_width, args.max_distance)
    write_histogram(args.out, counts, args.bin_width)

    print(
        f"Wrote {args.out} from {len(filtered_atoms)} atom(s); "
        f"{sum(counts)} pair(s) within {args.max_distance:g} A"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
