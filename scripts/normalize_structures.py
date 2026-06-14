#!/usr/bin/env python3
"""Create all-atom and heavy-atom normalized working PDBs from raw inputs."""

from __future__ import annotations

import csv
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from hexaplex_formation.pdb_utils import (  # noqa: E402
    dedupe_exact_atoms,
    heavy_atoms,
    is_hydrogen,
    load_pdb_atoms,
    write_pdb_atoms,
)


FIELDNAMES = [
    "structure_file",
    "original_atom_count",
    "allatom_deduped_atom_count",
    "heavy_atom_count",
    "heavy_deduped_atom_count",
    "hydrogens_detected",
    "duplicates_removed_allatom",
    "duplicates_removed_after_heavy_filter",
    "output_allatom_deduped_pdb",
    "output_heavy_pdb",
    "output_heavy_deduped_pdb",
]


def _relative(path: Path) -> str:
    return str(path.relative_to(REPO_ROOT))


def normalize_structure(input_path: Path, output_dir: Path) -> dict[str, str | int]:
    atoms = load_pdb_atoms(input_path)
    allatom_deduped = dedupe_exact_atoms(atoms)
    heavy = heavy_atoms(atoms)
    heavy_deduped = dedupe_exact_atoms(heavy)

    allatom_deduped_path = output_dir / f"{input_path.stem}_allatom_deduped.pdb"
    heavy_path = output_dir / f"{input_path.stem}_heavy.pdb"
    heavy_deduped_path = output_dir / f"{input_path.stem}_heavy_deduped.pdb"

    write_pdb_atoms(allatom_deduped, allatom_deduped_path)
    write_pdb_atoms(heavy, heavy_path)
    write_pdb_atoms(heavy_deduped, heavy_deduped_path)

    return {
        "structure_file": _relative(input_path),
        "original_atom_count": len(atoms),
        "allatom_deduped_atom_count": len(allatom_deduped),
        "heavy_atom_count": len(heavy),
        "heavy_deduped_atom_count": len(heavy_deduped),
        "hydrogens_detected": sum(1 for atom in atoms if is_hydrogen(atom)),
        "duplicates_removed_allatom": len(atoms) - len(allatom_deduped),
        "duplicates_removed_after_heavy_filter": len(heavy) - len(heavy_deduped),
        "output_allatom_deduped_pdb": _relative(allatom_deduped_path),
        "output_heavy_pdb": _relative(heavy_path),
        "output_heavy_deduped_pdb": _relative(heavy_deduped_path),
    }


def main() -> int:
    input_dir = REPO_ROOT / "inputs" / "structures"
    output_dir = REPO_ROOT / "outputs" / "intermediates" / "normalized_structures"
    summary_path = REPO_ROOT / "outputs" / "metrics" / "structure_normalization_summary.csv"

    output_dir.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)

    pdb_paths = sorted(input_dir.glob("*.pdb"))
    rows = [normalize_structure(path, output_dir) for path in pdb_paths]

    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Normalized {len(rows)} PDB file(s) from {_relative(input_dir)}")
    print(f"Wrote {_relative(summary_path)}")
    for row in rows:
        print(
            f"- {row['structure_file']}: {row['original_atom_count']} original, "
            f"{row['allatom_deduped_atom_count']} all-atom deduped, "
            f"{row['heavy_deduped_atom_count']} heavy deduped"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
