#!/usr/bin/env python3
"""Derive a no-CH2-COOH model from Nick's clarified ideal 16-mer PDB."""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from hexaplex_formation.pdb_utils import (  # noqa: E402
    dedupe_exact_atoms,
    heavy_atoms,
    load_pdb_atoms,
    write_pdb_atoms,
)


DEFAULT_INPUT = ROOT / "inputs" / "nick_ideal_models" / "Hexaplex_AntiParallel_30deg_Ideal.pdb"
DEFAULT_OUTPUT_DIR = ROOT / "inputs" / "nick_ideal_models" / "derived"
REMOVED_GLU_ATOMS = {"CG", "CD", "OE1", "OE2", "HG2", "HG3"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-pdb", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def element_counts(atoms) -> dict[str, int]:
    return dict(sorted(Counter(atom.element.upper() for atom in atoms).items()))


def should_remove(atom) -> bool:
    return atom.residue_name == "GLU" and atom.atom_name.strip().upper() in REMOVED_GLU_ATOMS


def write_xyz(path: Path, atoms) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f"{len(atoms)}\n")
        handle.write(
            "Heavy atoms only, exact deduplication applied after removing GLU "
            "CG/CD/OE1/OE2/HG2/HG3 no-CH2-COOH atoms.\n"
        )
        for atom in atoms:
            element = atom.element.upper()
            handle.write(f"{element:<6}{atom.x:12.6f}{atom.y:12.6f}{atom.z:12.6f}\n")


def main() -> int:
    args = parse_args()
    atoms = load_pdb_atoms(args.input_pdb)
    removed = [atom for atom in atoms if should_remove(atom)]
    kept = [atom for atom in atoms if not should_remove(atom)]
    heavy_deduped = dedupe_exact_atoms(heavy_atoms(kept))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    pdb_output = args.output_dir / "Hexaplex_AntiParallel_30deg_Ideal_no_CH2_COOH.pdb"
    xyz_output = args.output_dir / "Hexaplex_AntiParallel_30deg_Ideal_no_CH2_COOH.xyz"
    provenance_output = args.output_dir / "Hexaplex_AntiParallel_30deg_Ideal_no_CH2_COOH_provenance.json"
    readme_output = args.output_dir / "README.md"

    write_pdb_atoms(kept, pdb_output)
    write_xyz(xyz_output, heavy_deduped)

    metadata = {
        "source_pdb": str(args.input_pdb.relative_to(ROOT)),
        "derived_pdb": str(pdb_output.relative_to(ROOT)),
        "derived_xyz": str(xyz_output.relative_to(ROOT)),
        "atom_selection_rule": (
            "Remove atoms named CG, CD, OE1, OE2, HG2, and HG3 from every GLU residue. "
            "This removes the terminal GLU CH2-COOH group while preserving backbone atoms "
            "and CB. Hydrogens are excluded from the diffraction XYZ."
        ),
        "removed_glu_atom_names": sorted(REMOVED_GLU_ATOMS),
        "source_atom_count": len(atoms),
        "derived_pdb_atom_count": len(kept),
        "removed_atom_count": len(removed),
        "source_element_counts": element_counts(atoms),
        "removed_element_counts": element_counts(removed),
        "derived_pdb_element_counts": element_counts(kept),
        "derived_heavy_deduped_xyz_atom_count": len(heavy_deduped),
        "derived_heavy_deduped_xyz_element_counts": element_counts(heavy_deduped),
        "hydrogens_excluded_from_xyz": True,
        "exact_heavy_atom_deduplication_applied_to_xyz": True,
    }
    provenance_output.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    readme_output.write_text(
        "\n".join(
            [
                "# Nick Ideal 16-mer Derived Models",
                "",
                "`Hexaplex_AntiParallel_30deg_Ideal_no_CH2_COOH.pdb` is derived programmatically from "
                "`../Hexaplex_AntiParallel_30deg_Ideal.pdb`.",
                "",
                "Atom-selection rule: remove atoms named `CG`, `CD`, `OE1`, `OE2`, `HG2`, and `HG3` from every `GLU` residue. This removes the terminal GLU `CH2-COOH` group while preserving the peptide backbone and `CB` atom.",
                "",
                "The matching XYZ file uses the same corrected ideal-baseline convention: hydrogens excluded and exact heavy-atom deduplication applied.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(f"Wrote derived PDB: {pdb_output}")
    print(f"Wrote derived XYZ: {xyz_output}")
    print(f"Removed atoms: {len(removed)}")
    print(f"Derived heavy deduped XYZ atoms: {len(heavy_deduped)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
