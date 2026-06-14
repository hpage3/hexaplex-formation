#!/usr/bin/env python3
"""Decompose residue contacts by scaffold block and component class."""

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
from hexaplex_formation.pdb_utils import PDBAtom, heavy_atoms, load_pdb_atoms  # noqa: E402
from hexaplex_formation.strand_map import (  # noqa: E402
    StrandMap,
    block_for_atom,
    infer_component,
    load_strand_map,
    residue_label_from_atom,
)


FIELDNAMES = [
    "residue_i",
    "residue_j",
    "residue_name_i",
    "residue_name_j",
    "residue_number_i",
    "residue_number_j",
    "component_i",
    "component_j",
    "block_i",
    "block_j",
    "contact_category",
    "min_distance_A",
    "atom_i",
    "atom_j",
    "is_GLU_involved",
    "is_GLU_GLU",
    "is_backbone_O_to_GLU_sidechain_O",
    "is_GLU_sidechain_O_to_GLU_sidechain_O",
    "is_backbone_N_to_GLU_sidechain_O",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb", required=True, type=Path)
    parser.add_argument("--strand-map", type=Path, default=Path("inputs/metadata/strand_map_candidate.csv"))
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--cutoff", type=float, default=4.5)
    parser.add_argument("--heavy-only", action="store_true", default=True, help="Exclude hydrogens. This is the default.")
    parser.add_argument("--all-atoms", action="store_true", help="Include hydrogens.")
    args = parser.parse_args()
    if args.cutoff <= 0:
        parser.error("--cutoff must be greater than zero")
    return args


def classify_contact(component_i: str, component_j: str, block_i: str | None, block_j: str | None) -> str:
    if component_i == "scaffold" and component_j == "scaffold":
        if block_i == block_j:
            return "scaffold_within_block"
        return "scaffold_between_blocks"
    if component_i != "scaffold" and component_j != "scaffold":
        return "hexad_or_other_internal"
    return "scaffold_hexad_or_other"


def _is_glu_sidechain_o(residue_name: str, atom_name: str) -> bool:
    return residue_name.upper() == "GLU" and atom_name.upper() in {"OE1", "OE2"}


def motif_flags(atom_i: PDBAtom, atom_j: PDBAtom) -> dict[str, str]:
    res_i = atom_i.residue_name.upper()
    res_j = atom_j.residue_name.upper()
    atom_name_i = atom_i.atom_name.upper()
    atom_name_j = atom_j.atom_name.upper()
    i_side_o = _is_glu_sidechain_o(res_i, atom_name_i)
    j_side_o = _is_glu_sidechain_o(res_j, atom_name_j)
    return {
        "is_GLU_involved": "yes" if res_i == "GLU" or res_j == "GLU" else "no",
        "is_GLU_GLU": "yes" if res_i == "GLU" and res_j == "GLU" else "no",
        "is_backbone_O_to_GLU_sidechain_O": "yes"
        if (atom_name_i == "O" and j_side_o) or (atom_name_j == "O" and i_side_o)
        else "no",
        "is_GLU_sidechain_O_to_GLU_sidechain_O": "yes" if i_side_o and j_side_o else "no",
        "is_backbone_N_to_GLU_sidechain_O": "yes"
        if (atom_name_i == "N" and j_side_o) or (atom_name_j == "N" and i_side_o)
        else "no",
    }


def decomposition_rows(
    atoms: list[PDBAtom],
    strand_map: StrandMap,
    cutoff: float = 4.5,
    use_heavy_only: bool = True,
) -> list[dict[str, str]]:
    if use_heavy_only:
        atoms = heavy_atoms(atoms)

    residues = list(group_atoms_by_residue(atoms).items())
    rows: list[dict[str, str]] = []
    for index, (_key_i, atoms_i) in enumerate(residues):
        representative_i = atoms_i[0]
        for _key_j, atoms_j in residues[index + 1 :]:
            representative_j = atoms_j[0]
            minimum = min_inter_residue_distance(atoms_i, atoms_j, heavy_only=False)
            if minimum is None or minimum[0] > cutoff:
                continue

            min_distance, atom_i, atom_j = minimum
            component_i = infer_component(atom_i, strand_map)
            component_j = infer_component(atom_j, strand_map)
            block_i = block_for_atom(atom_i, strand_map) if component_i == "scaffold" else None
            block_j = block_for_atom(atom_j, strand_map) if component_j == "scaffold" else None
            row = {
                "residue_i": residue_label_from_atom(representative_i),
                "residue_j": residue_label_from_atom(representative_j),
                "residue_name_i": representative_i.residue_name,
                "residue_name_j": representative_j.residue_name,
                "residue_number_i": "" if representative_i.residue_number is None else str(representative_i.residue_number),
                "residue_number_j": "" if representative_j.residue_number is None else str(representative_j.residue_number),
                "component_i": component_i,
                "component_j": component_j,
                "block_i": block_i or "",
                "block_j": block_j or "",
                "contact_category": classify_contact(component_i, component_j, block_i, block_j),
                "min_distance_A": f"{min_distance:.6f}",
                "atom_i": atom_i.atom_name,
                "atom_j": atom_j.atom_name,
            }
            row.update(motif_flags(atom_i, atom_j))
            rows.append(row)
    return rows


def write_rows(rows: list[dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    rows = decomposition_rows(
        load_pdb_atoms(args.pdb),
        load_strand_map(args.strand_map),
        cutoff=args.cutoff,
        use_heavy_only=not args.all_atoms,
    )
    write_rows(rows, args.out)
    print(f"Wrote {args.out} with {len(rows)} block contact row(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
