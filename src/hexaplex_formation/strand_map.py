"""Candidate scaffold strand/block map helpers."""

from __future__ import annotations

import csv
from pathlib import Path

from hexaplex_formation.pdb_utils import PDBAtom


ResidueIdentity = tuple[str, str, int | None, str]
StrandMap = dict[ResidueIdentity, str]

HEXAD_OR_OTHER_ATOM_NAMES = {"N1", "C2", "N3", "C4", "N5", "C6", "OC2", "OC4", "OC6"}


def residue_identity_from_atom(atom: PDBAtom) -> ResidueIdentity:
    return atom.chain_id, atom.residue_name, atom.residue_number, atom.insertion_code


def _residue_identity_from_row(row: dict[str, str]) -> ResidueIdentity:
    residue_number = row["residue_number"].strip()
    return (
        row["chain_id"],
        row["residue_name"],
        int(residue_number) if residue_number else None,
        row["insertion_code"],
    )


def load_strand_map(path: str | Path) -> StrandMap:
    mapping: StrandMap = {}
    with Path(path).open("r", newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            mapping[_residue_identity_from_row(row)] = row["block_id"]
    return mapping


def block_for_atom(atom: PDBAtom, strand_map: StrandMap) -> str | None:
    return strand_map.get(residue_identity_from_atom(atom))


def residue_label_from_atom(atom: PDBAtom) -> str:
    number = "" if atom.residue_number is None else str(atom.residue_number)
    residue = f"{atom.residue_name}{number}{atom.insertion_code}"
    return f"{atom.chain_id}:{residue}" if atom.chain_id else residue


def infer_component(atom: PDBAtom, strand_map: StrandMap) -> str:
    if atom.atom_name.upper() in HEXAD_OR_OTHER_ATOM_NAMES:
        return "hexad_or_other"
    return "scaffold" if block_for_atom(atom, strand_map) is not None else "hexad_or_other"
