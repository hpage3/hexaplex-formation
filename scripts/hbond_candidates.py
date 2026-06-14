#!/usr/bin/env python3
"""Detect rough geometric hydrogen-bond candidates from all-atom PDBs."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from hexaplex_formation.geometry import distance, residue_key  # noqa: E402
from hexaplex_formation.pdb_utils import PDBAtom, is_hydrogen, load_pdb_atoms  # noqa: E402


FIELDNAMES = [
    "donor_residue",
    "donor_atom",
    "hydrogen_atom",
    "acceptor_residue",
    "acceptor_atom",
    "donor_acceptor_distance_A",
    "hydrogen_acceptor_distance_A",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--donor-acceptor-cutoff", type=float, default=3.5)
    parser.add_argument("--hydrogen-acceptor-cutoff", type=float, default=2.6)
    args = parser.parse_args()
    if args.donor_acceptor_cutoff <= 0 or args.hydrogen_acceptor_cutoff <= 0:
        parser.error("distance cutoffs must be greater than zero")
    return args


def residue_label(atom: PDBAtom) -> str:
    number = "" if atom.residue_number is None else str(atom.residue_number)
    residue = f"{atom.residue_name}{number}{atom.insertion_code}"
    return f"{atom.chain_id}:{residue}" if atom.chain_id else residue


def _starts_with_n_or_o(atom: PDBAtom) -> bool:
    name = atom.atom_name.strip().upper()
    element = atom.element.strip().upper()
    return element in {"N", "O"} or name.startswith(("N", "O"))


def donor_heavy_atoms(atoms: list[PDBAtom]) -> list[PDBAtom]:
    return [atom for atom in atoms if not is_hydrogen(atom) and _starts_with_n_or_o(atom)]


def acceptor_atoms(atoms: list[PDBAtom]) -> list[PDBAtom]:
    return [atom for atom in atoms if not is_hydrogen(atom) and _starts_with_n_or_o(atom)]


def infer_parent_donor(
    hydrogen: PDBAtom,
    candidate_donors: list[PDBAtom],
    max_distance: float = 1.2,
) -> PDBAtom | None:
    best: tuple[float, PDBAtom] | None = None
    h_key = residue_key(hydrogen)
    for donor in candidate_donors:
        if residue_key(donor) != h_key:
            continue
        donor_distance = distance(hydrogen, donor)
        if donor_distance <= max_distance and (best is None or donor_distance < best[0]):
            best = (donor_distance, donor)
    return None if best is None else best[1]


def hbond_candidate_rows(
    atoms: list[PDBAtom],
    donor_acceptor_cutoff: float = 3.5,
    hydrogen_acceptor_cutoff: float = 2.6,
) -> list[dict[str, str]]:
    donors = donor_heavy_atoms(atoms)
    acceptors = acceptor_atoms(atoms)
    hydrogens = [atom for atom in atoms if is_hydrogen(atom)]
    donor_hydrogen_pairs = [
        (donor, hydrogen)
        for hydrogen in hydrogens
        if (donor := infer_parent_donor(hydrogen, donors)) is not None
    ]

    rows: list[dict[str, str]] = []
    for donor, hydrogen in donor_hydrogen_pairs:
        donor_residue = residue_key(donor)
        for acceptor in acceptors:
            if residue_key(acceptor) == donor_residue:
                continue
            donor_acceptor_distance = distance(donor, acceptor)
            if donor_acceptor_distance > donor_acceptor_cutoff:
                continue
            hydrogen_acceptor_distance = distance(hydrogen, acceptor)
            if hydrogen_acceptor_distance > hydrogen_acceptor_cutoff:
                continue
            rows.append(
                {
                    "donor_residue": residue_label(donor),
                    "donor_atom": donor.atom_name,
                    "hydrogen_atom": hydrogen.atom_name,
                    "acceptor_residue": residue_label(acceptor),
                    "acceptor_atom": acceptor.atom_name,
                    "donor_acceptor_distance_A": f"{donor_acceptor_distance:.6f}",
                    "hydrogen_acceptor_distance_A": f"{hydrogen_acceptor_distance:.6f}",
                }
            )
    return rows


def main() -> int:
    args = parse_args()
    atoms = load_pdb_atoms(args.pdb)
    rows = hbond_candidate_rows(
        atoms,
        donor_acceptor_cutoff=args.donor_acceptor_cutoff,
        hydrogen_acceptor_cutoff=args.hydrogen_acceptor_cutoff,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {args.out} with {len(rows)} rough hydrogen-bond candidate(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
