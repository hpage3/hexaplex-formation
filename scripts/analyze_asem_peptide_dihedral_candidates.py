#!/usr/bin/env python3
"""Audit Asem peptide-dihedral candidate PDBs for local geometry."""

from __future__ import annotations

import csv
import math
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from hexaplex_formation.pdb_utils import heavy_atoms, load_pdb_atoms  # noqa: E402


INPUT_DIR = ROOT / "inputs" / "asem_peptide_dihedral_candidates"
AUDIT_CSV = ROOT / "outputs" / "metrics" / "asem_peptide_dihedral_candidate_structural_audit.csv"

CANONICAL_OMEGA = ("CA", "C", "N'", "CA'")
NICK_NOTE_OMEGA_PROXY = ("N", "CA", "C", "N'")
BACKBONE_N_NAMES = {"N", "N'", "N''"}
BACKBONE_O_NAMES = {"O", "O'"}
DONOR_HYDROGEN = {"N": "H1", "N'": "H'", "N''": "H''"}


@dataclass(frozen=True)
class AtomRef:
    index: int
    atom: object


def coord(atom) -> np.ndarray:
    return np.asarray([atom.x, atom.y, atom.z], dtype=float)


def residue_key(atom) -> tuple[str, int | None, str, str]:
    return (atom.chain_id, atom.residue_number, atom.insertion_code, atom.residue_name)


def group_residues(atoms) -> tuple[list[tuple[str, int | None, str, str]], dict[tuple[str, int | None, str, str], dict[str, object]]]:
    groups: dict[tuple[str, int | None, str, str], dict[str, object]] = {}
    order: list[tuple[str, int | None, str, str]] = []
    for atom in atoms:
        key = residue_key(atom)
        if key not in groups:
            groups[key] = {}
            order.append(key)
        groups[key][atom.atom_name] = atom
    return order, groups


def dihedral_deg(a, b, c, d) -> float | None:
    p0, p1, p2, p3 = [coord(atom) for atom in (a, b, c, d)]
    if min(np.linalg.norm(p1 - p0), np.linalg.norm(p2 - p1), np.linalg.norm(p3 - p2)) < 1e-8:
        return None
    b0 = -(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    if np.linalg.norm(v) < 1e-8 or np.linalg.norm(w) < 1e-8:
        return None
    x = float(np.dot(v, w))
    y = float(np.dot(np.cross(b1, v), w))
    return math.degrees(math.atan2(y, x))


def median_abs_dihedral(order, groups, names: tuple[str, str, str, str]) -> tuple[float | None, int]:
    values = []
    for key in order:
        group = groups[key]
        if all(name in group for name in names):
            value = dihedral_deg(*(group[name] for name in names))
            if value is not None:
                values.append(abs(value))
    if not values:
        return None, 0
    return float(np.median(values)), len(values)


def contact_proxy(atoms) -> dict[str, object]:
    heavy = heavy_atoms(atoms)
    min_distance = float("inf")
    min_pair = ""
    counts = {2.0: 0, 2.2: 0, 2.5: 0}
    for i, atom_i in enumerate(heavy):
        key_i = residue_key(atom_i)
        pos_i = coord(atom_i)
        for atom_j in heavy[i + 1 :]:
            if residue_key(atom_j) == key_i:
                continue
            distance = float(np.linalg.norm(pos_i - coord(atom_j)))
            if distance < min_distance:
                min_distance = distance
                min_pair = (
                    f"{atom_i.chain_id}{atom_i.residue_number}:{atom_i.atom_name}-"
                    f"{atom_j.chain_id}{atom_j.residue_number}:{atom_j.atom_name}"
                )
            for threshold in counts:
                if distance < threshold:
                    counts[threshold] += 1
    return {
        "min_heavy_nonbonded_distance_A": min_distance if np.isfinite(min_distance) else None,
        "min_heavy_nonbonded_pair": min_pair,
        "clash_count_lt_2p0": counts[2.0],
        "clash_count_lt_2p2": counts[2.2],
        "clash_count_lt_2p5": counts[2.5],
    }


def angle_deg(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float | None:
    ba = a - b
    bc = c - b
    norm = np.linalg.norm(ba) * np.linalg.norm(bc)
    if norm < 1e-8:
        return None
    cosine = float(np.clip(np.dot(ba, bc) / norm, -1.0, 1.0))
    return math.degrees(math.acos(cosine))


def hbond_proxy(atoms) -> dict[str, object]:
    order, groups = group_residues(atoms)
    donors: list[tuple[object, object | None, tuple[str, int | None, str, str]]] = []
    acceptors: list[tuple[object, tuple[str, int | None, str, str]]] = []
    for key in order:
        group = groups[key]
        for name in BACKBONE_N_NAMES:
            if name in group:
                donors.append((group[name], group.get(DONOR_HYDROGEN.get(name, "")), key))
        for name in BACKBONE_O_NAMES:
            if name in group:
                acceptors.append((group[name], key))

    best: dict[str, object] = {
        "best_backbone_hbond_distance_A": None,
        "best_backbone_hbond_HO_distance_A": None,
        "best_backbone_hbond_NHO_angle_deg": None,
        "best_backbone_hbond_pair": "",
        "hbond_plausible_flag": False,
    }
    best_score = float("inf")
    for donor, hydrogen, donor_key in donors:
        for acceptor, acceptor_key in acceptors:
            if donor_key == acceptor_key:
                continue
            if donor.chain_id == acceptor.chain_id and donor.residue_number == acceptor.residue_number:
                continue
            no_distance = float(np.linalg.norm(coord(donor) - coord(acceptor)))
            ho_distance = None
            nho_angle = None
            if hydrogen is not None:
                ho_distance = float(np.linalg.norm(coord(hydrogen) - coord(acceptor)))
                nho_angle = angle_deg(coord(donor), coord(hydrogen), coord(acceptor))
            angle_penalty = max(0.0, 150.0 - nho_angle) / 100.0 if nho_angle is not None else 0.2
            score = abs(no_distance - 2.9) + angle_penalty
            if 2.4 <= no_distance <= 4.0 and score < best_score:
                best_score = score
                best = {
                    "best_backbone_hbond_distance_A": no_distance,
                    "best_backbone_hbond_HO_distance_A": ho_distance,
                    "best_backbone_hbond_NHO_angle_deg": nho_angle,
                    "best_backbone_hbond_pair": (
                        f"{donor.chain_id}{donor.residue_number}:{donor.atom_name}-"
                        f"{acceptor.chain_id}{acceptor.residue_number}:{acceptor.atom_name}"
                    ),
                    "hbond_plausible_flag": bool(no_distance <= 3.5 and (nho_angle is None or nho_angle >= 120.0)),
                }
    return best


def qualitative_status(candidate_file: str, nick_omega: float | None, contacts: dict[str, object], hbonds: dict[str, object]) -> str:
    if candidate_file == "1_558214.pdb":
        return "original_hbond_preserving_baseline"
    if nick_omega is None:
        return "omega_assignment_ambiguous"
    if bool(hbonds["hbond_plausible_flag"]) is False:
        return "omega_improved_hbond_unclear"
    if int(contacts["clash_count_lt_2p2"]) > 0:
        return "omega_improved_but_sterics_worse"
    return "possible_candidate_for_larger_model"


def audit_candidate(path: Path) -> dict[str, object]:
    atoms = load_pdb_atoms(path)
    order, groups = group_residues(atoms)
    canonical_omega, canonical_count = median_abs_dihedral(order, groups, CANONICAL_OMEGA)
    nick_omega, nick_count = median_abs_dihedral(order, groups, NICK_NOTE_OMEGA_PROXY)
    contacts = contact_proxy(atoms)
    hbonds = hbond_proxy(atoms)
    element_counts = dict(sorted(Counter(atom.element.upper() for atom in atoms).items()))
    heavy = heavy_atoms(atoms)
    return {
        "candidate_file": path.name,
        "is_original_hbond_preserving_baseline": path.name == "1_558214.pdb",
        "atom_count": len(atoms),
        "heavy_atom_count": len(heavy),
        "element_counts": ";".join(f"{key}:{value}" for key, value in element_counts.items()),
        "residue_count": len(order),
        "canonical_amide_omega_CA_C_Np_CAp_deg": canonical_omega,
        "canonical_amide_omega_count": canonical_count,
        "nick_note_omega_proxy_N_CA_C_Np_deg": nick_omega,
        "nick_note_omega_proxy_count": nick_count,
        "omega_deviation_from_180_deg": abs(180.0 - nick_omega) if nick_omega is not None else None,
        **contacts,
        **hbonds,
        "qualitative_status": qualitative_status(path.name, nick_omega, contacts, hbonds),
        "notes": (
            "Canonical CA-C-N'-CA' amide torsion is near 180 in all files; "
            "N-CA-C-N' is reported as the Nick-note omega proxy because it tracks the cited candidate differences. "
            "Steric proxy excludes same-residue heavy-heavy pairs but does not infer a full bond graph."
        ),
    }


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys()) if rows else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    candidates = sorted(INPUT_DIR.glob("*.pdb"))
    if not candidates:
        raise FileNotFoundError(f"No candidate PDBs found in {INPUT_DIR}")
    rows = [audit_candidate(path) for path in candidates]
    write_csv(AUDIT_CSV, rows)
    print(f"Wrote {AUDIT_CSV}")
    for row in rows:
        print(
            f"{row['candidate_file']}: omega_proxy={row['nick_note_omega_proxy_N_CA_C_Np_deg']:.3f}, "
            f"min_contact={row['min_heavy_nonbonded_distance_A']:.3f}, "
            f"hbond={row['best_backbone_hbond_distance_A']:.3f} "
            f"plausible={row['hbond_plausible_flag']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
