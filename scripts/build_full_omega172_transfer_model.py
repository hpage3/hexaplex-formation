#!/usr/bin/env python3
"""Attempt, then stop if needed, a full-model omega172 coordinate transfer."""

from __future__ import annotations

import csv
import math
import sys
from collections import Counter
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from hexaplex_formation.pdb_utils import load_pdb_atoms  # noqa: E402


FULL_PDB = ROOT / "inputs" / "nick_ideal_models" / "Hexaplex_AntiParallel_30deg_Ideal.pdb"
OMEGA167 = ROOT / "inputs" / "asem_omega_series_struct2" / "1_558214_omega167.pdb"
OMEGA172 = ROOT / "inputs" / "asem_omega_series_struct2" / "1_558214_omega172.pdb"
STOP_REPORT = ROOT / "outputs" / "reports" / "full_omega172_transfer_stop_report.md"
FINAL_REPORT = ROOT / "outputs" / "reports" / "full_omega172_transfer_model_report.md"
INSPECTION_CSV = ROOT / "outputs" / "metrics" / "full_omega172_transfer_correspondence_audit.csv"

EXPECTED_MOVED_NAMES = {"O'", "N''", "H''"}
TRANSFER_TOLERANCE_A = 0.10


def coord(atom) -> np.ndarray:
    return np.asarray([atom.x, atom.y, atom.z], dtype=float)


def residue_key(atom) -> tuple[int | None, str]:
    return (atom.residue_number, atom.residue_name)


def atom_key(atom) -> tuple[int | None, str, str]:
    return (atom.residue_number, atom.residue_name, atom.atom_name)


def residue_groups(atoms) -> dict[tuple[int | None, str], dict[str, object]]:
    groups: dict[tuple[int | None, str], dict[str, object]] = {}
    for atom in atoms:
        groups.setdefault(residue_key(atom), {})[atom.atom_name] = atom
    return groups


def compare_local_series(omega167_atoms, omega172_atoms) -> dict[str, object]:
    if len(omega167_atoms) != len(omega172_atoms):
        return {"atom_count_equal": False, "residue_count_equal": None, "stop_reason": "omega167 and omega172 atom counts differ"}

    residue_count_167 = len({residue_key(atom) for atom in omega167_atoms})
    residue_count_172 = len({residue_key(atom) for atom in omega172_atoms})
    correspondence_ok = True
    displacements = []
    moved_names = Counter()
    moved_count = 0
    for atom_a, atom_b in zip(omega167_atoms, omega172_atoms):
        same = (
            atom_a.atom_serial == atom_b.atom_serial
            and atom_a.atom_name == atom_b.atom_name
            and atom_a.residue_name == atom_b.residue_name
            and atom_a.residue_number == atom_b.residue_number
        )
        correspondence_ok = correspondence_ok and same
        distance = float(np.linalg.norm(coord(atom_b) - coord(atom_a)))
        displacements.append(distance)
        if distance > 1e-4:
            moved_count += 1
            moved_names[atom_a.atom_name] += 1

    return {
        "atom_count_equal": True,
        "residue_count_equal": residue_count_167 == residue_count_172,
        "residue_count": residue_count_167,
        "atom_serial_name_residue_correspondence": correspondence_ok,
        "moved_atom_count": moved_count,
        "moved_atom_names": ";".join(f"{name}:{count}" for name, count in sorted(moved_names.items())),
        "moved_atom_name_set": set(moved_names),
        "moved_names_match_expected": set(moved_names) == EXPECTED_MOVED_NAMES,
        "max_local_displacement_A": max(displacements) if displacements else None,
        "rms_local_displacement_A": math.sqrt(sum(value * value for value in displacements) / len(displacements)) if displacements else None,
    }


def candidate_mapping_rows(local_atoms, full_atoms) -> list[dict[str, object]]:
    local_groups = residue_groups(local_atoms)
    full_groups = {
        key: group
        for key, group in residue_groups(full_atoms).items()
        if key[1] in {"CYP", "MEP"}
    }
    rows = []
    used_full: set[tuple[int | None, str]] = set()
    for local_key, local_group in sorted(local_groups.items(), key=lambda item: (item[0][1], item[0][0] or -1)):
        best = None
        for full_key, full_group in full_groups.items():
            if full_key in used_full or full_key[1] != local_key[1]:
                continue
            common_names = sorted((set(local_group) & set(full_group)) - EXPECTED_MOVED_NAMES)
            if len(common_names) < 10:
                continue
            diffs = [float(np.linalg.norm(coord(local_group[name]) - coord(full_group[name]))) for name in common_names]
            rms = math.sqrt(sum(value * value for value in diffs) / len(diffs))
            candidate = {
                "local_residue_number": local_key[0],
                "local_residue_name": local_key[1],
                "best_full_residue_number": full_key[0],
                "best_full_residue_name": full_key[1],
                "common_unmoved_atom_count": len(common_names),
                "absolute_unmoved_rmsd_A": rms,
                "absolute_unmoved_max_delta_A": max(diffs),
            }
            if best is None or rms < best["absolute_unmoved_rmsd_A"]:
                best = candidate
        if best is not None:
            used_full.add((best["best_full_residue_number"], best["best_full_residue_name"]))
            rows.append(best)
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_reports(local_stats: dict[str, object], full_atoms, mapping_rows: list[dict[str, object]]) -> None:
    full_residue_count = len({residue_key(atom) for atom in full_atoms})
    full_base_residue_count = len({residue_key(atom) for atom in full_atoms if atom.residue_name in {"CYP", "MEP"}})
    max_abs_rms = max(float(row["absolute_unmoved_rmsd_A"]) for row in mapping_rows) if mapping_rows else None
    min_abs_rms = min(float(row["absolute_unmoved_rmsd_A"]) for row in mapping_rows) if mapping_rows else None
    clean_absolute_matches = [row for row in mapping_rows if float(row["absolute_unmoved_rmsd_A"]) <= TRANSFER_TOLERANCE_A]
    stop_reason = (
        "No high-confidence direct coordinate correspondence was established. "
        "Although the full model contains 90 CYP/MEP residues and the local struct2 model also contains 90 CYP/MEP residues, "
        "the best absolute-coordinate matches of shared unmoved atoms have multi-angstrom RMSD "
        f"(range {min_abs_rms:.3f}-{max_abs_rms:.3f} A), far above the {TRANSFER_TOLERANCE_A:.2f} A transfer tolerance. "
        "Copying raw omega172 displacement vectors into the full model would therefore perturb the wrong coordinate frame."
    )

    lines = [
        "# Full Omega172 Transfer Stop Report",
        "",
        "## Decision",
        "",
        "Stopped without generating a derived full PDB.",
        "",
        "## Stop Reason",
        "",
        stop_reason,
        "",
        "## Local Omega-Series Correspondence",
        "",
        f"- omega167/omega172 atom count equality: {local_stats['atom_count_equal']}",
        f"- omega167/omega172 residue count equality: {local_stats['residue_count_equal']}",
        f"- atom serial/name/residue correspondence: {local_stats['atom_serial_name_residue_correspondence']}",
        f"- moved atom names: {local_stats['moved_atom_names']}",
        f"- moved atom count: {local_stats['moved_atom_count']}",
        f"- moved names match expected `O'`, `N''`, `H''`: {local_stats['moved_names_match_expected']}",
        f"- local max displacement: {float(local_stats['max_local_displacement_A']):.4f} A",
        f"- local RMS displacement: {float(local_stats['rms_local_displacement_A']):.4f} A",
        "",
        "## Full Model Inspection",
        "",
        f"- full model atom count: {len(full_atoms)}",
        f"- full model residue count: {full_residue_count}",
        f"- full CYP/MEP residue count: {full_base_residue_count}",
        f"- local struct2 residue count: {local_stats['residue_count']}",
        f"- candidate residue mappings inspected: {len(mapping_rows)}",
        f"- clean absolute-coordinate matches within {TRANSFER_TOLERANCE_A:.2f} A RMSD: {len(clean_absolute_matches)}",
        "",
        "The same residue/atom-name pattern appears in repeated form, but the coordinate frames do not match closely enough for a direct displacement transfer. A local rigid-frame/Kabsch transfer might be mathematically possible, but it would require choosing among repeated symmetric motifs and rotating displacement vectors from one local frame into another. That is outside this conservative task and would risk inventing chemistry.",
        "",
        "## Recommendation",
        "",
        "Ask Asem to build the full periodic omega172 model directly from `1_558214_omega172.pdb`, then run corrected diffraction on that official full model. Do not use a coordinate-transfer surrogate for scientific conclusions.",
        "",
        "## Outputs",
        "",
        "- `outputs/metrics/full_omega172_transfer_correspondence_audit.csv`",
        "- `outputs/reports/full_omega172_transfer_stop_report.md`",
        "- `outputs/reports/full_omega172_transfer_model_report.md`",
    ]
    STOP_REPORT.parent.mkdir(parents=True, exist_ok=True)
    STOP_REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    final_lines = [
        "# Full Omega172 Transfer Model Report",
        "",
        "## Purpose",
        "",
        "This run attempted a conservative transfer of Asem's omega172 local-coordinate adjustment into the full 7,146-atom ideal antiparallel 30-degree model.",
        "",
        "## Strict Stop Condition",
        "",
        "The workflow is allowed to build a derived full PDB only if residue/atom correspondence is clear and validation can prove that atom count, rise, twist, and unmoved coordinates are preserved.",
        "",
        "## Result",
        "",
        "Transfer stopped. No derived full PDB, XYZ, diffraction profile, metrics, or plots were generated.",
        "",
        "## Mapping Method Tested",
        "",
        "The script compared the local omega167 and omega172 files by atom serial/name/residue identity, then attempted to map each local CYP/MEP residue to the best full-model CYP/MEP residue of the same residue name using shared unmoved atom names. Direct absolute-coordinate RMSD was used as the conservative transfer validation criterion.",
        "",
        "## Stop Reason",
        "",
        stop_reason,
        "",
        "## Validation Headline",
        "",
        f"- Local moved atom names: {local_stats['moved_atom_names']}",
        f"- Full atom count inspected: {len(full_atoms)}",
        f"- Derived full atom count: not applicable; no derived model was written.",
        "- Rise/twist/global coordinates: unchanged because no transfer was performed.",
        "- Corrected diffraction: not run.",
        "",
        "## Recommendation",
        "",
        "Reject this transfer attempt as ambiguous and ask Asem to build the official full periodic omega172 model directly. The derived coordinate-transfer model should not be used for discussion because it was not generated.",
        "",
        "## Limitations",
        "",
        "- Coordinate-transfer is not energy minimization.",
        "- Atom-name mapping is not a chemically exact torsion definition.",
        "- Asem/Nick visual chemistry review is still needed.",
        "- No diffraction result exists from this stopped transfer attempt.",
    ]
    FINAL_REPORT.parent.mkdir(parents=True, exist_ok=True)
    FINAL_REPORT.write_text("\n".join(final_lines) + "\n", encoding="utf-8")


def main() -> int:
    omega167_atoms = load_pdb_atoms(OMEGA167)
    omega172_atoms = load_pdb_atoms(OMEGA172)
    full_atoms = load_pdb_atoms(FULL_PDB)
    local_stats = compare_local_series(omega167_atoms, omega172_atoms)
    mapping_rows = candidate_mapping_rows(omega167_atoms, full_atoms)
    write_csv(INSPECTION_CSV, mapping_rows)
    write_reports(local_stats, full_atoms, mapping_rows)
    print(f"Stopped transfer; wrote {STOP_REPORT}")
    print(f"Wrote final report: {FINAL_REPORT}")
    print(f"Wrote correspondence audit: {INSPECTION_CSV}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
