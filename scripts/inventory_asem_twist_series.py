#!/usr/bin/env python3
"""Inventory Asem's imported 29/30/31 twist-series folders."""

from __future__ import annotations

import csv
import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from hexaplex_formation.pdb_utils import load_pdb_atoms, residue_count, residue_names  # noqa: E402

INPUT_ROOT = ROOT / "inputs" / "asem_twist_series_29_30_31" / "raw"
BASELINE = ROOT / "inputs" / "nick_ideal_models" / "Hexaplex_AntiParallel_30deg_Ideal.pdb"
OUTPUT_CSV = ROOT / "outputs" / "metrics" / "asem_twist_series_inventory.csv"

PDB_SUFFIXES = (".pdb", ".pdb.txt", ".ent")
XYZ_SUFFIXES = (".xyz",)
YAML_SUFFIXES = (".yaml", ".yml")
PY_SUFFIXES = (".py",)
AMBER_NAMES = {"initial.in", "tleap.in", "leap.in"}
AMBER_SUFFIXES = (".in", ".mol2", ".frcmod", ".lib", ".off", ".prep", ".prmtop", ".inpcrd", ".rst7", ".parm7")
EXPECTED_RESIDUES = {"CYP", "MEP", "GLU"}


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT)).replace("\\", "/")
    except ValueError:
        return str(path)


def classify_kind(path: Path) -> str:
    name = path.name.lower()
    suffix = path.suffix.lower()
    if name.endswith(PDB_SUFFIXES):
        return "pdb"
    if suffix in XYZ_SUFFIXES:
        return "xyz"
    if suffix in YAML_SUFFIXES:
        return "yaml"
    if suffix in PY_SUFFIXES:
        return "python"
    if name in AMBER_NAMES or name.startswith("leaprc") or suffix in AMBER_SUFFIXES:
        return "amber_tleap"
    return "other"


def xyz_atom_count(path: Path) -> int | None:
    try:
        first = path.read_text(encoding="utf-8", errors="replace").splitlines()[0].strip()
    except (IndexError, OSError):
        return None
    return int(first) if re.fullmatch(r"\d+", first) else None


def pdb_summary(path: Path) -> dict[str, object]:
    atoms = load_pdb_atoms(path)
    residues = set(residue_names(atoms))
    atom_count = len(atoms)
    baseline_atoms = len(load_pdb_atoms(BASELINE)) if BASELINE.exists() else 0
    scale_ratio = atom_count / baseline_atoms if baseline_atoms else 0.0
    residue_counter = Counter(atom.residue_name for atom in atoms if atom.residue_name)
    ready = (
        residue_count(atoms) >= 0.8 * 180
        and residues.issuperset(EXPECTED_RESIDUES)
        and atom_count >= 0.45 * baseline_atoms
    )
    note = "full diffraction-ready candidate" if ready else "local/intermediate or ambiguous scale"
    return {
        "atom_count": atom_count,
        "residue_count": residue_count(atoms),
        "residue_names": ";".join(sorted(residues)),
        "residue_name_counts": ";".join(f"{k}:{v}" for k, v in sorted(residue_counter.items())),
        "full_model_scale_ratio": f"{scale_ratio:.6g}",
        "diffraction_readiness": "diffraction_ready" if ready else "excluded_or_ambiguous",
        "readiness_note": note,
    }


def build_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for folder in sorted(INPUT_ROOT.iterdir() if INPUT_ROOT.exists() else []):
        if not folder.is_dir():
            continue
        for path in sorted(folder.rglob("*")):
            if not path.is_file():
                continue
            kind = classify_kind(path)
            row: dict[str, object] = {
                "source_folder": folder.name,
                "relative_path": display_path(path),
                "filename": path.name,
                "file_kind": kind,
                "size_bytes": path.stat().st_size,
                "atom_count": "",
                "residue_count": "",
                "residue_names": "",
                "residue_name_counts": "",
                "full_model_scale_ratio": "",
                "diffraction_readiness": "not_applicable",
                "readiness_note": "",
            }
            if kind == "pdb":
                row.update(pdb_summary(path))
            elif kind == "xyz":
                count = xyz_atom_count(path)
                row["atom_count"] = "" if count is None else count
                row["readiness_note"] = "XYZ atom count parsed from first line" if count is not None else "XYZ atom count unavailable"
            elif kind in {"yaml", "python", "amber_tleap"}:
                row["readiness_note"] = "generation/provenance input; documented but not executed"
            rows.append(row)
    return rows


def main() -> int:
    rows = build_rows()
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "source_folder",
        "relative_path",
        "filename",
        "file_kind",
        "size_bytes",
        "atom_count",
        "residue_count",
        "residue_names",
        "residue_name_counts",
        "full_model_scale_ratio",
        "diffraction_readiness",
        "readiness_note",
    ]
    with OUTPUT_CSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {OUTPUT_CSV} ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
