#!/usr/bin/env python3
"""Create heavy-atom XYZ files with the distal GLU side chain removed."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


# Preserve the peptide backbone and the GLU recognition-core atoms through CB.
# Remove only the distal glutamate branch and its attached hydrogens.
REMOVED_GLU_ATOMS = {"CG", "HG2", "HG3", "CD", "OE1", "OE2"}
STRIPPING_RULE = "GLU atoms CG,HG2,HG3,CD,OE1,OE2 removed; all other residues/atoms preserved"

MANIFEST_FIELDS = [
    "model_id",
    "twist_deg",
    "rise_A",
    "source_pdb",
    "output_xyz",
    "source_atom_count",
    "source_heavy_atom_count",
    "removed_atom_count",
    "removed_heavy_atom_count",
    "skipped_retained_hydrogen_count",
    "duplicate_heavy_atom_count",
    "retained_atom_count",
    "stripping_rule",
    "status",
    "notes",
]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def infer_element(line: str, atom_name: str) -> str:
    element = line[76:78].strip() if len(line) >= 78 else ""
    if element:
        return element[0].upper() + element[1:].lower()
    letters = "".join(char for char in atom_name if char.isalpha())
    if not letters:
        raise ValueError(f"Cannot infer element from atom name {atom_name!r}")
    return letters[0].upper()


def parse_pdb_for_ablation(path: Path) -> tuple[list[tuple[str, float, float, float]], dict[str, int]]:
    xyz_atoms: list[tuple[str, float, float, float]] = []
    seen: set[tuple[str, float, float, float]] = set()
    counts = {
        "source_atom_count": 0,
        "source_heavy_atom_count": 0,
        "removed_atom_count": 0,
        "removed_heavy_atom_count": 0,
        "skipped_retained_hydrogen_count": 0,
        "duplicate_heavy_atom_count": 0,
    }

    for line_number, line in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), start=1):
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if len(line) < 54:
            raise ValueError(f"{path} line {line_number} is too short for PDB coordinates")

        counts["source_atom_count"] += 1
        atom_name = line[12:16].strip()
        residue_name = line[17:20].strip()
        element = infer_element(line, atom_name)
        is_hydrogen = element.upper() == "H"
        if not is_hydrogen:
            counts["source_heavy_atom_count"] += 1

        if residue_name == "GLU" and atom_name in REMOVED_GLU_ATOMS:
            counts["removed_atom_count"] += 1
            if not is_hydrogen:
                counts["removed_heavy_atom_count"] += 1
            continue
        if is_hydrogen:
            counts["skipped_retained_hydrogen_count"] += 1
            continue

        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError as exc:
            raise ValueError(f"{path} line {line_number} has malformed coordinates") from exc

        atom = (element, x, y, z)
        if atom in seen:
            counts["duplicate_heavy_atom_count"] += 1
            continue
        seen.add(atom)
        xyz_atoms.append(atom)

    if counts["source_atom_count"] == 0:
        raise ValueError(f"No ATOM/HETATM records found in {path}")
    if counts["removed_atom_count"] == 0:
        raise ValueError(f"No removable GLU side-chain atoms found in {path}")
    if not xyz_atoms:
        raise ValueError(f"No retained heavy atoms found in {path}")
    return xyz_atoms, counts


def write_xyz(path: Path, atoms: list[tuple[str, float, float, float]], comment: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [str(len(atoms)), comment]
    lines.extend(f"{element} {x:.6f} {y:.6f} {z:.6f}" for element, x, y, z in atoms)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def process_manifest(manifest: Path, output_dir: Path, output_manifest: Path) -> list[dict[str, object]]:
    source_rows = read_csv(manifest)
    if not source_rows:
        raise ValueError(f"No rows found in {manifest}")
    required = {"model_id", "twist_deg", "rise_A", "source_pdb"}
    missing = required - set(source_rows[0])
    if missing:
        raise ValueError(f"{manifest} is missing required columns: {sorted(missing)}")

    rows: list[dict[str, object]] = []
    for source in source_rows:
        model_id = source["model_id"].strip()
        source_pdb = Path(source["source_pdb"])
        output_xyz = output_dir / f"{model_id}.xyz"
        row: dict[str, object] = {
            "model_id": model_id,
            "twist_deg": source["twist_deg"],
            "rise_A": source["rise_A"],
            "source_pdb": str(source_pdb).replace("\\", "/"),
            "output_xyz": str(output_xyz).replace("\\", "/"),
            "source_atom_count": "",
            "source_heavy_atom_count": "",
            "removed_atom_count": "",
            "removed_heavy_atom_count": "",
            "skipped_retained_hydrogen_count": "",
            "duplicate_heavy_atom_count": "",
            "retained_atom_count": "",
            "stripping_rule": STRIPPING_RULE,
            "status": "pending",
            "notes": "",
        }
        try:
            atoms, counts = parse_pdb_for_ablation(source_pdb)
            write_xyz(
                output_xyz,
                atoms,
                f"{model_id}; {STRIPPING_RULE}; hydrogens excluded; exact duplicate heavy atoms removed",
            )
            row.update(counts)
            row["retained_atom_count"] = len(atoms)
            row["status"] = "complete"
            row["notes"] = (
                "retained_atom_count is the heavy-atom XYZ count after stripping, "
                "hydrogen exclusion, and exact-coordinate deduplication"
            )
        except Exception as exc:
            row["status"] = "failed"
            row["notes"] = f"{type(exc).__name__}: {exc}"
        rows.append(row)

    write_csv(output_manifest, rows)
    return rows


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path(
            "outputs/asem_sidechains_20260625/stage2_top5_refined/scoring/"
            "sidechain_stage2_top5_scored_manifest_ranked.csv"
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs/asem_sidechains_20260625/stage2_top5_no_glu_ablation/xyz"),
    )
    parser.add_argument(
        "--output-manifest",
        type=Path,
        default=Path(
            "outputs/asem_sidechains_20260625/stage2_top5_no_glu_ablation/"
            "sidechain_ablation_xyz_manifest.csv"
        ),
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    rows = process_manifest(args.manifest, args.output_dir, args.output_manifest)
    complete = sum(row["status"] == "complete" for row in rows)
    print(f"Wrote {complete}/{len(rows)} ablated XYZ files and {args.output_manifest}")
    return 1 if complete != len(rows) else 0


if __name__ == "__main__":
    raise SystemExit(main())
