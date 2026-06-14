#!/usr/bin/env python3
"""Compute simple z-axis residue centroid order metrics."""

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

from hexaplex_formation.geometry import group_atoms_by_residue, residue_centroid  # noqa: E402
from hexaplex_formation.pdb_utils import heavy_atoms, load_pdb_atoms  # noqa: E402


FIELDNAMES = [
    "residue",
    "chain",
    "residue_name",
    "residue_number",
    "centroid_x",
    "centroid_y",
    "centroid_z",
    "radius_xy",
    "angle_xy_rad",
]

SUMMARY_FIELDNAMES = [
    "pdb_file",
    "residue_count",
    "mean_radius_xy",
    "min_radius_xy",
    "max_radius_xy",
    "z_min",
    "z_max",
    "z_span",
    "angle_min",
    "angle_max",
    "angular_coverage_rad",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--residue-filter", default=None, help="Comma-separated residue names to include.")
    parser.add_argument("--heavy-only", action="store_true", help="Exclude hydrogens. This is the default.")
    parser.add_argument("--all-atoms", action="store_true", help="Include hydrogens.")
    args = parser.parse_args()
    if args.heavy_only and args.all_atoms:
        parser.error("--heavy-only and --all-atoms cannot be supplied together")
    return args


def _residue_filter(value: str | None) -> set[str] | None:
    if value is None:
        return None
    names = {part.strip().upper() for part in value.split(",") if part.strip()}
    return names or None


def residue_label(key: tuple[str, str, int | None, str]) -> str:
    chain_id, residue_name, residue_number, insertion_code = key
    number = "" if residue_number is None else str(residue_number)
    residue = f"{residue_name}{number}{insertion_code}"
    return f"{chain_id}:{residue}" if chain_id else residue


def helical_rows(atoms, use_heavy_only: bool = True, residue_filter: set[str] | None = None) -> list[dict[str, str]]:
    if use_heavy_only:
        atoms = heavy_atoms(atoms)
    if residue_filter is not None:
        atoms = [atom for atom in atoms if atom.residue_name.upper() in residue_filter]

    rows: list[dict[str, str]] = []
    for key, residue_atoms in group_atoms_by_residue(atoms).items():
        center = residue_centroid(residue_atoms)
        if center is None:
            continue
        x, y, z = center
        radius = math.sqrt(x * x + y * y)
        angle = math.atan2(y, x)
        chain_id, residue_name, residue_number, _icode = key
        rows.append(
            {
                "residue": residue_label(key),
                "chain": chain_id,
                "residue_name": residue_name,
                "residue_number": "" if residue_number is None else str(residue_number),
                "centroid_x": f"{x:.6f}",
                "centroid_y": f"{y:.6f}",
                "centroid_z": f"{z:.6f}",
                "radius_xy": f"{radius:.6f}",
                "angle_xy_rad": f"{angle:.6f}",
            }
        )
    return rows


def summary_row(pdb_path: Path, rows: list[dict[str, str]]) -> dict[str, str]:
    if not rows:
        return {
            "pdb_file": str(pdb_path),
            "residue_count": "0",
            "mean_radius_xy": "",
            "min_radius_xy": "",
            "max_radius_xy": "",
            "z_min": "",
            "z_max": "",
            "z_span": "",
            "angle_min": "",
            "angle_max": "",
            "angular_coverage_rad": "",
        }

    radii = [float(row["radius_xy"]) for row in rows]
    zs = [float(row["centroid_z"]) for row in rows]
    angles = [float(row["angle_xy_rad"]) for row in rows]
    radius_mean = sum(radii) / len(radii)
    z_min = min(zs)
    z_max = max(zs)
    angle_min = min(angles)
    angle_max = max(angles)
    return {
        "pdb_file": str(pdb_path),
        "residue_count": str(len(rows)),
        "mean_radius_xy": f"{radius_mean:.6f}",
        "min_radius_xy": f"{min(radii):.6f}",
        "max_radius_xy": f"{max(radii):.6f}",
        "z_min": f"{z_min:.6f}",
        "z_max": f"{z_max:.6f}",
        "z_span": f"{z_max - z_min:.6f}",
        "angle_min": f"{angle_min:.6f}",
        "angle_max": f"{angle_max:.6f}",
        "angular_coverage_rad": f"{angle_max - angle_min:.6f}",
    }


def main() -> int:
    args = parse_args()
    atoms = load_pdb_atoms(args.pdb)
    rows = helical_rows(atoms, use_heavy_only=not args.all_atoms, residue_filter=_residue_filter(args.residue_filter))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)

    summary_path = args.out.with_name(f"{args.out.stem}_summary.csv")
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDNAMES)
        writer.writeheader()
        writer.writerow(summary_row(args.pdb, rows))

    print(f"Wrote {args.out} and {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
