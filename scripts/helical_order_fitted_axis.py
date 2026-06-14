#!/usr/bin/env python3
"""Compute residue helical order metrics around a fitted principal axis."""

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

from hexaplex_formation.geometry import (  # noqa: E402
    build_perpendicular_basis,
    covariance_matrix_3d,
    dot,
    mean_point,
    norm,
    power_iteration_principal_axis,
    project_point_to_axis,
    residue_centroid,
    vector_sub,
)
from hexaplex_formation.geometry import group_atoms_by_residue  # noqa: E402
from hexaplex_formation.pdb_utils import heavy_atoms, load_pdb_atoms  # noqa: E402


FIELDNAMES = [
    "residue",
    "chain",
    "residue_name",
    "residue_number",
    "centroid_x",
    "centroid_y",
    "centroid_z",
    "axial_t",
    "radius_fitted",
    "angle_fitted_rad",
]

SUMMARY_FIELDNAMES = [
    "pdb_file",
    "residue_count",
    "axis_origin_x",
    "axis_origin_y",
    "axis_origin_z",
    "axis_x",
    "axis_y",
    "axis_z",
    "mean_radius_fitted",
    "min_radius_fitted",
    "max_radius_fitted",
    "axial_t_min",
    "axial_t_max",
    "axial_span",
    "angle_min",
    "angle_max",
    "angular_coverage_rad",
    "z_axis_angle_degrees",
    "unwrapped_angle_span_rad",
    "approximate_turns",
    "approximate_pitch_per_turn",
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


def unwrap_angles(angles: list[float]) -> list[float]:
    if not angles:
        return []
    unwrapped = [angles[0]]
    two_pi = 2.0 * math.pi
    for angle in angles[1:]:
        adjusted = angle
        previous = unwrapped[-1]
        while adjusted - previous > math.pi:
            adjusted -= two_pi
        while adjusted - previous < -math.pi:
            adjusted += two_pi
        unwrapped.append(adjusted)
    return unwrapped


def _empty_summary(pdb_path: Path) -> dict[str, str]:
    row = {field: "" for field in SUMMARY_FIELDNAMES}
    row["pdb_file"] = str(pdb_path)
    row["residue_count"] = "0"
    return row


def _format_point(prefix: str, point: tuple[float, float, float], row: dict[str, str]) -> None:
    row[f"{prefix}_x"] = f"{point[0]:.6f}"
    row[f"{prefix}_y"] = f"{point[1]:.6f}"
    row[f"{prefix}_z"] = f"{point[2]:.6f}"


def fitted_helical_rows(
    atoms,
    use_heavy_only: bool = True,
    residue_filter: set[str] | None = None,
) -> tuple[list[dict[str, str]], dict[str, tuple[float, float, float] | list[float]]]:
    if use_heavy_only:
        atoms = heavy_atoms(atoms)
    if residue_filter is not None:
        atoms = [atom for atom in atoms if atom.residue_name.upper() in residue_filter]

    residue_items = []
    centroids = []
    for key, residue_atoms in group_atoms_by_residue(atoms).items():
        center = residue_centroid(residue_atoms)
        if center is None:
            continue
        residue_items.append((key, center))
        centroids.append(center)

    if not centroids:
        return [], {}

    origin = mean_point(centroids)
    axis = power_iteration_principal_axis(covariance_matrix_3d(centroids))
    basis_u, basis_v = build_perpendicular_basis(axis)

    rows: list[dict[str, str]] = []
    for key, center in residue_items:
        axial_t, projected = project_point_to_axis(center, origin, axis)
        radial = vector_sub(center, projected)
        radius = norm(radial)
        angle = math.atan2(dot(radial, basis_v), dot(radial, basis_u))
        chain_id, residue_name, residue_number, _icode = key
        rows.append(
            {
                "residue": residue_label(key),
                "chain": chain_id,
                "residue_name": residue_name,
                "residue_number": "" if residue_number is None else str(residue_number),
                "centroid_x": f"{center[0]:.6f}",
                "centroid_y": f"{center[1]:.6f}",
                "centroid_z": f"{center[2]:.6f}",
                "axial_t": f"{axial_t:.6f}",
                "radius_fitted": f"{radius:.6f}",
                "angle_fitted_rad": f"{angle:.6f}",
            }
        )

    return rows, {"origin": origin, "axis": axis, "angles": [float(row["angle_fitted_rad"]) for row in rows]}


def summary_row(pdb_path: Path, rows: list[dict[str, str]], fit: dict[str, tuple[float, float, float] | list[float]]) -> dict[str, str]:
    if not rows:
        return _empty_summary(pdb_path)

    origin = fit["origin"]
    axis = fit["axis"]
    assert isinstance(origin, tuple)
    assert isinstance(axis, tuple)

    radii = [float(row["radius_fitted"]) for row in rows]
    axial_ts = [float(row["axial_t"]) for row in rows]
    angles = [float(row["angle_fitted_rad"]) for row in rows]
    unwrapped = unwrap_angles(angles)
    axial_min = min(axial_ts)
    axial_max = max(axial_ts)
    axial_span = axial_max - axial_min
    unwrapped_span = unwrapped[-1] - unwrapped[0] if len(unwrapped) > 1 else 0.0
    approximate_turns = abs(unwrapped_span) / (2.0 * math.pi)
    z_axis_angle = math.degrees(math.acos(max(-1.0, min(1.0, abs(axis[2])))))

    row = {
        "pdb_file": str(pdb_path),
        "residue_count": str(len(rows)),
        "axis_x": f"{axis[0]:.6f}",
        "axis_y": f"{axis[1]:.6f}",
        "axis_z": f"{axis[2]:.6f}",
        "mean_radius_fitted": f"{sum(radii) / len(radii):.6f}",
        "min_radius_fitted": f"{min(radii):.6f}",
        "max_radius_fitted": f"{max(radii):.6f}",
        "axial_t_min": f"{axial_min:.6f}",
        "axial_t_max": f"{axial_max:.6f}",
        "axial_span": f"{axial_span:.6f}",
        "angle_min": f"{min(angles):.6f}",
        "angle_max": f"{max(angles):.6f}",
        "angular_coverage_rad": f"{max(angles) - min(angles):.6f}",
        "z_axis_angle_degrees": f"{z_axis_angle:.6f}",
        "unwrapped_angle_span_rad": f"{unwrapped_span:.6f}",
        "approximate_turns": f"{approximate_turns:.6f}",
        "approximate_pitch_per_turn": f"{axial_span / approximate_turns:.6f}" if approximate_turns > 0 else "",
    }
    _format_point("axis_origin", origin, row)
    return row


def main() -> int:
    args = parse_args()
    rows, fit = fitted_helical_rows(
        load_pdb_atoms(args.pdb),
        use_heavy_only=not args.all_atoms,
        residue_filter=_residue_filter(args.residue_filter),
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    summary_path = args.out.with_name(f"{args.out.stem}_summary.csv")
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDNAMES, lineterminator="\n")
        writer.writeheader()
        writer.writerow(summary_row(args.pdb, rows, fit))

    print(f"Wrote {args.out} and {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
