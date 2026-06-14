#!/usr/bin/env python3
"""Generate local CYP/MEP base-length variants without running diffraction."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import defaultdict
from dataclasses import replace
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from hexaplex_formation.geometry import (  # noqa: E402
    covariance_matrix_3d,
    distance,
    group_atoms_by_residue,
    mean_point,
    power_iteration_principal_axis,
    project_point_to_axis,
)
from hexaplex_formation.pdb_utils import PDBAtom, heavy_atoms, load_pdb_atoms, write_pdb_atoms  # noqa: E402


BACKBONE_ATOM_NAMES = {"N", "CA", "C", "O", "OXT"}
TRANSFORMABLE_RESIDUES = {"CYP", "MEP"}
DEFAULT_SCALE_FACTORS = [0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20]

MANIFEST_FIELDNAMES = [
    "variant_id",
    "scale_factor",
    "input_pdb",
    "output_pdb",
    "total_atom_count",
    "transformed_atom_count",
    "fixed_atom_count",
    "transformed_residue_types",
    "transformed_atom_names",
    "anchor_rule",
    "transformation_rule",
    "notes",
    "warnings",
]

GEOMETRY_FIELDNAMES = [
    "variant_id",
    "scale_factor",
    "total_atom_count",
    "transformed_atom_count",
    "fixed_atom_count",
    "chain_count",
    "residue_count",
    "min_heavy_atom_distance_A",
    "suspicious_overlap_count",
    "radial_extent_min_A",
    "radial_extent_max_A",
    "radial_extent_mean_A",
    "mean_transformed_atom_displacement_A",
    "max_transformed_atom_displacement_A",
    "notes",
    "warnings",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pdb",
        type=Path,
        default=Path("outputs/intermediates/ai_candidate_inputs/full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain.pdb"),
    )
    parser.add_argument(
        "--selection-summary",
        type=Path,
        default=Path("outputs/metrics/hexad_arm_geometry_summary.csv"),
    )
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/base_length_sweep/structures"))
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("outputs/base_length_sweep/structures/base_length_variant_manifest.csv"),
    )
    parser.add_argument(
        "--geometry-out",
        type=Path,
        default=Path("outputs/metrics/base_length_variant_geometry.csv"),
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=Path("outputs/reports/base_length_variant_generation_report.md"),
    )
    parser.add_argument(
        "--scale-factors",
        default=",".join(f"{value:.2f}" for value in DEFAULT_SCALE_FACTORS),
        help="Comma-separated scale factors.",
    )
    parser.add_argument("--overlap-threshold", type=float, default=1.0)
    return parser.parse_args()


def parse_scale_factors(value: str) -> list[float]:
    factors = [float(part.strip()) for part in value.split(",") if part.strip()]
    if not factors:
        raise ValueError("at least one scale factor is required")
    if any(factor <= 0 for factor in factors):
        raise ValueError("scale factors must be greater than zero")
    return factors


def scale_token(scale_factor: float) -> str:
    return f"{scale_factor:.2f}".replace(".", "p")


def variant_id_for_scale(scale_factor: float) -> str:
    return f"hexaplex_base_length_scale_{scale_token(scale_factor)}"


def read_transformable_atom_names(path: Path) -> dict[str, set[str]]:
    if not path.exists():
        raise FileNotFoundError(f"selection summary is missing: {path}")
    selected: dict[str, set[str]] = defaultdict(set)
    with path.open("r", newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            residue_name = row.get("residue_name", "").strip().upper()
            atom_name = row.get("atom_name", "").strip().upper()
            if (
                residue_name in TRANSFORMABLE_RESIDUES
                and atom_name
                and row.get("selection_recommendation") == "transformable_candidate"
            ):
                selected[residue_name].add(atom_name)
    return dict(selected)


def atom_key(atom: PDBAtom) -> tuple[str, str, int | None, str, str]:
    return atom.chain_id, atom.residue_name, atom.residue_number, atom.insertion_code, atom.atom_name


def atom_xyz(atom: PDBAtom) -> tuple[float, float, float]:
    return atom.x, atom.y, atom.z


def residue_anchor(residue_atoms: list[PDBAtom]) -> tuple[tuple[float, float, float], str, str]:
    ca_atoms = [atom for atom in residue_atoms if atom.atom_name.strip().upper() == "CA"]
    if ca_atoms:
        return atom_xyz(ca_atoms[0]), "CA", ""
    backbone_atoms = [atom for atom in residue_atoms if atom.atom_name.strip().upper() in BACKBONE_ATOM_NAMES]
    if backbone_atoms:
        return mean_point(atom_xyz(atom) for atom in backbone_atoms), "backbone_centroid", ""
    return mean_point(atom_xyz(atom) for atom in residue_atoms), "residue_centroid", "missing backbone-like anchor atoms"


def transformed_atom_indices(
    atoms: list[PDBAtom],
    transformable_names: dict[str, set[str]],
) -> set[int]:
    indices: set[int] = set()
    for index, atom in enumerate(atoms):
        residue_name = atom.residue_name.upper()
        atom_name = atom.atom_name.strip().upper()
        if residue_name in transformable_names and atom_name in transformable_names[residue_name]:
            indices.add(index)
    return indices


def generate_variant_atoms(
    atoms: list[PDBAtom],
    transformable_names: dict[str, set[str]],
    scale_factor: float,
) -> tuple[list[PDBAtom], set[int], list[str]]:
    output = list(atoms)
    transformed_indices = transformed_atom_indices(atoms, transformable_names)
    warnings: list[str] = []
    start_by_identity = {id(atom): index for index, atom in enumerate(atoms)}

    for _residue_key, residue_atoms in group_atoms_by_residue(atoms).items():
        anchor, _anchor_rule, anchor_warning = residue_anchor(residue_atoms)
        if anchor_warning:
            warnings.append(anchor_warning)
        for atom in residue_atoms:
            index = start_by_identity[id(atom)]
            if index not in transformed_indices:
                continue
            dx = atom.x - anchor[0]
            dy = atom.y - anchor[1]
            dz = atom.z - anchor[2]
            output[index] = replace(
                atom,
                x=anchor[0] + scale_factor * dx,
                y=anchor[1] + scale_factor * dy,
                z=anchor[2] + scale_factor * dz,
            )
    return output, transformed_indices, sorted(set(warnings))


def write_csv(rows: list[dict[str, str]], path: Path, fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def infer_axis(atoms: list[PDBAtom]) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    heavy = heavy_atoms(atoms)
    backbone = [atom for atom in heavy if atom.atom_name.strip().upper() in BACKBONE_ATOM_NAMES]
    axis_atoms = backbone if len(backbone) >= 3 else heavy
    origin = mean_point(atom_xyz(atom) for atom in axis_atoms)
    axis = power_iteration_principal_axis(covariance_matrix_3d(atom_xyz(atom) for atom in axis_atoms))
    return origin, axis


def radial_extents(atoms: list[PDBAtom]) -> tuple[float, float, float]:
    origin, axis = infer_axis(atoms)
    values: list[float] = []
    for atom in atoms:
        _t, projected = project_point_to_axis(atom_xyz(atom), origin, axis)
        values.append(distance(atom_xyz(atom), projected))
    return min(values), max(values), sum(values) / len(values)


def min_heavy_distance_and_overlaps(
    atoms: list[PDBAtom],
    overlap_threshold: float,
) -> tuple[float, int]:
    heavy = heavy_atoms(atoms)
    min_distance = math.inf
    overlap_count = 0
    for i, atom_i in enumerate(heavy):
        for atom_j in heavy[i + 1 :]:
            atom_distance = distance(atom_i, atom_j)
            if atom_distance < min_distance:
                min_distance = atom_distance
            if atom_distance < overlap_threshold:
                overlap_count += 1
    return min_distance, overlap_count


def chain_count(atoms: list[PDBAtom]) -> int:
    return len({atom.chain_id for atom in atoms})


def residue_count(atoms: list[PDBAtom]) -> int:
    return len(group_atoms_by_residue(atoms))


def displacement_summary(
    baseline_atoms: list[PDBAtom],
    variant_atoms: list[PDBAtom],
    transformed_indices: set[int],
) -> tuple[float, float]:
    if not transformed_indices:
        return 0.0, 0.0
    values = [distance(baseline_atoms[index], variant_atoms[index]) for index in transformed_indices]
    return sum(values) / len(values), max(values)


def coordinates_match(
    atoms_a: list[PDBAtom],
    atoms_b: list[PDBAtom],
    tolerance: float = 1e-9,
) -> bool:
    if len(atoms_a) != len(atoms_b):
        return False
    for atom_a, atom_b in zip(atoms_a, atoms_b):
        if atom_key(atom_a) != atom_key(atom_b):
            return False
        if distance(atom_a, atom_b) > tolerance:
            return False
    return True


def format_float(value: float) -> str:
    return f"{value:.6f}"


def build_variant_rows(
    input_pdb: Path,
    output_dir: Path,
    scale_factors: list[float],
    atoms: list[PDBAtom],
    transformable_names: dict[str, set[str]],
    overlap_threshold: float,
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, object]]]:
    manifest_rows: list[dict[str, str]] = []
    geometry_rows: list[dict[str, str]] = []
    variant_summaries: list[dict[str, object]] = []
    output_dir.mkdir(parents=True, exist_ok=True)
    transformed_residue_types = sorted(transformable_names)
    transformed_atom_names = sorted({atom_name for names in transformable_names.values() for atom_name in names})
    anchor_rule = "CA if present; otherwise centroid of backbone-like N/CA/C/O/OXT atoms; otherwise residue centroid"
    transformation_rule = "new_xyz = anchor_xyz + scale_factor * (old_xyz - anchor_xyz)"

    for scale_factor in scale_factors:
        variant_id = variant_id_for_scale(scale_factor)
        output_pdb = output_dir / f"{variant_id}.pdb"
        variant_atoms, transformed_indices, generation_warnings = generate_variant_atoms(
            atoms,
            transformable_names,
            scale_factor,
        )
        write_pdb_atoms(variant_atoms, output_pdb)
        mean_displacement, max_displacement = displacement_summary(atoms, variant_atoms, transformed_indices)
        radial_min, radial_max, radial_mean = radial_extents(variant_atoms)
        min_distance, overlap_count = min_heavy_distance_and_overlaps(variant_atoms, overlap_threshold)
        scale_one_match = scale_factor == 1.0 and coordinates_match(atoms, variant_atoms)
        notes = []
        warnings = list(generation_warnings)
        if scale_factor == 1.0:
            notes.append("scale 1.00 baseline coordinate match" if scale_one_match else "scale 1.00 baseline mismatch")
            if not scale_one_match:
                warnings.append("scale 1.00 did not match baseline coordinates before PDB formatting")
        if overlap_count:
            warnings.append(f"{overlap_count} heavy-atom pair(s) below {overlap_threshold:.2f} A")

        fixed_count = len(atoms) - len(transformed_indices)
        manifest_rows.append(
            {
                "variant_id": variant_id,
                "scale_factor": f"{scale_factor:.2f}",
                "input_pdb": str(input_pdb),
                "output_pdb": str(output_pdb),
                "total_atom_count": str(len(atoms)),
                "transformed_atom_count": str(len(transformed_indices)),
                "fixed_atom_count": str(fixed_count),
                "transformed_residue_types": ";".join(transformed_residue_types),
                "transformed_atom_names": ";".join(transformed_atom_names),
                "anchor_rule": anchor_rule,
                "transformation_rule": transformation_rule,
                "notes": "; ".join(notes),
                "warnings": "; ".join(sorted(set(warnings))),
            }
        )
        geometry_rows.append(
            {
                "variant_id": variant_id,
                "scale_factor": f"{scale_factor:.2f}",
                "total_atom_count": str(len(atoms)),
                "transformed_atom_count": str(len(transformed_indices)),
                "fixed_atom_count": str(fixed_count),
                "chain_count": str(chain_count(variant_atoms)),
                "residue_count": str(residue_count(variant_atoms)),
                "min_heavy_atom_distance_A": format_float(min_distance),
                "suspicious_overlap_count": str(overlap_count),
                "radial_extent_min_A": format_float(radial_min),
                "radial_extent_max_A": format_float(radial_max),
                "radial_extent_mean_A": format_float(radial_mean),
                "mean_transformed_atom_displacement_A": format_float(mean_displacement),
                "max_transformed_atom_displacement_A": format_float(max_displacement),
                "notes": "; ".join(notes),
                "warnings": "; ".join(sorted(set(warnings))),
            }
        )
        variant_summaries.append(
            {
                "variant_id": variant_id,
                "scale_factor": scale_factor,
                "output_pdb": output_pdb,
                "transformed_atom_count": len(transformed_indices),
                "scale_one_match": scale_one_match if scale_factor == 1.0 else None,
                "suspicious_overlap_count": overlap_count,
                "mean_displacement": mean_displacement,
                "max_displacement": max_displacement,
            }
        )
    return manifest_rows, geometry_rows, variant_summaries


def markdown_table(rows: list[dict[str, str]], columns: list[str]) -> str:
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(row.get(column, "") for column in columns) + " |")
    return "\n".join(lines)


def write_report(
    path: Path,
    input_pdb: Path,
    manifest_rows: list[dict[str, str]],
    geometry_rows: list[dict[str, str]],
    transformable_names: dict[str, set[str]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    scale_values = [row["scale_factor"] for row in manifest_rows]
    transformed_counts = sorted({row["transformed_atom_count"] for row in manifest_rows})
    scale_one = next((row for row in geometry_rows if row["scale_factor"] == "1.00"), None)
    overlap_total = sum(int(row["suspicious_overlap_count"]) for row in geometry_rows)
    selected_atoms = "; ".join(
        f"{residue}:{','.join(sorted(names))}" for residue, names in sorted(transformable_names.items())
    )
    lines = [
        "# Base-Length Variant Generation Report",
        "",
        "## Scope and cautions",
        "",
        "- This is a geometry sensitivity study, not a structural determination.",
        "- CYP/MEP candidate arm atoms were axis-facing in the fitted-axis inspection; variants scale local anchor-to-atom vectors, not global outward radial vectors.",
        "- The diffraction sweep is not run by this workflow.",
        "- Native PDB residue names, chain IDs, atom order, and fixed atom coordinates are preserved by construction before PDB formatting.",
        "",
        "## Baseline",
        "",
        f"- Input PDB: {input_pdb}",
        "",
        "## Operational definition",
        "",
        "Base/hexad-arm length is defined as the local distance from each CYP/MEP residue backbone anchor to its selected non-backbone candidate arm atoms.",
        "",
        "## Atom-selection rule",
        "",
        "- Fixed atoms: all GLU atoms and CYP/MEP N, CA, C, O, OXT atoms.",
        f"- Transformable atoms from geometry inspection: {selected_atoms}",
        "",
        "## Anchor and transformation rules",
        "",
        "- Anchor rule: CA if present; otherwise centroid of available backbone-like N/CA/C/O/OXT atoms; otherwise residue centroid.",
        "- Transformation rule: `new_xyz = anchor_xyz + scale_factor * (old_xyz - anchor_xyz)`.",
        "",
        "## Generated scale factors",
        "",
        "- " + ", ".join(scale_values),
        "",
        "## Geometry sanity-check summary",
        "",
        f"- Variants generated: {len(manifest_rows)}",
        f"- Transformed atom count(s): {', '.join(transformed_counts)}",
        f"- Scale 1.00 matched baseline coordinates before formatting: {scale_one['notes'] if scale_one else 'not generated'}",
        f"- Total suspicious overlap counts across variants: {overlap_total}",
        "",
        markdown_table(
            geometry_rows,
            [
                "variant_id",
                "scale_factor",
                "transformed_atom_count",
                "min_heavy_atom_distance_A",
                "suspicious_overlap_count",
                "mean_transformed_atom_displacement_A",
                "max_transformed_atom_displacement_A",
                "notes",
                "warnings",
            ],
        ),
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    scale_factors = parse_scale_factors(args.scale_factors)
    atoms = load_pdb_atoms(args.pdb)
    transformable_names = read_transformable_atom_names(args.selection_summary)
    if not transformable_names:
        raise SystemExit(f"No transformable atoms found in {args.selection_summary}")
    manifest_rows, geometry_rows, _variant_summaries = build_variant_rows(
        args.pdb,
        args.out_dir,
        scale_factors,
        atoms,
        transformable_names,
        args.overlap_threshold,
    )
    write_csv(manifest_rows, args.manifest, MANIFEST_FIELDNAMES)
    write_csv(geometry_rows, args.geometry_out, GEOMETRY_FIELDNAMES)
    write_report(args.report, args.pdb, manifest_rows, geometry_rows, transformable_names)
    print(f"Wrote {len(manifest_rows)} variant PDB(s) to {args.out_dir}")
    print(f"Wrote {args.manifest}")
    print(f"Wrote {args.geometry_out}")
    print(f"Wrote {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
