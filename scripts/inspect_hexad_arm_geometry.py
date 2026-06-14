#!/usr/bin/env python3
"""Inspect nonstandard CYP/MEP atom geometry for future base-length scaling."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from hexaplex_formation.geometry import (  # noqa: E402
    covariance_matrix_3d,
    group_atoms_by_residue,
    mean_point,
    power_iteration_principal_axis,
    project_point_to_axis,
)
from hexaplex_formation.pdb_utils import PDBAtom, heavy_atoms, load_pdb_atoms  # noqa: E402


BACKBONE_ATOM_NAMES = {"N", "CA", "C", "O", "OXT"}
TARGET_RESIDUES = {"CYP", "MEP", "GLU"}
ARM_RESIDUES = {"CYP", "MEP"}

INVENTORY_FIELDNAMES = [
    "residue_name",
    "atom_name",
    "classification",
    "atom_count",
    "residue_occurrence_count",
    "chain_count",
    "chains",
    "elements",
]

GEOMETRY_FIELDNAMES = [
    "residue_name",
    "atom_name",
    "classification",
    "atom_count",
    "residue_occurrence_count",
    "mean_radial_distance_A",
    "min_radial_distance_A",
    "max_radial_distance_A",
    "stddev_radial_distance_A",
    "mean_delta_from_residue_backbone_A",
    "outward_candidate",
    "axis_facing_candidate",
    "selection_recommendation",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pdb",
        type=Path,
        default=Path("outputs/intermediates/ai_candidate_inputs/full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain.pdb"),
    )
    parser.add_argument("--out-inventory", type=Path, default=Path("outputs/metrics/hexad_arm_atom_inventory.csv"))
    parser.add_argument("--out-summary", type=Path, default=Path("outputs/metrics/hexad_arm_geometry_summary.csv"))
    parser.add_argument("--out-md", type=Path, default=Path("outputs/reports/hexad_arm_geometry_report.md"))
    parser.add_argument("--outward-delta", type=float, default=1.0)
    parser.add_argument("--min-occurrence-fraction", type=float, default=0.8)
    return parser.parse_args()


def atom_classification(atom: PDBAtom) -> str:
    atom_name = atom.atom_name.strip().upper()
    residue_name = atom.residue_name.strip().upper()
    if atom_name in BACKBONE_ATOM_NAMES:
        return "backbone_like"
    if residue_name == "GLU":
        return "glu_scaffold"
    if residue_name in ARM_RESIDUES:
        return "candidate_sidechain_or_arm"
    return "other"


def _coord(atom: PDBAtom) -> tuple[float, float, float]:
    return atom.x, atom.y, atom.z


def infer_axis(atoms: list[PDBAtom]) -> tuple[tuple[float, float, float], tuple[float, float, float], str]:
    heavy = heavy_atoms(atoms)
    backbone = [atom for atom in heavy if atom.atom_name.strip().upper() in BACKBONE_ATOM_NAMES]
    axis_atoms = backbone if len(backbone) >= 3 else heavy
    origin = mean_point(_coord(atom) for atom in axis_atoms)
    axis = power_iteration_principal_axis(covariance_matrix_3d(_coord(atom) for atom in axis_atoms))
    source = "backbone_like_heavy_atoms" if axis_atoms is backbone else "all_heavy_atoms"
    return origin, axis, source


def residue_label(atom: PDBAtom) -> str:
    number = "" if atom.residue_number is None else str(atom.residue_number)
    label = f"{atom.residue_name}{number}{atom.insertion_code}"
    return f"{atom.chain_id}:{label}" if atom.chain_id else label


def radial_distances(
    atoms: list[PDBAtom],
    origin: tuple[float, float, float],
    axis: tuple[float, float, float],
) -> dict[int, float]:
    distances: dict[int, float] = {}
    for index, atom in enumerate(atoms):
        _axial_t, projected = project_point_to_axis(_coord(atom), origin, axis)
        dx = atom.x - projected[0]
        dy = atom.y - projected[1]
        dz = atom.z - projected[2]
        distances[index] = math.sqrt(dx * dx + dy * dy + dz * dz)
    return distances


def _format_float(value: float | None) -> str:
    return "" if value is None else f"{value:.6f}"


def _mean(values: list[float]) -> float | None:
    return sum(values) / len(values) if values else None


def _stddev(values: list[float]) -> float | None:
    if not values:
        return None
    mean_value = sum(values) / len(values)
    return math.sqrt(sum((value - mean_value) ** 2 for value in values) / len(values))


def build_inventory_rows(atoms: list[PDBAtom]) -> list[dict[str, str]]:
    grouped: dict[tuple[str, str], list[PDBAtom]] = defaultdict(list)
    for atom in atoms:
        if atom.residue_name.upper() in TARGET_RESIDUES:
            grouped[(atom.residue_name.upper(), atom.atom_name.strip().upper())].append(atom)

    rows: list[dict[str, str]] = []
    for (residue_name, atom_name), atom_group in sorted(grouped.items()):
        residue_keys = {
            (atom.chain_id, atom.residue_name, atom.residue_number, atom.insertion_code)
            for atom in atom_group
        }
        chains = sorted({atom.chain_id or "(blank)" for atom in atom_group})
        elements = sorted({atom.element.upper() for atom in atom_group if atom.element})
        rows.append(
            {
                "residue_name": residue_name,
                "atom_name": atom_name,
                "classification": atom_classification(atom_group[0]),
                "atom_count": str(len(atom_group)),
                "residue_occurrence_count": str(len(residue_keys)),
                "chain_count": str(len(chains)),
                "chains": ";".join(chains),
                "elements": ";".join(elements),
            }
        )
    return rows


def _residue_backbone_means(atoms: list[PDBAtom], radial_by_index: dict[int, float]) -> dict[str, float]:
    values: dict[str, list[float]] = defaultdict(list)
    for index, atom in enumerate(atoms):
        residue_name = atom.residue_name.upper()
        if residue_name in TARGET_RESIDUES and atom.atom_name.strip().upper() in BACKBONE_ATOM_NAMES:
            values[residue_name].append(radial_by_index[index])
    return {residue_name: sum(distances) / len(distances) for residue_name, distances in values.items() if distances}


def build_geometry_rows(
    atoms: list[PDBAtom],
    radial_by_index: dict[int, float],
    outward_delta: float = 1.0,
    min_occurrence_fraction: float = 0.8,
) -> list[dict[str, str]]:
    grouped_indices: dict[tuple[str, str], list[int]] = defaultdict(list)
    residue_counts: dict[str, int] = defaultdict(int)
    for _key, residue_atoms in group_atoms_by_residue(atoms).items():
        residue_name = residue_atoms[0].residue_name.upper()
        if residue_name in TARGET_RESIDUES:
            residue_counts[residue_name] += 1
    for index, atom in enumerate(atoms):
        residue_name = atom.residue_name.upper()
        if residue_name in TARGET_RESIDUES:
            grouped_indices[(residue_name, atom.atom_name.strip().upper())].append(index)

    backbone_means = _residue_backbone_means(atoms, radial_by_index)
    rows: list[dict[str, str]] = []
    for (residue_name, atom_name), indices in sorted(grouped_indices.items()):
        atom_group = [atoms[index] for index in indices]
        values = [radial_by_index[index] for index in indices]
        mean_value = _mean(values)
        backbone_mean = backbone_means.get(residue_name)
        delta = mean_value - backbone_mean if mean_value is not None and backbone_mean is not None else None
        residue_occurrences = {
            (atom.chain_id, atom.residue_name, atom.residue_number, atom.insertion_code)
            for atom in atom_group
        }
        occurrence_fraction = len(residue_occurrences) / residue_counts[residue_name] if residue_counts[residue_name] else 0.0
        classification = atom_classification(atom_group[0])
        is_consistent_arm = (
            residue_name in ARM_RESIDUES
            and classification == "candidate_sidechain_or_arm"
            and delta is not None
            and abs(delta) >= outward_delta
            and occurrence_fraction >= min_occurrence_fraction
        )
        is_outward = is_consistent_arm and delta is not None and delta > 0
        is_axis_facing = is_consistent_arm and delta is not None and delta < 0
        if is_consistent_arm:
            recommendation = "transformable_candidate"
        elif classification in {"backbone_like", "glu_scaffold"}:
            recommendation = "fixed_candidate"
        else:
            recommendation = "review"
        rows.append(
            {
                "residue_name": residue_name,
                "atom_name": atom_name,
                "classification": classification,
                "atom_count": str(len(indices)),
                "residue_occurrence_count": str(len(residue_occurrences)),
                "mean_radial_distance_A": _format_float(mean_value),
                "min_radial_distance_A": _format_float(min(values)),
                "max_radial_distance_A": _format_float(max(values)),
                "stddev_radial_distance_A": _format_float(_stddev(values)),
                "mean_delta_from_residue_backbone_A": _format_float(delta),
                "outward_candidate": "yes" if is_outward else "no",
                "axis_facing_candidate": "yes" if is_axis_facing else "no",
                "selection_recommendation": recommendation,
            }
        )
    return rows


def chain_inventory(atoms: list[PDBAtom]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    grouped: dict[str, list[tuple[str, int | None, str]]] = defaultdict(list)
    for key, residue_atoms in group_atoms_by_residue(atoms).items():
        chain_id, residue_name, residue_number, insertion_code = key
        grouped[chain_id].append((residue_name, residue_number, insertion_code))
    for chain_id in sorted(grouped):
        residues = grouped[chain_id]
        rows.append(
            {
                "chain_id": chain_id or "(blank)",
                "residue_count": str(len(residues)),
                "first_residue": f"{residues[0][0]}{residues[0][1] or ''}{residues[0][2]}",
                "last_residue": f"{residues[-1][0]}{residues[-1][1] or ''}{residues[-1][2]}",
                "residue_pattern": "-".join(residue_name for residue_name, _number, _icode in residues),
            }
        )
    return rows


def write_csv(rows: list[dict[str, str]], path: Path, fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _markdown_table(rows: list[dict[str, str]], columns: list[str]) -> str:
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(row.get(column, "") for column in columns) + " |")
    return "\n".join(lines)


def write_report(
    atoms: list[PDBAtom],
    inventory_rows: list[dict[str, str]],
    geometry_rows: list[dict[str, str]],
    pdb_path: Path,
    out_path: Path,
    axis_origin: tuple[float, float, float],
    axis: tuple[float, float, float],
    axis_source: str,
) -> None:
    chains = chain_inventory(atoms)
    fixed_atoms = [
        row
        for row in geometry_rows
        if row["selection_recommendation"] == "fixed_candidate" and row["residue_name"] in TARGET_RESIDUES
    ]
    transformable_atoms = [row for row in geometry_rows if row["selection_recommendation"] == "transformable_candidate"]
    review_atoms = [
        row
        for row in geometry_rows
        if row["residue_name"] in ARM_RESIDUES and row["selection_recommendation"] == "review"
    ]
    lines = [
        "# Hexad Arm Atom Geometry Inspection",
        "",
        "## Scope and cautions",
        "",
        "- This report inspects atom-level geometry only; it does not assign biological meaning to CYP or MEP.",
        "- The central axis is an approximate fitted axis used for radial summaries.",
        "- Proposed fixed/transformable atoms are operational candidates for a later base-length scaling workflow, not chemistry claims.",
        "- The diffraction sweep is not implemented here.",
        "",
        "## Input and inferred axis",
        "",
        f"- PDB: {pdb_path}",
        f"- Atom count: {len(atoms)}",
        f"- Residue count: {len(group_atoms_by_residue(atoms))}",
        f"- Axis source: {axis_source}",
        f"- Axis origin: ({axis_origin[0]:.6f}, {axis_origin[1]:.6f}, {axis_origin[2]:.6f})",
        f"- Axis direction: ({axis[0]:.6f}, {axis[1]:.6f}, {axis[2]:.6f})",
        "",
        "## Inferred chain structure",
        "",
        _markdown_table(chains, ["chain_id", "residue_count", "first_residue", "last_residue"]),
        "",
        "## Atom-name inventory",
        "",
        _markdown_table(
            inventory_rows,
            ["residue_name", "atom_name", "classification", "atom_count", "residue_occurrence_count", "chains"],
        ),
        "",
        "## Proposed fixed atoms",
        "",
        _markdown_table(
            fixed_atoms,
            ["residue_name", "atom_name", "classification", "mean_radial_distance_A", "selection_recommendation"],
        ),
        "",
        "## Proposed transformable atoms",
        "",
        _markdown_table(
            transformable_atoms,
            [
                "residue_name",
                "atom_name",
                "mean_radial_distance_A",
                "mean_delta_from_residue_backbone_A",
                "axis_facing_candidate",
                "outward_candidate",
                "selection_recommendation",
            ],
        ),
        "",
        "## Ambiguous/review atoms",
        "",
        _markdown_table(
            review_atoms,
            ["residue_name", "atom_name", "classification", "mean_radial_distance_A", "mean_delta_from_residue_backbone_A"],
        ),
        "",
        "## Recommended operational base-length definition",
        "",
        "For the later base-length sweep, keep all GLU atoms and all N/CA/C/O/OXT backbone-like atoms fixed. For CYP and MEP, treat non-backbone atoms flagged as `transformable_candidate` as candidate arm atoms whose radial component from the fitted central axis may be scaled. In the current six-chain geometry these candidates are axis-facing rather than outward relative to the CYP/MEP backbone-like radial mean. Leave atoms marked `review` fixed until a visual or chemistry-aware inspection justifies including them.",
        "",
        "The base length should be measured as the radial distance from the fitted central axis to the selected transformable CYP/MEP atom group, summarized per residue and chain. Scaling should preserve residue identity and should be reported as a geometric perturbation only.",
        "",
        "## Output files",
        "",
        "- `outputs/metrics/hexad_arm_atom_inventory.csv`",
        "- `outputs/metrics/hexad_arm_geometry_summary.csv`",
        "",
    ]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    atoms = load_pdb_atoms(args.pdb)
    axis_origin, axis, axis_source = infer_axis(atoms)
    radial_by_index = radial_distances(atoms, axis_origin, axis)
    inventory_rows = build_inventory_rows(atoms)
    geometry_rows = build_geometry_rows(
        atoms,
        radial_by_index,
        outward_delta=args.outward_delta,
        min_occurrence_fraction=args.min_occurrence_fraction,
    )
    write_csv(inventory_rows, args.out_inventory, INVENTORY_FIELDNAMES)
    write_csv(geometry_rows, args.out_summary, GEOMETRY_FIELDNAMES)
    write_report(atoms, inventory_rows, geometry_rows, args.pdb, args.out_md, axis_origin, axis, axis_source)
    print(f"Wrote {args.out_inventory}")
    print(f"Wrote {args.out_summary}")
    print(f"Wrote {args.out_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
