#!/usr/bin/env python3
"""Generate rigid-body seed-formation ensembles and order parameters."""

from __future__ import annotations

import argparse
import csv
import math
import random
import sys
from collections import OrderedDict
from dataclasses import dataclass, replace
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from hexaplex_formation.geometry import (  # noqa: E402
    build_perpendicular_basis,
    covariance_matrix_3d,
    distance,
    dot,
    mean_point,
    norm,
    normalize_vector,
    power_iteration_principal_axis,
    project_point_to_axis,
    vector_sub,
)
from hexaplex_formation.pdb_utils import (  # noqa: E402
    PDBAtom,
    chain_ids,
    dedupe_exact_atoms,
    heavy_atoms,
    is_hydrogen,
    load_pdb_atoms,
    write_pdb_atoms,
)


BASE_RESIDUES = {"CYP", "MEP"}
BACKBONE_LIKE_ATOMS = {"N", "CA", "C", "O", "OXT"}
DEFAULT_UNIT_COUNTS = [4, 5, 6, 7, 8]
CONTACT_CUTOFF_A = 4.5
SUMMARY_COLUMNS = [
    "unit_count",
    "sample_id",
    "ensemble_type",
    "chain_count",
    "total_atom_count",
    "radius_of_gyration_A",
    "interchain_centroid_mean_distance_A",
    "interchain_centroid_min_distance_A",
    "RMSD_to_formed_seed_A",
    "contact_fraction_vs_target",
    "number_of_target_contacts",
    "number_of_contacts_present",
    "CYP_MEP_contact_fraction_vs_target",
    "backbone_contact_fraction_vs_target",
    "axial_extent_A",
    "radial_extent_A",
    "helical_axis_alignment_score",
    "axial_register_score",
    "angular_phase_order_score",
    "refined_angular_phase_score",
    "compactness_score",
    "seed_formation_score",
    "notes",
    "warnings",
]

ENDPOINT_MEAN_COLUMNS = [
    "unit_count",
    "ensemble_type",
    "sample_count",
    "radius_of_gyration_A_mean",
    "compactness_score_mean",
    "contact_fraction_vs_target_mean",
    "CYP_MEP_contact_fraction_vs_target_mean",
    "axial_register_score_mean",
    "angular_phase_order_score_mean",
    "refined_angular_phase_score_mean",
    "RMSD_to_formed_seed_A_mean",
    "seed_formation_score_mean",
]


Coordinate = tuple[float, float, float]
Matrix3 = tuple[Coordinate, Coordinate, Coordinate]


@dataclass(frozen=True)
class SeedReference:
    unit_count: int
    path: Path
    atoms: tuple[PDBAtom, ...]
    chain_atoms: OrderedDict[str, tuple[PDBAtom, ...]]
    target_contacts: frozenset[tuple[int, int]]
    target_cyp_mep_contacts: frozenset[tuple[int, int]]
    target_backbone_contacts: frozenset[tuple[int, int]]
    axis_origin: Coordinate
    axis: Coordinate
    chain_centroids: dict[str, Coordinate]
    chain_axial_positions: dict[str, float]
    chain_angles: dict[str, float]
    radius_of_gyration: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--structures-dir", type=Path, default=Path("outputs/mini_hexaplex/structures"))
    parser.add_argument("--unit-counts", default=",".join(str(value) for value in DEFAULT_UNIT_COUNTS))
    parser.add_argument("--samples-per-ensemble", type=int, default=25)
    parser.add_argument("--loose-mode", choices=["preserved_angular", "angular_randomized", "radially_separated", "all"], default="all")
    parser.add_argument("--random-seed", type=int, default=20260608)
    parser.add_argument("--contact-cutoff", type=float, default=CONTACT_CUTOFF_A)
    parser.add_argument("--ensemble-dir", type=Path, default=Path("outputs/seed_formation/ensembles"))
    parser.add_argument("--plot-dir", type=Path, default=Path("outputs/seed_formation/plots"))
    parser.add_argument("--out-csv", type=Path, default=Path("outputs/metrics/seed_formation_order_parameters.csv"))
    parser.add_argument("--endpoint-means-csv", type=Path, default=Path("outputs/metrics/seed_endpoint_metric_means_by_ensemble.csv"))
    parser.add_argument("--out-report", type=Path, default=Path("outputs/reports/seed_formation_order_parameters_report.md"))
    parser.add_argument("--save-examples", type=int, default=3)
    return parser.parse_args()


def atom_xyz(atom: PDBAtom) -> Coordinate:
    return atom.x, atom.y, atom.z


def format_float(value: float | None, digits: int = 6) -> str:
    if value is None or not math.isfinite(value):
        return ""
    return f"{value:.{digits}f}"


def parse_unit_counts(value: str) -> list[int]:
    counts: list[int] = []
    for token in value.split(","):
        stripped = token.strip()
        if not stripped:
            continue
        counts.append(int(stripped))
    if not counts:
        raise ValueError("At least one unit count is required")
    return counts


def group_by_chain(atoms: list[PDBAtom]) -> OrderedDict[str, tuple[PDBAtom, ...]]:
    grouped: OrderedDict[str, list[PDBAtom]] = OrderedDict()
    for atom in atoms:
        grouped.setdefault(atom.chain_id, []).append(atom)
    return OrderedDict((chain_id, tuple(chain_atoms)) for chain_id, chain_atoms in grouped.items())


def residue_keys_for_chain(atoms: tuple[PDBAtom, ...]) -> list[tuple[str, int | None, str]]:
    seen: OrderedDict[tuple[str, int | None, str], None] = OrderedDict()
    for atom in atoms:
        seen.setdefault((atom.residue_name, atom.residue_number, atom.insertion_code), None)
    return list(seen.keys())


def validate_seed_atoms(atoms: list[PDBAtom], unit_count: int, path: Path) -> OrderedDict[str, tuple[PDBAtom, ...]]:
    if not atoms:
        raise ValueError(f"No ATOM/HETATM records found in {path}")
    chain_order = chain_ids(atoms)
    if len(chain_order) != 6:
        raise ValueError(f"{path} has {len(chain_order)} chains ({','.join(chain_order)}), expected six")
    chain_atoms = group_by_chain(atoms)
    for chain_id, atoms_for_chain in chain_atoms.items():
        residues = residue_keys_for_chain(atoms_for_chain)
        observed_units = len(residues) / 2.0
        if len(residues) != unit_count * 2:
            raise ValueError(
                f"{path} chain {chain_id} has {len(residues)} residues ({observed_units:g} CYP/MEP-GLU units), "
                f"expected {unit_count * 2} residues for {unit_count} units"
            )
    return chain_atoms


def load_seed_reference(path: Path, unit_count: int, contact_cutoff: float) -> tuple[SeedReference, list[str]]:
    if not path.exists():
        raise FileNotFoundError(f"Mini-hexaplex PDB not found: {path}")
    raw_atoms = load_pdb_atoms(path)
    atoms = dedupe_exact_atoms(raw_atoms)
    warnings: list[str] = []
    if len(atoms) != len(raw_atoms):
        warnings.append(f"removed {len(raw_atoms) - len(atoms)} exact duplicate atoms from {path.name}")
    chain_atoms = validate_seed_atoms(atoms, unit_count, path)
    target_contacts, target_cyp_mep, target_backbone = find_interchain_contacts(atoms, contact_cutoff)
    axis_origin, axis = infer_axis(atoms)
    basis_u, basis_v = build_perpendicular_basis(axis)
    chain_centroids = {chain_id: mean_point(atom_xyz(atom) for atom in atoms_for_chain) for chain_id, atoms_for_chain in chain_atoms.items()}
    chain_axial = {chain_id: project_point_to_axis(point, axis_origin, axis)[0] for chain_id, point in chain_centroids.items()}
    chain_angles = {
        chain_id: math.atan2(
            dot(vector_sub(point, project_point_to_axis(point, axis_origin, axis)[1]), basis_v),
            dot(vector_sub(point, project_point_to_axis(point, axis_origin, axis)[1]), basis_u),
        )
        for chain_id, point in chain_centroids.items()
    }
    reference = SeedReference(
        unit_count=unit_count,
        path=path,
        atoms=tuple(atoms),
        chain_atoms=chain_atoms,
        target_contacts=frozenset(target_contacts),
        target_cyp_mep_contacts=frozenset(target_cyp_mep),
        target_backbone_contacts=frozenset(target_backbone),
        axis_origin=axis_origin,
        axis=axis,
        chain_centroids=chain_centroids,
        chain_axial_positions=chain_axial,
        chain_angles=chain_angles,
        radius_of_gyration=radius_of_gyration(atoms),
    )
    return reference, warnings


def infer_axis(atoms: list[PDBAtom] | tuple[PDBAtom, ...]) -> tuple[Coordinate, Coordinate]:
    heavy = heavy_atoms(atoms)
    axis_atoms = [atom for atom in heavy if atom.atom_name.strip().upper() in BACKBONE_LIKE_ATOMS]
    if len(axis_atoms) < 3:
        axis_atoms = heavy
    if len(axis_atoms) < 3:
        raise ValueError("At least three heavy atoms are required to infer a seed axis")
    origin = mean_point(atom_xyz(atom) for atom in axis_atoms)
    axis = power_iteration_principal_axis(covariance_matrix_3d(atom_xyz(atom) for atom in axis_atoms))
    return origin, axis


def radius_of_gyration(atoms: list[PDBAtom] | tuple[PDBAtom, ...]) -> float:
    points = [atom_xyz(atom) for atom in atoms]
    center = mean_point(points)
    return math.sqrt(sum(distance(point, center) ** 2 for point in points) / len(points))


def mat_vec_mul(matrix: Matrix3, vector: Coordinate) -> Coordinate:
    return (
        matrix[0][0] * vector[0] + matrix[0][1] * vector[1] + matrix[0][2] * vector[2],
        matrix[1][0] * vector[0] + matrix[1][1] * vector[1] + matrix[1][2] * vector[2],
        matrix[2][0] * vector[0] + matrix[2][1] * vector[1] + matrix[2][2] * vector[2],
    )


def rotation_matrix_from_axis_angle(axis: Coordinate, angle: float) -> Matrix3:
    ux, uy, uz = normalize_vector(axis)
    c = math.cos(angle)
    s = math.sin(angle)
    one_c = 1.0 - c
    return (
        (c + ux * ux * one_c, ux * uy * one_c - uz * s, ux * uz * one_c + uy * s),
        (uy * ux * one_c + uz * s, c + uy * uy * one_c, uy * uz * one_c - ux * s),
        (uz * ux * one_c - uy * s, uz * uy * one_c + ux * s, c + uz * uz * one_c),
    )


def random_unit_vector(rng: random.Random) -> Coordinate:
    z = rng.uniform(-1.0, 1.0)
    theta = rng.uniform(0.0, 2.0 * math.pi)
    r = math.sqrt(max(0.0, 1.0 - z * z))
    return r * math.cos(theta), r * math.sin(theta), z


def random_rotation_matrix(rng: random.Random, max_angle_deg: float) -> Matrix3:
    return rotation_matrix_from_axis_angle(random_unit_vector(rng), math.radians(rng.uniform(-max_angle_deg, max_angle_deg)))


def transform_chain(
    atoms: tuple[PDBAtom, ...],
    rotation: Matrix3,
    translation: Coordinate,
    center: Coordinate | None = None,
) -> list[PDBAtom]:
    if center is None:
        center = mean_point(atom_xyz(atom) for atom in atoms)
    transformed: list[PDBAtom] = []
    for atom in atoms:
        relative = vector_sub(atom_xyz(atom), center)
        rotated = mat_vec_mul(rotation, relative)
        x = rotated[0] + center[0] + translation[0]
        y = rotated[1] + center[1] + translation[1]
        z = rotated[2] + center[2] + translation[2]
        transformed.append(replace(atom, x=x, y=y, z=z))
    return transformed


def transform_chain_about_point(
    atoms: tuple[PDBAtom, ...],
    rotation: Matrix3,
    origin: Coordinate,
    translation: Coordinate = (0.0, 0.0, 0.0),
) -> list[PDBAtom]:
    transformed: list[PDBAtom] = []
    for atom in atoms:
        relative = vector_sub(atom_xyz(atom), origin)
        rotated = mat_vec_mul(rotation, relative)
        transformed.append(
            replace(
                atom,
                x=rotated[0] + origin[0] + translation[0],
                y=rotated[1] + origin[1] + translation[1],
                z=rotated[2] + origin[2] + translation[2],
            )
        )
    return transformed


def generate_formed_perturbed(reference: SeedReference, sample_index: int, rng: random.Random) -> list[PDBAtom]:
    atoms: list[PDBAtom] = []
    for chain_atoms in reference.chain_atoms.values():
        rotation = random_rotation_matrix(rng, max_angle_deg=4.0)
        translation = (
            rng.gauss(0.0, 0.35),
            rng.gauss(0.0, 0.35),
            rng.gauss(0.0, 0.35),
        )
        atoms.extend(transform_chain(chain_atoms, rotation, translation))
    return atoms


def generate_loose_initial(reference: SeedReference, sample_index: int, rng: random.Random) -> list[PDBAtom]:
    atoms: list[PDBAtom] = []
    axis = reference.axis
    basis_u, basis_v = build_perpendicular_basis(axis)
    for chain_index, (chain_id, chain_atoms) in enumerate(reference.chain_atoms.items()):
        centroid = reference.chain_centroids[chain_id]
        _, projected = project_point_to_axis(centroid, reference.axis_origin, axis)
        radial = vector_sub(centroid, projected)
        if norm(radial) < 1e-6:
            angle = chain_index * 2.0 * math.pi / 6.0
            radial_dir = (
                math.cos(angle) * basis_u[0] + math.sin(angle) * basis_v[0],
                math.cos(angle) * basis_u[1] + math.sin(angle) * basis_v[1],
                math.cos(angle) * basis_u[2] + math.sin(angle) * basis_v[2],
            )
        else:
            radial_dir = normalize_vector(radial)
        radial_separation = rng.uniform(8.0, 18.0)
        axial_separation = rng.uniform(-5.0, 5.0) + (chain_index - 2.5) * rng.uniform(0.6, 1.8)
        jitter = random_unit_vector(rng)
        translation = (
            radial_dir[0] * radial_separation + axis[0] * axial_separation + jitter[0] * rng.uniform(0.0, 2.0),
            radial_dir[1] * radial_separation + axis[1] * axial_separation + jitter[1] * rng.uniform(0.0, 2.0),
            radial_dir[2] * radial_separation + axis[2] * axial_separation + jitter[2] * rng.uniform(0.0, 2.0),
        )
        rotation = random_rotation_matrix(rng, max_angle_deg=75.0)
        atoms.extend(transform_chain(chain_atoms, rotation, translation))
    return atoms


def generate_angular_randomized_loose_initial(reference: SeedReference, sample_index: int, rng: random.Random) -> list[PDBAtom]:
    atoms: list[PDBAtom] = []
    axis = normalize_vector(reference.axis)
    basis_u, basis_v = build_perpendicular_basis(axis)
    used_angles: list[float] = []
    for chain_index, (chain_id, chain_atoms) in enumerate(reference.chain_atoms.items()):
        centroid = reference.chain_centroids[chain_id]
        axial_t, projected = project_point_to_axis(centroid, reference.axis_origin, axis)
        radial_distance = max(distance(centroid, projected), 1.0)
        formed_angle = reference.chain_angles[chain_id]
        # Keep trying until the chain label is clearly not near its formed phase.
        random_angle = rng.uniform(0.0, 2.0 * math.pi)
        for _ in range(40):
            too_close_to_formed = abs(math.atan2(math.sin(random_angle - formed_angle), math.cos(random_angle - formed_angle))) < math.radians(35.0)
            too_close_to_used = any(
                abs(math.atan2(math.sin(random_angle - used), math.cos(random_angle - used))) < math.radians(18.0)
                for used in used_angles
            )
            if not too_close_to_formed and not too_close_to_used:
                break
            random_angle = rng.uniform(0.0, 2.0 * math.pi)
        used_angles.append(random_angle)
        target_radial_distance = radial_distance + rng.uniform(8.0, 18.0)
        target_axial = axial_t + rng.uniform(-6.0, 6.0) + (chain_index - 2.5) * rng.uniform(0.4, 1.6)
        target_centroid = (
            reference.axis_origin[0]
            + axis[0] * target_axial
            + target_radial_distance * (math.cos(random_angle) * basis_u[0] + math.sin(random_angle) * basis_v[0]),
            reference.axis_origin[1]
            + axis[1] * target_axial
            + target_radial_distance * (math.cos(random_angle) * basis_u[1] + math.sin(random_angle) * basis_v[1]),
            reference.axis_origin[2]
            + axis[2] * target_axial
            + target_radial_distance * (math.cos(random_angle) * basis_u[2] + math.sin(random_angle) * basis_v[2]),
        )
        phase_rotation = rotation_matrix_from_axis_angle(axis, random_angle - formed_angle)
        rotated_chain = transform_chain_about_point(chain_atoms, phase_rotation, reference.axis_origin)
        rotated_centroid = mean_point(atom_xyz(atom) for atom in rotated_chain)
        translation = vector_sub(target_centroid, rotated_centroid)
        tumble = random_rotation_matrix(rng, max_angle_deg=45.0)
        atoms.extend(transform_chain(tuple(rotated_chain), tumble, translation, center=rotated_centroid))
    return atoms


def generate_radially_separated(reference: SeedReference, sample_index: int, rng: random.Random) -> list[PDBAtom]:
    """Dry-down-motivated loose class that mostly isolates radial compaction cost."""

    atoms: list[PDBAtom] = []
    axis = reference.axis
    basis_u, basis_v = build_perpendicular_basis(axis)
    for chain_index, (chain_id, chain_atoms) in enumerate(reference.chain_atoms.items()):
        centroid = reference.chain_centroids[chain_id]
        _, projected = project_point_to_axis(centroid, reference.axis_origin, axis)
        radial = vector_sub(centroid, projected)
        if norm(radial) < 1e-6:
            angle = chain_index * 2.0 * math.pi / 6.0
            radial_dir = (
                math.cos(angle) * basis_u[0] + math.sin(angle) * basis_v[0],
                math.cos(angle) * basis_u[1] + math.sin(angle) * basis_v[1],
                math.cos(angle) * basis_u[2] + math.sin(angle) * basis_v[2],
            )
        else:
            radial_dir = normalize_vector(radial)
        radial_separation = rng.uniform(18.0, 32.0)
        axial_separation = rng.uniform(-1.2, 1.2)
        tangential_jitter = rng.uniform(-1.0, 1.0)
        tangent = (
            -radial_dir[1] * axis[2] + radial_dir[2] * axis[1],
            -radial_dir[2] * axis[0] + radial_dir[0] * axis[2],
            -radial_dir[0] * axis[1] + radial_dir[1] * axis[0],
        )
        if norm(tangent) < 1e-6:
            tangent = basis_v
        tangent = normalize_vector(tangent)
        translation = (
            radial_dir[0] * radial_separation + axis[0] * axial_separation + tangent[0] * tangential_jitter,
            radial_dir[1] * radial_separation + axis[1] * axial_separation + tangent[1] * tangential_jitter,
            radial_dir[2] * radial_separation + axis[2] * axial_separation + tangent[2] * tangential_jitter,
        )
        rotation = random_rotation_matrix(rng, max_angle_deg=12.0)
        atoms.extend(transform_chain(chain_atoms, rotation, translation))
    return atoms


def is_cyp_mep_contact(atom_a: PDBAtom, atom_b: PDBAtom) -> bool:
    return atom_a.residue_name in BASE_RESIDUES or atom_b.residue_name in BASE_RESIDUES


def is_backbone_contact(atom_a: PDBAtom, atom_b: PDBAtom) -> bool:
    return atom_a.atom_name.strip().upper() in BACKBONE_LIKE_ATOMS and atom_b.atom_name.strip().upper() in BACKBONE_LIKE_ATOMS


def find_interchain_contacts(
    atoms: list[PDBAtom] | tuple[PDBAtom, ...],
    cutoff: float,
) -> tuple[set[tuple[int, int]], set[tuple[int, int]], set[tuple[int, int]]]:
    heavy_indexed = [(index, atom) for index, atom in enumerate(atoms) if not is_hydrogen(atom)]
    contacts: set[tuple[int, int]] = set()
    cyp_mep_contacts: set[tuple[int, int]] = set()
    backbone_contacts: set[tuple[int, int]] = set()
    for left_pos, (i, atom_i) in enumerate(heavy_indexed):
        for j, atom_j in heavy_indexed[left_pos + 1 :]:
            if atom_i.chain_id == atom_j.chain_id:
                continue
            if distance(atom_i, atom_j) > cutoff:
                continue
            key = (i, j)
            contacts.add(key)
            if is_cyp_mep_contact(atom_i, atom_j):
                cyp_mep_contacts.add(key)
            if is_backbone_contact(atom_i, atom_j):
                backbone_contacts.add(key)
    return contacts, cyp_mep_contacts, backbone_contacts


def contact_fraction(target_contacts: frozenset[tuple[int, int]], sample_contacts: set[tuple[int, int]]) -> tuple[float | None, int]:
    if not target_contacts:
        return None, 0
    present = len(target_contacts.intersection(sample_contacts))
    return present / len(target_contacts), present


def largest_eigenvector_symmetric_4(matrix: tuple[tuple[float, float, float, float], ...]) -> tuple[float, float, float, float]:
    vector = (1.0, 0.0, 0.0, 0.0)
    for _ in range(80):
        next_vector = tuple(sum(matrix[i][j] * vector[j] for j in range(4)) for i in range(4))
        length = math.sqrt(sum(value * value for value in next_vector))
        if length == 0:
            return vector
        next_vector = tuple(value / length for value in next_vector)
        if math.sqrt(sum((next_vector[i] - vector[i]) ** 2 for i in range(4))) < 1e-12:
            return next_vector  # type: ignore[return-value]
        vector = next_vector  # type: ignore[assignment]
    return vector


def quaternion_to_rotation_matrix(q: tuple[float, float, float, float]) -> Matrix3:
    w, x, y, z = q
    return (
        (1.0 - 2.0 * (y * y + z * z), 2.0 * (x * y - z * w), 2.0 * (x * z + y * w)),
        (2.0 * (x * y + z * w), 1.0 - 2.0 * (x * x + z * z), 2.0 * (y * z - x * w)),
        (2.0 * (x * z - y * w), 2.0 * (y * z + x * w), 1.0 - 2.0 * (x * x + y * y)),
    )


def optimal_alignment_rotation(moving_points: list[Coordinate], reference_points: list[Coordinate]) -> Matrix3:
    moving_center = mean_point(moving_points)
    reference_center = mean_point(reference_points)
    sxx = sxy = sxz = syx = syy = syz = szx = szy = szz = 0.0
    for moving, reference in zip(moving_points, reference_points):
        mx, my, mz = vector_sub(moving, moving_center)
        rx, ry, rz = vector_sub(reference, reference_center)
        sxx += mx * rx
        sxy += mx * ry
        sxz += mx * rz
        syx += my * rx
        syy += my * ry
        syz += my * rz
        szx += mz * rx
        szy += mz * ry
        szz += mz * rz
    trace = sxx + syy + szz
    matrix = (
        (trace, syz - szy, szx - sxz, sxy - syx),
        (syz - szy, sxx - syy - szz, sxy + syx, szx + sxz),
        (szx - sxz, sxy + syx, -sxx + syy - szz, syz + szy),
        (sxy - syx, szx + sxz, syz + szy, -sxx - syy + szz),
    )
    return quaternion_to_rotation_matrix(largest_eigenvector_symmetric_4(matrix))


def aligned_points(moving_points: list[Coordinate], reference_points: list[Coordinate]) -> list[Coordinate]:
    if len(moving_points) != len(reference_points) or not moving_points:
        raise ValueError("Alignment requires equal non-empty point sets")
    moving_center = mean_point(moving_points)
    reference_center = mean_point(reference_points)
    rotation = optimal_alignment_rotation(moving_points, reference_points)
    aligned: list[Coordinate] = []
    for point in moving_points:
        rotated = mat_vec_mul(rotation, vector_sub(point, moving_center))
        aligned.append((rotated[0] + reference_center[0], rotated[1] + reference_center[1], rotated[2] + reference_center[2]))
    return aligned


def rmsd_to_reference(sample_atoms: list[PDBAtom], reference_atoms: tuple[PDBAtom, ...]) -> float:
    moving = [atom_xyz(atom) for atom in sample_atoms]
    reference = [atom_xyz(atom) for atom in reference_atoms]
    aligned = aligned_points(moving, reference)
    return math.sqrt(sum(distance(a, b) ** 2 for a, b in zip(aligned, reference)) / len(reference))


def chain_centroid_metrics(chain_atoms: OrderedDict[str, tuple[PDBAtom, ...]]) -> tuple[float | None, float | None, dict[str, Coordinate]]:
    centroids = {chain_id: mean_point(atom_xyz(atom) for atom in atoms) for chain_id, atoms in chain_atoms.items()}
    distances: list[float] = []
    centroid_items = list(centroids.items())
    for index, (_, left) in enumerate(centroid_items):
        for _, right in centroid_items[index + 1 :]:
            distances.append(distance(left, right))
    if not distances:
        return None, None, centroids
    return sum(distances) / len(distances), min(distances), centroids


def axial_radial_extent(atoms: list[PDBAtom], axis_origin: Coordinate, axis: Coordinate) -> tuple[float, float]:
    axial_values: list[float] = []
    radial_values: list[float] = []
    for atom in atoms:
        z, projected = project_point_to_axis(atom_xyz(atom), axis_origin, axis)
        axial_values.append(z)
        radial_values.append(distance(atom_xyz(atom), projected))
    return max(axial_values) - min(axial_values), max(radial_values) - min(radial_values)


def axis_alignment_score(sample_atoms: list[PDBAtom], reference_axis: Coordinate) -> float | None:
    try:
        _, sample_axis = infer_axis(sample_atoms)
    except ValueError:
        return None
    return abs(dot(normalize_vector(sample_axis), normalize_vector(reference_axis)))


def axial_register_score(reference: SeedReference, centroids: dict[str, Coordinate]) -> float | None:
    deviations: list[float] = []
    sample_positions = {
        chain_id: project_point_to_axis(point, reference.axis_origin, reference.axis)[0]
        for chain_id, point in centroids.items()
        if chain_id in reference.chain_axial_positions
    }
    if len(sample_positions) != len(reference.chain_axial_positions):
        return None
    mean_reference = sum(reference.chain_axial_positions.values()) / len(reference.chain_axial_positions)
    mean_sample = sum(sample_positions.values()) / len(sample_positions)
    for chain_id, sample_z in sample_positions.items():
        deviations.append(abs((sample_z - mean_sample) - (reference.chain_axial_positions[chain_id] - mean_reference)))
    mean_deviation = sum(deviations) / len(deviations)
    return 1.0 / (1.0 + mean_deviation / 2.0)


def circular_difference(a: float, b: float) -> float:
    return abs(math.atan2(math.sin(a - b), math.cos(a - b)))


def signed_circular_difference(a: float, b: float) -> float:
    return math.atan2(math.sin(a - b), math.cos(a - b))


def angular_phase_order_score(reference: SeedReference, centroids: dict[str, Coordinate]) -> float | None:
    basis_u, basis_v = build_perpendicular_basis(reference.axis)
    differences: list[float] = []
    sample_angles: dict[str, float] = {}
    for chain_id, point in centroids.items():
        if chain_id not in reference.chain_angles:
            continue
        _, projected = project_point_to_axis(point, reference.axis_origin, reference.axis)
        radial = vector_sub(point, projected)
        sample_angles[chain_id] = math.atan2(dot(radial, basis_v), dot(radial, basis_u))
    if len(sample_angles) != len(reference.chain_angles):
        return None
    offsets = [circular_difference(sample_angles[chain_id], reference.chain_angles[chain_id]) for chain_id in sample_angles]
    mean_offset = sum(offsets) / len(offsets)
    for chain_id, angle in sample_angles.items():
        differences.append(circular_difference(angle - mean_offset, reference.chain_angles[chain_id]))
    mean_difference = sum(differences) / len(differences)
    return max(0.0, 1.0 - mean_difference / math.pi)


def refined_angular_phase_score(reference: SeedReference, centroids: dict[str, Coordinate]) -> float | None:
    sample_angles: dict[str, float] = {}
    basis_u, basis_v = build_perpendicular_basis(reference.axis)
    for chain_id, point in centroids.items():
        if chain_id not in reference.chain_angles:
            continue
        _, projected = project_point_to_axis(point, reference.axis_origin, reference.axis)
        radial = vector_sub(point, projected)
        if norm(radial) < 1e-9:
            continue
        sample_angles[chain_id] = math.atan2(dot(radial, basis_v), dot(radial, basis_u))
    shared_chains = sorted(set(reference.chain_angles).intersection(sample_angles))
    if len(shared_chains) < 3:
        return None
    deltas = [signed_circular_difference(sample_angles[chain_id], reference.chain_angles[chain_id]) for chain_id in shared_chains]
    offset = math.atan2(sum(math.sin(delta) for delta in deltas), sum(math.cos(delta) for delta in deltas))
    residuals = [signed_circular_difference(delta, offset) for delta in deltas]
    mean_cos = sum(math.cos(residual) for residual in residuals) / len(residuals)
    return max(0.0, min(1.0, (1.0 + mean_cos) / 2.0))


def compactness_score(sample_rg: float, reference_rg: float) -> float:
    if sample_rg <= 0:
        return 0.0
    return max(0.0, min(1.0, reference_rg / sample_rg))


def seed_formation_score(
    contact_fraction_value: float | None,
    rmsd: float | None,
    axial_score: float | None,
    angular_score: float | None,
    compactness: float | None,
) -> float | None:
    components = [
        contact_fraction_value,
        1.0 / (1.0 + rmsd / 5.0) if rmsd is not None else None,
        axial_score,
        angular_score,
        compactness,
    ]
    available = [value for value in components if value is not None and math.isfinite(value)]
    if not available:
        return None
    return sum(available) / len(available)


def analyze_sample(
    reference: SeedReference,
    sample_atoms: list[PDBAtom],
    sample_id: str,
    ensemble_type: str,
    contact_cutoff: float,
    reference_warnings: list[str] | None = None,
) -> dict[str, str]:
    chain_atoms = group_by_chain(sample_atoms)
    warnings = list(reference_warnings or [])
    if len(chain_atoms) != 6:
        warnings.append(f"sample has {len(chain_atoms)} chains, expected six")
    sample_contacts, sample_cyp_mep, sample_backbone = find_interchain_contacts(sample_atoms, contact_cutoff)
    contact_fraction_value, contacts_present = contact_fraction(reference.target_contacts, sample_contacts)
    cyp_mep_fraction, _ = contact_fraction(reference.target_cyp_mep_contacts, sample_cyp_mep)
    backbone_fraction, _ = contact_fraction(reference.target_backbone_contacts, sample_backbone)
    rg = radius_of_gyration(sample_atoms)
    centroid_mean, centroid_min, centroids = chain_centroid_metrics(chain_atoms)
    axial_extent, radial_extent = axial_radial_extent(sample_atoms, reference.axis_origin, reference.axis)
    try:
        rmsd = rmsd_to_reference(sample_atoms, reference.atoms)
    except ValueError as exc:
        rmsd = None
        warnings.append(str(exc))
    axis_score = axis_alignment_score(sample_atoms, reference.axis)
    register_score = axial_register_score(reference, centroids)
    phase_score = angular_phase_order_score(reference, centroids)
    refined_phase_score = refined_angular_phase_score(reference, centroids)
    compactness = compactness_score(rg, reference.radius_of_gyration)
    formation_score = seed_formation_score(contact_fraction_value, rmsd, register_score, phase_score, compactness)
    return {
        "unit_count": str(reference.unit_count),
        "sample_id": sample_id,
        "ensemble_type": ensemble_type,
        "chain_count": str(len(chain_atoms)),
        "total_atom_count": str(len(sample_atoms)),
        "radius_of_gyration_A": format_float(rg),
        "interchain_centroid_mean_distance_A": format_float(centroid_mean),
        "interchain_centroid_min_distance_A": format_float(centroid_min),
        "RMSD_to_formed_seed_A": format_float(rmsd),
        "contact_fraction_vs_target": format_float(contact_fraction_value),
        "number_of_target_contacts": str(len(reference.target_contacts)),
        "number_of_contacts_present": str(contacts_present),
        "CYP_MEP_contact_fraction_vs_target": format_float(cyp_mep_fraction),
        "backbone_contact_fraction_vs_target": format_float(backbone_fraction),
        "axial_extent_A": format_float(axial_extent),
        "radial_extent_A": format_float(radial_extent),
        "helical_axis_alignment_score": format_float(axis_score),
        "axial_register_score": format_float(register_score),
        "angular_phase_order_score": format_float(phase_score),
        "refined_angular_phase_score": format_float(refined_phase_score),
        "compactness_score": format_float(compactness),
        "seed_formation_score": format_float(formation_score),
        "notes": "contacts are evaluated in the sample frame using preserved atom indices; RMSD is separately rigid-aligned to the formed seed",
        "warnings": "; ".join(warnings),
    }


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def safe_float(value: str) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def write_plots(rows: list[dict[str, str]], plot_dir: Path) -> list[Path]:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover - optional dependency
        print(f"WARNING: matplotlib unavailable; skipping plots: {exc}", file=sys.stderr)
        return []

    plot_dir.mkdir(parents=True, exist_ok=True)
    plot_paths: list[Path] = []

    def grouped_xy(metric: str) -> dict[str, tuple[list[float], list[float]]]:
        grouped: dict[str, dict[int, list[float]]] = {}
        for row in rows:
            value = safe_float(row.get(metric, ""))
            unit_count = safe_float(row.get("unit_count", ""))
            if value is None or unit_count is None:
                continue
            grouped.setdefault(row["ensemble_type"], {}).setdefault(int(unit_count), []).append(value)
        out: dict[str, tuple[list[float], list[float]]] = {}
        for ensemble_type, by_unit in grouped.items():
            units = sorted(by_unit)
            means = [sum(by_unit[unit]) / len(by_unit[unit]) for unit in units]
            out[ensemble_type] = ([float(unit) for unit in units], means)
        return out

    def line_plot(metric: str, ylabel: str, filename: str) -> None:
        fig, ax = plt.subplots(figsize=(7.0, 4.5))
        for ensemble_type, (x_values, y_values) in grouped_xy(metric).items():
            ax.plot(x_values, y_values, marker="o", label=ensemble_type)
        ax.set_xlabel("Unit count")
        ax.set_ylabel(ylabel)
        ax.legend()
        ax.grid(True, alpha=0.25)
        fig.tight_layout()
        path = plot_dir / filename
        fig.savefig(path, dpi=160)
        plt.close(fig)
        plot_paths.append(path)

    def scatter_plot(x_metric: str, y_metric: str, xlabel: str, ylabel: str, filename: str) -> None:
        fig, ax = plt.subplots(figsize=(6.0, 5.0))
        for ensemble_type in sorted({row["ensemble_type"] for row in rows}):
            x_values = [safe_float(row.get(x_metric, "")) for row in rows if row["ensemble_type"] == ensemble_type]
            y_values = [safe_float(row.get(y_metric, "")) for row in rows if row["ensemble_type"] == ensemble_type]
            pairs = [(x, y) for x, y in zip(x_values, y_values) if x is not None and y is not None]
            if not pairs:
                continue
            ax.scatter([pair[0] for pair in pairs], [pair[1] for pair in pairs], s=18, alpha=0.7, label=ensemble_type)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.legend()
        ax.grid(True, alpha=0.25)
        fig.tight_layout()
        path = plot_dir / filename
        fig.savefig(path, dpi=160)
        plt.close(fig)
        plot_paths.append(path)

    line_plot("seed_formation_score", "Exploratory seed formation score", "unit_count_vs_seed_formation_score.png")
    line_plot("contact_fraction_vs_target", "Target contact fraction", "unit_count_vs_contact_fraction_vs_target.png")
    line_plot("RMSD_to_formed_seed_A", "RMSD to formed seed (A)", "unit_count_vs_RMSD_to_formed_seed_A.png")
    line_plot("CYP_MEP_contact_fraction_vs_target", "CYP/MEP target contact fraction", "unit_count_vs_CYP_MEP_contact_fraction_vs_target.png")
    scatter_plot(
        "RMSD_to_formed_seed_A",
        "contact_fraction_vs_target",
        "RMSD to formed seed (A)",
        "Target contact fraction",
        "RMSD_to_formed_seed_vs_contact_fraction.png",
    )
    scatter_plot(
        "radius_of_gyration_A",
        "contact_fraction_vs_target",
        "Radius of gyration (A)",
        "Target contact fraction",
        "radius_of_gyration_vs_contact_fraction.png",
    )
    return plot_paths


def mean_by_unit_and_ensemble(rows: list[dict[str, str]], metric: str) -> dict[tuple[int, str], float]:
    grouped: dict[tuple[int, str], list[float]] = {}
    for row in rows:
        value = safe_float(row.get(metric, ""))
        if value is None:
            continue
        grouped.setdefault((int(row["unit_count"]), row["ensemble_type"]), []).append(value)
    return {key: sum(values) / len(values) for key, values in grouped.items()}


def write_endpoint_metric_means(path: Path, rows: list[dict[str, str]]) -> None:
    metrics = [
        "radius_of_gyration_A",
        "compactness_score",
        "contact_fraction_vs_target",
        "CYP_MEP_contact_fraction_vs_target",
        "axial_register_score",
        "angular_phase_order_score",
        "refined_angular_phase_score",
        "RMSD_to_formed_seed_A",
        "seed_formation_score",
    ]
    grouped: dict[tuple[int, str], list[dict[str, str]]] = {}
    for row in rows:
        grouped.setdefault((int(row["unit_count"]), row["ensemble_type"]), []).append(row)
    mean_rows: list[dict[str, str]] = []
    for (unit_count, ensemble_type), group_rows in sorted(grouped.items()):
        out = {"unit_count": str(unit_count), "ensemble_type": ensemble_type, "sample_count": str(len(group_rows))}
        for metric in metrics:
            values = [safe_float(row.get(metric, "")) for row in group_rows]
            finite = [value for value in values if value is not None]
            out[f"{metric}_mean"] = format_float(sum(finite) / len(finite) if finite else None)
        mean_rows.append(out)
    write_csv(path, mean_rows, ENDPOINT_MEAN_COLUMNS)


def existing_rows_for_preserved_loose(path: Path, unit_counts: list[int]) -> list[dict[str, str]]:
    if not path.exists():
        return []
    unit_set = {str(unit_count) for unit_count in unit_counts}
    with path.open("r", newline="", encoding="utf-8") as handle:
        return [
            row
            for row in csv.DictReader(handle)
            if row.get("unit_count") in unit_set and row.get("ensemble_type") == "loose_initial"
        ]


def loose_modes_to_generate(loose_mode: str) -> list[str]:
    if loose_mode == "all":
        return ["preserved_angular", "angular_randomized", "radially_separated"]
    return [loose_mode]


def generate_loose_atoms_for_mode(
    reference: SeedReference,
    sample_index: int,
    rng: random.Random,
    mode: str,
) -> tuple[str, str, list[PDBAtom]]:
    if mode == "angular_randomized":
        return (
            "angular_randomized_loose_initial",
            f"central{reference.unit_count}_angular_randomized_loose_initial_{sample_index:04d}",
            generate_angular_randomized_loose_initial(reference, sample_index, rng),
        )
    if mode == "radially_separated":
        return (
            "radially_separated",
            f"central{reference.unit_count}_radially_separated_{sample_index:04d}",
            generate_radially_separated(reference, sample_index, rng),
        )
    return (
        "loose_initial",
        f"central{reference.unit_count}_loose_initial_{sample_index:04d}",
        generate_loose_initial(reference, sample_index, rng),
    )


def reagent_class_summary(reference: SeedReference) -> str:
    classes: OrderedDict[str, list[str]] = OrderedDict()
    for chain_id, chain_atoms in reference.chain_atoms.items():
        base_residue = next((atom.residue_name for atom in chain_atoms if atom.residue_name in BASE_RESIDUES), "unknown")
        classes.setdefault(base_residue, []).append(chain_id)
    if set(classes) >= BASE_RESIDUES and all(len(classes.get(residue, [])) == 3 for residue in BASE_RESIDUES):
        return "; ".join(f"{residue}: chains {','.join(chains)}" for residue, chains in classes.items())
    return "reagent-to-chain mapping is not explicitly encoded beyond observed residue identities; composition-aware costs require a future mapping step"


def ensemble_types_in_order(rows: list[dict[str, str]]) -> list[str]:
    preferred = ["formed_perturbed", "loose_initial", "angular_randomized_loose_initial", "radially_separated"]
    observed = {row["ensemble_type"] for row in rows}
    return [name for name in preferred if name in observed] + sorted(observed.difference(preferred))


def write_report(
    rows: list[dict[str, str]],
    plot_paths: list[Path],
    report_path: Path,
    contact_cutoff: float,
    samples_per_ensemble: int,
    loose_mode: str,
    reagent_summary: str,
) -> None:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    units = sorted({int(row["unit_count"]) for row in rows})
    score_means = mean_by_unit_and_ensemble(rows, "seed_formation_score")
    contact_means = mean_by_unit_and_ensemble(rows, "contact_fraction_vs_target")
    rmsd_means = mean_by_unit_and_ensemble(rows, "RMSD_to_formed_seed_A")
    ensemble_types = ensemble_types_in_order(rows)
    summary_columns = ["unit_count"]
    for ensemble_type in ensemble_types:
        summary_columns.append(f"{ensemble_type} score")
    for ensemble_type in ensemble_types:
        summary_columns.append(f"{ensemble_type} contacts")
    for ensemble_type in ensemble_types:
        summary_columns.append(f"{ensemble_type} RMSD A")
    lines = [
        "# Seed-Formation Order Parameters",
        "",
        "## Purpose",
        "",
        "This workflow prepares order parameters for asking how a short six-strand Hexaplex seed could form from initially separate or loosely associated strand fragments. The formed target is the central mini-hexaplex coordinate segment for each requested unit count.",
        "",
        "This is distinct from diffraction analysis. The metrics here are real-space structural descriptors for synthetic rigid-body ensembles, not reciprocal-space scattering features or assignments of diffraction peaks.",
        "",
        "The workflow is a preparatory layer for later Schrodinger bridge or stochastic pathway modeling: it defines target and initial ensembles plus measurable coordinates that can become bridge constraints or diagnostics. It does not implement a bridge.",
        "",
        "## Formed Seed Target",
        "",
        f"Targets were loaded from `outputs/mini_hexaplex/structures/mini_hexaplex_centralN_units.pdb` for unit counts {', '.join(str(unit) for unit in units)}. Each target was deduplicated by exact atom identity/coordinate, checked for six chains, and checked for N CYP/MEP-GLU units per chain.",
        "",
        "## Ensemble Generation",
        "",
        f"For each unit count, `{samples_per_ensemble}` `formed_perturbed` samples were generated by applying small independent rigid-body rotations and translations to each chain. The requested loose mode was `{loose_mode}`. `preserved_angular` writes `loose_initial`; `angular_randomized` writes `angular_randomized_loose_initial` with chain-label angular phases randomized around the formed seed axis before loose radial/axial separation; `radially_separated` writes a dry-down-motivated class that preserves rough angular phase while increasing radial separation and limiting axial offsets. Intrachain geometry is preserved by construction.",
        "",
        "The lab preparation context motivates treating dry-down as a potentially relevant concentration/compaction phase. The `radially_separated` class is only a synthetic geometric proxy for that idea; it does not model solvent evaporation, crystallization, molecular dynamics, or a full atomistic Schrodinger bridge.",
        "",
        "Severe overlaps are reduced operationally by separating chains radially and axially, but no energy model or full steric optimizer is used.",
        "",
        f"Observed residue-identity composition in the first requested target: {reagent_summary}.",
        "",
        "## Order Parameters",
        "",
        f"Target contacts use heavy interchain atom pairs within {contact_cutoff:.2f} A in the formed seed. Contact fractions are evaluated in the generated sample frame using preserved atom indices because the fragments are generated directly from the target atom ordering. RMSD is computed separately after optimal rigid alignment to the formed seed.",
        "",
        "Computed metrics include radius of gyration, interchain centroid distances, rigid-aligned RMSD, target contact recovery, CYP/MEP-involving contact recovery, backbone-like contact recovery, axial and radial extents, fitted-axis alignment, axial register, angular phase order, compactness, and the exploratory score.",
        "",
        "The CYP/MEP contact subset includes target contacts where either residue is CYP or MEP. The backbone-like subset includes target contacts where both atom names are in N/CA/C/O/OXT.",
        "",
        "## Seed Formation Score",
        "",
        "`seed_formation_score` is an exploratory unweighted average of visible component scores: target contact fraction, `1 / (1 + RMSD_to_formed_seed_A / 5)`, axial register score, angular phase order score, and compactness score. The weights are arbitrary and should be treated as a diagnostic convenience, not a validated reaction coordinate.",
        "",
        "## Summary",
        "",
        "| " + " | ".join(summary_columns) + " |",
        "|" + "|".join("---:" for _ in summary_columns) + "|",
    ]
    for unit in units:
        values = [str(unit)]
        values.extend(format_float(score_means.get((unit, ensemble_type)), 4) for ensemble_type in ensemble_types)
        values.extend(format_float(contact_means.get((unit, ensemble_type)), 4) for ensemble_type in ensemble_types)
        values.extend(format_float(rmsd_means.get((unit, ensemble_type)), 3) for ensemble_type in ensemble_types)
        lines.append("| " + " | ".join(values) + " |")
    lines.extend(["", "## Summary Plots", ""])
    for path in plot_paths:
        relative = path.relative_to(REPO_ROOT) if path.is_absolute() and path.is_relative_to(REPO_ROOT) else path
        lines.append(f"- `{relative}`")
    lines.extend(
        [
            "",
            "## Conservative Interpretation",
            "",
            "The formed-perturbed ensemble is expected to remain compact, contact-rich, and low-RMSD. The loose-start ensembles are expected to be expanded, contact-poor, and high-RMSD. The `radially_separated` class is intended to isolate radial compaction more than angular scrambling. A clear separation between those distributions means the order parameters can distinguish formed-like seeds from dispersed rigid-chain arrangements; it does not show that loose strands spontaneously assemble.",
            "",
            "Metrics that are especially relevant for later bridge modeling are contact recovery, RMSD to formed seed, compactness/radius of gyration, axial register, and angular phase order because they describe complementary aspects of closure, shape, register, and six-strand phasing.",
            "",
            "## Limitations",
            "",
            "- Rigid-body strand fragments only.",
            "- No molecular dynamics.",
            "- No solvent or counterions.",
            "- No energy model yet.",
            "- Ensemble generation is synthetic.",
            "- Order parameters are exploratory.",
            "- This is not evidence of spontaneous stability or formation.",
            "",
        ]
    )
    report_path.write_text("\n".join(lines), encoding="utf-8")


def run(args: argparse.Namespace) -> list[dict[str, str]]:
    if args.samples_per_ensemble < 1:
        raise ValueError("--samples-per-ensemble must be positive")
    unit_counts = parse_unit_counts(args.unit_counts)
    rng = random.Random(args.random_seed)
    rows: list[dict[str, str]] = []
    reagent_summary = "not evaluated"
    args.ensemble_dir.mkdir(parents=True, exist_ok=True)
    for unit_count in unit_counts:
        path = args.structures_dir / f"mini_hexaplex_central{unit_count}_units.pdb"
        reference, reference_warnings = load_seed_reference(path, unit_count, args.contact_cutoff)
        if reagent_summary == "not evaluated":
            reagent_summary = reagent_class_summary(reference)
        for sample_index in range(args.samples_per_ensemble):
            sample_id = f"central{unit_count}_formed_perturbed_{sample_index:04d}"
            atoms = generate_formed_perturbed(reference, sample_index, rng)
            rows.append(analyze_sample(reference, atoms, sample_id, "formed_perturbed", args.contact_cutoff, reference_warnings))
            if sample_index < args.save_examples:
                write_pdb_atoms(atoms, args.ensemble_dir / f"{sample_id}.pdb")
        for mode in loose_modes_to_generate(args.loose_mode):
            for sample_index in range(args.samples_per_ensemble):
                ensemble_type, sample_id, atoms = generate_loose_atoms_for_mode(reference, sample_index, rng, mode)
                rows.append(analyze_sample(reference, atoms, sample_id, ensemble_type, args.contact_cutoff, reference_warnings))
                if sample_index < args.save_examples:
                    write_pdb_atoms(atoms, args.ensemble_dir / f"{sample_id}.pdb")
    output_rows = rows
    if args.loose_mode == "angular_randomized":
        output_rows = existing_rows_for_preserved_loose(args.out_csv, unit_counts) + rows
    write_csv(args.out_csv, output_rows, SUMMARY_COLUMNS)
    write_endpoint_metric_means(args.endpoint_means_csv, output_rows)
    plot_paths = write_plots(output_rows, args.plot_dir)
    write_report(output_rows, plot_paths, args.out_report, args.contact_cutoff, args.samples_per_ensemble, args.loose_mode, reagent_summary)
    return output_rows


def main() -> None:
    try:
        rows = run(parse_args())
    except Exception as exc:
        raise SystemExit(f"ERROR: {exc}") from exc
    print(f"Wrote {len(rows)} seed-formation order-parameter rows")


if __name__ == "__main__":
    main()
