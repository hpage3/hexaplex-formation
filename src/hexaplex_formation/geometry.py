"""Small geometry helpers for formation-oriented structure metrics."""

from __future__ import annotations

import math
from collections import OrderedDict
from typing import Iterable

from hexaplex_formation.pdb_utils import PDBAtom, heavy_atoms

Coordinate = tuple[float, float, float]
ResidueKey = tuple[str, str, int | None, str]


def _coord(value: PDBAtom | Coordinate) -> Coordinate:
    if isinstance(value, PDBAtom):
        return value.x, value.y, value.z
    return value


def distance(a: PDBAtom | Coordinate, b: PDBAtom | Coordinate) -> float:
    ax, ay, az = _coord(a)
    bx, by, bz = _coord(b)
    dx = ax - bx
    dy = ay - by
    dz = az - bz
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def vector_sub(a: PDBAtom | Coordinate, b: PDBAtom | Coordinate) -> Coordinate:
    ax, ay, az = _coord(a)
    bx, by, bz = _coord(b)
    return ax - bx, ay - by, az - bz


def vector_add(u: Coordinate, v: Coordinate) -> Coordinate:
    return u[0] + v[0], u[1] + v[1], u[2] + v[2]


def scalar_mul(value: float, u: Coordinate) -> Coordinate:
    return value * u[0], value * u[1], value * u[2]


def dot(u: Coordinate, v: Coordinate) -> float:
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]


def norm(u: Coordinate) -> float:
    return math.sqrt(dot(u, u))


def normalize_vector(v: Coordinate) -> Coordinate:
    length = norm(v)
    if length == 0:
        raise ValueError("Cannot normalize zero-length vector")
    return v[0] / length, v[1] / length, v[2] / length


def cross(u: Coordinate, v: Coordinate) -> Coordinate:
    return (
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    )


def mean_point(points: Iterable[Coordinate]) -> Coordinate:
    point_list = list(points)
    if not point_list:
        raise ValueError("mean_point requires at least one point")
    count = len(point_list)
    return (
        sum(point[0] for point in point_list) / count,
        sum(point[1] for point in point_list) / count,
        sum(point[2] for point in point_list) / count,
    )


def covariance_matrix_3d(points: Iterable[Coordinate]) -> tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]:
    point_list = list(points)
    if not point_list:
        raise ValueError("covariance_matrix_3d requires at least one point")
    center = mean_point(point_list)
    accum = [[0.0, 0.0, 0.0] for _ in range(3)]
    for point in point_list:
        delta = vector_sub(point, center)
        for i in range(3):
            for j in range(3):
                accum[i][j] += delta[i] * delta[j]
    count = len(point_list)
    return tuple(tuple(value / count for value in row) for row in accum)  # type: ignore[return-value]


def _matrix_vector_mul(
    matrix: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]],
    vector: Coordinate,
) -> Coordinate:
    return (
        matrix[0][0] * vector[0] + matrix[0][1] * vector[1] + matrix[0][2] * vector[2],
        matrix[1][0] * vector[0] + matrix[1][1] * vector[1] + matrix[1][2] * vector[2],
        matrix[2][0] * vector[0] + matrix[2][1] * vector[1] + matrix[2][2] * vector[2],
    )


def _apply_axis_sign_convention(axis: Coordinate) -> Coordinate:
    eps = 1e-12
    if axis[2] < -eps or (abs(axis[2]) <= eps and axis[0] < 0):
        return -axis[0], -axis[1], -axis[2]
    return axis


def power_iteration_principal_axis(
    cov_matrix: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]],
    iterations: int = 100,
    tolerance: float = 1e-12,
) -> Coordinate:
    axis: Coordinate = normalize_vector((1.0, 1.0, 1.0))
    for _ in range(iterations):
        next_axis = _matrix_vector_mul(cov_matrix, axis)
        if norm(next_axis) == 0:
            return 0.0, 0.0, 1.0
        next_axis = normalize_vector(next_axis)
        delta = min(norm(vector_sub(next_axis, axis)), norm(vector_add(next_axis, axis)))
        axis = next_axis
        if delta < tolerance:
            break
    return _apply_axis_sign_convention(axis)


def project_point_to_axis(point: Coordinate, origin: Coordinate, axis: Coordinate) -> tuple[float, Coordinate]:
    normalized_axis = normalize_vector(axis)
    offset = vector_sub(point, origin)
    axial_t = dot(offset, normalized_axis)
    projected = vector_add(origin, scalar_mul(axial_t, normalized_axis))
    return axial_t, projected


def build_perpendicular_basis(axis: Coordinate) -> tuple[Coordinate, Coordinate]:
    normalized_axis = normalize_vector(axis)
    if abs(normalized_axis[2]) < 0.9:
        reference: Coordinate = (0.0, 0.0, 1.0)
    else:
        reference = (1.0, 0.0, 0.0)
    basis_u = normalize_vector(cross(reference, normalized_axis))
    basis_v = normalize_vector(cross(normalized_axis, basis_u))
    if basis_u[0] < -1e-12 or (abs(basis_u[0]) <= 1e-12 and basis_u[1] < -1e-12):
        basis_u = -basis_u[0], -basis_u[1], -basis_u[2]
        basis_v = -basis_v[0], -basis_v[1], -basis_v[2]
    return basis_u, basis_v


def centroid_of_atoms(atoms: Iterable[PDBAtom]) -> Coordinate | None:
    atom_list = list(atoms)
    if not atom_list:
        return None
    count = len(atom_list)
    return (
        sum(atom.x for atom in atom_list) / count,
        sum(atom.y for atom in atom_list) / count,
        sum(atom.z for atom in atom_list) / count,
    )


def residue_key(atom: PDBAtom) -> ResidueKey:
    return atom.chain_id, atom.residue_name, atom.residue_number, atom.insertion_code


def group_atoms_by_residue(atoms: Iterable[PDBAtom]) -> OrderedDict[ResidueKey, list[PDBAtom]]:
    grouped: OrderedDict[ResidueKey, list[PDBAtom]] = OrderedDict()
    for atom in atoms:
        key = residue_key(atom)
        grouped.setdefault(key, []).append(atom)
    return grouped


def residue_centroid(atoms_for_residue: Iterable[PDBAtom]) -> Coordinate | None:
    return centroid_of_atoms(atoms_for_residue)


def min_inter_residue_distance(
    res_atoms_a: Iterable[PDBAtom],
    res_atoms_b: Iterable[PDBAtom],
    heavy_only: bool = False,
) -> tuple[float, PDBAtom, PDBAtom] | None:
    atoms_a = list(res_atoms_a)
    atoms_b = list(res_atoms_b)
    if heavy_only:
        atoms_a = heavy_atoms(atoms_a)
        atoms_b = heavy_atoms(atoms_b)
    if not atoms_a or not atoms_b:
        return None

    best: tuple[float, PDBAtom, PDBAtom] | None = None
    for atom_a in atoms_a:
        for atom_b in atoms_b:
            atom_distance = distance(atom_a, atom_b)
            if best is None or atom_distance < best[0]:
                best = (atom_distance, atom_a, atom_b)
    return best
