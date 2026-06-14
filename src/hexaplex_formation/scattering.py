"""Simplified scattering helpers for comparative radial profiles."""

from __future__ import annotations

import math
from typing import Callable, Iterable, TypedDict

from hexaplex_formation.pdb_utils import PDBAtom


class DebyeDistanceHistogram(TypedDict):
    self_term: float
    bins: dict[float, float]


def q_from_d(d_A: float) -> float:
    if d_A <= 0:
        raise ValueError("d-spacing must be greater than zero")
    return 2.0 * math.pi / d_A


def d_from_q(q_Ainv: float) -> float:
    if q_Ainv <= 0:
        raise ValueError("q must be greater than zero")
    return 2.0 * math.pi / q_Ainv


def make_q_grid(q_min: float, q_max: float, q_step: float) -> list[float]:
    if q_min <= 0 or q_max <= 0 or q_step <= 0:
        raise ValueError("q_min, q_max, and q_step must be greater than zero")
    if q_max < q_min:
        raise ValueError("q_max must be greater than or equal to q_min")
    values: list[float] = []
    index = 0
    while True:
        value = q_min + index * q_step
        if value > q_max + q_step * 0.5:
            break
        values.append(round(value, 10))
        index += 1
    return values


def _atom_weight(atom: PDBAtom, element_weights: dict[str, float] | None) -> float:
    if element_weights is None:
        return 1.0
    return element_weights.get(atom.element.upper(), 1.0)


def _pair_distances_and_weights(
    atoms: list[PDBAtom],
    element_weights: dict[str, float] | None,
) -> tuple[float, list[tuple[float, float]]]:
    weights = [_atom_weight(atom, element_weights) for atom in atoms]
    self_intensity = sum(weight * weight for weight in weights)
    pairs: list[tuple[float, float]] = []
    for i, atom_i in enumerate(atoms):
        weight_i = weights[i]
        for j in range(i + 1, len(atoms)):
            atom_j = atoms[j]
            dx = atom_i.x - atom_j.x
            dy = atom_i.y - atom_j.y
            dz = atom_i.z - atom_j.z
            distance = math.sqrt(dx * dx + dy * dy + dz * dz)
            pairs.append((distance, 2.0 * weight_i * weights[j]))
    return self_intensity, pairs


def pair_distance_histogram_for_debye(
    atoms: Iterable[PDBAtom],
    bin_width: float = 0.05,
    max_distance: float | None = None,
    element_weights: dict[str, float] | None = None,
) -> DebyeDistanceHistogram:
    if bin_width <= 0:
        raise ValueError("bin_width must be greater than zero")
    if max_distance is not None and max_distance <= 0:
        raise ValueError("max_distance must be greater than zero")

    atom_list = list(atoms)
    weights = [_atom_weight(atom, element_weights) for atom in atom_list]
    bins: dict[float, float] = {}
    self_term = sum(weight * weight for weight in weights)
    decimals = max(0, int(math.ceil(-math.log10(bin_width))) + 3)

    for i, atom_i in enumerate(atom_list):
        weight_i = weights[i]
        for j in range(i + 1, len(atom_list)):
            atom_j = atom_list[j]
            dx = atom_i.x - atom_j.x
            dy = atom_i.y - atom_j.y
            dz = atom_i.z - atom_j.z
            distance = math.sqrt(dx * dx + dy * dy + dz * dz)
            if max_distance is not None and distance > max_distance:
                continue
            bin_center = round(round(distance / bin_width) * bin_width, decimals)
            bins[bin_center] = bins.get(bin_center, 0.0) + 2.0 * weight_i * weights[j]

    return {"self_term": self_term, "bins": bins}


def _sinc(value: float) -> float:
    if abs(value) < 1e-12:
        return 1.0
    return math.sin(value) / value


def debye_intensity_from_distance_histogram(
    histogram: DebyeDistanceHistogram,
    q_values: Iterable[float],
) -> list[float]:
    bins = histogram["bins"]
    self_term = histogram["self_term"]
    intensities: list[float] = []
    for q_value in q_values:
        intensity = self_term
        for distance, pair_weight in bins.items():
            intensity += pair_weight * _sinc(q_value * distance)
        intensities.append(intensity)
    return intensities


def _default_sampling_key(atom: PDBAtom) -> tuple[str, int | None, str, str]:
    return (atom.chain_id, atom.residue_number, atom.insertion_code, atom.residue_name)


def _even_positions(length: int, count: int) -> list[int]:
    if count >= length:
        return list(range(length))
    return [min(length - 1, int((index + 0.5) * length / count)) for index in range(count)]


def stratified_sample_atoms(
    atoms: Iterable[PDBAtom],
    max_atoms: int,
    key_fn: Callable[[PDBAtom], object] | None = None,
) -> list[PDBAtom]:
    if max_atoms <= 0:
        raise ValueError("max_atoms must be greater than zero")

    atom_list = list(atoms)
    if len(atom_list) <= max_atoms:
        return atom_list

    key_fn = key_fn or _default_sampling_key
    grouped: dict[object, list[tuple[int, PDBAtom]]] = {}
    group_keys: list[object] = []
    for atom_index, atom in enumerate(atom_list):
        key = key_fn(atom)
        if key not in grouped:
            grouped[key] = []
            group_keys.append(key)
        grouped[key].append((atom_index, atom))

    if max_atoms < len(group_keys):
        selected_group_positions = set(_even_positions(len(group_keys), max_atoms))
        selected: list[tuple[int, PDBAtom]] = []
        for group_position, key in enumerate(group_keys):
            if group_position not in selected_group_positions:
                continue
            group_atoms = grouped[key]
            selected.append(group_atoms[_even_positions(len(group_atoms), 1)[0]])
        return [atom for _index, atom in sorted(selected, key=lambda item: item[0])]

    quotas: dict[object, int] = {key: 1 for key in group_keys}
    remaining = max_atoms - len(group_keys)
    fractional: list[tuple[float, int, object]] = []
    for group_position, key in enumerate(group_keys):
        group_size = len(grouped[key])
        extra_exact = remaining * group_size / len(atom_list)
        extra_floor = min(group_size - 1, int(math.floor(extra_exact)))
        quotas[key] += extra_floor
        fractional.append((extra_exact - extra_floor, group_position, key))

    assigned = sum(quotas.values())
    for _fraction, _group_position, key in sorted(fractional, reverse=True):
        if assigned >= max_atoms:
            break
        if quotas[key] < len(grouped[key]):
            quotas[key] += 1
            assigned += 1

    selected = []
    for key in group_keys:
        group_atoms = grouped[key]
        for position in _even_positions(len(group_atoms), quotas[key]):
            selected.append(group_atoms[position])
    return [atom for _index, atom in sorted(selected, key=lambda item: item[0])]


def debye_intensity(
    atoms: Iterable[PDBAtom],
    q_values: Iterable[float],
    element_weights: dict[str, float] | None = None,
) -> list[float]:
    atom_list = list(atoms)
    self_intensity, pairs = _pair_distances_and_weights(atom_list, element_weights)
    intensities: list[float] = []
    for q_value in q_values:
        if q_value == 0:
            total_weight = sum(_atom_weight(atom, element_weights) for atom in atom_list)
            intensities.append(total_weight * total_weight)
            continue
        intensity = self_intensity
        for distance, pair_weight in pairs:
            qr = q_value * distance
            intensity += pair_weight * _sinc(qr)
        intensities.append(intensity)
    return intensities


def integrate_window(
    profile_rows: Iterable[dict[str, str | float]],
    d_center: float | None = None,
    d_min: float | None = None,
    d_max: float | None = None,
) -> dict[str, str]:
    if d_center is not None and (d_min is None or d_max is None):
        d_min = d_center
        d_max = d_center
    if d_min is None or d_max is None:
        raise ValueError("provide d_center or d_min/d_max")
    if d_max < d_min:
        raise ValueError("d_max must be greater than or equal to d_min")

    selected: list[float] = []
    for row in profile_rows:
        d_value = float(row["d_A"])
        if d_min <= d_value <= d_max:
            selected.append(float(row["intensity"]))

    if not selected:
        return {
            "point_count": "0",
            "mean_intensity": "",
            "max_intensity": "",
            "integrated_intensity": "0.000000",
        }
    integrated = sum(selected)
    return {
        "point_count": str(len(selected)),
        "mean_intensity": f"{integrated / len(selected):.6f}",
        "max_intensity": f"{max(selected):.6f}",
        "integrated_intensity": f"{integrated:.6f}",
    }
