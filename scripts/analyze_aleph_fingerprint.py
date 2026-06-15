#!/usr/bin/env python3
"""Prototype Aleph geometric fingerprint for hexaplex structures."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from hexaplex_formation.geometry import (  # noqa: E402
    build_perpendicular_basis,
    covariance_matrix_3d,
    cross,
    distance,
    dot,
    group_atoms_by_residue,
    mean_point,
    normalize_vector,
    power_iteration_principal_axis,
    project_point_to_axis,
    vector_sub,
)
from hexaplex_formation.pdb_utils import PDBAtom, chain_ids, dedupe_exact_atoms, heavy_atoms, load_pdb_atoms  # noqa: E402


BASE_RESIDUES = {"CYP", "MEP"}
BACKBONE_LIKE_ATOMS = {"N", "CA", "C", "O", "OXT"}
DEFAULT_STRUCTURES = [
    ("full", Path("outputs/intermediates/ai_candidate_inputs/full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain.pdb")),
    ("central6", Path("outputs/mini_hexaplex/structures/mini_hexaplex_central6_units.pdb")),
    ("central7", Path("outputs/mini_hexaplex/structures/mini_hexaplex_central7_units.pdb")),
]

PER_UNIT_COLUMNS = [
    "structure_id",
    "source_pdb",
    "unit_index",
    "chain_count",
    "axial_position_A",
    "layer_centroid_x",
    "layer_centroid_y",
    "layer_centroid_z",
    "base_mean_radius_A",
    "base_radial_spread_A",
    "aleph_phase_raw_deg",
    "aleph_phase_deg",
    "local_twist_deg",
    "local_abs_twist_deg",
    "local_rise_A",
    "base_plane_normal_x",
    "base_plane_normal_y",
    "base_plane_normal_z",
    "aleph_base_plane_bend_deg",
    "scaffold_plane_normal_x",
    "scaffold_plane_normal_y",
    "scaffold_plane_normal_z",
    "aleph_scaffold_plane_bend_deg",
    "chain_angular_dispersion_deg",
    "chain_mean_resultant_length",
    "axis_flipped",
    "phase_unwrap_ok",
    "missing_chain_count",
    "local_twist_warning",
    "local_rise_warning",
    "plane_fit_warning",
    "warnings",
]

SUMMARY_COLUMNS = [
    "structure_id",
    "source_pdb",
    "unit_count",
    "chain_count",
    "mean_local_twist_deg",
    "std_local_twist_deg",
    "mean_abs_local_twist_deg",
    "std_abs_local_twist_deg",
    "mean_local_rise_A",
    "std_local_rise_A",
    "mean_base_radial_spread_A",
    "std_base_radial_spread_A",
    "mean_base_plane_bend_deg",
    "mean_scaffold_plane_bend_deg",
    "mean_chain_angular_dispersion_deg",
    "mean_chain_resultant_length",
    "aleph_regular_score",
    "axis_flipped",
    "phase_unwrap_ok",
    "warnings",
]

FFT_COLUMNS = [
    "structure_id",
    "signal_name",
    "sample_count",
    "dominant_frequency_index",
    "dominant_amplitude",
    "normalized_dominant_amplitude",
    "too_short_for_reliable_interpretation",
    "warnings",
]

QC_COLUMNS = [
    "structure_id",
    "source_pdb",
    "axis_flipped",
    "phase_unwrap_ok",
    "expected_unit_count",
    "observed_unit_count",
    "expected_chain_count",
    "observed_chain_count",
    "units_with_missing_chains",
    "missing_chain_count_by_unit",
    "mean_abs_local_twist_deg",
    "std_abs_local_twist_deg",
    "negative_rise_count",
    "local_twist_warning",
    "local_rise_warning",
    "plane_fit_warning",
    "warnings",
]


@dataclass(frozen=True)
class ResidueRecord:
    chain_id: str
    residue_name: str
    residue_number: int | None
    insertion_code: str
    atoms: tuple[PDBAtom, ...]


@dataclass(frozen=True)
class RepeatUnit:
    chain_id: str
    unit_index: int
    base_residue: ResidueRecord
    scaffold_residue: ResidueRecord | None


@dataclass(frozen=True)
class StructureFingerprint:
    structure_id: str
    path: Path
    raw_atom_count: int
    deduped_atom_count: int
    heavy_atom_count: int
    chain_count: int
    per_unit_rows: list[dict[str, str]]
    summary_row: dict[str, str]
    qc_row: dict[str, str]
    warnings: tuple[str, ...]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--per-unit-csv", type=Path, default=Path("outputs/metrics/aleph_fingerprint_per_unit.csv"))
    parser.add_argument("--summary-csv", type=Path, default=Path("outputs/metrics/aleph_fingerprint_summary.csv"))
    parser.add_argument("--fft-csv", type=Path, default=Path("outputs/metrics/aleph_fingerprint_fft_summary.csv"))
    parser.add_argument("--qc-csv", type=Path, default=Path("outputs/metrics/aleph_fingerprint_qc.csv"))
    parser.add_argument("--report", type=Path, default=Path("outputs/reports/aleph_fingerprint_report.md"))
    parser.add_argument("--plot-dir", type=Path, default=Path("outputs/plots/aleph_fingerprint"))
    parser.add_argument("--include-optional-twists", action="store_true", default=True)
    return parser.parse_args()


def atom_xyz(atom: PDBAtom) -> tuple[float, float, float]:
    return atom.x, atom.y, atom.z


def format_float(value: float | None, digits: int = 6) -> str:
    if value is None or not math.isfinite(value):
        return ""
    return f"{value:.{digits}f}"


def safe_float(value: object) -> float | None:
    try:
        parsed = float(str(value))
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def mean(values: list[float]) -> float | None:
    return sum(values) / len(values) if values else None


def sample_std(values: list[float]) -> float | None:
    if not values:
        return None
    if len(values) == 1:
        return 0.0
    avg = sum(values) / len(values)
    return math.sqrt(sum((value - avg) ** 2 for value in values) / (len(values) - 1))


def atom_class(atom: PDBAtom) -> str:
    if atom.atom_name.strip().upper() in BACKBONE_LIKE_ATOMS:
        return "backbone_like"
    if atom.residue_name in BASE_RESIDUES:
        return "base_like"
    if atom.residue_name == "GLU":
        return "scaffold_linker"
    return "other"


def infer_axis(points: list[tuple[float, float, float]]) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    if len(points) < 3:
        raise ValueError("At least three points are required for Aleph axis fitting")
    origin = mean_point(points)
    axis = power_iteration_principal_axis(covariance_matrix_3d(points))
    return origin, axis


def residue_records(atoms: list[PDBAtom]) -> OrderedDict[str, list[ResidueRecord]]:
    grouped: OrderedDict[str, list[ResidueRecord]] = OrderedDict()
    for (chain_id, residue_name, residue_number, insertion_code), residue_atoms in group_atoms_by_residue(atoms).items():
        grouped.setdefault(chain_id, []).append(
            ResidueRecord(
                chain_id=chain_id,
                residue_name=residue_name,
                residue_number=residue_number,
                insertion_code=insertion_code,
                atoms=tuple(residue_atoms),
            )
        )
    return grouped


def repeat_units_by_chain(atoms: list[PDBAtom]) -> tuple[dict[str, list[RepeatUnit]], list[str]]:
    warnings: list[str] = []
    units_by_chain: dict[str, list[RepeatUnit]] = {}
    for chain_id, residues in residue_records(atoms).items():
        units: list[RepeatUnit] = []
        base_records = [(index, residue) for index, residue in enumerate(residues) if residue.residue_name in BASE_RESIDUES]
        for unit_index, (residue_index, base_residue) in enumerate(base_records, start=1):
            scaffold = residues[residue_index + 1] if residue_index + 1 < len(residues) and residues[residue_index + 1].residue_name == "GLU" else None
            if scaffold is None:
                warnings.append(f"chain {chain_id} unit {unit_index} has no following GLU scaffold residue")
            units.append(RepeatUnit(chain_id=chain_id, unit_index=unit_index, base_residue=base_residue, scaffold_residue=scaffold))
        units_by_chain[chain_id] = units
    return units_by_chain, warnings


def centroid(atoms: list[PDBAtom] | tuple[PDBAtom, ...]) -> tuple[float, float, float] | None:
    selected = list(atoms)
    if not selected:
        return None
    return mean_point(atom_xyz(atom) for atom in selected)


def radius_from_axis(point: tuple[float, float, float], origin: tuple[float, float, float], axis: tuple[float, float, float]) -> float:
    _, projected = project_point_to_axis(point, origin, axis)
    return distance(point, projected)


def angle_about_axis(
    point: tuple[float, float, float],
    origin: tuple[float, float, float],
    axis: tuple[float, float, float],
    basis_u: tuple[float, float, float],
    basis_v: tuple[float, float, float],
) -> float:
    _, projected = project_point_to_axis(point, origin, axis)
    radial = vector_sub(point, projected)
    return math.atan2(dot(radial, basis_v), dot(radial, basis_u))


def unwrap_radians(values: list[float | None]) -> list[float | None]:
    unwrapped: list[float | None] = []
    previous: float | None = None
    offset = 0.0
    for value in values:
        if value is None:
            unwrapped.append(None)
            continue
        adjusted = value + offset
        if previous is not None:
            while adjusted - previous > math.pi:
                offset -= 2.0 * math.pi
                adjusted = value + offset
            while adjusted - previous < -math.pi:
                offset += 2.0 * math.pi
                adjusted = value + offset
        unwrapped.append(adjusted)
        previous = adjusted
    return unwrapped


def phase_unwrap_ok(raw_values: list[float | None], unwrapped_values: list[float | None]) -> bool:
    raw_count = sum(value is not None for value in raw_values)
    unwrapped_count = sum(value is not None for value in unwrapped_values)
    if raw_count != unwrapped_count:
        return False
    finite = [value for value in unwrapped_values if value is not None and math.isfinite(value)]
    return len(finite) == unwrapped_count


def wrap_degrees(value: float) -> float:
    wrapped = (value + 180.0) % 360.0 - 180.0
    if wrapped == -180.0:
        return 180.0
    return wrapped


def circular_mean_resultant_length(angles: list[float]) -> float | None:
    if not angles:
        return None
    sin_mean = sum(math.sin(value) for value in angles) / len(angles)
    cos_mean = sum(math.cos(value) for value in angles) / len(angles)
    return min(1.0, max(0.0, math.sqrt(sin_mean * sin_mean + cos_mean * cos_mean)))


def circular_dispersion_deg(angles: list[float]) -> float | None:
    resultant = circular_mean_resultant_length(angles)
    if resultant is None:
        return None
    if resultant <= 0:
        return 180.0
    return min(180.0, math.degrees(math.sqrt(max(0.0, -2.0 * math.log(resultant)))))


def orient_axis_by_unit_order(
    axis: tuple[float, float, float],
    origin: tuple[float, float, float],
    anchor_points: list[tuple[float, float, float] | None],
) -> tuple[tuple[float, float, float], bool]:
    positions = [project_point_to_axis(point, origin, axis)[0] for point in anchor_points if point is not None]
    deltas = [positions[index + 1] - positions[index] for index in range(len(positions) - 1)]
    valid = [delta for delta in deltas if abs(delta) > 1e-9]
    if valid and sum(valid) / len(valid) < 0.0:
        return (-axis[0], -axis[1], -axis[2]), True
    return axis, False


def plane_normal(points: list[tuple[float, float, float]]) -> tuple[float, float, float] | None:
    if len(points) < 3:
        return None
    center = mean_point(points)
    centered = [vector_sub(point, center) for point in points]
    best = None
    best_length = 0.0
    for i, first in enumerate(centered):
        for second in centered[i + 1 :]:
            candidate = cross(first, second)
            length = math.sqrt(dot(candidate, candidate))
            if length > best_length:
                best = candidate
                best_length = length
    if best is None or best_length == 0:
        return None
    return normalize_vector(best)


def normal_angle_deg(a: tuple[float, float, float] | None, b: tuple[float, float, float] | None) -> float | None:
    if a is None or b is None:
        return None
    value = max(-1.0, min(1.0, abs(dot(a, b))))
    return math.degrees(math.acos(value))


def find_optional_twists() -> list[tuple[str, Path]]:
    patterns = {
        "compact24": ["outputs/**/compact_hexaplex_twist_24.pdb", "outputs/**/*twist_24*.pdb"],
        "compact30": ["outputs/**/compact_hexaplex_twist_30.pdb", "outputs/**/*twist_30*.pdb"],
        "compact36": ["outputs/**/compact_hexaplex_twist_36.pdb", "outputs/**/*twist_36*.pdb"],
    }
    found: list[tuple[str, Path]] = []
    for structure_id, globs in patterns.items():
        match = None
        for pattern in globs:
            matches = sorted(REPO_ROOT.glob(pattern))
            if matches:
                match = matches[0]
                break
        if match is not None:
            found.append((structure_id, match.relative_to(REPO_ROOT)))
    return found


def load_structure_inputs(include_optional_twists: bool = True) -> tuple[list[tuple[str, Path]], list[str]]:
    inputs = list(DEFAULT_STRUCTURES)
    warnings: list[str] = []
    optional = find_optional_twists() if include_optional_twists else []
    if optional:
        inputs.extend(optional)
    else:
        warnings.append("optional compact 24/30/36 variants were not found in this repo")
    return inputs, warnings


def analyze_structure(structure_id: str, path: Path) -> StructureFingerprint:
    raw_atoms = load_pdb_atoms(path)
    atoms = dedupe_exact_atoms(raw_atoms)
    selected_heavy = heavy_atoms(atoms)
    warnings: list[str] = []
    if len(raw_atoms) != len(atoms):
        warnings.append(f"removed {len(raw_atoms) - len(atoms)} exact duplicate atoms")
    units_by_chain, unit_warnings = repeat_units_by_chain(selected_heavy)
    warnings.extend(unit_warnings)
    unit_count = max((len(units) for units in units_by_chain.values()), default=0)
    chains = sorted(units_by_chain)
    expected_chain_count = len(chains)
    expected_unit_count = unit_count
    base_centroids = []
    for units in units_by_chain.values():
        for unit in units:
            atoms_for_base = [atom for atom in heavy_atoms(unit.base_residue.atoms) if atom_class(atom) == "base_like"]
            point = centroid(atoms_for_base) or centroid(heavy_atoms(unit.base_residue.atoms))
            if point is not None:
                base_centroids.append(point)
    preliminary_layer_centroids = []
    anchor_points_for_orientation: list[tuple[float, float, float] | None] = []
    for unit_index in range(1, unit_count + 1):
        layer_points = []
        for chain in chains:
            if len(units_by_chain[chain]) < unit_index:
                continue
            unit = units_by_chain[chain][unit_index - 1]
            atoms_for_base = [atom for atom in heavy_atoms(unit.base_residue.atoms) if atom_class(atom) == "base_like"]
            point = centroid(atoms_for_base) or centroid(heavy_atoms(unit.base_residue.atoms))
            if point is not None:
                layer_points.append(point)
        anchor_points_for_orientation.append(layer_points[0] if layer_points else None)
        if layer_points:
            preliminary_layer_centroids.append(mean_point(layer_points))
    axis_points = preliminary_layer_centroids if len(preliminary_layer_centroids) >= 3 else base_centroids
    origin, axis = infer_axis(axis_points if len(axis_points) >= 3 else [atom_xyz(atom) for atom in selected_heavy])
    axis, axis_flipped = orient_axis_by_unit_order(axis, origin, anchor_points_for_orientation)
    basis_u, basis_v = build_perpendicular_basis(axis)

    raw_layers = []
    for unit_index in range(1, unit_count + 1):
        chain_units = [units_by_chain[chain][unit_index - 1] for chain in chains if len(units_by_chain[chain]) >= unit_index]
        base_points = []
        scaffold_points = []
        chain_angles = []
        for unit in chain_units:
            base_atoms = [atom for atom in heavy_atoms(unit.base_residue.atoms) if atom_class(atom) == "base_like"]
            base_point = centroid(base_atoms) or centroid(heavy_atoms(unit.base_residue.atoms))
            if base_point is not None:
                base_points.append(base_point)
                chain_angles.append(angle_about_axis(base_point, origin, axis, basis_u, basis_v))
            if unit.scaffold_residue is not None:
                scaffold_atoms = [atom for atom in heavy_atoms(unit.scaffold_residue.atoms) if atom_class(atom) == "scaffold_linker"]
                scaffold_point = centroid(scaffold_atoms) or centroid(heavy_atoms(unit.scaffold_residue.atoms))
                if scaffold_point is not None:
                    scaffold_points.append(scaffold_point)
        layer_centroid = mean_point(base_points) if base_points else None
        anchor_point = base_points[0] if base_points else None
        phase = angle_about_axis(anchor_point, origin, axis, basis_u, basis_v) if anchor_point is not None else None
        axial_position = project_point_to_axis(anchor_point, origin, axis)[0] if anchor_point is not None else None
        radii = [radius_from_axis(point, origin, axis) for point in base_points]
        raw_layers.append(
            {
                "unit_index": unit_index,
                "chain_count": len(chain_units),
                "missing_chain_count": max(0, expected_chain_count - len(chain_units)),
                "layer_centroid": layer_centroid,
                "axial_position": axial_position,
                "base_mean_radius": mean(radii),
                "base_radial_spread": sample_std(radii),
                "phase": phase,
                "base_plane_normal": plane_normal(base_points),
                "scaffold_plane_normal": plane_normal(scaffold_points),
                "chain_angular_dispersion": circular_dispersion_deg(chain_angles),
                "chain_mean_resultant_length": circular_mean_resultant_length(chain_angles),
                "warnings": "",
            }
        )

    raw_phases = [layer["phase"] for layer in raw_layers]
    unwrapped_phase = unwrap_radians(raw_phases)
    unwrap_ok = phase_unwrap_ok(raw_phases, unwrapped_phase)
    missing_chain_count_by_unit = "; ".join(f"{layer['unit_index']}:{layer['missing_chain_count']}" for layer in raw_layers if layer["missing_chain_count"])
    units_with_missing_chains = sum(1 for layer in raw_layers if layer["missing_chain_count"])
    per_unit_rows: list[dict[str, str]] = []
    twist_values: list[float] = []
    abs_twist_values: list[float] = []
    rise_values: list[float] = []
    spread_values: list[float] = []
    base_bends: list[float] = []
    scaffold_bends: list[float] = []
    dispersions: list[float] = []
    resultants: list[float] = []
    negative_rise_count = 0
    local_twist_warnings: list[str] = []
    local_rise_warnings: list[str] = []
    plane_fit_warnings: list[str] = []
    for index, layer in enumerate(raw_layers):
        next_layer = raw_layers[index + 1] if index + 1 < len(raw_layers) else None
        local_twist = None
        local_abs_twist = None
        if next_layer is not None and unwrapped_phase[index] is not None and unwrapped_phase[index + 1] is not None:
            local_twist = math.degrees(unwrapped_phase[index + 1] - unwrapped_phase[index])
            wrapped_twist = wrap_degrees(local_twist)
            local_abs_twist = abs(wrapped_twist)
            twist_values.append(local_twist)
            abs_twist_values.append(local_abs_twist)
            if abs(local_twist) > 180.0:
                local_twist_warnings.append(f"unit {layer['unit_index']} to {next_layer['unit_index']} has large unwrapped twist {local_twist:.2f} deg")
        local_rise = None
        if next_layer is not None and layer["axial_position"] is not None and next_layer["axial_position"] is not None:
            local_rise = next_layer["axial_position"] - layer["axial_position"]
            rise_values.append(local_rise)
            if local_rise < -1e-6:
                negative_rise_count += 1
                local_rise_warnings.append(f"unit {layer['unit_index']} to {next_layer['unit_index']} has negative rise {local_rise:.2f} A")
        base_bend = normal_angle_deg(layer["base_plane_normal"], next_layer["base_plane_normal"]) if next_layer else None
        if base_bend is not None:
            base_bends.append(base_bend)
        elif next_layer is not None:
            plane_fit_warnings.append(f"unit {layer['unit_index']} to {next_layer['unit_index']} missing base plane fit")
        scaffold_bend = normal_angle_deg(layer["scaffold_plane_normal"], next_layer["scaffold_plane_normal"]) if next_layer else None
        if scaffold_bend is not None:
            scaffold_bends.append(scaffold_bend)
        elif next_layer is not None:
            plane_fit_warnings.append(f"unit {layer['unit_index']} to {next_layer['unit_index']} missing scaffold plane fit")
        if layer["base_radial_spread"] is not None:
            spread_values.append(layer["base_radial_spread"])
        if layer["chain_angular_dispersion"] is not None:
            dispersions.append(layer["chain_angular_dispersion"])
        if layer["chain_mean_resultant_length"] is not None:
            resultants.append(layer["chain_mean_resultant_length"])
        centroid_point = layer["layer_centroid"]
        base_normal = layer["base_plane_normal"]
        scaffold_normal = layer["scaffold_plane_normal"]
        per_unit_rows.append(
            {
                "structure_id": structure_id,
                "source_pdb": str(path),
                "unit_index": str(layer["unit_index"]),
                "chain_count": str(layer["chain_count"]),
                "axial_position_A": format_float(layer["axial_position"]),
                "layer_centroid_x": format_float(centroid_point[0] if centroid_point else None),
                "layer_centroid_y": format_float(centroid_point[1] if centroid_point else None),
                "layer_centroid_z": format_float(centroid_point[2] if centroid_point else None),
                "base_mean_radius_A": format_float(layer["base_mean_radius"]),
                "base_radial_spread_A": format_float(layer["base_radial_spread"]),
                "aleph_phase_raw_deg": format_float(math.degrees(layer["phase"]) if layer["phase"] is not None else None),
                "aleph_phase_deg": format_float(math.degrees(unwrapped_phase[index]) if unwrapped_phase[index] is not None else None),
                "local_twist_deg": format_float(local_twist),
                "local_abs_twist_deg": format_float(local_abs_twist),
                "local_rise_A": format_float(local_rise),
                "base_plane_normal_x": format_float(base_normal[0] if base_normal else None),
                "base_plane_normal_y": format_float(base_normal[1] if base_normal else None),
                "base_plane_normal_z": format_float(base_normal[2] if base_normal else None),
                "aleph_base_plane_bend_deg": format_float(base_bend),
                "scaffold_plane_normal_x": format_float(scaffold_normal[0] if scaffold_normal else None),
                "scaffold_plane_normal_y": format_float(scaffold_normal[1] if scaffold_normal else None),
                "scaffold_plane_normal_z": format_float(scaffold_normal[2] if scaffold_normal else None),
                "aleph_scaffold_plane_bend_deg": format_float(scaffold_bend),
                "chain_angular_dispersion_deg": format_float(layer["chain_angular_dispersion"]),
                "chain_mean_resultant_length": format_float(layer["chain_mean_resultant_length"]),
                "axis_flipped": "true" if axis_flipped else "false",
                "phase_unwrap_ok": "true" if unwrap_ok else "false",
                "missing_chain_count": str(layer["missing_chain_count"]),
                "local_twist_warning": "; ".join(warning for warning in local_twist_warnings if f"unit {layer['unit_index']} " in warning),
                "local_rise_warning": "; ".join(warning for warning in local_rise_warnings if f"unit {layer['unit_index']} " in warning),
                "plane_fit_warning": "; ".join(warning for warning in plane_fit_warnings if f"unit {layer['unit_index']} " in warning),
                "warnings": layer["warnings"],
            }
        )

    twist_std = sample_std(twist_values)
    abs_twist_std = sample_std(abs_twist_values)
    rise_std = sample_std(rise_values)
    spread_std = sample_std(spread_values)
    local_twist_warning = "; ".join(local_twist_warnings)
    if twist_values and mean(abs_twist_values) is not None and abs((mean(abs_twist_values) or 0.0) - 30.0) > 10.0:
        local_twist_warning = f"{local_twist_warning}; mean absolute local twist differs from nominal 30 deg by >10 deg".strip("; ")
    local_rise_warning = "; ".join(local_rise_warnings)
    if negative_rise_count:
        local_rise_warning = f"{local_rise_warning}; {negative_rise_count} negative local rise values after axis orientation".strip("; ")
    plane_fit_warning = "; ".join(plane_fit_warnings)
    regular_score = None
    if twist_std is not None and rise_std is not None and spread_std is not None:
        regular_score = 1.0 / (1.0 + abs(twist_std) / 30.0 + abs(rise_std) / 3.4 + abs(spread_std))
    summary_row = {
        "structure_id": structure_id,
        "source_pdb": str(path),
        "unit_count": str(unit_count),
        "chain_count": str(len(chains)),
        "mean_local_twist_deg": format_float(mean(twist_values)),
        "std_local_twist_deg": format_float(twist_std),
        "mean_abs_local_twist_deg": format_float(mean(abs_twist_values)),
        "std_abs_local_twist_deg": format_float(abs_twist_std),
        "mean_local_rise_A": format_float(mean(rise_values)),
        "std_local_rise_A": format_float(rise_std),
        "mean_base_radial_spread_A": format_float(mean(spread_values)),
        "std_base_radial_spread_A": format_float(spread_std),
        "mean_base_plane_bend_deg": format_float(mean(base_bends)),
        "mean_scaffold_plane_bend_deg": format_float(mean(scaffold_bends)),
        "mean_chain_angular_dispersion_deg": format_float(mean(dispersions)),
        "mean_chain_resultant_length": format_float(mean(resultants)),
        "aleph_regular_score": format_float(regular_score),
        "axis_flipped": "true" if axis_flipped else "false",
        "phase_unwrap_ok": "true" if unwrap_ok else "false",
        "warnings": "; ".join(warnings),
    }
    qc_row = {
        "structure_id": structure_id,
        "source_pdb": str(path),
        "axis_flipped": "true" if axis_flipped else "false",
        "phase_unwrap_ok": "true" if unwrap_ok else "false",
        "expected_unit_count": str(expected_unit_count),
        "observed_unit_count": str(len(raw_layers)),
        "expected_chain_count": str(expected_chain_count),
        "observed_chain_count": str(len(chains)),
        "units_with_missing_chains": str(units_with_missing_chains),
        "missing_chain_count_by_unit": missing_chain_count_by_unit,
        "mean_abs_local_twist_deg": format_float(mean(abs_twist_values)),
        "std_abs_local_twist_deg": format_float(abs_twist_std),
        "negative_rise_count": str(negative_rise_count),
        "local_twist_warning": local_twist_warning,
        "local_rise_warning": local_rise_warning,
        "plane_fit_warning": plane_fit_warning,
        "warnings": "; ".join(warnings),
    }
    return StructureFingerprint(
        structure_id=structure_id,
        path=path,
        raw_atom_count=len(raw_atoms),
        deduped_atom_count=len(atoms),
        heavy_atom_count=len(selected_heavy),
        chain_count=len(chains),
        per_unit_rows=per_unit_rows,
        summary_row=summary_row,
        qc_row=qc_row,
        warnings=tuple(warnings),
    )


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def fft_summary(structure_id: str, signal_name: str, values: list[float]) -> dict[str, str]:
    warning = ""
    if len(values) < 4:
        return {
            "structure_id": structure_id,
            "signal_name": signal_name,
            "sample_count": str(len(values)),
            "dominant_frequency_index": "",
            "dominant_amplitude": "",
            "normalized_dominant_amplitude": "",
            "too_short_for_reliable_interpretation": "true",
            "warnings": "signal is too short for discrete spectral interpretation",
        }
    centered = [value - (sum(values) / len(values)) for value in values]
    try:
        import numpy as np

        arr = np.asarray(centered, dtype=float)
        spectrum = np.fft.rfft(arr)
        amplitudes = [float(value) for value in np.abs(spectrum)]
    except Exception:
        warning = "numpy unavailable; used standard-library DFT fallback"
        amplitudes = []
        count = len(centered)
        for frequency_index in range(count // 2 + 1):
            real = 0.0
            imag = 0.0
            for sample_index, value in enumerate(centered):
                angle = -2.0 * math.pi * frequency_index * sample_index / count
                real += value * math.cos(angle)
                imag += value * math.sin(angle)
            amplitudes.append(math.sqrt(real * real + imag * imag))
    if len(amplitudes) <= 1:
        dominant_index = 0
        dominant_amplitude = 0.0
    else:
        search = amplitudes[1:]
        dominant_index = max(range(1, len(amplitudes)), key=lambda index: amplitudes[index])
        dominant_amplitude = float(amplitudes[dominant_index])
    total = float(sum(amplitudes[1:]))
    normalized = dominant_amplitude / total if total > 0 else 0.0
    if len(values) < 10:
        warning = f"{warning}; short signal; spectral interpretation is limited" if warning else "short signal; spectral interpretation is limited"
    return {
        "structure_id": structure_id,
        "signal_name": signal_name,
        "sample_count": str(len(values)),
        "dominant_frequency_index": str(dominant_index),
        "dominant_amplitude": format_float(dominant_amplitude),
        "normalized_dominant_amplitude": format_float(normalized),
        "too_short_for_reliable_interpretation": "true" if len(values) < 10 else "false",
        "warnings": warning,
    }


def collect_signal(rows: list[dict[str, str]], structure_id: str, column: str) -> list[float]:
    values = []
    for row in rows:
        if row["structure_id"] != structure_id:
            continue
        value = safe_float(row.get(column))
        if value is not None:
            values.append(value)
    return values


FINGERPRINT_FEATURES = [
    ("local_abs_twist_deg", "abs twist", 0.0, 60.0, False),
    ("local_twist_deg", "signed twist", -60.0, 60.0, False),
    ("local_rise_A", "rise", -4.0, 6.0, False),
    ("base_radial_spread_A", "radial spread", 0.0, 4.0, False),
    ("aleph_phase_deg", "phase", None, None, True),
    ("aleph_base_plane_bend_deg", "base bend", 0.0, 90.0, False),
    ("aleph_scaffold_plane_bend_deg", "scaffold bend", 0.0, 90.0, False),
    ("chain_mean_resultant_length", "chain coherence", 0.0, 1.0, False),
]


def normalize_feature_value(value: object, minimum: float | None, maximum: float | None, cyclic: bool = False) -> float | None:
    parsed = safe_float(value)
    if parsed is None:
        return None
    if cyclic:
        parsed = parsed % 360.0
        minimum = 0.0
        maximum = 360.0
    if minimum is None or maximum is None or maximum == minimum:
        return 0.5
    return min(1.0, max(0.0, (parsed - minimum) / (maximum - minimum)))


def value_color(norm_value: float | None, warning: bool = False) -> str:
    if warning:
        return "#d62728"
    if norm_value is None:
        return "#e6e6e6"
    # Blue to white to orange, deliberately simple for SVG portability.
    if norm_value <= 0.5:
        t = norm_value / 0.5
        r = int(76 + (245 - 76) * t)
        g = int(120 + (245 - 120) * t)
        b = int(168 + (245 - 168) * t)
    else:
        t = (norm_value - 0.5) / 0.5
        r = int(245 + (245 - 245) * t)
        g = int(245 + (133 - 245) * t)
        b = int(245 + (24 - 245) * t)
    return f"#{r:02x}{g:02x}{b:02x}"


def row_has_qc_warning(row: dict[str, str]) -> bool:
    warning_fields = ["local_twist_warning", "local_rise_warning", "plane_fit_warning", "warnings"]
    if any(row.get(field, "").strip() for field in warning_fields):
        return True
    return safe_float(row.get("missing_chain_count")) not in (None, 0.0)


def qc_warning_count(row: dict[str, str]) -> int:
    count = 0
    for field in ["local_twist_warning", "local_rise_warning", "plane_fit_warning", "warnings"]:
        if row.get(field, "").strip():
            count += 1
    negative = safe_float(row.get("negative_rise_count")) or 0.0
    if negative > 0:
        count += int(negative)
    if row.get("phase_unwrap_ok") == "false":
        count += 1
    return count


def series_points(rows: list[dict[str, str]], structure_id: str, column: str) -> list[tuple[float, float, bool]]:
    points: list[tuple[float, float, bool]] = []
    for row in sorted([item for item in rows if item["structure_id"] == structure_id], key=lambda item: int(item["unit_index"])):
        value = safe_float(row.get(column))
        transition = safe_float(row.get("unit_index"))
        if value is None or transition is None:
            continue
        points.append((transition, value, row_has_qc_warning(row)))
    return points


def svg_series_plot(
    title: str,
    points_by_label: OrderedDict[str, list[tuple[float, float, bool]]],
    y_label: str,
    path: Path,
    reference_y: float | None = None,
    stacked: bool = False,
) -> None:
    if not points_by_label or not any(points for points in points_by_label.values()):
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    width = 920
    panel_h = 230 if stacked else 330
    top = 48
    left = 78
    right = 34
    gap = 36
    panel_count = len(points_by_label) if stacked else 1
    height = top + panel_count * panel_h + (panel_count - 1) * gap + 58
    colors = ["#4c78a8", "#f58518", "#54a24b", "#b279a2", "#e45756"]
    all_points = [point for points in points_by_label.values() for point in points]
    min_x = min(point[0] for point in all_points)
    max_x = max(point[0] for point in all_points)
    if max_x == min_x:
        max_x += 1.0
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width / 2:.0f}" y="26" text-anchor="middle" font-family="Arial" font-size="17">{title}</text>',
    ]
    plot_w = width - left - right

    def sx(value: float) -> float:
        return left + (value - min_x) / (max_x - min_x) * plot_w

    panels = list(points_by_label.items()) if stacked else [("combined", all_points)]
    for panel_index, (label, panel_points) in enumerate(panels):
        y0 = top + panel_index * (panel_h + gap)
        panel_values = [point[1] for point in panel_points]
        min_y = min(panel_values + ([reference_y] if reference_y is not None else []))
        max_y = max(panel_values + ([reference_y] if reference_y is not None else []))
        pad = max(1.0, (max_y - min_y) * 0.12)
        min_y -= pad
        max_y += pad
        if max_y == min_y:
            max_y += 1.0

        def sy(value: float) -> float:
            return y0 + 12 + (max_y - value) / (max_y - min_y) * (panel_h - 34)

        axis_bottom = y0 + panel_h - 22
        parts.append(f'<line x1="{left}" y1="{y0 + 12}" x2="{left}" y2="{axis_bottom}" stroke="#333"/>')
        parts.append(f'<line x1="{left}" y1="{axis_bottom}" x2="{width - right}" y2="{axis_bottom}" stroke="#333"/>')
        if stacked:
            parts.append(f'<text x="18" y="{y0 + 26}" font-family="Arial" font-size="12" font-weight="bold">{label}</text>')
        if reference_y is not None and min_y <= reference_y <= max_y:
            ref_y = sy(reference_y)
            parts.append(f'<line x1="{left}" y1="{ref_y:.2f}" x2="{width - right}" y2="{ref_y:.2f}" stroke="#777" stroke-dasharray="4 4"/>')
            parts.append(f'<text x="{width - right - 4}" y="{ref_y - 4:.2f}" text-anchor="end" font-family="Arial" font-size="10">30 deg reference</text>')
        if stacked:
            drawable = [(label, panel_points, colors[panel_index % len(colors)])]
        else:
            drawable = [(series_label, series_points, colors[index % len(colors)]) for index, (series_label, series_points) in enumerate(points_by_label.items())]
        for series_label, series_points_for_label, color in drawable:
            series_points_for_label = sorted(series_points_for_label)
            coords = " ".join(f"{sx(x):.2f},{sy(y):.2f}" for x, y, _ in series_points_for_label)
            if coords:
                parts.append(f'<polyline points="{coords}" fill="none" stroke="{color}" stroke-width="2"/>')
            for x, y, warning in series_points_for_label:
                parts.append(f'<circle cx="{sx(x):.2f}" cy="{sy(y):.2f}" r="4" fill="{"#d62728" if warning else color}" stroke="white" stroke-width="1"/>')
        if not stacked:
            for index, (series_label, _) in enumerate(points_by_label.items()):
                legend_y = y0 + 18 + index * 18
                parts.append(f'<rect x="735" y="{legend_y}" width="11" height="11" fill="{colors[index % len(colors)]}"/>')
                parts.append(f'<text x="752" y="{legend_y + 10}" font-family="Arial" font-size="11">{series_label}</text>')
        parts.append(f'<text x="18" y="{y0 + panel_h / 2:.0f}" text-anchor="middle" font-family="Arial" font-size="11" transform="rotate(-90 18,{y0 + panel_h / 2:.0f})">{y_label}</text>')
    parts.append(f'<text x="{width / 2:.0f}" y="{height - 16}" text-anchor="middle" font-family="Arial" font-size="12">Aleph unit-transition index</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def svg_series_fingerprint(structure_id: str, rows: list[dict[str, str]], path: Path) -> None:
    points = series_points(rows, structure_id, "local_twist_deg")
    svg_series_plot(
        f"{structure_id} Aleph series fingerprint",
        OrderedDict([(structure_id, points)]),
        "signed local twist (deg)",
        path,
        reference_y=30.0,
    )


def svg_series_comparison(rows: list[dict[str, str]], path: Path) -> None:
    points_by_label = OrderedDict(
        (structure_id, series_points(rows, structure_id, "local_twist_deg"))
        for structure_id in ["full", "central6", "central7"]
        if series_points(rows, structure_id, "local_twist_deg")
    )
    svg_series_plot(
        "Aleph series fingerprint comparison",
        points_by_label,
        "signed local twist (deg)",
        path,
        reference_y=30.0,
        stacked=True,
    )


def svg_companion_traces(structure_id: str, rows: list[dict[str, str]], path: Path) -> None:
    traces = OrderedDict()
    for column, label in [
        ("local_twist_deg", "signed twist"),
        ("local_rise_A", "local rise"),
        ("aleph_base_plane_bend_deg", "base-plane bend"),
        ("aleph_scaffold_plane_bend_deg", "scaffold-plane bend"),
    ]:
        points = series_points(rows, structure_id, column)
        if points:
            traces[label] = points
    svg_series_plot(
        f"{structure_id} Aleph series companion traces",
        traces,
        "feature value",
        path,
        reference_y=None,
        stacked=True,
    )


def svg_line_plot(title: str, rows: list[dict[str, str]], y_column: str, y_label: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    width = 900
    height = 430
    left = 70
    right = 30
    top = 45
    bottom = 60
    series = OrderedDict()
    for row in rows:
        value = safe_float(row.get(y_column))
        unit = safe_float(row.get("unit_index"))
        if value is None or unit is None:
            continue
        series.setdefault(row["structure_id"], []).append((unit, value))
    all_points = [point for points in series.values() for point in points]
    if not all_points:
        return
    min_x = min(point[0] for point in all_points)
    max_x = max(point[0] for point in all_points)
    min_y = min(point[1] for point in all_points)
    max_y = max(point[1] for point in all_points)
    if max_x == min_x:
        max_x += 1
    if max_y == min_y:
        max_y += 1
    plot_w = width - left - right
    plot_h = height - top - bottom
    colors = ["#4c78a8", "#f58518", "#54a24b", "#b279a2", "#e45756", "#72b7b2"]

    def sx(value: float) -> float:
        return left + (value - min_x) / (max_x - min_x) * plot_w

    def sy(value: float) -> float:
        return top + (max_y - value) / (max_y - min_y) * plot_h

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width / 2:.0f}" y="26" text-anchor="middle" font-family="Arial" font-size="17">{title}</text>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" stroke="#333"/>',
        f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}" stroke="#333"/>',
        f'<text x="{width / 2:.0f}" y="{height - 15}" text-anchor="middle" font-family="Arial" font-size="12">Aleph unit index</text>',
        f'<text x="18" y="{top + plot_h / 2:.0f}" text-anchor="middle" font-family="Arial" font-size="12" transform="rotate(-90 18,{top + plot_h / 2:.0f})">{y_label}</text>',
    ]
    for index, (structure_id, points) in enumerate(series.items()):
        points = sorted(points)
        color = colors[index % len(colors)]
        coords = " ".join(f"{sx(x):.2f},{sy(y):.2f}" for x, y in points)
        parts.append(f'<polyline points="{coords}" fill="none" stroke="{color}" stroke-width="2"/>')
        for x, y in points:
            parts.append(f'<circle cx="{sx(x):.2f}" cy="{sy(y):.2f}" r="2.5" fill="{color}"/>')
        legend_y = 50 + index * 20
        parts.append(f'<rect x="700" y="{legend_y}" width="12" height="12" fill="{color}"/>')
        parts.append(f'<text x="718" y="{legend_y + 11}" font-family="Arial" font-size="12">{structure_id}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def svg_fft_plot(rows: list[dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    width = 900
    height = 430
    left = 80
    bottom = 340
    filtered = [row for row in rows if row["dominant_frequency_index"]]
    if not filtered:
        return
    labels = [f"{row['structure_id']} {row['signal_name']}" for row in filtered]
    values = [safe_float(row["normalized_dominant_amplitude"]) or 0.0 for row in filtered]
    max_value = max(values + [1.0])
    bar_w = max(10, min(28, int(650 / max(1, len(values)))))
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="450" y="28" text-anchor="middle" font-family="Arial" font-size="17">Aleph FFT normalized dominant amplitudes</text>',
        f'<line x1="{left}" y1="55" x2="{left}" y2="{bottom}" stroke="#333"/>',
        f'<line x1="{left}" y1="{bottom}" x2="{width - 30}" y2="{bottom}" stroke="#333"/>',
    ]
    for index, value in enumerate(values):
        h = value / max_value * 250 if max_value else 0
        x = left + 20 + index * (bar_w + 8)
        y = bottom - h
        parts.append(f'<rect x="{x}" y="{y:.2f}" width="{bar_w}" height="{h:.2f}" fill="#4c78a8"/>')
        parts.append(f'<text x="{x + bar_w / 2:.0f}" y="{bottom + 14}" text-anchor="end" font-family="Arial" font-size="8" transform="rotate(-55 {x + bar_w / 2:.0f},{bottom + 14})">{labels[index]}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def svg_fingerprint_strip(structure_id: str, rows: list[dict[str, str]], path: Path, title: str | None = None) -> None:
    selected = sorted([row for row in rows if row["structure_id"] == structure_id], key=lambda row: int(row["unit_index"]))
    if not selected:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    cell_w = 38
    cell_h = 24
    left = 160
    top = 54
    row_gap = 5
    width = left + cell_w * len(selected) + 45
    height = top + (len(FINGERPRINT_FEATURES) + 1) * (cell_h + row_gap) + 58
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width / 2:.0f}" y="26" text-anchor="middle" font-family="Arial" font-size="17">{title or structure_id + " Aleph fingerprint strip"}</text>',
        f'<text x="{left + (cell_w * len(selected)) / 2:.0f}" y="{height - 18}" text-anchor="middle" font-family="Arial" font-size="12">Aleph unit index</text>',
    ]
    for index, row in enumerate(selected):
        x = left + index * cell_w
        parts.append(f'<text x="{x + cell_w / 2:.0f}" y="46" text-anchor="middle" font-family="Arial" font-size="10">{row["unit_index"]}</text>')
    for feature_index, (column, label, minimum, maximum, cyclic) in enumerate(FINGERPRINT_FEATURES):
        y = top + feature_index * (cell_h + row_gap)
        parts.append(f'<text x="{left - 8}" y="{y + 17}" text-anchor="end" font-family="Arial" font-size="11">{label}</text>')
        for unit_index, row in enumerate(selected):
            x = left + unit_index * cell_w
            color = value_color(normalize_feature_value(row.get(column), minimum, maximum, cyclic))
            value = row.get(column, "")
            parts.append(f'<rect x="{x}" y="{y}" width="{cell_w - 2}" height="{cell_h}" fill="{color}" stroke="#ffffff"/>')
            if value:
                parts.append(f'<title>{structure_id} unit {row["unit_index"]} {label}: {value}</title>')
    y = top + len(FINGERPRINT_FEATURES) * (cell_h + row_gap)
    parts.append(f'<text x="{left - 8}" y="{y + 17}" text-anchor="end" font-family="Arial" font-size="11">QC flag</text>')
    for unit_index, row in enumerate(selected):
        x = left + unit_index * cell_w
        flagged = row_has_qc_warning(row)
        parts.append(f'<rect x="{x}" y="{y}" width="{cell_w - 2}" height="{cell_h}" fill="{value_color(None, flagged)}" stroke="#ffffff"/>')
        parts.append(f'<text x="{x + cell_w / 2:.0f}" y="{y + 16}" text-anchor="middle" font-family="Arial" font-size="12" fill="{"white" if flagged else "#666"}">{"!" if flagged else "-"}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def svg_fingerprint_comparison(rows: list[dict[str, str]], path: Path) -> None:
    structures = [structure for structure in ["full", "central6", "central7"] if any(row["structure_id"] == structure for row in rows)]
    if not structures:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    max_units = max(int(row["unit_index"]) for row in rows)
    cell_w = 27
    cell_h = 14
    left = 176
    top = 56
    row_gap = 3
    panel_gap = 18
    feature_count = len(FINGERPRINT_FEATURES) + 1
    panel_h = feature_count * (cell_h + row_gap) + 18
    width = left + max_units * cell_w + 50
    height = top + len(structures) * (panel_h + panel_gap) + 36
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width / 2:.0f}" y="26" text-anchor="middle" font-family="Arial" font-size="17">Aleph fingerprint comparison</text>',
    ]
    for s_index, structure_id in enumerate(structures):
        selected = {int(row["unit_index"]): row for row in rows if row["structure_id"] == structure_id}
        panel_y = top + s_index * (panel_h + panel_gap)
        parts.append(f'<text x="18" y="{panel_y + 14}" font-family="Arial" font-size="13" font-weight="bold">{structure_id}</text>')
        for unit_index in range(1, max_units + 1):
            x = left + (unit_index - 1) * cell_w
            if s_index == 0:
                parts.append(f'<text x="{x + cell_w / 2:.0f}" y="48" text-anchor="middle" font-family="Arial" font-size="9">{unit_index}</text>')
        for feature_index, (column, label, minimum, maximum, cyclic) in enumerate(FINGERPRINT_FEATURES):
            y = panel_y + feature_index * (cell_h + row_gap)
            parts.append(f'<text x="{left - 8}" y="{y + 11}" text-anchor="end" font-family="Arial" font-size="9">{label}</text>')
            for unit_index in range(1, max_units + 1):
                row = selected.get(unit_index)
                color = value_color(normalize_feature_value(row.get(column), minimum, maximum, cyclic) if row else None)
                x = left + (unit_index - 1) * cell_w
                parts.append(f'<rect x="{x}" y="{y}" width="{cell_w - 1}" height="{cell_h}" fill="{color}" stroke="#ffffff"/>')
        y = panel_y + len(FINGERPRINT_FEATURES) * (cell_h + row_gap)
        parts.append(f'<text x="{left - 8}" y="{y + 11}" text-anchor="end" font-family="Arial" font-size="9">QC flag</text>')
        for unit_index in range(1, max_units + 1):
            row = selected.get(unit_index)
            flagged = bool(row and row_has_qc_warning(row))
            x = left + (unit_index - 1) * cell_w
            parts.append(f'<rect x="{x}" y="{y}" width="{cell_w - 1}" height="{cell_h}" fill="{value_color(None, flagged)}" stroke="#ffffff"/>')
    parts.append("</svg>")
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def svg_qc_summary_plot(summary_rows: list[dict[str, str]], qc_rows: list[dict[str, str]], path: Path) -> None:
    if not summary_rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    qc_by_structure = {row["structure_id"]: row for row in qc_rows}
    metrics = [
        ("mean_abs_local_twist_deg", "mean abs twist", "summary"),
        ("std_abs_local_twist_deg", "twist std", "summary"),
        ("negative_rise_count", "negative rises", "qc"),
        ("aleph_regular_score", "regular score", "summary"),
        ("qc_warning_count", "QC warnings", "derived"),
    ]
    width = 920
    height = 430
    left = 85
    top = 52
    bottom = 338
    group_w = 150
    bar_w = 20
    colors = ["#4c78a8", "#f58518", "#54a24b", "#b279a2", "#e45756"]
    structures = [row["structure_id"] for row in summary_rows]
    values_by_metric: list[list[float]] = []
    for column, _, source in metrics:
        values = []
        for summary in summary_rows:
            qc = qc_by_structure.get(summary["structure_id"], {})
            if source == "summary":
                value = safe_float(summary.get(column)) or 0.0
            elif source == "qc":
                value = safe_float(qc.get(column)) or 0.0
            else:
                value = float(qc_warning_count(qc))
            values.append(value)
        max_value = max(values + [1.0])
        values_by_metric.append([value / max_value if max_value else 0.0 for value in values])
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="460" y="28" text-anchor="middle" font-family="Arial" font-size="17">Aleph QC summary</text>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{bottom}" stroke="#333"/>',
        f'<line x1="{left}" y1="{bottom}" x2="{width - 30}" y2="{bottom}" stroke="#333"/>',
    ]
    for s_index, structure_id in enumerate(structures):
        group_x = left + 30 + s_index * group_w
        parts.append(f'<text x="{group_x + 45}" y="{bottom + 20}" text-anchor="middle" font-family="Arial" font-size="12">{structure_id}</text>')
        for m_index, (_, label, _) in enumerate(metrics):
            value = values_by_metric[m_index][s_index]
            h = value * 245
            x = group_x + m_index * (bar_w + 5)
            y = bottom - h
            parts.append(f'<rect x="{x}" y="{y:.2f}" width="{bar_w}" height="{h:.2f}" fill="{colors[m_index]}"/>')
    for m_index, (_, label, _) in enumerate(metrics):
        legend_x = 575
        legend_y = 62 + m_index * 21
        parts.append(f'<rect x="{legend_x}" y="{legend_y}" width="12" height="12" fill="{colors[m_index]}"/>')
        parts.append(f'<text x="{legend_x + 18}" y="{legend_y + 11}" font-family="Arial" font-size="11">{label} normalized across structures</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def write_plots(
    per_unit_rows: list[dict[str, str]],
    summary_rows: list[dict[str, str]],
    qc_rows: list[dict[str, str]],
    fft_rows: list[dict[str, str]],
    plot_dir: Path,
) -> list[Path]:
    specs = [
        ("aleph_local_twist_vs_unit.svg", "Aleph local twist after phase unwrapping", "local_twist_deg", "local twist (deg)"),
        ("aleph_local_rise_vs_unit.svg", "Aleph local rise after axis orientation", "local_rise_A", "local rise (A)"),
        ("aleph_radial_spread_vs_unit.svg", "Aleph radial spread vs unit index", "base_radial_spread_A", "base radial spread (A)"),
        ("aleph_phase_raw_vs_unit.svg", "Aleph raw angular phase progression", "aleph_phase_raw_deg", "raw Aleph phase (deg)"),
        ("aleph_phase_progression_vs_unit.svg", "Aleph unwrapped angular phase progression", "aleph_phase_deg", "unwrapped Aleph phase (deg)"),
        ("aleph_chain_resultant_vs_unit.svg", "Aleph chain angular mean resultant length", "chain_mean_resultant_length", "mean resultant length"),
    ]
    paths = []
    for filename, title, column, label in specs:
        path = plot_dir / filename
        svg_line_plot(title, per_unit_rows, column, label, path)
        if path.exists():
            paths.append(path)
    for structure_id in ["full", "central6", "central7"]:
        path = plot_dir / f"aleph_series_fingerprint_{structure_id}.svg"
        svg_series_fingerprint(structure_id, per_unit_rows, path)
        if path.exists():
            paths.append(path)
    series_comparison_path = plot_dir / "aleph_series_fingerprint_comparison.svg"
    svg_series_comparison(per_unit_rows, series_comparison_path)
    if series_comparison_path.exists():
        paths.append(series_comparison_path)
    for structure_id in ["full", "central7"]:
        path = plot_dir / f"aleph_series_companion_traces_{structure_id}.svg"
        svg_companion_traces(structure_id, per_unit_rows, path)
        if path.exists():
            paths.append(path)
    for structure_id in ["full", "central6", "central7"]:
        path = plot_dir / f"aleph_fingerprint_strip_{structure_id}.svg"
        svg_fingerprint_strip(structure_id, per_unit_rows, path)
        if path.exists():
            paths.append(path)
    comparison_path = plot_dir / "aleph_fingerprint_comparison.svg"
    svg_fingerprint_comparison(per_unit_rows, comparison_path)
    if comparison_path.exists():
        paths.append(comparison_path)
    qc_path = plot_dir / "aleph_fingerprint_qc_summary.svg"
    svg_qc_summary_plot(summary_rows, qc_rows, qc_path)
    if qc_path.exists():
        paths.append(qc_path)
    fft_path = plot_dir / "aleph_fft_dominant_amplitudes.svg"
    svg_fft_plot(fft_rows, fft_path)
    if fft_path.exists():
        paths.append(fft_path)
    return paths


def report_text(
    fingerprints: list[StructureFingerprint],
    fft_rows: list[dict[str, str]],
    qc_rows: list[dict[str, str]],
    plot_paths: list[Path],
    optional_warning: str,
    args: argparse.Namespace,
) -> str:
    lines = [
        "# Aleph Geometric Fingerprint Prototype",
        "",
        "Aleph is an exploratory one-dimensional geometric fingerprint along the fitted hexaplex axis. It is not a diffraction simulator and does not prove formation, stability, or experimental correctness.",
        "",
        "Aleph asks whether each model has an ordered repeating geometric signature along its axis that can be plotted as per-unit signals and summarized with a discrete FFT.",
        "",
        "## Inputs",
        "",
    ]
    for fingerprint in fingerprints:
        lines.append(
            f"- `{fingerprint.structure_id}`: `{fingerprint.path}` ({fingerprint.chain_count} chains, {fingerprint.summary_row['unit_count']} units, {fingerprint.heavy_atom_count} heavy atoms)"
        )
    if optional_warning:
        lines.append(f"- Optional compact twist variants: {optional_warning}")

    lines.extend(
        [
            "",
            "## Aleph Summary",
            "",
            "| Structure | Units | Mean signed twist | Mean abs twist | Twist std | Mean rise | Rise std | Regular score |",
            "|---|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for fingerprint in fingerprints:
        row = fingerprint.summary_row
        lines.append(
            f"| {row['structure_id']} | {row['unit_count']} | {row['mean_local_twist_deg']} | {row['mean_abs_local_twist_deg']} | {row['std_local_twist_deg']} | "
            f"{row['mean_local_rise_A']} | {row['std_local_rise_A']} | {row['aleph_regular_score']} |"
        )

    lines.extend(
        [
            "",
            "## Geometry QC",
            "",
            "| Structure | Axis flipped | Phase unwrap ok | Expected units | Observed units | Missing-chain units | Negative rises | Twist warning | Rise warning | Plane warning |",
            "|---|---|---|---:|---:|---:|---:|---|---|---|",
        ]
    )
    for row in qc_rows:
        lines.append(
            f"| {row['structure_id']} | {row['axis_flipped']} | {row['phase_unwrap_ok']} | {row['expected_unit_count']} | "
            f"{row['observed_unit_count']} | {row['units_with_missing_chains']} | {row['negative_rise_count']} | "
            f"{row['local_twist_warning']} | {row['local_rise_warning']} | {row['plane_fit_warning']} |"
        )

    lines.extend(
        [
            "",
            "## Series-style Aleph fingerprint",
            "",
            "The Aleph series fingerprint is the primary visual fingerprint in this prototype. It is closer to a true structural fingerprint because it represents the hexaplex as an ordered one-dimensional trace rather than a multi-feature dashboard.",
            "",
            "The primary Aleph series is signed local twist between adjacent units, in degrees. The x-axis is Aleph unit-transition index, so transition 1 corresponds to the step from unit 1 to unit 2. The y-axis is signed local twist after phase unwrapping. A 30 deg reference line is shown where useful, and markers with QC warnings are highlighted.",
            "",
            "Under the current Aleph definitions, central7 currently looks like the cleanest 30 deg-like Aleph series fingerprint because its mean absolute local twist is near 30 deg, its rise is positive, and it has no QC warnings. central6 is shorter and has positive rise, but its signed twist trace and mean absolute twist deviate from the nominal 30 deg value. The full model remains a geometry-definition diagnostic case because its current warnings indicate that layer assignment and antiparallel ordering require further inspection.",
            "",
            "This representation provides a more natural ordered signal for future DFT/FFT exploration than the feature comparison panel. Whether spectral analysis adds value remains an open question.",
            "",
            "Series outputs:",
            "",
            "- `outputs\\plots\\aleph_fingerprint\\aleph_series_fingerprint_full.svg`",
            "- `outputs\\plots\\aleph_fingerprint\\aleph_series_fingerprint_central6.svg`",
            "- `outputs\\plots\\aleph_fingerprint\\aleph_series_fingerprint_central7.svg`",
            "- `outputs\\plots\\aleph_fingerprint\\aleph_series_fingerprint_comparison.svg`",
            "- `outputs\\plots\\aleph_fingerprint\\aleph_series_companion_traces_full.svg`",
            "- `outputs\\plots\\aleph_fingerprint\\aleph_series_companion_traces_central7.svg`",
            "",
            "## Visual Aleph feature comparison panel",
            "",
            "Aleph converts the hexaplex into ordered per-unit geometric traces. The feature comparison panel places unit index on the x-axis and stacks Aleph features as rows, so changes in twist, rise, radial spread, phase, plane bend, chain coherence, and QC flags can be scanned compactly.",
            "",
            "Rows labeled `abs twist` and `signed twist` show the local angular step after phase unwrapping; `rise` shows the axis-oriented local axial step; `radial spread` tracks base-like radial variability; `phase` shows the unwrapped angular progression; `base bend` and `scaffold bend` compare adjacent fitted plane normals; `chain coherence` reports the circular mean resultant length; and the QC row marks per-unit warnings.",
            "",
            "The visualization is a structural fingerprint, not a diffraction simulation. It reveals axial ordering, local geometric irregularity, phase progression, and QC behavior that pair-distance counts do not show directly.",
            "",
            "Feature comparison panel outputs:",
            "",
            "- `outputs\\plots\\aleph_fingerprint\\aleph_fingerprint_strip_full.svg`",
            "- `outputs\\plots\\aleph_fingerprint\\aleph_fingerprint_strip_central6.svg`",
            "- `outputs\\plots\\aleph_fingerprint\\aleph_fingerprint_strip_central7.svg`",
            "- `outputs\\plots\\aleph_fingerprint\\aleph_fingerprint_comparison.svg`",
            "- `outputs\\plots\\aleph_fingerprint\\aleph_fingerprint_qc_summary.svg`",
        ]
    )

    lines.extend(
        [
            "",
            "## FFT Diagnostic Summary",
            "",
            "| Structure | Signal | Samples | Dominant index | Normalized amplitude | Warning |",
            "|---|---|---:|---:|---:|---|",
        ]
    )
    for row in fft_rows:
        lines.append(
            f"| {row['structure_id']} | {row['signal_name']} | {row['sample_count']} | {row['dominant_frequency_index']} | "
            f"{row['normalized_dominant_amplitude']} | {row['warnings']} |"
        )

    full = next((fingerprint for fingerprint in fingerprints if fingerprint.structure_id == "full"), None)
    central6 = next((fingerprint for fingerprint in fingerprints if fingerprint.structure_id == "central6"), None)
    central7 = next((fingerprint for fingerprint in fingerprints if fingerprint.structure_id == "central7"), None)
    lines.extend(["", "## Interpretation", ""])
    if full is not None:
        lines.append(
            "The full model remains the best first Aleph diagnostic target because it has the longest ordered per-unit signal; shorter central6 and central7 fragments are useful for local comparison but have limited spectral resolution."
        )
    if central6 is not None and central7 is not None:
        lines.append(
            "central6 and central7 can be compared to the full model by local twist, local rise, radial spread, and phase progression rather than by pair-distance counts alone."
        )
    lines.append(
        "Backbone/scaffold Aleph bend and base-like Aleph bend are reported separately where per-layer plane normals are computable; disagreement between them would indicate different scaffold and base-like geometric signals."
    )
    lines.append(
        "This QC pass should be read before expanding FFT interpretation: stable axis orientation, phase unwrapping, positive rise convention, and bounded circular dispersion are prerequisites for treating Aleph signals as repeat fingerprints."
    )
    lines.extend(
        [
            "",
                "## Assumptions And Cautions",
                "",
            "- Unit assignment uses base-like CYP/MEP residues as repeat anchors by residue order within each chain.",
            "- The Aleph phase and axial position use the first available chain-specific base-like centroid as an angular/axial anchor for each unit, while the layer centroid is still reported separately. This avoids symmetry cancellation of the six-strand layer centroid.",
            "- The fitted axis is flipped when needed so the chain-specific anchor axial coordinate generally increases with unit index.",
            "- Chain angular dispersion is computed with bounded circular statistics and paired with mean resultant length.",
            "- FFT summaries for 6- and 7-unit models are explicitly marked as short-signal diagnostics.",
            "- Aleph is a geometric fingerprint, not a diffraction simulator or structural mechanism.",
            "",
            "## Outputs",
            "",
            f"- Per-unit CSV: `{args.per_unit_csv}`",
            f"- Summary CSV: `{args.summary_csv}`",
            f"- FFT CSV: `{args.fft_csv}`",
            f"- QC CSV: `{args.qc_csv}`",
            f"- Plot directory: `{args.plot_dir}`",
        ]
    )
    for path in plot_paths:
        lines.append(f"- Plot: `{path}`")
    warnings = [warning for fingerprint in fingerprints for warning in fingerprint.warnings]
    if warnings:
        lines.extend(["", "## Warnings", ""])
        for warning in warnings:
            lines.append(f"- {warning}")
    return "\n".join(lines) + "\n"


def run(args: argparse.Namespace) -> dict[str, object]:
    inputs, input_warnings = load_structure_inputs(args.include_optional_twists)
    fingerprints: list[StructureFingerprint] = []
    for structure_id, path in inputs:
        if path.exists():
            fingerprints.append(analyze_structure(structure_id, path))
    per_unit_rows = [row for fingerprint in fingerprints for row in fingerprint.per_unit_rows]
    summary_rows = [fingerprint.summary_row for fingerprint in fingerprints]
    qc_rows = [fingerprint.qc_row for fingerprint in fingerprints]
    fft_rows: list[dict[str, str]] = []
    for fingerprint in fingerprints:
        for signal_name, column in [
            ("local_twist_deg", "local_twist_deg"),
            ("local_rise_A", "local_rise_A"),
            ("radial_spread_A", "base_radial_spread_A"),
        ]:
            fft_rows.append(fft_summary(fingerprint.structure_id, signal_name, collect_signal(per_unit_rows, fingerprint.structure_id, column)))

    write_csv(args.per_unit_csv, per_unit_rows, PER_UNIT_COLUMNS)
    write_csv(args.summary_csv, summary_rows, SUMMARY_COLUMNS)
    write_csv(args.fft_csv, fft_rows, FFT_COLUMNS)
    write_csv(args.qc_csv, qc_rows, QC_COLUMNS)
    plot_paths = write_plots(per_unit_rows, summary_rows, qc_rows, fft_rows, args.plot_dir)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(report_text(fingerprints, fft_rows, qc_rows, plot_paths, "; ".join(input_warnings), args), encoding="utf-8")
    return {
        "per_unit_rows": len(per_unit_rows),
        "summary_rows": len(summary_rows),
        "fft_rows": len(fft_rows),
        "qc_rows": len(qc_rows),
        "plots": plot_paths,
        "fingerprints": fingerprints,
        "input_warnings": input_warnings,
    }


def main() -> None:
    result = run(parse_args())
    print(f"Wrote {result['per_unit_rows']} Aleph per-unit rows")
    print(f"Wrote {result['summary_rows']} Aleph summary rows")
    print(f"Wrote {result['fft_rows']} Aleph FFT rows")
    print(f"Wrote {result['qc_rows']} Aleph QC rows")
    print(f"Wrote {len(result['plots'])} plots")
    if result["input_warnings"]:
        print("; ".join(result["input_warnings"]))


if __name__ == "__main__":
    main()
