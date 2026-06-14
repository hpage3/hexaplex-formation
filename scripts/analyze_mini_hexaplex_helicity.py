#!/usr/bin/env python3
"""Analyze fitted-axis helicity of generated mini-hexaplex variants."""

from __future__ import annotations

import argparse
import csv
import math
import re
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
    distance,
    dot,
    group_atoms_by_residue,
    mean_point,
    power_iteration_principal_axis,
    project_point_to_axis,
    vector_sub,
)
from hexaplex_formation.pdb_utils import PDBAtom, chain_ids, heavy_atoms, load_pdb_atoms, residue_count  # noqa: E402


BASE_RESIDUES = {"CYP", "MEP"}
BACKBONE_ATOM_NAMES = {"N", "CA", "C", "O", "OXT"}
FULL_BASELINE_ID = "full_length_baseline"
DEFAULT_BASELINE_PDB = Path(
    "outputs/intermediates/ai_candidate_inputs/full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain.pdb"
)

SUMMARY_COLUMNS = [
    "variant_id",
    "truncation_rule",
    "units_per_chain",
    "chain_count",
    "residues_per_chain",
    "total_residue_count",
    "total_atom_count",
    "representative_point_rule",
    "axis_reference",
    "mean_helical_r2",
    "median_helical_r2",
    "mean_angular_residual_deg",
    "mean_twist_per_unit_deg",
    "std_twist_per_unit_deg",
    "mean_pitch_A",
    "axial_extent_A",
    "normalized_axial_extent_vs_full",
    "coherent_helical_turns",
    "normalized_coherent_helical_turns_vs_full",
    "helical_coherence_score",
    "six_strand_phase_coherence_score",
    "circular_std_phase_deg",
    "structural_coherence_flag",
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
    glu_residue: ResidueRecord


@dataclass(frozen=True)
class UnitPoint:
    chain_id: str
    unit_index: int
    z: float
    theta: float
    radius: float
    point_rule: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--structures-dir", type=Path, default=Path("outputs/mini_hexaplex/structures"))
    parser.add_argument("--baseline-pdb", type=Path, default=DEFAULT_BASELINE_PDB)
    parser.add_argument("--manifest", type=Path, default=Path("outputs/mini_hexaplex/mini_hexaplex_variant_manifest.csv"))
    parser.add_argument("--geometry", type=Path, default=Path("outputs/metrics/mini_hexaplex_geometry_summary.csv"))
    parser.add_argument("--feature-summary", type=Path, default=Path("outputs/metrics/mini_hexaplex_feature_summary.csv"))
    parser.add_argument("--out-csv", type=Path, default=Path("outputs/metrics/mini_hexaplex_helicity_summary.csv"))
    parser.add_argument("--plot-dir", type=Path, default=Path("outputs/mini_hexaplex/plots"))
    parser.add_argument("--out-report", type=Path, default=Path("outputs/reports/mini_hexaplex_helicity_report.md"))
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


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def infer_axis(atoms: list[PDBAtom]) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    heavy = heavy_atoms(atoms)
    backbone = [atom for atom in heavy if atom.atom_name.strip().upper() in BACKBONE_ATOM_NAMES]
    axis_atoms = backbone if len(backbone) >= 3 else heavy
    if len(axis_atoms) < 3:
        raise ValueError("At least three atoms are required to infer a fitted axis")
    origin = mean_point(atom_xyz(atom) for atom in axis_atoms)
    axis = power_iteration_principal_axis(covariance_matrix_3d(atom_xyz(atom) for atom in axis_atoms))
    return origin, axis


def axial_extent_A(
    atoms: list[PDBAtom],
    axis_origin: tuple[float, float, float],
    axis: tuple[float, float, float],
) -> float | None:
    axial_values = [project_point_to_axis(atom_xyz(atom), axis_origin, axis)[0] for atom in atoms]
    if not axial_values:
        return None
    return max(axial_values) - min(axial_values)


def build_residues(atoms: list[PDBAtom]) -> OrderedDict[str, list[ResidueRecord]]:
    by_chain: OrderedDict[str, list[ResidueRecord]] = OrderedDict()
    for (chain_id, residue_name, residue_number, insertion_code), residue_atoms in group_atoms_by_residue(atoms).items():
        by_chain.setdefault(chain_id, []).append(
            ResidueRecord(
                chain_id=chain_id,
                residue_name=residue_name,
                residue_number=residue_number,
                insertion_code=insertion_code,
                atoms=tuple(residue_atoms),
            )
        )
    return by_chain


def group_repeat_units(by_chain: OrderedDict[str, list[ResidueRecord]]) -> tuple[dict[str, list[RepeatUnit]], list[str]]:
    units_by_chain: dict[str, list[RepeatUnit]] = {}
    warnings: list[str] = []
    for chain_id, residues in by_chain.items():
        if len(residues) % 2:
            warnings.append(f"chain {chain_id} has odd residue count {len(residues)}")
        units: list[RepeatUnit] = []
        for unit_index, start in enumerate(range(0, len(residues) - 1, 2), start=1):
            base = residues[start]
            glu = residues[start + 1]
            if base.residue_name not in BASE_RESIDUES or glu.residue_name != "GLU":
                warnings.append(
                    f"chain {chain_id} unit {unit_index} is {base.residue_name}/{glu.residue_name}, expected CYP|MEP/GLU"
                )
            units.append(RepeatUnit(chain_id=chain_id, unit_index=unit_index, base_residue=base, glu_residue=glu))
        units_by_chain[chain_id] = units
    return units_by_chain, warnings


def representative_point(unit: RepeatUnit) -> tuple[tuple[float, float, float], str]:
    ca_atoms = [atom for atom in unit.base_residue.atoms if atom.atom_name.strip().upper() == "CA"]
    if ca_atoms:
        return atom_xyz(ca_atoms[0]), "base residue CA atom"
    heavy = heavy_atoms(unit.base_residue.atoms)
    if heavy:
        return mean_point(atom_xyz(atom) for atom in heavy), "base residue heavy-atom centroid"
    return mean_point(atom_xyz(atom) for atom in unit.base_residue.atoms), "base residue all-atom centroid"


def cylindrical_coordinates(
    point: tuple[float, float, float],
    origin: tuple[float, float, float],
    axis: tuple[float, float, float],
    basis_u: tuple[float, float, float],
    basis_v: tuple[float, float, float],
) -> tuple[float, float, float]:
    z, projected = project_point_to_axis(point, origin, axis)
    radial_vector = vector_sub(point, projected)
    theta = math.atan2(dot(radial_vector, basis_v), dot(radial_vector, basis_u))
    radius = distance(point, projected)
    return z, theta, radius


def unwrap_angles(angles: list[float]) -> list[float]:
    if not angles:
        return []
    unwrapped = [angles[0]]
    offset = 0.0
    previous = angles[0]
    for angle in angles[1:]:
        delta = angle - previous
        if delta > math.pi:
            offset -= 2.0 * math.pi
        elif delta < -math.pi:
            offset += 2.0 * math.pi
        unwrapped.append(angle + offset)
        previous = angle
    return unwrapped


def linear_fit(x_values: list[float], y_values: list[float]) -> tuple[float, float, float, float]:
    if len(x_values) != len(y_values) or len(x_values) < 2:
        raise ValueError("linear_fit requires at least two paired values")
    x_mean = sum(x_values) / len(x_values)
    y_mean = sum(y_values) / len(y_values)
    ss_xx = sum((x - x_mean) ** 2 for x in x_values)
    if ss_xx == 0:
        raise ValueError("linear_fit requires non-identical x values")
    slope = sum((x - x_mean) * (y - y_mean) for x, y in zip(x_values, y_values)) / ss_xx
    intercept = y_mean - slope * x_mean
    residuals = [y - (slope * x + intercept) for x, y in zip(x_values, y_values)]
    ss_res = sum(residual * residual for residual in residuals)
    ss_tot = sum((y - y_mean) ** 2 for y in y_values)
    r2 = 1.0 if ss_tot == 0 and ss_res == 0 else 1.0 - ss_res / ss_tot if ss_tot else 0.0
    rms_residual = math.sqrt(ss_res / len(residuals))
    return slope, intercept, max(0.0, min(1.0, r2)), rms_residual


def coefficient_of_determination(x_values: list[float], y_values: list[float]) -> float:
    return linear_fit(x_values, y_values)[2]


def chain_helicity_metrics(points: list[UnitPoint]) -> dict[str, float | None]:
    ordered = sorted(points, key=lambda point: (point.z, point.unit_index))
    if len(ordered) < 2:
        return {
            "helical_r2": None,
            "angular_residual_deg": None,
            "mean_twist_per_unit_deg": None,
            "std_twist_per_unit_deg": None,
            "pitch_A": None,
        }
    z_values = [point.z for point in ordered]
    theta_values = unwrap_angles([point.theta for point in ordered])
    slope, _intercept, r2, rms_residual = linear_fit(z_values, theta_values)
    twists = [math.degrees(theta_values[index + 1] - theta_values[index]) for index in range(len(theta_values) - 1)]
    mean_twist = sum(twists) / len(twists) if twists else None
    if twists:
        twist_mean = mean_twist or 0.0
        std_twist = math.sqrt(sum((twist - twist_mean) ** 2 for twist in twists) / len(twists))
    else:
        std_twist = None
    pitch = abs((2.0 * math.pi) / slope) if abs(slope) > 1e-12 else None
    return {
        "helical_r2": r2,
        "angular_residual_deg": math.degrees(rms_residual),
        "mean_twist_per_unit_deg": mean_twist,
        "std_twist_per_unit_deg": std_twist,
        "pitch_A": pitch,
    }


def median(values: list[float]) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    midpoint = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[midpoint]
    return (ordered[midpoint - 1] + ordered[midpoint]) / 2.0


def mean(values: list[float]) -> float | None:
    return sum(values) / len(values) if values else None


def circular_std_deg(angles: list[float]) -> float | None:
    if not angles:
        return None
    c = sum(math.cos(angle) for angle in angles) / len(angles)
    s = sum(math.sin(angle) for angle in angles) / len(angles)
    resultant = min(1.0, max(0.0, math.hypot(c, s)))
    if resultant <= 0:
        return 180.0
    return math.degrees(math.sqrt(-2.0 * math.log(resultant)))


def six_strand_phase_coherence(unit_points_by_chain: dict[str, list[UnitPoint]]) -> tuple[float | None, float | None]:
    if len(unit_points_by_chain) < 2:
        return None, None
    sorted_by_chain = {
        chain_id: sorted(points, key=lambda point: (point.z, point.unit_index))
        for chain_id, points in unit_points_by_chain.items()
    }
    comparable_count = min((len(points) for points in sorted_by_chain.values()), default=0)
    if comparable_count == 0:
        return None, None
    per_level_std: list[float] = []
    for index in range(comparable_count):
        angles = [points[index].theta for points in sorted_by_chain.values()]
        ordered_angles = sorted(angles)
        gaps = [
            (ordered_angles[(gap_index + 1) % len(ordered_angles)] - ordered_angles[gap_index]) % (2.0 * math.pi)
            for gap_index in range(len(ordered_angles))
        ]
        ideal_gap = (2.0 * math.pi) / len(ordered_angles)
        gap_errors = [((gap - ideal_gap + math.pi) % (2.0 * math.pi)) - math.pi for gap in gaps]
        per_level_std_value = math.sqrt(sum(error * error for error in gap_errors) / len(gap_errors))
        per_level_std.append(math.degrees(per_level_std_value))
    mean_std = mean(per_level_std)
    if mean_std is None:
        return None, None
    score = max(0.0, min(1.0, 1.0 - mean_std / 60.0))
    return score, mean_std


def unit_count_from_variant_id(variant_id: str) -> str:
    if variant_id == FULL_BASELINE_ID:
        return "15"
    match = re.search(r"(?:first|central)([0-9]+)_units$", variant_id)
    return match.group(1) if match else ""


def residues_per_chain_text(units_by_chain: dict[str, list[RepeatUnit]]) -> str:
    return ";".join(f"{chain_id}:{len(units) * 2}" for chain_id, units in sorted(units_by_chain.items()))


def analyze_structure(
    variant_id: str,
    pdb_path: Path,
    axis_origin: tuple[float, float, float],
    axis: tuple[float, float, float],
    manifest_by_variant: dict[str, dict[str, str]],
    geometry_by_variant: dict[str, dict[str, str]],
) -> dict[str, str]:
    atoms = load_pdb_atoms(pdb_path)
    residues = build_residues(atoms)
    units_by_chain, grouping_warnings = group_repeat_units(residues)
    basis_u, basis_v = build_perpendicular_basis(axis)
    unit_points_by_chain: dict[str, list[UnitPoint]] = {}
    point_rules: set[str] = set()
    for chain_id, units in units_by_chain.items():
        chain_points: list[UnitPoint] = []
        for unit in units:
            point, rule = representative_point(unit)
            point_rules.add(rule)
            z, theta, radius = cylindrical_coordinates(point, axis_origin, axis, basis_u, basis_v)
            chain_points.append(UnitPoint(chain_id, unit.unit_index, z, theta, radius, rule))
        unit_points_by_chain[chain_id] = chain_points

    chain_metrics = [chain_helicity_metrics(points) for points in unit_points_by_chain.values()]
    r2_values = [value for value in (safe_float(metric["helical_r2"]) for metric in chain_metrics) if value is not None]
    residual_values = [value for value in (safe_float(metric["angular_residual_deg"]) for metric in chain_metrics) if value is not None]
    twist_values = [value for value in (safe_float(metric["mean_twist_per_unit_deg"]) for metric in chain_metrics) if value is not None]
    std_twist_values = [value for value in (safe_float(metric["std_twist_per_unit_deg"]) for metric in chain_metrics) if value is not None]
    pitch_values = [value for value in (safe_float(metric["pitch_A"]) for metric in chain_metrics) if value is not None]
    phase_score, phase_std = six_strand_phase_coherence(unit_points_by_chain)

    manifest_row = manifest_by_variant.get(variant_id, {})
    geometry_row = geometry_by_variant.get(variant_id, {})
    units_per_chain = manifest_row.get("units_per_chain") or geometry_row.get("units_per_chain") or unit_count_from_variant_id(variant_id)
    structural_flag = geometry_row.get("structural_coherence_flag") or ("coherent" if variant_id == FULL_BASELINE_ID else "")
    axial_extent = safe_float(geometry_row.get("axial_extent_A")) or axial_extent_A(atoms, axis_origin, axis)
    truncation_rule = manifest_row.get("truncation_rule") or ("full-length baseline reference" if variant_id == FULL_BASELINE_ID else "")
    warnings = list(grouping_warnings)
    chain_count = len(chain_ids(atoms))
    expected_count = safe_float(units_per_chain)
    if expected_count is not None:
        for chain_id, units in sorted(units_by_chain.items()):
            if len(units) != int(expected_count):
                warnings.append(f"chain {chain_id} has {len(units)} repeat units, expected {int(expected_count)}")
    if chain_count != 6:
        warnings.append(f"expected six chains, found {chain_count}")

    mean_r2 = mean(r2_values)
    score = max(0.0, min(1.0, mean_r2 if mean_r2 is not None else 0.0))
    return {
        "variant_id": variant_id,
        "truncation_rule": truncation_rule,
        "units_per_chain": units_per_chain,
        "chain_count": str(chain_count),
        "residues_per_chain": manifest_row.get("residues_per_chain") or geometry_row.get("residues_per_chain") or residues_per_chain_text(units_by_chain),
        "total_residue_count": manifest_row.get("total_residue_count") or geometry_row.get("total_residue_count") or str(residue_count(atoms)),
        "total_atom_count": manifest_row.get("total_atom_count") or geometry_row.get("total_atom_count") or str(len(atoms)),
        "representative_point_rule": "; ".join(sorted(point_rules)),
        "axis_reference": str(DEFAULT_BASELINE_PDB),
        "mean_helical_r2": format_float(mean_r2),
        "median_helical_r2": format_float(median(r2_values)),
        "mean_angular_residual_deg": format_float(mean(residual_values)),
        "mean_twist_per_unit_deg": format_float(mean(twist_values)),
        "std_twist_per_unit_deg": format_float(mean(std_twist_values)),
        "mean_pitch_A": format_float(mean(pitch_values)),
        "axial_extent_A": format_float(axial_extent),
        "normalized_axial_extent_vs_full": "",
        "coherent_helical_turns": "",
        "normalized_coherent_helical_turns_vs_full": "",
        "helical_coherence_score": format_float(score),
        "six_strand_phase_coherence_score": format_float(phase_score),
        "circular_std_phase_deg": format_float(phase_std),
        "structural_coherence_flag": structural_flag,
        "warnings": "; ".join(warnings),
    }


def add_normalized_extent_metrics(rows: list[dict[str, str]]) -> None:
    baseline_row = next((row for row in rows if row.get("variant_id") == FULL_BASELINE_ID), None)
    if baseline_row is None:
        return
    baseline_axial_extent = safe_float(baseline_row.get("axial_extent_A"))
    baseline_pitch = safe_float(baseline_row.get("mean_pitch_A"))
    baseline_turns = (
        baseline_axial_extent / baseline_pitch
        if baseline_axial_extent is not None and baseline_pitch is not None and baseline_pitch > 0
        else None
    )
    for row in rows:
        axial_extent = safe_float(row.get("axial_extent_A"))
        mean_pitch = safe_float(row.get("mean_pitch_A"))
        coherent_turns = axial_extent / mean_pitch if axial_extent is not None and mean_pitch is not None and mean_pitch > 0 else None
        row["coherent_helical_turns"] = format_float(coherent_turns)
        row["normalized_axial_extent_vs_full"] = format_float(
            axial_extent / baseline_axial_extent
            if axial_extent is not None and baseline_axial_extent is not None and baseline_axial_extent > 0
            else None
        )
        row["normalized_coherent_helical_turns_vs_full"] = format_float(
            coherent_turns / baseline_turns
            if coherent_turns is not None and baseline_turns is not None and baseline_turns > 0
            else None
        )


def structure_paths(structures_dir: Path, baseline_pdb: Path) -> list[tuple[str, Path]]:
    paths = [(FULL_BASELINE_ID, baseline_pdb)]
    for path in sorted(structures_dir.glob("mini_hexaplex_*.pdb")):
        variant_id = path.stem.removeprefix("mini_hexaplex_")
        if variant_id.startswith("literal_first"):
            continue
        paths.append((variant_id, path))
    return paths


def write_plots(rows: list[dict[str, str]], feature_rows: list[dict[str, str]], plot_dir: Path) -> list[Path]:
    plot_dir.mkdir(parents=True, exist_ok=True)
    expected_paths = [
        plot_dir / "mini_hexaplex_units_vs_helical_coherence_score.png",
        plot_dir / "mini_hexaplex_units_vs_mean_angular_residual_deg.png",
        plot_dir / "mini_hexaplex_units_vs_std_twist_per_unit_deg.png",
        plot_dir / "mini_hexaplex_units_vs_coherent_helical_turns.png",
        plot_dir / "mini_hexaplex_units_vs_normalized_axial_extent_vs_full.png",
        plot_dir / "mini_hexaplex_units_vs_ratio_to_full_4p5_5p0.png",
        plot_dir / "mini_hexaplex_helicity_vs_4p5_5p0_response.png",
        plot_dir / "mini_hexaplex_coherent_turns_and_4p5_5p0_response.png",
    ]
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover
        print(f"WARNING: matplotlib unavailable; skipping plots: {exc}", file=sys.stderr)
        return [path for path in expected_paths if path.exists()]

    outputs: list[Path] = []
    feature_by_variant = {row["variant_id"]: row for row in feature_rows}
    families = {
        "central": [row for row in rows if row["variant_id"].startswith("central") or row["variant_id"] == FULL_BASELINE_ID],
        "lower_end": [row for row in rows if row["variant_id"].startswith("lower_end") or row["variant_id"] == FULL_BASELINE_ID],
    }

    def plot_metric(path: Path, y_key: str, y_label: str, title: str) -> None:
        fig, ax = plt.subplots(figsize=(8.0, 4.8), dpi=180)
        for family, family_rows in families.items():
            points = sorted(
                (
                    safe_float(row.get("units_per_chain")),
                    safe_float(row.get(y_key)),
                )
                for row in family_rows
            )
            x_values = [point[0] for point in points if point[0] is not None and point[1] is not None]
            y_values = [point[1] for point in points if point[0] is not None and point[1] is not None]
            if x_values:
                ax.plot(x_values, y_values, marker="o", linewidth=1.2, label=family)
        ax.set_xlabel("base/GLU units per chain")
        ax.set_ylabel(y_label)
        ax.set_title(title)
        ax.grid(True, alpha=0.25, linewidth=0.6)
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(path)
        plt.close(fig)
        outputs.append(path)

    def plot_feature_metric(path: Path, y_key: str, y_label: str, title: str) -> None:
        fig, ax = plt.subplots(figsize=(8.0, 4.8), dpi=180)
        for family, family_rows in families.items():
            points = []
            for row in family_rows:
                feature_row = feature_by_variant.get(row["variant_id"], {})
                points.append((safe_float(row.get("units_per_chain")), safe_float(feature_row.get(y_key))))
            points = sorted(points)
            x_values = [point[0] for point in points if point[0] is not None and point[1] is not None]
            y_values = [point[1] for point in points if point[0] is not None and point[1] is not None]
            if x_values:
                ax.plot(x_values, y_values, marker="o", linewidth=1.2, label=family)
        ax.set_xlabel("base/GLU units per chain")
        ax.set_ylabel(y_label)
        ax.set_title(title)
        ax.grid(True, alpha=0.25, linewidth=0.6)
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(path)
        plt.close(fig)
        outputs.append(path)

    plot_metric(
        plot_dir / "mini_hexaplex_units_vs_helical_coherence_score.png",
        "helical_coherence_score",
        "helical coherence score",
        "Mini-hexaplex length vs helical coherence",
    )
    plot_metric(
        plot_dir / "mini_hexaplex_units_vs_mean_angular_residual_deg.png",
        "mean_angular_residual_deg",
        "mean angular residual (deg)",
        "Mini-hexaplex length vs angular residual",
    )
    plot_metric(
        plot_dir / "mini_hexaplex_units_vs_std_twist_per_unit_deg.png",
        "std_twist_per_unit_deg",
        "std twist per unit (deg)",
        "Mini-hexaplex length vs twist variability",
    )
    plot_metric(
        plot_dir / "mini_hexaplex_units_vs_coherent_helical_turns.png",
        "coherent_helical_turns",
        "coherent helical turns",
        "Mini-hexaplex length vs coherent helical extent",
    )
    plot_metric(
        plot_dir / "mini_hexaplex_units_vs_normalized_axial_extent_vs_full.png",
        "normalized_axial_extent_vs_full",
        "axial extent / full axial extent",
        "Mini-hexaplex length vs normalized axial extent",
    )
    plot_feature_metric(
        plot_dir / "mini_hexaplex_units_vs_ratio_to_full_4p5_5p0.png",
        "ratio_to_full_length_4p5_5A",
        "4.5-5.0 A ratio to full length",
        "Mini-hexaplex length vs 4.5-5.0 A integrated intensity ratio",
    )

    path = plot_dir / "mini_hexaplex_helicity_vs_4p5_5p0_response.png"
    fig, axes = plt.subplots(2, 1, figsize=(8.0, 7.0), dpi=180, sharex=True)
    for family, family_rows in families.items():
        helicity_points = sorted(
            (
                safe_float(row.get("units_per_chain")),
                safe_float(row.get("helical_coherence_score")),
                row["variant_id"],
            )
            for row in family_rows
        )
        x_values = [point[0] for point in helicity_points if point[0] is not None and point[1] is not None]
        y_values = [point[1] for point in helicity_points if point[0] is not None and point[1] is not None]
        if x_values:
            axes[0].plot(x_values, y_values, marker="o", linewidth=1.2, label=family)
        response_points = []
        for unit_count, _score, variant_id in helicity_points:
            feature_row = feature_by_variant.get(variant_id)
            response = safe_float(feature_row.get("ratio_to_full_length_4p5_5A")) if feature_row else None
            if unit_count is not None and response is not None:
                response_points.append((unit_count, response))
        if response_points:
            axes[1].plot(
                [point[0] for point in response_points],
                [point[1] for point in response_points],
                marker="o",
                linewidth=1.2,
                label=family,
            )
    axes[0].set_ylabel("helical coherence score")
    axes[0].set_title("Helicity and 4.5-5.0 A response by mini-hexaplex length")
    axes[1].set_xlabel("base/GLU units per chain")
    axes[1].set_ylabel("4.5-5.0 A ratio to full length")
    for ax in axes:
        ax.grid(True, alpha=0.25, linewidth=0.6)
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    outputs.append(path)

    path = plot_dir / "mini_hexaplex_coherent_turns_and_4p5_5p0_response.png"
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.6), dpi=180, sharex=True)
    for family, family_rows in families.items():
        extent_points = sorted(
            (
                safe_float(row.get("units_per_chain")),
                safe_float(row.get("coherent_helical_turns")),
                row["variant_id"],
            )
            for row in family_rows
        )
        x_values = [point[0] for point in extent_points if point[0] is not None and point[1] is not None]
        y_values = [point[1] for point in extent_points if point[0] is not None and point[1] is not None]
        if x_values:
            axes[0].plot(x_values, y_values, marker="o", linewidth=1.2, label=family)
        response_points = []
        for unit_count, _extent, variant_id in extent_points:
            feature_row = feature_by_variant.get(variant_id)
            response = safe_float(feature_row.get("ratio_to_full_length_4p5_5A")) if feature_row else None
            if unit_count is not None and response is not None:
                response_points.append((unit_count, response))
        if response_points:
            axes[1].plot(
                [point[0] for point in response_points],
                [point[1] for point in response_points],
                marker="o",
                linewidth=1.2,
                label=family,
            )
    axes[0].set_xlabel("base/GLU units per chain")
    axes[0].set_ylabel("coherent helical turns")
    axes[0].set_title("Panel A: coherent helical extent")
    axes[1].set_xlabel("base/GLU units per chain")
    axes[1].set_ylabel("4.5-5.0 A ratio to full")
    axes[1].set_title("Panel B: 4.5-5.0 A response")
    for ax in axes:
        ax.grid(True, alpha=0.25, linewidth=0.6)
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    outputs.append(path)
    return outputs


def _markdown_table(rows: list[dict[str, str]], columns: list[str]) -> str:
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join("---" for _ in columns) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(row.get(column, "") for column in columns) + " |")
    return "\n".join(lines)


def trend_note(rows: list[dict[str, str]], key: str) -> str:
    families = {
        "central": [row for row in rows if row["variant_id"].startswith("central")],
        "lower-end": [row for row in rows if row["variant_id"].startswith("lower_end")],
    }
    notes: list[str] = []
    for family, family_rows in families.items():
        points = sorted(
            (
                safe_float(row.get("units_per_chain")),
                safe_float(row.get(key)),
            )
            for row in family_rows
        )
        values = [point for point in points if point[0] is not None and point[1] is not None]
        if len(values) < 2:
            notes.append(f"{family}: insufficient data")
            continue
        y_values = [value for _unit_count, value in values]
        if max(y_values) - min(y_values) < 1e-4:
            notes.append(f"{family}: saturated or effectively flat across sampled counts")
            continue
        increases = sum(1 for index in range(1, len(values)) if values[index][1] >= values[index - 1][1])
        if increases == len(values) - 1:
            notes.append(f"{family}: monotonic nondecreasing")
        elif values[-1][1] > values[0][1]:
            notes.append(f"{family}: overall increase with local variation")
        else:
            notes.append(f"{family}: no clear increase")
    return "; ".join(notes)


def write_report(rows: list[dict[str, str]], feature_rows: list[dict[str, str]], plot_paths: list[Path], path: Path) -> None:
    ordered_rows = sorted(rows, key=lambda row: (safe_float(row.get("units_per_chain")) or 0.0, row["variant_id"]))
    feature_by_variant = {row["variant_id"]: row for row in feature_rows}
    combined_rows = []
    for row in ordered_rows:
        copied = dict(row)
        feature_row = feature_by_variant.get(row["variant_id"], {})
        copied["ratio_to_full_length_4p5_5A"] = feature_row.get("ratio_to_full_length_4p5_5A", "")
        combined_rows.append(copied)
    coherent_counts = sorted(
        {
            int(safe_float(row.get("units_per_chain")) or 0)
            for row in ordered_rows
            if row.get("structural_coherence_flag") == "coherent" and row["variant_id"] != FULL_BASELINE_ID
        }
    )
    shortest_coherent = min(coherent_counts) if coherent_counts else "not established"
    table_columns = [
        "variant_id",
        "units_per_chain",
        "structural_coherence_flag",
        "axial_extent_A",
        "normalized_axial_extent_vs_full",
        "coherent_helical_turns",
        "normalized_coherent_helical_turns_vs_full",
        "helical_coherence_score",
        "mean_helical_r2",
        "mean_angular_residual_deg",
        "std_twist_per_unit_deg",
        "mean_pitch_A",
        "ratio_to_full_length_4p5_5A",
    ]
    lines = [
        "# Mini-Hexaplex Helicity Report",
        "",
        "## Purpose",
        "",
        "This geometry-only analysis asks how helix-like each coordinate-truncated mini-hexaplex is as a function of base/GLU unit count, without rerunning diffraction.",
        "",
        "## Metric Definition",
        "",
        "- The common axis is fitted once from the full cleaned six-chain baseline using the same principal-axis logic as the mini truncation workflow.",
        "- Each repeat unit is represented by the base-like CYP/MEP residue CA atom when present, otherwise by the CYP/MEP heavy-atom centroid.",
        "- Representative points are converted to cylindrical coordinates around the full-baseline fitted axis.",
        "- For each chain, theta is unwrapped after sorting by axial coordinate, then theta = omega*z + phi is fit by least squares.",
        "- The raw helical score is mean_helical_r2 across chains, clipped to 0-1 as helical_coherence_score; it is retained as a local helix-fit diagnostic, not the main length-response metric.",
        "- Because these mini-hexaplexes are coordinate truncations from an ideal helix, theta-vs-z linearity saturates at 1.0 and cannot distinguish structural length across the sampled truncations.",
        "- The main structural-length metric is coherent_helical_turns = axial_extent_A / mean_pitch_A, with normalized_coherent_helical_turns_vs_full comparing each truncation against the full-length baseline.",
        "- normalized_axial_extent_vs_full uses the existing geometry axial extent where available and the full-baseline PDB axial extent as the denominator.",
        "- Six-strand phase coherence is estimated from angular spacing consistency at comparable axial ranks; it is secondary to the raw chain R2 score.",
        "",
        "## Helicity Summary",
        "",
        _markdown_table(combined_rows, table_columns),
        "",
        "## Plots",
        "",
    ]
    lines.extend(f"- {plot_path}" for plot_path in plot_paths)
    lines.extend(
        [
            "",
            "## Conservative Interpretation",
            "",
            f"- Raw helix-linearity trend: {trend_note(ordered_rows, 'helical_coherence_score')}.",
            f"- Coherent helical turns trend: {trend_note(ordered_rows, 'coherent_helical_turns')}.",
            f"- Normalized axial extent trend: {trend_note(ordered_rows, 'normalized_axial_extent_vs_full')}.",
            f"- Shortest unit count still flagged geometrically coherent by the existing structural summary: {shortest_coherent}.",
            "- In this coordinate-truncation set, the raw theta-vs-z R2 score is saturated at 1.0 for all sampled lengths, so it does not by itself define the minimum meaningful mini-hexaplex length.",
            "- The useful structural metric is coherent helical extent: axial length expressed as fitted helical turns, interpreted alongside normalized axial extent and the existing structural_coherence_flag.",
            "- The 4.5-5.0 A intensity response rises with length while the raw helix-linearity score is already saturated; the diffraction response therefore appears to track structural extent more than local helix fit quality.",
            "- The metric is most useful for comparing central and lower-end truncation families across length, not for proving independent stability.",
            "",
            "## Limitations",
            "",
            "- These are coordinate truncations only; no relaxation or minimization was performed.",
            "- The metric depends on representative atom choice.",
            "- Anti-parallel chains require careful angle unwrapping, and short chains can have inflated linear-fit scores.",
            "- This analysis does not prove independent stability or a formation mechanism.",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    baseline_atoms = load_pdb_atoms(args.baseline_pdb)
    axis_origin, axis = infer_axis(baseline_atoms)
    manifest_by_variant = {row["variant_id"]: row for row in read_csv_rows(args.manifest)}
    geometry_by_variant = {row["variant_id"]: row for row in read_csv_rows(args.geometry)}
    feature_rows = read_csv_rows(args.feature_summary)
    rows = [
        analyze_structure(variant_id, path, axis_origin, axis, manifest_by_variant, geometry_by_variant)
        for variant_id, path in structure_paths(args.structures_dir, args.baseline_pdb)
    ]
    add_normalized_extent_metrics(rows)
    write_csv(args.out_csv, rows, SUMMARY_COLUMNS)
    plot_paths = write_plots(rows, feature_rows, args.plot_dir)
    write_report(rows, feature_rows, plot_paths, args.out_report)
    print(f"Wrote {args.out_csv}")
    for plot_path in plot_paths:
        print(f"Wrote {plot_path}")
    print(f"Wrote {args.out_report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
