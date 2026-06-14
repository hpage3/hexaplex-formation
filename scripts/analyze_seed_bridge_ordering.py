#!/usr/bin/env python3
"""Estimate SB-inspired feature ordering along loose-to-formed seed bridges."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import Counter, defaultdict
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
    mean_point,
    power_iteration_principal_axis,
    project_point_to_axis,
    vector_sub,
)
from hexaplex_formation.pdb_utils import PDBAtom, chain_ids, dedupe_exact_atoms, heavy_atoms, load_pdb_atoms  # noqa: E402


BASE_FEATURE_SPECS = [
    ("compactness_score", "compactness_score", "larger compactness_score is more formed-like"),
    ("contact_fraction_vs_target", "contact_fraction_vs_target", "larger target-contact recovery is more formed-like"),
    (
        "CYP_MEP_contact_fraction_vs_target",
        "CYP_MEP_contact_fraction_vs_target",
        "larger CYP/MEP-involving target-contact recovery is more formed-like",
    ),
    ("axial_register_score", "axial_register_score", "larger axial register score is more formed-like"),
    ("rmsd_formedness_score", "RMSD_to_formed_seed_A", "RMSD converted to 1 / (1 + RMSD_to_formed_seed_A / 5)"),
    ("seed_formation_score", "seed_formation_score", "larger exploratory seed formation score is more formed-like"),
]

LEGACY_PHASE_SPEC = (
    "angular_phase_order_score",
    "angular_phase_order_score",
    "legacy angular phase metric from seed-formation CSV; retained only for comparison, not a physical ordering claim",
)
REFINED_PHASE_SPEC = (
    "refined_angular_phase_score",
    "refined_angular_phase_score",
    "chain-label-aware angular agreement with the formed seed after optimizing one global rotation around the formed fitted axis",
)

FEATURE_SPECS = BASE_FEATURE_SPECS + [LEGACY_PHASE_SPEC]

PAIR_COLUMNS = [
    "unit_count",
    "pair_id",
    "loose_sample_id",
    "formed_sample_id",
    "coupling_weight",
    "matching_cost",
    "endpoint_distance",
    "method_used",
    "notes",
]

PATH_COLUMNS = [
    "unit_count",
    "pair_id",
    "t",
    "compactness_score",
    "contact_fraction_vs_target",
    "CYP_MEP_contact_fraction_vs_target",
    "axial_register_score",
    "angular_phase_order_score",
    "refined_angular_phase_score",
    "rmsd_formedness_score",
    "seed_formation_score",
    "notes",
]

ACTIVATION_COLUMNS = [
    "unit_count",
    "feature",
    "threshold_mode",
    "threshold_fraction",
    "loose_mean",
    "formed_mean",
    "threshold_value",
    "mean_activation_time",
    "median_activation_time",
    "std_activation_time",
    "fraction_crossed",
    "mean_rank",
    "median_rank",
    "first_feature_frequency",
    "second_feature_frequency",
    "interpretation_note",
    "warnings",
]

ORDERING_COLUMNS = [
    "unit_count",
    "first_feature",
    "second_feature",
    "third_feature",
    "ordering_summary",
    "ordering_confidence",
    "interpretation_note",
    "warnings",
]


@dataclass(frozen=True)
class Sample:
    sample_id: str
    raw: dict[str, str]
    features: dict[str, float]


@dataclass(frozen=True)
class EndpointPair:
    pair_id: int
    loose: Sample
    formed: Sample
    coupling_weight: float | None
    matching_cost: float
    endpoint_distance: float
    method_used: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-csv", type=Path, default=Path("outputs/metrics/seed_formation_order_parameters.csv"))
    parser.add_argument("--contact-network-summary", type=Path, default=Path("outputs/metrics/seed_contact_network_summary.csv"))
    parser.add_argument("--structures-dir", type=Path, default=Path("outputs/mini_hexaplex/structures"))
    parser.add_argument("--ensemble-dir", type=Path, default=Path("outputs/seed_formation/ensembles"))
    parser.add_argument("--unit-counts", default="4,5,6,7")
    parser.add_argument("--loose-ensemble-type", default="loose_initial")
    parser.add_argument("--angular-phase-mode", choices=["legacy", "refined", "both"], default="refined")
    parser.add_argument("--coordinate-backed-only", action="store_true")
    parser.add_argument("--n-time-points", type=int, default=21)
    parser.add_argument("--time-points", default="")
    parser.add_argument("--threshold-mode", default="formed_fraction", choices=["formed_fraction"])
    parser.add_argument("--formed-fraction", type=float, default=0.75)
    parser.add_argument("--formed-fractions", default="")
    parser.add_argument("--sinkhorn-epsilon", type=float, default=0.75)
    parser.add_argument("--sinkhorn-iterations", type=int, default=500)
    parser.add_argument("--max-pairs", type=int, default=0)
    parser.add_argument("--pairs-csv", type=Path, default=Path("outputs/metrics/seed_bridge_endpoint_pairs.csv"))
    parser.add_argument("--paths-csv", type=Path, default=Path("outputs/metrics/seed_bridge_order_parameter_paths.csv"))
    parser.add_argument("--activation-csv", type=Path, default=Path("outputs/metrics/seed_bridge_activation_summary.csv"))
    parser.add_argument("--ordering-csv", type=Path, default=Path("outputs/metrics/seed_bridge_feature_ordering.csv"))
    parser.add_argument("--plot-dir", type=Path, default=Path("outputs/seed_formation/plots"))
    parser.add_argument("--out-report", type=Path, default=Path("outputs/reports/seed_bridge_ordering_report.md"))
    return parser.parse_args()


def selected_feature_specs(angular_phase_mode: str) -> list[tuple[str, str, str]]:
    specs = list(BASE_FEATURE_SPECS)
    if angular_phase_mode in {"legacy", "both"}:
        specs.append(LEGACY_PHASE_SPEC)
    if angular_phase_mode in {"refined", "both"}:
        specs.append(REFINED_PHASE_SPEC)
    return specs


def parse_int_list(value: str) -> list[int]:
    values = [int(token.strip()) for token in value.split(",") if token.strip()]
    if not values:
        raise ValueError("At least one unit count is required")
    return values


def parse_float_list(value: str) -> list[float]:
    values = [float(token.strip()) for token in value.split(",") if token.strip()]
    if not values:
        raise ValueError("At least one numeric value is required")
    return values


def time_grid(args: argparse.Namespace) -> list[float]:
    if args.time_points.strip():
        values = parse_float_list(args.time_points)
    else:
        if args.n_time_points < 2:
            raise ValueError("--n-time-points must be at least 2")
        values = [index / (args.n_time_points - 1) for index in range(args.n_time_points)]
    if values[0] != 0.0 or values[-1] != 1.0:
        raise ValueError("Bridge time points must start at 0.0 and end at 1.0")
    return values


def threshold_fractions(args: argparse.Namespace) -> list[float]:
    if args.formed_fractions.strip():
        values = parse_float_list(args.formed_fractions)
    else:
        values = [args.formed_fraction]
    for value in values:
        if value < 0.0 or value > 1.0:
            raise ValueError("formed fractions must be between 0 and 1")
    return values


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
        raise FileNotFoundError(f"Input CSV not found: {path}")
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def atom_xyz(atom: PDBAtom) -> tuple[float, float, float]:
    return atom.x, atom.y, atom.z


def infer_axis(atoms: list[PDBAtom]) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    axis_atoms = [atom for atom in heavy_atoms(atoms) if atom.atom_name.strip().upper() in {"N", "CA", "C", "O", "OXT"}]
    if len(axis_atoms) < 3:
        axis_atoms = heavy_atoms(atoms)
    if len(axis_atoms) < 3:
        raise ValueError("At least three heavy atoms are required to infer phase axis")
    origin = mean_point(atom_xyz(atom) for atom in axis_atoms)
    axis = power_iteration_principal_axis(covariance_matrix_3d(atom_xyz(atom) for atom in axis_atoms))
    return origin, axis


def group_atoms_by_chain(atoms: list[PDBAtom]) -> dict[str, list[PDBAtom]]:
    grouped: dict[str, list[PDBAtom]] = defaultdict(list)
    for atom in atoms:
        grouped[atom.chain_id].append(atom)
    return grouped


def circular_difference(a: float, b: float) -> float:
    return math.atan2(math.sin(a - b), math.cos(a - b))


def chain_angles_about_axis(
    atoms: list[PDBAtom],
    axis_origin: tuple[float, float, float],
    axis: tuple[float, float, float],
) -> dict[str, float]:
    basis_u, basis_v = build_perpendicular_basis(axis)
    angles: dict[str, float] = {}
    for chain_id, chain_atoms in group_atoms_by_chain(atoms).items():
        centroid = mean_point(atom_xyz(atom) for atom in heavy_atoms(chain_atoms))
        _, projected = project_point_to_axis(centroid, axis_origin, axis)
        radial = vector_sub(centroid, projected)
        if distance(centroid, projected) < 1e-9:
            continue
        angles[chain_id] = math.atan2(dot(radial, basis_v), dot(radial, basis_u))
    return angles


def refined_angular_phase_score(sample_atoms: list[PDBAtom], reference_atoms: list[PDBAtom]) -> float:
    """Chain-label-aware angular agreement after one global rotation around the reference axis."""

    reference_origin, reference_axis = infer_axis(reference_atoms)
    reference_angles = chain_angles_about_axis(reference_atoms, reference_origin, reference_axis)
    sample_angles = chain_angles_about_axis(sample_atoms, reference_origin, reference_axis)
    shared_chains = sorted(set(reference_angles).intersection(sample_angles))
    if len(shared_chains) < 3:
        raise ValueError("Need at least three shared chains for refined angular phase score")
    deltas = [circular_difference(sample_angles[chain_id], reference_angles[chain_id]) for chain_id in shared_chains]
    offset = math.atan2(sum(math.sin(delta) for delta in deltas), sum(math.cos(delta) for delta in deltas))
    residuals = [circular_difference(delta, offset) for delta in deltas]
    mean_cos = sum(math.cos(residual) for residual in residuals) / len(residuals)
    return max(0.0, min(1.0, (1.0 + mean_cos) / 2.0))


def sample_pdb_path(sample_id: str, ensemble_dir: Path) -> Path:
    return ensemble_dir / f"{sample_id}.pdb"


def add_refined_phase_values(
    rows: list[dict[str, str]],
    unit_counts: list[int],
    structures_dir: Path,
    ensemble_dir: Path,
) -> tuple[list[dict[str, str]], dict[int, str]]:
    reference_atoms: dict[int, list[PDBAtom]] = {}
    for unit in unit_counts:
        path = structures_dir / f"mini_hexaplex_central{unit}_units.pdb"
        if not path.exists():
            raise FileNotFoundError(f"Formed seed reference PDB not found: {path}")
        reference_atoms[unit] = dedupe_exact_atoms(load_pdb_atoms(path))
    warnings: dict[int, str] = {}
    augmented: list[dict[str, str]] = []
    missing_counts: Counter[int] = Counter()
    scored_counts: Counter[int] = Counter()
    for row in rows:
        copied = dict(row)
        unit = int(row["unit_count"])
        if unit in reference_atoms:
            path = sample_pdb_path(row["sample_id"], ensemble_dir)
            if path.exists():
                sample_atoms = dedupe_exact_atoms(load_pdb_atoms(path))
                copied["refined_angular_phase_score"] = format_float(
                    refined_angular_phase_score(sample_atoms, reference_atoms[unit])
                )
                scored_counts[unit] += 1
            else:
                if safe_float(copied.get("refined_angular_phase_score", "")) is not None:
                    scored_counts[unit] += 1
                else:
                    copied["refined_angular_phase_score"] = ""
                    missing_counts[unit] += 1
        augmented.append(copied)
    for unit in unit_counts:
        if missing_counts[unit]:
            warnings[unit] = (
                f"refined_angular_phase_score available for {scored_counts[unit]} coordinate-backed samples; "
                f"{missing_counts[unit]} rows lacked saved ensemble PDBs and were excluded when refined phase was selected"
            )
    return augmented, warnings


def rmsd_formedness_score(rmsd: float) -> float:
    return 1.0 / (1.0 + rmsd / 5.0)


def oriented_feature_value(row: dict[str, str], feature: str, source_column: str) -> float | None:
    raw = safe_float(row.get(source_column, ""))
    if raw is None:
        return None
    if feature == "rmsd_formedness_score":
        return rmsd_formedness_score(raw)
    return raw


def sample_from_row(row: dict[str, str], feature_specs: list[tuple[str, str, str]]) -> Sample | None:
    features: dict[str, float] = {}
    for feature, source_column, _ in feature_specs:
        value = oriented_feature_value(row, feature, source_column)
        if value is None:
            return None
        features[feature] = value
    return Sample(sample_id=row["sample_id"], raw=row, features=features)


def empirically_orient_unit_samples(
    loose: list[Sample],
    formed: list[Sample],
    features: list[str],
) -> tuple[list[Sample], list[Sample], list[str]]:
    """Validate that larger feature values point from loose toward formed endpoints."""

    warnings: list[str] = []
    for feature in features:
        loose_values = [sample.features[feature] for sample in loose]
        formed_values = [sample.features[feature] for sample in formed]
        if mean(formed_values) <= mean(loose_values):
            warnings.append(
                f"{feature} formed mean is not higher than loose mean; treated as non-activating for formed-like ordering"
            )
    return loose, formed, warnings


def samples_by_unit(
    rows: list[dict[str, str]],
    unit_counts: list[int],
    ensemble_type: str,
    feature_specs: list[tuple[str, str, str]],
    coordinate_backed_only: bool,
    ensemble_dir: Path,
) -> dict[int, list[Sample]]:
    grouped: dict[int, list[Sample]] = {unit: [] for unit in unit_counts}
    for row in rows:
        if row.get("ensemble_type") != ensemble_type:
            continue
        unit = int(row["unit_count"])
        if unit not in grouped:
            continue
        if coordinate_backed_only and not sample_pdb_path(row["sample_id"], ensemble_dir).exists():
            continue
        sample = sample_from_row(row, feature_specs)
        if sample is not None:
            grouped[unit].append(sample)
    return grouped


def mean(values: list[float]) -> float:
    return sum(values) / len(values)


def std(values: list[float]) -> float:
    if not values:
        return float("nan")
    m = mean(values)
    return math.sqrt(sum((value - m) ** 2 for value in values) / len(values))


def median(values: list[float]) -> float:
    ordered = sorted(values)
    if not ordered:
        return float("nan")
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[mid]
    return (ordered[mid - 1] + ordered[mid]) / 2.0


def pooled_feature_stats(loose: list[Sample], formed: list[Sample], features: list[str]) -> dict[str, tuple[float, float]]:
    stats: dict[str, tuple[float, float]] = {}
    for feature in features:
        values = [sample.features[feature] for sample in loose + formed]
        m = mean(values)
        s = std(values)
        stats[feature] = (m, s if s > 1e-12 else 1.0)
    return stats


def standardized_matrix(samples: list[Sample], features: list[str], stats: dict[str, tuple[float, float]]):
    import numpy as np

    matrix = np.zeros((len(samples), len(features)), dtype=float)
    for row_index, sample in enumerate(samples):
        for col_index, feature in enumerate(features):
            m, s = stats[feature]
            matrix[row_index, col_index] = (sample.features[feature] - m) / s
    return matrix


def cost_matrix(loose: list[Sample], formed: list[Sample], features: list[str], stats: dict[str, tuple[float, float]]):
    import numpy as np

    loose_matrix = standardized_matrix(loose, features, stats)
    formed_matrix = standardized_matrix(formed, features, stats)
    diff = loose_matrix[:, None, :] - formed_matrix[None, :, :]
    return np.sqrt(np.sum(diff * diff, axis=2))


def greedy_pairs_from_cost(cost, max_pairs: int | None = None) -> list[tuple[int, int, float, float | None]]:
    candidates = [(float(cost[i, j]), i, j) for i in range(cost.shape[0]) for j in range(cost.shape[1])]
    candidates.sort()
    used_loose: set[int] = set()
    used_formed: set[int] = set()
    pairs: list[tuple[int, int, float, float | None]] = []
    target_count = max_pairs or min(cost.shape[0], cost.shape[1])
    for value, i, j in candidates:
        if i in used_loose or j in used_formed:
            continue
        pairs.append((i, j, value, None))
        used_loose.add(i)
        used_formed.add(j)
        if len(pairs) >= target_count:
            break
    return pairs


def sinkhorn_coupling(cost, epsilon: float, iterations: int):
    import numpy as np

    if epsilon <= 0:
        raise ValueError("sinkhorn epsilon must be positive")
    n, m = cost.shape
    a = np.full(n, 1.0 / n)
    b = np.full(m, 1.0 / m)
    scaled = cost / max(epsilon, 1e-12)
    scaled -= scaled.min()
    kernel = np.exp(-scaled)
    kernel = np.maximum(kernel, 1e-300)
    u = np.ones(n)
    v = np.ones(m)
    for _ in range(iterations):
        u = a / (kernel @ v)
        v = b / (kernel.T @ u)
        if not (np.all(np.isfinite(u)) and np.all(np.isfinite(v))):
            raise FloatingPointError("non-finite Sinkhorn scaling")
    coupling = (u[:, None] * kernel) * v[None, :]
    if not np.all(np.isfinite(coupling)):
        raise FloatingPointError("non-finite Sinkhorn coupling")
    return coupling


def pairs_from_coupling(cost, coupling, max_pairs: int | None = None) -> list[tuple[int, int, float, float | None]]:
    candidates = [
        (float(coupling[i, j]), float(cost[i, j]), i, j)
        for i in range(coupling.shape[0])
        for j in range(coupling.shape[1])
    ]
    candidates.sort(reverse=True)
    used_loose: set[int] = set()
    used_formed: set[int] = set()
    pairs: list[tuple[int, int, float, float | None]] = []
    target_count = max_pairs or min(coupling.shape[0], coupling.shape[1])
    for weight, value, i, j in candidates:
        if i in used_loose or j in used_formed:
            continue
        pairs.append((i, j, value, weight))
        used_loose.add(i)
        used_formed.add(j)
        if len(pairs) >= target_count:
            break
    return pairs


def endpoint_matching(
    loose: list[Sample],
    formed: list[Sample],
    features: list[str],
    epsilon: float = 0.75,
    iterations: int = 500,
    max_pairs: int | None = None,
) -> tuple[list[EndpointPair], str, str]:
    stats = pooled_feature_stats(loose, formed, features)
    cost = cost_matrix(loose, formed, features, stats)
    try:
        coupling = sinkhorn_coupling(cost, epsilon, iterations)
        raw_pairs = pairs_from_coupling(cost, coupling, max_pairs)
        method = "sinkhorn_greedy_unique"
        note = f"entropic OT coupling epsilon={epsilon:g}; unique endpoint pairs extracted by descending coupling weight"
    except Exception as exc:
        raw_pairs = greedy_pairs_from_cost(cost, max_pairs)
        method = "greedy_minimum_cost_fallback"
        note = f"Sinkhorn failed ({exc}); used greedy minimum-cost matching"
    endpoint_pairs = [
        EndpointPair(
            pair_id=pair_index,
            loose=loose[i],
            formed=formed[j],
            matching_cost=value,
            endpoint_distance=value,
            coupling_weight=weight,
            method_used=method,
        )
        for pair_index, (i, j, value, weight) in enumerate(raw_pairs)
    ]
    return endpoint_pairs, method, note


def interpolate_feature(loose_value: float, formed_value: float, t: float) -> float:
    return loose_value + t * (formed_value - loose_value)


def threshold_value(loose_mean: float, formed_mean: float, fraction: float) -> float:
    return loose_mean + fraction * (formed_mean - loose_mean)


def activation_time(times: list[float], values: list[float], threshold: float) -> float:
    for t, value in zip(times, values):
        if value >= threshold:
            return t
    return float("nan")


def average_ranks(values: dict[str, float]) -> dict[str, float]:
    finite = [(feature, value) for feature, value in values.items() if math.isfinite(value)]
    finite.sort(key=lambda item: item[1])
    ranks: dict[str, float] = {feature: float("nan") for feature in values}
    index = 0
    while index < len(finite):
        end = index + 1
        while end < len(finite) and abs(finite[end][1] - finite[index][1]) < 1e-12:
            end += 1
        avg_rank = (index + 1 + end) / 2.0
        for tie_index in range(index, end):
            ranks[finite[tie_index][0]] = avg_rank
        index = end
    return ranks


def feature_order_from_medians(medians: dict[str, float]) -> list[str]:
    return [
        feature
        for feature, value in sorted(
            medians.items(),
            key=lambda item: (math.inf if not math.isfinite(item[1]) else item[1], item[0]),
        )
    ]


def build_paths(unit: int, pairs: list[EndpointPair], times: list[float], features: list[str]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for pair in pairs:
        for t in times:
            row = {
                "unit_count": str(unit),
                "pair_id": str(pair.pair_id),
                "t": format_float(t, 3),
                "notes": "linear interpolation in oriented order-parameter space; not atomistic dynamics",
            }
            for feature in features:
                row[feature] = format_float(interpolate_feature(pair.loose.features[feature], pair.formed.features[feature], t))
            rows.append(row)
    return rows


def build_pair_rows(unit: int, pairs: list[EndpointPair], note: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for pair in pairs:
        rows.append(
            {
                "unit_count": str(unit),
                "pair_id": str(pair.pair_id),
                "loose_sample_id": pair.loose.sample_id,
                "formed_sample_id": pair.formed.sample_id,
                "coupling_weight": format_float(pair.coupling_weight, 10),
                "matching_cost": format_float(pair.matching_cost),
                "endpoint_distance": format_float(pair.endpoint_distance),
                "method_used": pair.method_used,
                "notes": note,
            }
        )
    return rows


def activation_records_for_unit(
    unit: int,
    loose: list[Sample],
    formed: list[Sample],
    pairs: list[EndpointPair],
    times: list[float],
    features: list[str],
    fractions: list[float],
    warnings: list[str] | None = None,
) -> tuple[list[dict[str, str]], dict[float, dict[str, float]], dict[float, dict[int, dict[str, float]]]]:
    activation_rows: list[dict[str, str]] = []
    median_by_fraction: dict[float, dict[str, float]] = {}
    ranks_by_fraction: dict[float, dict[int, dict[str, float]]] = {}
    loose_means = {feature: mean([sample.features[feature] for sample in loose]) for feature in features}
    formed_means = {feature: mean([sample.features[feature] for sample in formed]) for feature in features}
    non_activating = {feature for feature in features if formed_means[feature] <= loose_means[feature]}
    for fraction in fractions:
        thresholds = {
            feature: (
                float("nan")
                if feature in non_activating
                else threshold_value(loose_means[feature], formed_means[feature], fraction)
            )
            for feature in features
        }
        pair_times: dict[int, dict[str, float]] = defaultdict(dict)
        for pair in pairs:
            for feature in features:
                if feature in non_activating:
                    pair_times[pair.pair_id][feature] = float("nan")
                else:
                    values = [interpolate_feature(pair.loose.features[feature], pair.formed.features[feature], t) for t in times]
                    pair_times[pair.pair_id][feature] = activation_time(times, values, thresholds[feature])
        pair_ranks = {pair_id: average_ranks(values) for pair_id, values in pair_times.items()}
        ranks_by_fraction[fraction] = pair_ranks
        median_by_fraction[fraction] = {}
        first_counts: Counter[str] = Counter()
        second_counts: Counter[str] = Counter()
        for pair_id, values in pair_times.items():
            ordered = feature_order_from_medians(values)
            if ordered:
                first_counts[ordered[0]] += 1
            if len(ordered) > 1:
                second_counts[ordered[1]] += 1
        for feature in features:
            values = [pair_times[pair.pair_id][feature] for pair in pairs]
            crossed = [value for value in values if math.isfinite(value)]
            ranks = [pair_ranks[pair.pair_id][feature] for pair in pairs if math.isfinite(pair_ranks[pair.pair_id][feature])]
            median_time = median(crossed) if crossed else float("nan")
            median_by_fraction[fraction][feature] = median_time
            activation_rows.append(
                {
                    "unit_count": str(unit),
                    "feature": feature,
                    "threshold_mode": "formed_fraction",
                    "threshold_fraction": format_float(fraction, 2),
                    "loose_mean": format_float(loose_means[feature]),
                    "formed_mean": format_float(formed_means[feature]),
                    "threshold_value": format_float(thresholds[feature]),
                    "mean_activation_time": format_float(mean(crossed) if crossed else None),
                    "median_activation_time": format_float(median_time),
                    "std_activation_time": format_float(std(crossed) if crossed else None),
                    "fraction_crossed": format_float(len(crossed) / len(values) if values else None),
                    "mean_rank": format_float(mean(ranks) if ranks else None),
                    "median_rank": format_float(median(ranks) if ranks else None),
                    "first_feature_frequency": format_float(first_counts[feature] / len(pairs) if pairs else None),
                    "second_feature_frequency": format_float(second_counts[feature] / len(pairs) if pairs else None),
                    "interpretation_note": "earlier activation means this oriented order parameter reaches the formed-fraction threshold sooner along paired linear bridge paths",
                    "warnings": "; ".join(warnings or []),
                }
            )
    return activation_rows, median_by_fraction, ranks_by_fraction


def ordering_confidence(first_frequency: float, first_median: float, second_median: float) -> str:
    if first_frequency >= 0.6 and (not math.isfinite(second_median) or second_median - first_median >= 0.05):
        return "moderate"
    if first_frequency >= 0.4:
        return "low_to_moderate"
    return "ambiguous"


def ordering_row(unit: int, activation_rows: list[dict[str, str]], primary_fraction: float) -> dict[str, str]:
    relevant = [
        row
        for row in activation_rows
        if int(row["unit_count"]) == unit and abs(float(row["threshold_fraction"]) - primary_fraction) < 1e-9
    ]
    medians = {row["feature"]: safe_float(row["median_activation_time"]) or float("nan") for row in relevant}
    ordered = feature_order_from_medians(medians)
    first = ordered[0] if ordered else ""
    second = ordered[1] if len(ordered) > 1 else ""
    third = ordered[2] if len(ordered) > 2 else ""
    first_row = next((row for row in relevant if row["feature"] == first), None)
    first_frequency = safe_float(first_row["first_feature_frequency"]) if first_row else None
    first_median = medians.get(first, float("nan"))
    second_median = medians.get(second, float("nan"))
    confidence = ordering_confidence(first_frequency or 0.0, first_median, second_median)
    summary = " -> ".join(ordered[:5])
    return {
        "unit_count": str(unit),
        "first_feature": first,
        "second_feature": second,
        "third_feature": third,
        "ordering_summary": summary,
        "ordering_confidence": confidence,
        "interpretation_note": "ordering is based on median activation time at the primary formed-fraction threshold",
        "warnings": "" if confidence != "ambiguous" else "first-feature frequencies and activation times are close; treat ordering as ambiguous",
    }


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_plots(
    path_rows: list[dict[str, str]],
    activation_rows: list[dict[str, str]],
    ordering_rows: list[dict[str, str]],
    features: list[str],
    units: list[int],
    primary_fraction: float,
    plot_dir: Path,
) -> list[Path]:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except Exception as exc:  # pragma: no cover - optional dependency
        print(f"WARNING: matplotlib unavailable; skipping plots: {exc}", file=sys.stderr)
        return []

    plot_dir.mkdir(parents=True, exist_ok=True)
    plot_paths: list[Path] = []
    feature_labels = {
        "compactness_score": "compactness",
        "contact_fraction_vs_target": "total contacts",
        "CYP_MEP_contact_fraction_vs_target": "CYP/MEP contacts",
        "axial_register_score": "axial register",
        "angular_phase_order_score": "angular phase",
        "refined_angular_phase_score": "refined angular phase",
        "rmsd_formedness_score": "RMSD formedness",
        "seed_formation_score": "seed score",
    }
    for unit in units:
        unit_rows = [row for row in path_rows if int(row["unit_count"]) == unit]
        times = sorted({float(row["t"]) for row in unit_rows})
        fig, ax = plt.subplots(figsize=(8.0, 5.0))
        for feature in features:
            means = []
            start_values = [safe_float(row[feature]) for row in unit_rows if float(row["t"]) == 0.0]
            end_values = [safe_float(row[feature]) for row in unit_rows if float(row["t"]) == 1.0]
            start = mean([value for value in start_values if value is not None])
            end = mean([value for value in end_values if value is not None])
            denom = end - start
            for t in times:
                values = [safe_float(row[feature]) for row in unit_rows if abs(float(row["t"]) - t) < 1e-9]
                value_mean = mean([value for value in values if value is not None])
                means.append((value_mean - start) / denom if abs(denom) > 1e-12 else 0.0)
            ax.plot(times, means, marker="o", linewidth=1.5, label=feature_labels.get(feature, feature))
        ax.set_xlabel("Bridge time t")
        ax.set_ylabel("Normalized formedness")
        ax.set_title(f"{unit}-unit bridge mean trajectories")
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=8)
        fig.tight_layout()
        path = plot_dir / f"seed_bridge_mean_trajectories_unit{unit}.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        plot_paths.append(path)

        relevant = [
            row
            for row in activation_rows
            if int(row["unit_count"]) == unit and abs(float(row["threshold_fraction"]) - primary_fraction) < 1e-9
        ]
        fig, ax = plt.subplots(figsize=(8.0, 4.8))
        values = [safe_float(row["median_activation_time"]) for row in relevant]
        labels = [feature_labels.get(row["feature"], row["feature"]) for row in relevant]
        ax.bar(range(len(labels)), [value if value is not None else math.nan for value in values])
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=35, ha="right")
        ax.set_ylabel("Median activation time")
        ax.set_title(f"{unit}-unit activation times")
        fig.tight_layout()
        path = plot_dir / f"seed_bridge_activation_times_unit{unit}.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        plot_paths.append(path)

    heatmap = np.full((len(units), len(features)), np.nan)
    for row in activation_rows:
        fraction = float(row["threshold_fraction"])
        if abs(fraction - primary_fraction) > 1e-9:
            continue
        unit = int(row["unit_count"])
        if unit not in units or row["feature"] not in features:
            continue
        heatmap[units.index(unit), features.index(row["feature"])] = safe_float(row["median_activation_time"]) or math.nan
    fig, ax = plt.subplots(figsize=(9.5, 4.5))
    image = ax.imshow(heatmap, vmin=0.0, vmax=1.0, cmap="viridis_r")
    ax.set_yticks(range(len(units)))
    ax.set_yticklabels([str(unit) for unit in units])
    ax.set_xticks(range(len(features)))
    ax.set_xticklabels([feature_labels.get(feature, feature) for feature in features], rotation=35, ha="right")
    ax.set_ylabel("Unit count")
    fig.colorbar(image, ax=ax, label="Median activation time")
    fig.tight_layout()
    path = plot_dir / "seed_bridge_feature_ordering_heatmap.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    plot_paths.append(path)

    fig, ax = plt.subplots(figsize=(9.5, 5.0))
    x = np.arange(len(features))
    width = 0.8 / len(units)
    for index, unit in enumerate(units):
        values = []
        for feature in features:
            row = next(
                (
                    row
                    for row in activation_rows
                    if int(row["unit_count"]) == unit
                    and row["feature"] == feature
                    and abs(float(row["threshold_fraction"]) - primary_fraction) < 1e-9
                ),
                None,
            )
            values.append(safe_float(row["first_feature_frequency"]) if row else 0.0)
        ax.bar(x + index * width - 0.4 + width / 2.0, values, width=width, label=str(unit))
    ax.set_xticks(x)
    ax.set_xticklabels([feature_labels.get(feature, feature) for feature in features], rotation=35, ha="right")
    ax.set_ylabel("First-feature frequency")
    ax.legend(title="units")
    fig.tight_layout()
    path = plot_dir / "seed_bridge_first_feature_frequency.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    plot_paths.append(path)

    compare_units = [unit for unit in [4, 7] if unit in units]
    if len(compare_units) == 2:
        fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.5), sharey=True)
        for ax, unit in zip(axes, compare_units):
            unit_rows = [row for row in path_rows if int(row["unit_count"]) == unit]
            times = sorted({float(row["t"]) for row in unit_rows})
            for feature in features:
                start = mean([float(row[feature]) for row in unit_rows if float(row["t"]) == 0.0])
                end = mean([float(row[feature]) for row in unit_rows if float(row["t"]) == 1.0])
                denom = end - start
                means = []
                for t in times:
                    value_mean = mean([float(row[feature]) for row in unit_rows if abs(float(row["t"]) - t) < 1e-9])
                    means.append((value_mean - start) / denom if abs(denom) > 1e-12 else 0.0)
                ax.plot(times, means, label=feature_labels.get(feature, feature))
            ax.set_title(f"{unit} units")
            ax.set_xlabel("Bridge time t")
            ax.grid(True, alpha=0.25)
        axes[0].set_ylabel("Normalized formedness")
        axes[1].legend(fontsize=8)
        fig.tight_layout()
        path = plot_dir / "seed_bridge_mean_trajectories_4_vs_7.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        plot_paths.append(path)
    return plot_paths


def write_report(
    units: list[int],
    features: list[str],
    pair_rows: list[dict[str, str]],
    activation_rows: list[dict[str, str]],
    ordering_rows: list[dict[str, str]],
    plot_paths: list[Path],
    threshold_fractions_used: list[float],
    primary_fraction: float,
    orientation_warnings: dict[int, list[str]],
    feature_specs: list[tuple[str, str, str]],
    coordinate_warnings: dict[int, str],
    angular_phase_mode: str,
    loose_ensemble_type: str,
    report_path: Path,
) -> None:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    method_counts = Counter(row["method_used"] for row in pair_rows)
    lines = [
        "# Seed Bridge Feature-Ordering Analysis",
        "",
        "## Purpose",
        "",
        f"This is a first-pass, SB-inspired order-parameter bridge analysis for the loose-to-formed mini-hexaplex seed transition. It estimates which order parameters become formed-like earliest when `{loose_ensemble_type}` endpoints are paired to `formed_perturbed` endpoints by a low-cost entropic transport coupling in feature space.",
        "",
        "This is not a full atomistic Schrödinger bridge. It does not include molecular dynamics, a force field, solvent, counterions, or physical transition probabilities.",
        "",
        "## Inputs",
        "",
        "- `outputs/metrics/seed_formation_order_parameters.csv`",
        f"- Existing `{loose_ensemble_type}` and `formed_perturbed` ensembles summarized in that CSV.",
        "- Contact-network outputs are used as prior context only, not as path constraints in this first bridge script.",
        "",
        "## Endpoint Pairing",
        "",
        f"Endpoint pairs were built per unit count using methods: {dict(method_counts)}. The primary method is Sinkhorn entropic optimal transport in standardized oriented-feature space, followed by unique pair extraction from highest coupling weights. Greedy minimum-cost matching is recorded if Sinkhorn fails.",
        "",
        "## Feature Transformations",
        "",
        "All selected features are oriented so larger means more formed-like:",
    ]
    for feature, source, note in feature_specs:
        if feature in features:
            lines.append(f"- `{feature}` from `{source}`: {note}.")
    lines.extend(
        [
            "",
            f"Angular phase mode: `{angular_phase_mode}`. The legacy `angular_phase_order_score` is retained only as a comparison coordinate and should not be interpreted as a physical ordering claim. The refined coordinate measures chain-label-aware angular agreement to the formed seed after optimizing one global rotation around the formed fitted axis.",
        ]
    )
    if coordinate_warnings:
        lines.extend(["", "Coordinate-backed sample availability:"])
        for unit in sorted(coordinate_warnings):
            lines.append(f"- unit {unit}: {coordinate_warnings[unit]}")
    if any(orientation_warnings.values()):
        lines.extend(
            [
                "",
                "Endpoint distribution checks found feature(s) where the formed_perturbed mean was not higher than the loose_initial mean. Those features were retained for endpoint matching but treated as non-activating in formed-like ordering:",
            ]
        )
        for unit in sorted(orientation_warnings):
            for warning in orientation_warnings[unit]:
                lines.append(f"- unit {unit}: {warning}")
    lines.extend(
        [
            "",
            "## Thresholds And Activation Times",
            "",
            "For each feature, the formed-like threshold is `loose_mean + fraction * (formed_mean - loose_mean)`. Activation time is the first bridge time where the paired linear order-parameter path reaches that threshold. Threshold fractions evaluated: "
            + ", ".join(format_float(value, 2) for value in threshold_fractions_used)
            + ".",
            "",
            "## Feature Ordering",
            "",
            "| unit_count | first | second | third | confidence | summary |",
            "|---:|---|---|---|---|---|",
        ]
    )
    for row in ordering_rows:
        lines.append(
            f"| {row['unit_count']} | {row['first_feature']} | {row['second_feature']} | {row['third_feature']} | {row['ordering_confidence']} | {row['ordering_summary']} |"
        )
    lines.extend(["", f"Ordering table uses the primary threshold fraction `{primary_fraction:.2f}`.", "", "## Plots", ""])
    for path in plot_paths:
        lines.append(f"- `{path}`")
    lines.extend(
        [
            "",
            "## Conservative Interpretation",
            "",
            "Earliest activation in this workflow means earliest crossing along a paired linear interpolation in order-parameter space. It should be read as a candidate ordering for later SB collective-variable design, not as proof of physical time ordering.",
            "",
            "If compactness or RMSD formedness activates early, closure-like coordinates may be good bridge progress variables. If CYP/MEP contact recovery activates early or frequently, base/contact specificity may need to be represented explicitly. Axial register and angular phase order are especially useful if they activate later or distinguish 7-unit behavior from 4-6 units, because they encode register and six-strand phasing beyond compactness.",
            "",
            (
                "For `angular_randomized_loose_initial`, refined angular phase should be interpreted as a test coordinate for deliberately broken chain-label angular placement. If it becomes early, angular phasing is informative under randomized loose conditions; if contacts, compactness, register, or RMSD formedness remain first, those coordinates are still the stronger bridge-ordering candidates."
                if loose_ensemble_type == "angular_randomized_loose_initial"
                else "For preserved-angular loose ensembles, angular phase coordinates may already be close to formed-like and should not be interpreted as appearing early unless endpoint means show genuine loose-to-formed separation."
            ),
            "",
            "The known contact-network result that 4 units is already six-chain connected means this bridge analysis should emphasize ordering among redundancy, contact recovery, register, and phase-order metrics rather than binary connectivity.",
            "",
            "## Limitations",
            "",
            "- Order-parameter interpolation only.",
            "- No atomistic dynamics.",
            "- No solvent or counterions.",
            "- No energy model.",
            "- Endpoint ensembles are synthetic.",
            "- Activation thresholds are exploratory.",
            "- Does not prove spontaneous formation mechanism.",
            "",
        ]
    )
    report_path.write_text("\n".join(lines), encoding="utf-8")


def run(args: argparse.Namespace) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    units = parse_int_list(args.unit_counts)
    times = time_grid(args)
    fractions = threshold_fractions(args)
    primary_fraction = args.formed_fraction if args.formed_fraction in fractions else fractions[0]
    if args.loose_ensemble_type == "angular_randomized_loose_initial":
        if args.pairs_csv == Path("outputs/metrics/seed_bridge_endpoint_pairs.csv"):
            args.pairs_csv = Path("outputs/metrics/seed_bridge_endpoint_pairs_angular_randomized.csv")
        if args.paths_csv == Path("outputs/metrics/seed_bridge_order_parameter_paths.csv"):
            args.paths_csv = Path("outputs/metrics/seed_bridge_order_parameter_paths_angular_randomized.csv")
        if args.activation_csv == Path("outputs/metrics/seed_bridge_activation_summary.csv"):
            args.activation_csv = Path("outputs/metrics/seed_bridge_activation_summary_angular_randomized.csv")
        if args.ordering_csv == Path("outputs/metrics/seed_bridge_feature_ordering.csv"):
            args.ordering_csv = Path("outputs/metrics/seed_bridge_feature_ordering_angular_randomized.csv")
        if args.out_report == Path("outputs/reports/seed_bridge_ordering_report.md"):
            args.out_report = Path("outputs/reports/seed_bridge_angular_randomized_ordering_report.md")
    rows = read_csv_rows(args.input_csv)
    feature_specs = selected_feature_specs(args.angular_phase_mode)
    coordinate_warnings: dict[int, str] = {}
    if args.angular_phase_mode in {"refined", "both"}:
        rows, coordinate_warnings = add_refined_phase_values(rows, units, args.structures_dir, args.ensemble_dir)
    features = [feature for feature, _, _ in feature_specs]
    coordinate_backed_only = args.coordinate_backed_only
    loose_by_unit = samples_by_unit(
        rows,
        units,
        args.loose_ensemble_type,
        feature_specs,
        coordinate_backed_only,
        args.ensemble_dir,
    )
    formed_by_unit = samples_by_unit(
        rows,
        units,
        "formed_perturbed",
        feature_specs,
        coordinate_backed_only,
        args.ensemble_dir,
    )
    pair_rows: list[dict[str, str]] = []
    path_rows: list[dict[str, str]] = []
    activation_rows: list[dict[str, str]] = []
    ordering_rows: list[dict[str, str]] = []
    orientation_warnings: dict[int, list[str]] = {}
    for unit in units:
        loose = loose_by_unit.get(unit, [])
        formed = formed_by_unit.get(unit, [])
        if not loose or not formed:
            raise ValueError(f"Missing {args.loose_ensemble_type} or formed_perturbed samples for unit count {unit}")
        loose, formed, unit_orientation_warnings = empirically_orient_unit_samples(loose, formed, features)
        orientation_warnings[unit] = unit_orientation_warnings
        max_pairs = args.max_pairs if args.max_pairs > 0 else None
        pairs, _, note = endpoint_matching(
            loose,
            formed,
            features,
            epsilon=args.sinkhorn_epsilon,
            iterations=args.sinkhorn_iterations,
            max_pairs=max_pairs,
        )
        pair_rows.extend(build_pair_rows(unit, pairs, note))
        path_rows.extend(build_paths(unit, pairs, times, features))
        unit_activation_rows, _, _ = activation_records_for_unit(
            unit,
            loose,
            formed,
            pairs,
            times,
            features,
            fractions,
            unit_orientation_warnings,
        )
        activation_rows.extend(unit_activation_rows)
        ordering_rows.append(ordering_row(unit, unit_activation_rows, primary_fraction))
    write_csv(args.pairs_csv, pair_rows, PAIR_COLUMNS)
    write_csv(args.paths_csv, path_rows, PATH_COLUMNS)
    write_csv(args.activation_csv, activation_rows, ACTIVATION_COLUMNS)
    write_csv(args.ordering_csv, ordering_rows, ORDERING_COLUMNS)
    plot_paths = write_plots(path_rows, activation_rows, ordering_rows, features, units, primary_fraction, args.plot_dir)
    write_report(
        units,
        features,
        pair_rows,
        activation_rows,
        ordering_rows,
        plot_paths,
        fractions,
        primary_fraction,
        orientation_warnings,
        feature_specs,
        coordinate_warnings,
        args.angular_phase_mode,
        args.loose_ensemble_type,
        args.out_report,
    )
    return pair_rows, path_rows, activation_rows, ordering_rows


def main() -> None:
    try:
        pair_rows, path_rows, activation_rows, ordering_rows = run(parse_args())
    except Exception as exc:
        raise SystemExit(f"ERROR: {exc}") from exc
    print(
        "Wrote "
        f"{len(pair_rows)} endpoint pairs, {len(path_rows)} path rows, "
        f"{len(activation_rows)} activation rows, and {len(ordering_rows)} ordering rows"
    )


if __name__ == "__main__":
    main()
