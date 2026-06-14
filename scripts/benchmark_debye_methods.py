#!/usr/bin/env python3
"""Sanity benchmark for direct and histogram Debye approximations."""

from __future__ import annotations

import argparse
import math
import sys
import time
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from hexaplex_formation.pdb_utils import heavy_atoms, load_pdb_atoms  # noqa: E402
from hexaplex_formation.scattering import (  # noqa: E402
    debye_intensity,
    debye_intensity_from_distance_histogram,
    make_q_grid,
    pair_distance_histogram_for_debye,
    stratified_sample_atoms,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pdb",
        type=Path,
        default=Path("outputs/intermediates/normalized_structures/full_hexaplex_anti_parallel_30deg_ideal_heavy_deduped.pdb"),
    )
    parser.add_argument("--sample-atoms", type=int, default=400)
    parser.add_argument("--q-min", type=float, default=0.2)
    parser.add_argument("--q-max", type=float, default=2.5)
    parser.add_argument("--q-step", type=float, default=0.02)
    parser.add_argument("--distance-bin-width", type=float, default=0.05)
    args = parser.parse_args()
    if args.sample_atoms <= 0:
        parser.error("--sample-atoms must be greater than zero")
    if args.distance_bin_width <= 0:
        parser.error("--distance-bin-width must be greater than zero")
    return args


def pearson(xs: list[float], ys: list[float]) -> float | None:
    if len(xs) < 3 or len(xs) != len(ys):
        return None
    x_mean = sum(xs) / len(xs)
    y_mean = sum(ys) / len(ys)
    numerator = sum((x - x_mean) * (y - y_mean) for x, y in zip(xs, ys))
    x_var = sum((x - x_mean) ** 2 for x in xs)
    y_var = sum((y - y_mean) ** 2 for y in ys)
    if x_var == 0.0 or y_var == 0.0:
        return None
    return numerator / math.sqrt(x_var * y_var)


def main() -> int:
    args = parse_args()
    atoms = heavy_atoms(load_pdb_atoms(args.pdb))
    sampled_atoms = stratified_sample_atoms(atoms, args.sample_atoms)
    q_values = make_q_grid(args.q_min, args.q_max, args.q_step)

    direct_start = time.perf_counter()
    direct = debye_intensity(sampled_atoms, q_values)
    direct_seconds = time.perf_counter() - direct_start

    histogram_start = time.perf_counter()
    histogram = pair_distance_histogram_for_debye(sampled_atoms, bin_width=args.distance_bin_width)
    histogram_profile = debye_intensity_from_distance_histogram(histogram, q_values)
    histogram_seconds = time.perf_counter() - histogram_start

    differences = [abs(a - b) for a, b in zip(direct, histogram_profile)]
    correlation = pearson(direct, histogram_profile)

    print(f"pdb: {args.pdb}")
    print(f"input_atoms: {len(atoms)}")
    print(f"sampled_atoms: {len(sampled_atoms)}")
    print(f"q_points: {len(q_values)}")
    print(f"direct_seconds: {direct_seconds:.6f}")
    print(f"histogram_seconds: {histogram_seconds:.6f}")
    print(f"distance_bin_width: {args.distance_bin_width}")
    print(f"max_absolute_difference: {max(differences):.6f}")
    print(f"mean_absolute_difference: {sum(differences) / len(differences):.6f}")
    print(f"pearson_correlation: {correlation:.6f}" if correlation is not None else "pearson_correlation: n/a")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
