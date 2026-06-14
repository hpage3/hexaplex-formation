#!/usr/bin/env python3
"""Compute a simplified Debye-style isotropic radial scattering profile."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from hexaplex_formation.pdb_utils import heavy_atoms, load_pdb_atoms  # noqa: E402
from hexaplex_formation.scattering import (  # noqa: E402
    debye_intensity,
    debye_intensity_from_distance_histogram,
    d_from_q,
    make_q_grid,
    pair_distance_histogram_for_debye,
    stratified_sample_atoms,
)


FIELDNAMES = ["q_Ainv", "d_A", "intensity"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--q-min", type=float, default=0.2)
    parser.add_argument("--q-max", type=float, default=2.5)
    parser.add_argument("--q-step", type=float, default=0.01)
    parser.add_argument("--method", choices=["direct", "histogram"], default="histogram")
    parser.add_argument("--distance-bin-width", type=float, default=0.05)
    parser.add_argument("--heavy-only", action="store_true", default=True, help="Exclude hydrogens. This is the default.")
    parser.add_argument("--all-atoms", action="store_true", help="Include hydrogens.")
    parser.add_argument("--sample-atoms", type=int, default=None)
    parser.add_argument("--max-atoms", type=int, default=None)
    parser.add_argument("--unit-weights", action="store_true", default=True)
    args = parser.parse_args()
    if args.distance_bin_width <= 0:
        parser.error("--distance-bin-width must be greater than zero")
    if args.sample_atoms is not None and args.sample_atoms <= 0:
        parser.error("--sample-atoms must be greater than zero")
    if args.max_atoms is not None and args.max_atoms <= 0:
        parser.error("--max-atoms must be greater than zero")
    return args


def _metadata_path(out_path: Path) -> Path:
    return out_path.with_name(f"{out_path.stem}.metadata.json")


def write_profile(
    out_path: Path,
    q_values: list[float],
    intensities: list[float],
    metadata: dict[str, object],
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, lineterminator="\n")
        writer.writeheader()
        for q_value, intensity in zip(q_values, intensities):
            writer.writerow(
                {
                    "q_Ainv": f"{q_value:.6f}",
                    "d_A": f"{d_from_q(q_value):.6f}",
                    "intensity": f"{intensity:.6f}",
                }
            )
    _metadata_path(out_path).write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    atoms = load_pdb_atoms(args.pdb)
    if not args.all_atoms:
        atoms = heavy_atoms(atoms)
    source_atom_count = len(atoms)
    sampling_used = False
    first_n_truncated = False
    if args.sample_atoms is not None and len(atoms) > args.sample_atoms:
        atoms = stratified_sample_atoms(atoms, args.sample_atoms)
        sampling_used = True
    elif args.max_atoms is not None and len(atoms) > args.max_atoms:
        print(
            "WARNING: --max-atoms uses first-N atom truncation and may bias hexads-plus-scaffold structures. "
            "Prefer --sample-atoms for deterministic stratified sampling.",
            file=sys.stderr,
        )
        atoms = atoms[: args.max_atoms]
        first_n_truncated = True

    q_values = make_q_grid(args.q_min, args.q_max, args.q_step)
    print(
        f"Computing {args.method} Debye profile for {args.pdb}: {len(atoms)} atoms, {len(q_values)} q value(s)",
        file=sys.stderr,
    )
    if args.method == "direct":
        intensities = debye_intensity(atoms, q_values, element_weights=None)
    else:
        histogram = pair_distance_histogram_for_debye(atoms, bin_width=args.distance_bin_width, element_weights=None)
        intensities = debye_intensity_from_distance_histogram(histogram, q_values)

    metadata = {
        "pdb": str(args.pdb),
        "atom_count_input": source_atom_count,
        "atom_count_used": len(atoms),
        "method": args.method,
        "distance_bin_width": args.distance_bin_width,
        "sample_atoms": args.sample_atoms,
        "max_atoms": args.max_atoms,
        "sampling_used": sampling_used,
        "first_n_truncated": first_n_truncated,
        "truncated": first_n_truncated,
        "q_min": args.q_min,
        "q_max": args.q_max,
        "q_step": args.q_step,
        "caution": (
            "Simplified Debye-style comparative approximation; not a full fiber-diffraction simulation. "
            "The d ~= 4.5 A feature is reciprocal-space-like, not a literal atom contact. "
            "Component profiles are comparative controls, not additive decompositions."
        ),
    }
    write_profile(args.out, q_values, intensities, metadata)
    print(f"Wrote {args.out}")
    print(f"Wrote {_metadata_path(args.out)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
