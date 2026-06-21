#!/usr/bin/env python3
"""Generate rigid-layer rise variants of Nick's ideal antiparallel model."""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from dataclasses import replace
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from hexaplex_formation.pdb_utils import (  # noqa: E402
    dedupe_exact_atoms,
    heavy_atoms,
    load_pdb_atoms,
    write_pdb_atoms,
)


DEFAULT_SOURCE = ROOT / "inputs" / "nick_ideal_models" / "Hexaplex_AntiParallel_30deg_Ideal.pdb"
DEFAULT_OUTPUT_DIR = ROOT / "inputs" / "nick_ideal_models" / "rise_variants"
RISES_A = [3.36, 3.37, 3.38, 3.39, 3.40]
BASE_RESIDUES = {"CYP", "MEP"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-pdb", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--rises", type=float, nargs="+", default=RISES_A)
    return parser.parse_args()


def token(value: float) -> str:
    return f"{value:.2f}".replace(".", "p")


def residue_key(atom) -> tuple[str, int | None, str, str]:
    return (atom.chain_id, atom.residue_number, atom.insertion_code, atom.residue_name)


def residue_centroids(atoms) -> dict[tuple[str, int | None, str, str], tuple[str, float]]:
    groups: dict[tuple[str, int | None, str, str], list[object]] = {}
    for atom in atoms:
        groups.setdefault(residue_key(atom), []).append(atom)
    return {
        key: (members[0].residue_name, float(np.mean([atom.z for atom in members])))
        for key, members in groups.items()
    }


def infer_hexad_layers(atoms) -> list[dict[str, object]]:
    """Infer six-residue CYP/MEP base planes from z-centroids.

    The PDB has 90 CYP/MEP residues arranged as 15 six-residue base planes.
    Each plane contains three CYP and three MEP residues. This is the rigid
    layer used for the controlled rise perturbation.
    """

    base_rows = [
        (z, key, residue_name)
        for key, (residue_name, z) in residue_centroids(atoms).items()
        if residue_name in BASE_RESIDUES
    ]
    base_rows.sort(key=lambda item: item[0])
    layers: list[list[tuple[float, tuple[str, int | None, str, str], str]]] = []
    for row in base_rows:
        if not layers or row[0] - layers[-1][-1][0] > 0.6:
            layers.append([row])
        else:
            layers[-1].append(row)

    if len(layers) != 15:
        raise ValueError(f"Expected 15 six-residue base planes; inferred {len(layers)}")
    for index, layer in enumerate(layers):
        if len(layer) != 6:
            raise ValueError(f"Layer {index} has {len(layer)} base residues, expected 6")
        counts = Counter(row[2] for row in layer)
        if counts != {"CYP": 3, "MEP": 3}:
            raise ValueError(f"Layer {index} residue mix is {dict(counts)}, expected 3 CYP and 3 MEP")

    return [
        {
            "layer_index": index,
            "original_center_z_A": float(np.mean([row[0] for row in layer])),
            "base_residue_numbers": ";".join(str(row[1][1]) for row in layer),
        }
        for index, layer in enumerate(layers)
    ]


def assign_residues_to_layers(atoms, layers: list[dict[str, object]]) -> dict[tuple[str, int | None, str, str], int]:
    centers = np.asarray([float(layer["original_center_z_A"]) for layer in layers], dtype=float)
    assignments = {}
    for key, (_residue_name, z) in residue_centroids(atoms).items():
        assignments[key] = int(np.argmin(np.abs(centers - z)))
    return assignments


def transform_atoms(atoms, layers: list[dict[str, object]], assignments, target_rise: float):
    original_centers = np.asarray([float(layer["original_center_z_A"]) for layer in layers], dtype=float)
    measured_spacing = float(np.mean(np.diff(original_centers)))
    anchor_index = len(layers) // 2
    anchor_z = original_centers[anchor_index]
    target_centers = anchor_z + (np.arange(len(layers)) - anchor_index) * target_rise
    deltas = target_centers - original_centers
    transformed = []
    for atom in atoms:
        layer_index = assignments[residue_key(atom)]
        transformed.append(replace(atom, z=atom.z + float(deltas[layer_index])))
    return transformed, original_centers, target_centers, measured_spacing, anchor_index


def write_xyz(pdb_atoms, xyz_path: Path) -> dict[str, object]:
    heavy = heavy_atoms(pdb_atoms)
    deduped = dedupe_exact_atoms(heavy)
    xyz_path.parent.mkdir(parents=True, exist_ok=True)
    with xyz_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f"{len(deduped)}\n")
        handle.write("Heavy atoms only, exact deduplication applied from ideal rise variant PDB.\n")
        for atom in deduped:
            handle.write(f"{atom.element.upper():<6}{atom.x:12.6f}{atom.y:12.6f}{atom.z:12.6f}\n")
    return {
        "heavy_atom_count_before_dedup": len(heavy),
        "xyz_atom_count": len(deduped),
        "duplicate_heavy_atoms_removed": len(heavy) - len(deduped),
    }


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_readme(output_dir: Path, source: Path) -> None:
    lines = [
        "# Ideal Antiparallel Rise Variants",
        "",
        f"Generated from `{source.relative_to(ROOT)}`.",
        "",
        "These are rigid-layer z-position variants for Nick's clarified ideal antiparallel 30-degree model.",
        "The generator preserves x/y coordinates, twist geometry, atom ordering, and intralayer geometry.",
        "",
        "Layer convention:",
        "",
        "- The source PDB does not store explicit chain/layer IDs.",
        "- The script infers 15 six-residue CYP/MEP base planes from residue z-centroids.",
        "- Each base plane contains three CYP and three MEP residues.",
        "- Non-base GLU residues are assigned to the nearest inferred base plane and moved rigidly with that plane.",
        "- The central inferred base plane is fixed as the z anchor.",
        "",
        "This is a controlled rise-sensitivity perturbation only. It does not adjust peptide omega, twist angle, side-chain removal logic, or perform minimization.",
    ]
    (output_dir / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    atoms = load_pdb_atoms(args.source_pdb)
    layers = infer_hexad_layers(atoms)
    assignments = assign_residues_to_layers(atoms, layers)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "xyz").mkdir(parents=True, exist_ok=True)
    write_readme(args.output_dir, args.source_pdb)

    validation_rows = []
    source_atom_count = len(atoms)
    for rise in args.rises:
        transformed, original_centers, target_centers, original_rise, anchor_index = transform_atoms(
            atoms,
            layers,
            assignments,
            rise,
        )
        pdb_path = args.output_dir / f"Hexaplex_AntiParallel_30deg_Ideal_rise_{token(rise)}.pdb"
        xyz_path = args.output_dir / "xyz" / f"Hexaplex_AntiParallel_30deg_Ideal_rise_{token(rise)}.xyz"
        write_pdb_atoms(transformed, pdb_path)
        xyz_info = write_xyz(transformed, xyz_path)
        measured = float(np.mean(np.diff(target_centers)))
        max_error = float(np.max(np.abs(np.diff(target_centers) - rise)))
        max_abs_z_delta = float(np.max(np.abs(target_centers - original_centers)))
        validation_rows.append(
            {
                "rise_A": f"{rise:.2f}",
                "pdb_path": str(pdb_path.relative_to(ROOT)),
                "xyz_path": str(xyz_path.relative_to(ROOT)),
                "source_atom_count": source_atom_count,
                "variant_atom_count": len(transformed),
                "atom_count_preserved": len(transformed) == source_atom_count,
                "inferred_base_layer_count": len(layers),
                "base_residues_per_layer": 6,
                "central_anchor_layer_index": anchor_index,
                "source_mean_base_layer_rise_A": original_rise,
                "measured_variant_base_layer_rise_A": measured,
                "max_adjacent_rise_error_A": max_error,
                "max_abs_layer_z_delta_A": max_abs_z_delta,
                **xyz_info,
            }
        )
    write_csv(args.output_dir / "rise_variant_validation.csv", validation_rows)
    print(f"Wrote {len(args.rises)} rise variants to {args.output_dir}")
    print(f"Wrote validation table: {args.output_dir / 'rise_variant_validation.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
