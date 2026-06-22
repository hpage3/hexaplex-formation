#!/usr/bin/env python3
"""Generate focused 30-degree rise variants for geometric sensitivity testing."""

from __future__ import annotations

import argparse
import csv
import statistics
import sys
from collections import Counter
from dataclasses import replace
from pathlib import Path


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
DEFAULT_OUTPUT_DIR = ROOT / "outputs" / "rise_twist_variants"
DEFAULT_RISES = [3.40, 3.38]
BASE_RESIDUES = {"CYP", "MEP"}


def rise_token(value: float) -> str:
    return f"{value:.2f}".replace(".", "p")


def model_id_for(twist_deg: float, rise_A: float) -> str:
    twist_token = f"{twist_deg:.0f}" if float(twist_deg).is_integer() else str(twist_deg).replace(".", "p")
    return f"twist{twist_token}_rise{rise_token(rise_A)}"


def parse_rise_values(values: list[str]) -> list[float]:
    rises = [float(value) for value in values]
    if any(value <= 0 for value in rises):
        raise ValueError("Rise values must be positive")
    return rises


def residue_key(atom) -> tuple[str, int | None, str, str]:
    return (atom.chain_id, atom.residue_number, atom.insertion_code, atom.residue_name)


def residue_centroids(atoms) -> dict[tuple[str, int | None, str, str], tuple[str, float]]:
    grouped = {}
    for atom in atoms:
        grouped.setdefault(residue_key(atom), []).append(atom)
    return {
        key: (members[0].residue_name, statistics.fmean(atom.z for atom in members))
        for key, members in grouped.items()
    }


def infer_base_layers(atoms) -> list[dict[str, object]]:
    base_rows = [
        (z, key, residue_name)
        for key, (residue_name, z) in residue_centroids(atoms).items()
        if residue_name in BASE_RESIDUES
    ]
    base_rows.sort(key=lambda item: item[0])
    layers = []
    for row in base_rows:
        if not layers or row[0] - layers[-1][-1][0] > 0.6:
            layers.append([row])
        else:
            layers[-1].append(row)

    if len(layers) != 15:
        raise ValueError(f"Expected 15 base layers; inferred {len(layers)}")
    for index, layer in enumerate(layers):
        counts = Counter(row[2] for row in layer)
        if len(layer) != 6 or counts != {"CYP": 3, "MEP": 3}:
            raise ValueError(f"Layer {index} has {len(layer)} residues and mix {dict(counts)}")

    return [
        {
            "layer_index": index,
            "original_center_z_A": statistics.fmean(row[0] for row in layer),
            "base_residue_numbers": ";".join(str(row[1][1]) for row in layer),
        }
        for index, layer in enumerate(layers)
    ]


def nearest_layer_index(z_value: float, layer_centers: list[float]) -> int:
    return min(range(len(layer_centers)), key=lambda index: abs(layer_centers[index] - z_value))


def assign_residues_to_layers(atoms, layers: list[dict[str, object]]) -> dict[tuple[str, int | None, str, str], int]:
    centers = [float(layer["original_center_z_A"]) for layer in layers]
    return {
        key: nearest_layer_index(z_value, centers)
        for key, (_residue_name, z_value) in residue_centroids(atoms).items()
    }


def transform_atoms(atoms, layers: list[dict[str, object]], assignments, target_rise_A: float):
    original_centers = [float(layer["original_center_z_A"]) for layer in layers]
    original_steps = [
        original_centers[index + 1] - original_centers[index]
        for index in range(len(original_centers) - 1)
    ]
    source_mean_rise = statistics.fmean(original_steps)
    anchor_index = len(original_centers) // 2
    anchor_z = original_centers[anchor_index]
    target_centers = [
        anchor_z + (index - anchor_index) * target_rise_A
        for index in range(len(original_centers))
    ]
    deltas = [
        target_centers[index] - original_centers[index]
        for index in range(len(original_centers))
    ]
    transformed = [
        replace(atom, z=atom.z + deltas[assignments[residue_key(atom)]])
        for atom in atoms
    ]
    return transformed, original_centers, target_centers, source_mean_rise, anchor_index


def write_xyz(pdb_atoms, xyz_path: Path) -> dict[str, object]:
    heavy = heavy_atoms(pdb_atoms)
    deduped = dedupe_exact_atoms(heavy)
    xyz_path.parent.mkdir(parents=True, exist_ok=True)
    with xyz_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f"{len(deduped)}\n")
        handle.write("Heavy atoms only, exact deduplication applied from rise-variant PDB.\n")
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
    text = "\n".join(
        [
            "# Rise Twist Variants",
            "",
            f"Source PDB: `{source.relative_to(ROOT)}`",
            "",
            "These PDBs are rigid-layer z-position variants of the ideal antiparallel 30-degree model.",
            "They preserve atom count, atom identity, residue identity, x/y coordinates, and intralayer geometry.",
            "Only inferred layer z positions are changed. These are geometric sensitivity-test models, not relaxed physical models.",
            "",
        ]
    )
    (output_dir / "README.md").write_text(text, encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-pdb", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--twist-deg", type=float, default=30.0)
    parser.add_argument("--rises", nargs="+", default=[str(value) for value in DEFAULT_RISES])
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rises = parse_rise_values(args.rises)
    atoms = load_pdb_atoms(args.source_pdb)
    layers = infer_base_layers(atoms)
    assignments = assign_residues_to_layers(atoms, layers)
    output_dir = args.output_dir
    structure_dir = output_dir / "structures"
    xyz_dir = output_dir / "xyz"
    structure_dir.mkdir(parents=True, exist_ok=True)
    xyz_dir.mkdir(parents=True, exist_ok=True)
    write_readme(output_dir, args.source_pdb)

    validation_rows = []
    source_atom_count = len(atoms)
    for rise_A in rises:
        model_id = model_id_for(args.twist_deg, rise_A)
        transformed, original_centers, target_centers, source_mean_rise, anchor_index = transform_atoms(
            atoms,
            layers,
            assignments,
            rise_A,
        )
        pdb_path = structure_dir / f"{model_id}.pdb"
        xyz_path = xyz_dir / f"{model_id}.xyz"
        write_pdb_atoms(transformed, pdb_path)
        xyz_info = write_xyz(transformed, xyz_path)
        target_steps = [
            target_centers[index + 1] - target_centers[index]
            for index in range(len(target_centers) - 1)
        ]
        validation_rows.append(
            {
                "model_id": model_id,
                "twist_deg": f"{args.twist_deg:.6g}",
                "rise_A": f"{rise_A:.2f}",
                "pdb_path": str(pdb_path.relative_to(ROOT)),
                "xyz_path": str(xyz_path.relative_to(ROOT)),
                "source_atom_count": source_atom_count,
                "variant_atom_count": len(transformed),
                "atom_count_preserved": len(transformed) == source_atom_count,
                "inferred_base_layer_count": len(layers),
                "central_anchor_layer_index": anchor_index,
                "source_mean_base_layer_rise_A": f"{source_mean_rise:.8g}",
                "measured_variant_base_layer_rise_A": f"{statistics.fmean(target_steps):.8g}",
                "max_adjacent_rise_error_A": f"{max(abs(step - rise_A) for step in target_steps):.8g}",
                "max_abs_layer_z_delta_A": f"{max(abs(target_centers[index] - original_centers[index]) for index in range(len(layers))):.8g}",
                **xyz_info,
            }
        )

    write_csv(output_dir / "rise_variant_validation.csv", validation_rows)
    print(f"Wrote {len(validation_rows)} rise variants to {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
