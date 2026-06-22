#!/usr/bin/env python3
"""Build synthetic layer-equivalent extensions of Asem's twist-series models."""

from __future__ import annotations

import csv
import math
import sys
from collections import defaultdict
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

INPUTS = {
    "29": ROOT / "inputs" / "asem_twist_series_29_30_31" / "raw" / "29" / "29.pdb.txt",
    "30": ROOT / "inputs" / "asem_twist_series_29_30_31" / "raw" / "30" / "30.pdb.txt",
    "31": ROOT / "inputs" / "asem_twist_series_29_30_31" / "raw" / "31" / "31.pdb.txt",
}
TARGETS = [32, 64, 100]
OUTPUT_DIR = ROOT / "inputs" / "asem_twist_series_29_30_31" / "extended_stacks"
XYZ_DIR = OUTPUT_DIR / "xyz"
AUDIT_CSV = ROOT / "outputs" / "metrics" / "asem_twist_extended_stack_build_audit.csv"
REPORT = ROOT / "outputs" / "reports" / "asem_twist_extended_stack_build_report.md"


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT)).replace("\\", "/")
    except ValueError:
        return str(path)


def fmt(value: float | None) -> str:
    return "" if value is None or not math.isfinite(value) else f"{value:.12g}"


def residue_key(atom) -> tuple[str, int | None, str, str]:
    return (atom.chain_id, atom.residue_number, atom.insertion_code, atom.residue_name)


def layer_groups(atoms: list) -> tuple[list[tuple[tuple[str, int | None, str, str], list]], list[list]]:
    residues: dict[tuple[str, int | None, str, str], list] = defaultdict(list)
    for atom in atoms:
        residues[residue_key(atom)].append(atom)
    ordered_keys = sorted(residues, key=lambda key: (key[1] if key[1] is not None else -1, key[3], key[0], key[2]))
    if len(ordered_keys) % 6 != 0:
        raise ValueError(f"Residue count {len(ordered_keys)} is not divisible by six.")
    ordered_residues = [(key, residues[key]) for key in ordered_keys]
    layers = []
    for start in range(0, len(ordered_residues), 6):
        atoms_in_layer = []
        for _key, group in ordered_residues[start : start + 6]:
            atoms_in_layer.extend(group)
        layers.append(atoms_in_layer)
    return ordered_residues, layers


def coords(atoms: list) -> np.ndarray:
    return np.asarray([[atom.x, atom.y, atom.z] for atom in atoms], dtype=float)


def infer_axis_and_transform(layers: list[list]) -> dict[str, object]:
    if len(layers) < 3:
        raise ValueError("Need at least three layers to infer axis/twist.")
    centers = np.asarray([coords(layer).mean(axis=0) for layer in layers], dtype=float)
    center_mean = centers.mean(axis=0)
    _u, _s, vh = np.linalg.svd(centers - center_mean, full_matrices=False)
    axis = vh[0]
    if axis[2] < 0:
        axis = -axis
    axis = axis / np.linalg.norm(axis)
    ref = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(ref, axis)) > 0.9:
        ref = np.array([0.0, 1.0, 0.0])
    basis_x = ref - axis * np.dot(ref, axis)
    basis_x = basis_x / np.linalg.norm(basis_x)
    basis_y = np.cross(axis, basis_x)
    rel = centers - center_mean
    axial = rel @ axis
    if np.mean(np.diff(axial)) < 0:
        axis = -axis
        basis_y = np.cross(axis, basis_x)
        axial = rel @ axis
    perp_x = rel @ basis_x
    perp_y = rel @ basis_y
    angles = np.unwrap(np.arctan2(perp_y, perp_x))
    rises = np.diff(axial)
    twists = np.diff(angles)
    if np.any(rises <= 0):
        # The centroid axis can be slightly noisy; use absolute rise for the synthetic repeat.
        rises = np.abs(rises)
    rise_mean = float(np.mean(rises))
    twist_mean = float(np.mean(twists))
    if not math.isfinite(rise_mean) or rise_mean <= 0:
        raise ValueError("Nonpositive or invalid mean rise.")
    if not math.isfinite(twist_mean) or abs(twist_mean) < math.radians(1):
        raise ValueError("Invalid mean twist.")
    return {
        "axis_point": center_mean,
        "axis": axis,
        "rise_mean": rise_mean,
        "rise_sd": float(np.std(rises, ddof=1)) if len(rises) > 1 else 0.0,
        "twist_mean_rad": twist_mean,
        "twist_sd_rad": float(np.std(twists, ddof=1)) if len(twists) > 1 else 0.0,
        "axis_method": "PCA of six-residue layer centroids; complete layers repeated by mean screw transform",
    }


def rotate_about_axis(points: np.ndarray, axis_point: np.ndarray, axis: np.ndarray, angle: float) -> np.ndarray:
    rel = points - axis_point
    cos_a = math.cos(angle)
    sin_a = math.sin(angle)
    cross = np.cross(axis, rel)
    dot = rel @ axis
    rotated = rel * cos_a + cross * sin_a + np.outer(dot, axis) * (1.0 - cos_a)
    return rotated + axis_point


def build_extended_atoms(atoms: list, layers: list[list], transform: dict[str, object], target_layers: int) -> list:
    source_layers = len(layers)
    axis_point = transform["axis_point"]
    axis = transform["axis"]
    rise = float(transform["rise_mean"])
    twist = float(transform["twist_mean_rad"])
    extended = []
    for target_layer in range(target_layers):
        source_layer = target_layer % source_layers
        cycle = target_layer // source_layers
        angle = cycle * source_layers * twist
        shift = cycle * source_layers * rise * axis
        layer_atoms = layers[source_layer]
        layer_coords = coords(layer_atoms)
        transformed_coords = rotate_about_axis(layer_coords, axis_point, axis, angle) + shift
        residue_order = []
        seen = set()
        for atom in layer_atoms:
            key = residue_key(atom)
            if key not in seen:
                seen.add(key)
                residue_order.append(key)
        residue_index = {key: i for i, key in enumerate(residue_order)}
        for atom, xyz in zip(layer_atoms, transformed_coords):
            new_residue_number = target_layer * 6 + residue_index[residue_key(atom)] + 1
            extended.append(
                replace(
                    atom,
                    chain_id="A",
                    residue_number=new_residue_number,
                    insertion_code="",
                    x=float(xyz[0]),
                    y=float(xyz[1]),
                    z=float(xyz[2]),
                )
            )
    return extended


def write_xyz(path: Path, atoms: list, source_note: str) -> int:
    heavy = heavy_atoms(atoms)
    deduped = dedupe_exact_atoms(heavy)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f"{len(deduped)}\n")
        handle.write(source_note + "\n")
        for atom in deduped:
            handle.write(f"{atom.element.upper():<6}{atom.x:12.6f}{atom.y:12.6f}{atom.z:12.6f}\n")
    return len(heavy) - len(deduped)


def build_all() -> list[dict[str, object]]:
    rows = []
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    XYZ_DIR.mkdir(parents=True, exist_ok=True)
    for twist, path in INPUTS.items():
        atoms = load_pdb_atoms(path)
        try:
            _residues, layers = layer_groups(atoms)
            transform = infer_axis_and_transform(layers)
            stop_reason = ""
        except Exception as exc:  # noqa: BLE001 - recorded as audit stop reason
            layers = []
            transform = {}
            stop_reason = str(exc)
        for target in TARGETS:
            base_row = {
                "source_twist_deg": twist,
                "source_file": display_path(path),
                "source_atom_count": len(atoms),
                "source_heavy_atom_count": len(heavy_atoms(atoms)),
                "source_layer_count": len(layers),
                "target_repeat_count": target,
                "built_repeat_count": "",
                "output_pdb": "",
                "output_xyz": "",
                "output_atom_count": "",
                "output_heavy_atom_count": "",
                "dedup_removed_count": "",
                "axis_method": transform.get("axis_method", ""),
                "rise_mean_A": fmt(transform.get("rise_mean") if transform else None),
                "rise_sd_A": fmt(transform.get("rise_sd") if transform else None),
                "twist_mean_deg": fmt(math.degrees(transform.get("twist_mean_rad")) if transform else None),
                "twist_sd_deg": fmt(math.degrees(transform.get("twist_sd_rad")) if transform else None),
                "build_status": "stopped" if stop_reason else "built",
                "stop_reason": stop_reason,
                "note": "Synthetic layer-equivalent extension; not minimized or independently generated.",
            }
            if stop_reason:
                rows.append(base_row)
                continue
            try:
                extended = build_extended_atoms(atoms, layers, transform, target)
                pdb_path = OUTPUT_DIR / f"asem_{twist}deg_extended_{target}layer_equiv.pdb"
                xyz_path = XYZ_DIR / f"asem_{twist}deg_extended_{target}layer_equiv_heavy_dedup.xyz"
                write_pdb_atoms(extended, pdb_path)
                dedup_removed = write_xyz(
                    xyz_path,
                    extended,
                    f"Heavy atoms only, exact deduplication applied from synthetic {twist} degree {target} layer-equivalent extension.",
                )
                base_row.update(
                    {
                        "built_repeat_count": target,
                        "output_pdb": display_path(pdb_path),
                        "output_xyz": display_path(xyz_path),
                        "output_atom_count": len(extended),
                        "output_heavy_atom_count": len(heavy_atoms(extended)),
                        "dedup_removed_count": dedup_removed,
                    }
                )
            except Exception as exc:  # noqa: BLE001
                base_row["build_status"] = "stopped"
                base_row["stop_reason"] = str(exc)
            rows.append(base_row)
    return rows


def write_report(rows: list[dict[str, object]]) -> None:
    built = [row for row in rows if row["build_status"] == "built"]
    stopped = [row for row in rows if row["build_status"] != "built"]
    lines = [
        "# Asem Twist Extended Stack Build Report",
        "",
        "## Purpose",
        "",
        "Build synthetic layer-equivalent extensions of Asem's 29, 30, and 31 degree twist models for finite-length diffraction sensitivity testing.",
        "",
        "## Build Method",
        "",
        "Each source PDB was grouped into consecutive six-residue layer-equivalent units. A helical axis was estimated from layer centroids by PCA. The mean adjacent-layer rise and twist were estimated from centroid projections around that axis. Complete source layers were then repeated by the model's own mean screw transform. New residue IDs and chain `A` were assigned to avoid duplicate PDB identifiers.",
        "",
        "These are synthetic periodic/stack extensions for diffraction sensitivity only. They are not chemically searched, independently generated, minimized, or refined structures.",
        "",
        "## Build Summary",
        "",
        "| twist | source layers | target | status | atoms | heavy atoms | rise mean A | twist mean deg | output |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in rows:
        lines.append(
            f"| {row['source_twist_deg']} | {row['source_layer_count']} | {row['target_repeat_count']} | {row['build_status']} | {row['output_atom_count']} | {row['output_heavy_atom_count']} | {row['rise_mean_A']} | {row['twist_mean_deg']} | {row['output_pdb']} |"
        )
    if stopped:
        lines.extend(["", "## Stop Reasons", ""])
        for row in stopped:
            lines.append(f"- {row['source_twist_deg']} target {row['target_repeat_count']}: {row['stop_reason']}")
    lines.extend(
        [
            "",
            "## Validation Notes",
            "",
            f"Built variants: {len(built)} of {len(rows)} requested. Layer grouping was accepted only when the residue count was divisible by six. Heavy-atom XYZ files exclude hydrogens and apply exact heavy-atom deduplication.",
            "",
            "## Limitations",
            "",
            "- Synthetic repeated stacks may exaggerate coherent diffraction relative to disordered real fibers.",
            "- No AmberTools/tleap, minimization, pNAB, notebooks, or new candidate search was run.",
            "- Layer-equivalent count is used instead of claiming independently generated hexads.",
        ]
    )
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    rows = build_all()
    fields = [
        "source_twist_deg",
        "source_file",
        "source_atom_count",
        "source_heavy_atom_count",
        "source_layer_count",
        "target_repeat_count",
        "built_repeat_count",
        "output_pdb",
        "output_xyz",
        "output_atom_count",
        "output_heavy_atom_count",
        "dedup_removed_count",
        "axis_method",
        "rise_mean_A",
        "rise_sd_A",
        "twist_mean_deg",
        "twist_sd_deg",
        "build_status",
        "stop_reason",
        "note",
    ]
    AUDIT_CSV.parent.mkdir(parents=True, exist_ok=True)
    write_csv(AUDIT_CSV, rows, fields)
    write_report(rows)
    print(f"Wrote {AUDIT_CSV}")
    print(f"Wrote {REPORT}")
    return 0


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    raise SystemExit(main())
