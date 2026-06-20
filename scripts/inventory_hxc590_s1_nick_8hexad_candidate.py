#!/usr/bin/env python3
"""Inventory and profile Nick's included HXC590 8-hexad XYZ candidate."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from hexaplex_formation.pdb_utils import PDBAtom  # noqa: E402
from hexaplex_formation.scattering import (  # noqa: E402
    debye_intensity_from_distance_histogram,
    d_from_q,
    make_q_grid,
    pair_distance_histogram_for_debye,
)


DEFAULT_SOURCE = Path(r"C:\Users\hpage3\research\nick-diffraction-reference\working\2026-06-19\Hexaplex_8Hexads.xyz")
DEFAULT_COPIED = Path("inputs/candidates/nick_hexaplex_8hexads.xyz")
DEFAULT_INVENTORY = Path("outputs/metrics/hxc590_s1_nick_8hexad_candidate_inventory.csv")
DEFAULT_REPORT = Path("outputs/reports/hxc590_s1_nick_8hexad_candidate_inventory_report.md")
DEFAULT_PROFILE = Path("outputs/metrics/hxc590_s1_nick_hexaplex_8hexads_profile.csv")

INVENTORY_COLUMNS = [
    "candidate_label",
    "source_path",
    "copied_path",
    "filename",
    "sha256",
    "file_size_bytes",
    "xyz_atom_count_header",
    "xyz_atom_count_observed",
    "element_counts",
    "coordinate_span_x_a",
    "coordinate_span_y_a",
    "coordinate_span_z_a",
    "hydrogens_present",
    "duplicate_coordinate_atom_line_count",
    "profile_path",
    "profile_method",
    "provenance_note",
    "not_16mer_note",
]

PROFILE_COLUMNS = ["q_Ainv", "d_A", "intensity"]


@dataclass(frozen=True)
class XYZAtom:
    element: str
    x: float
    y: float
    z: float
    raw_line: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-xyz", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--copied-xyz", type=Path, default=DEFAULT_COPIED)
    parser.add_argument("--inventory-csv", type=Path, default=DEFAULT_INVENTORY)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--profile-csv", type=Path, default=DEFAULT_PROFILE)
    parser.add_argument("--q-min", type=float, default=0.2)
    parser.add_argument("--q-max", type=float, default=3.2)
    parser.add_argument("--q-step", type=float, default=0.01)
    parser.add_argument("--distance-bin-width", type=float, default=0.05)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_xyz(path: Path) -> tuple[int | None, str, list[XYZAtom]]:
    if not path.exists():
        raise FileNotFoundError(f"XYZ file not found: {path}")
    lines = path.read_text(encoding="utf-8").splitlines()
    header_count: int | None = None
    comment = ""
    start = 0
    if lines:
        try:
            header_count = int(lines[0].strip())
            comment = lines[1] if len(lines) > 1 else ""
            start = 2
        except ValueError:
            start = 0
    atoms: list[XYZAtom] = []
    for line_number, line in enumerate(lines[start:], start=start + 1):
        if not line.strip():
            continue
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(f"Invalid XYZ atom row at line {line_number}: {line}")
        atoms.append(
            XYZAtom(
                element=parts[0],
                x=float(parts[1]),
                y=float(parts[2]),
                z=float(parts[3]),
                raw_line="\t".join([parts[0], f"{float(parts[1]):.5f}", f"{float(parts[2]):.5f}", f"{float(parts[3]):.5f}"]),
            )
        )
    return header_count, comment, atoms


def to_pdb_atoms(atoms: list[XYZAtom]) -> list[PDBAtom]:
    return [
        PDBAtom(
            record_type="HETATM",
            atom_serial=index,
            atom_name=atom.element,
            alt_loc="",
            residue_name="N8H",
            chain_id="",
            residue_number=1,
            insertion_code="",
            x=atom.x,
            y=atom.y,
            z=atom.z,
            occupancy=None,
            temp_factor=None,
            element=atom.element,
        )
        for index, atom in enumerate(atoms, start=1)
    ]


def span(values: list[float]) -> float:
    return max(values) - min(values) if values else 0.0


def format_float(value: float, digits: int = 6) -> str:
    if not math.isfinite(value):
        return ""
    return f"{value:.{digits}f}"


def write_profile(path: Path, atoms: list[XYZAtom], q_min: float, q_max: float, q_step: float, bin_width: float) -> dict[str, object]:
    q_values = make_q_grid(q_min, q_max, q_step)
    histogram = pair_distance_histogram_for_debye(to_pdb_atoms(atoms), bin_width=bin_width, element_weights=None)
    intensities = debye_intensity_from_distance_histogram(histogram, q_values)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=PROFILE_COLUMNS, lineterminator="\n")
        writer.writeheader()
        for q_value, intensity in zip(q_values, intensities):
            writer.writerow(
                {
                    "q_Ainv": f"{q_value:.6f}",
                    "d_A": f"{d_from_q(q_value):.6f}",
                    "intensity": f"{intensity:.6f}",
                }
            )
    metadata = {
        "candidate_label": "nick_hexaplex_8hexads",
        "source_xyz": str(DEFAULT_SOURCE),
        "copied_xyz": str(DEFAULT_COPIED),
        "atom_count": len(atoms),
        "method": "Debye-style histogram profile from Nick-provided XYZ",
        "distance_bin_width": bin_width,
        "q_min": q_min,
        "q_max": q_max,
        "q_step": q_step,
        "caution": "Simplified comparative radial profile for scoring windows; not a full powder refinement or unique phase assignment.",
    }
    path.with_name(f"{path.stem}.metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return metadata


def build_inventory_row(source: Path, copied: Path, profile: Path, atoms: list[XYZAtom], header_count: int | None) -> dict[str, str]:
    element_counts = Counter(atom.element for atom in atoms)
    duplicate_count = len(atoms) - len({atom.raw_line for atom in atoms})
    return {
        "candidate_label": "nick_hexaplex_8hexads",
        "source_path": str(source),
        "copied_path": str(copied),
        "filename": copied.name,
        "sha256": sha256(copied),
        "file_size_bytes": str(copied.stat().st_size),
        "xyz_atom_count_header": "" if header_count is None else str(header_count),
        "xyz_atom_count_observed": str(len(atoms)),
        "element_counts": ";".join(f"{element}:{count}" for element, count in sorted(element_counts.items())),
        "coordinate_span_x_a": format_float(span([atom.x for atom in atoms])),
        "coordinate_span_y_a": format_float(span([atom.y for atom in atoms])),
        "coordinate_span_z_a": format_float(span([atom.z for atom in atoms])),
        "hydrogens_present": "yes" if element_counts.get("H", 0) else "no",
        "duplicate_coordinate_atom_line_count": str(duplicate_count),
        "profile_path": str(profile),
        "profile_method": "Debye-style histogram profile from copied XYZ using q/d/intensity scorer-compatible columns",
        "provenance_note": "Nick-provided Hexaplex_8Hexads.xyz, copied unchanged into the HXC590 scoring repo for this screen.",
        "not_16mer_note": "This candidate is treated as Nick's included 8-hexad candidate, not a 16-mer.",
    }


def write_inventory(path: Path, row: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=INVENTORY_COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)


def write_report(path: Path, row: dict[str, str], profile_metadata: dict[str, object]) -> None:
    lines = [
        "# HXC590 S1 Nick 8-Hexad Candidate Inventory",
        "",
        "## Purpose",
        "",
        "This note inventories Nick's included `Hexaplex_8Hexads.xyz` candidate for the corrected HXC590 S1 powder compatibility screen.",
        "",
        "The file is treated as Nick's included 8-hexad candidate, not a 16-mer. No pNAB generation, structural model generation, or coordinate modification was performed.",
        "",
        "## Inventory",
        "",
        f"- Candidate label: `{row['candidate_label']}`",
        f"- Source path: `{row['source_path']}`",
        f"- Copied path: `{row['copied_path']}`",
        f"- SHA256: `{row['sha256']}`",
        f"- File size: {row['file_size_bytes']} bytes",
        f"- XYZ atom count: header {row['xyz_atom_count_header']}, observed {row['xyz_atom_count_observed']}",
        f"- Element counts: {row['element_counts']}",
        f"- Coordinate spans: x {row['coordinate_span_x_a']} A, y {row['coordinate_span_y_a']} A, z {row['coordinate_span_z_a']} A",
        f"- Hydrogens present: {row['hydrogens_present']}",
        f"- Duplicate coordinate/atom lines detected: {row['duplicate_coordinate_atom_line_count']}",
        "",
        "## Profile",
        "",
        f"- Profile path: `{row['profile_path']}`",
        f"- Method: {row['profile_method']}",
        f"- q range: {profile_metadata['q_min']} to {profile_metadata['q_max']} A^-1, step {profile_metadata['q_step']}",
        f"- Distance-bin width: {profile_metadata['distance_bin_width']} A",
        "",
        "This profile uses the scorer-compatible `q_Ainv,d_A,intensity` columns for a falsification-style comparison. It is a simplified comparative radial profile, not a calibrated experimental fit and not a unique refined phase assignment.",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(args: argparse.Namespace) -> dict[str, object]:
    source_header, _source_comment, source_atoms = read_xyz(args.source_xyz)
    copied_header, _copied_comment, copied_atoms = read_xyz(args.copied_xyz)
    if source_header != copied_header or len(source_atoms) != len(copied_atoms):
        raise ValueError("Copied Nick 8-hexad XYZ does not match source atom count")
    if sha256(args.source_xyz) != sha256(args.copied_xyz):
        raise ValueError("Copied Nick 8-hexad XYZ checksum differs from source")
    profile_metadata = write_profile(args.profile_csv, copied_atoms, args.q_min, args.q_max, args.q_step, args.distance_bin_width)
    row = build_inventory_row(args.source_xyz, args.copied_xyz, args.profile_csv, copied_atoms, copied_header)
    write_inventory(args.inventory_csv, row)
    write_report(args.report, row, profile_metadata)
    return {
        "candidate_label": row["candidate_label"],
        "atom_count": row["xyz_atom_count_observed"],
        "element_counts": row["element_counts"],
        "duplicate_count": row["duplicate_coordinate_atom_line_count"],
        "profile_path": row["profile_path"],
    }


def main() -> None:
    result = run(parse_args())
    print(f"Inventoried {result['candidate_label']}: {result['atom_count']} atoms")
    print(f"Element counts: {result['element_counts']}")
    print(f"Duplicate coordinate/atom lines: {result['duplicate_count']}")
    print(f"Wrote profile: {result['profile_path']}")


if __name__ == "__main__":
    main()
