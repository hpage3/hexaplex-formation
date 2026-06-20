#!/usr/bin/env python3
"""Inventory the ideal antiparallel 30-degree HXC590 hexaplex model."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import sys
from collections import Counter
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from hexaplex_formation.pdb_utils import (  # noqa: E402
    atom_identity_key,
    bounding_box,
    dedupe_exact_atoms,
    load_pdb_atoms,
)


DEFAULT_SOURCE = Path(r"C:\Users\Public\hexaplex-tools\repos\fiber-diffraction-original\inputs\Hexaplex_AntiParallel_30deg_Ideal.pdb")
DEFAULT_COPIED = Path("inputs/candidates/ideal_antiparallel_30deg_hexaplex.pdb")
DEFAULT_EXISTING_PROFILE = Path("outputs/mini_hexaplex/radial_profiles/full_length_baseline_radial.csv")
DEFAULT_PROFILE = Path("outputs/metrics/hxc590_s1_ideal_antiparallel_30deg_hexaplex_profile.csv")
DEFAULT_INVENTORY = Path("outputs/metrics/hxc590_s1_ideal_antiparallel_30deg_hexaplex_inventory.csv")
DEFAULT_REPORT = Path("outputs/reports/hxc590_s1_ideal_antiparallel_30deg_hexaplex_inventory_report.md")
DEFAULT_EXISTING_DEDUPED = Path("outputs/intermediates/ai_candidate_inputs/full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain.pdb")

INVENTORY_COLUMNS = [
    "candidate_label",
    "source_path",
    "copied_path",
    "filename",
    "sha256",
    "file_size_bytes",
    "total_pdb_atom_count",
    "element_counts",
    "residue_counts",
    "chain_ids",
    "chain_count",
    "coordinate_span_x_a",
    "coordinate_span_y_a",
    "coordinate_span_z_a",
    "hydrogens_present",
    "duplicate_atom_coordinate_line_count",
    "duplicate_atom_coordinate_lines_detected",
    "deduped_atom_count",
    "deduped_reference_path",
    "deduped_reference_sha256",
    "profile_path",
    "profile_source_path",
    "profile_method",
    "provenance_note",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-pdb", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--copied-pdb", type=Path, default=DEFAULT_COPIED)
    parser.add_argument("--existing-profile", type=Path, default=DEFAULT_EXISTING_PROFILE)
    parser.add_argument("--profile-csv", type=Path, default=DEFAULT_PROFILE)
    parser.add_argument("--inventory-csv", type=Path, default=DEFAULT_INVENTORY)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--deduped-reference", type=Path, default=DEFAULT_EXISTING_DEDUPED)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def format_float(value: float) -> str:
    return f"{value:.6f}"


def duplicate_count(atoms) -> int:
    seen = set()
    duplicates = 0
    for atom in atoms:
        key = atom_identity_key(atom)
        if key in seen:
            duplicates += 1
        else:
            seen.add(key)
    return duplicates


def counts_text(counter: Counter[str]) -> str:
    return ";".join(f"{key}:{value}" for key, value in sorted(counter.items()))


def build_inventory_row(args: argparse.Namespace) -> dict[str, str]:
    if not args.source_pdb.exists():
        raise FileNotFoundError(f"Ideal antiparallel 30-degree source PDB not found: {args.source_pdb}")
    if not args.copied_pdb.exists():
        raise FileNotFoundError(f"Copied candidate PDB not found: {args.copied_pdb}")
    if sha256(args.source_pdb) != sha256(args.copied_pdb):
        raise ValueError("Copied candidate PDB checksum differs from source")
    if not args.existing_profile.exists():
        raise FileNotFoundError(f"Existing full-length baseline profile not found: {args.existing_profile}")

    atoms = load_pdb_atoms(args.copied_pdb)
    element_counts = Counter((atom.element or "").upper() for atom in atoms)
    residue_counts = Counter(atom.residue_name for atom in atoms)
    chains = sorted({atom.chain_id or "(blank)" for atom in atoms})
    bbox = bounding_box(atoms)
    if bbox is None:
        raise ValueError("No atoms found in copied candidate PDB")
    min_x, max_x, min_y, max_y, min_z, max_z = bbox
    duplicates = duplicate_count(atoms)
    deduped_atoms = dedupe_exact_atoms(atoms)

    args.profile_csv.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(args.existing_profile, args.profile_csv)
    metadata = {
        "candidate_label": "ideal_antiparallel_30deg_hexaplex",
        "source_pdb": str(args.source_pdb),
        "copied_pdb": str(args.copied_pdb),
        "copied_pdb_sha256": sha256(args.copied_pdb),
        "profile_source_path": str(args.existing_profile),
        "profile_source_note": "Reused existing full-length baseline radial profile for the same original PDB after the repo's established deduped/full-length profile workflow.",
        "deduped_reference": str(args.deduped_reference),
        "caution": "Nick clarified that his 16-mer simulation benchmark shorthand refers to Hexaplex_AntiParallel_30deg_Ideal.pdb. The scoring repo labels this file by ideal antiparallel 30-degree hexaplex model provenance instead.",
    }
    args.profile_csv.with_name(f"{args.profile_csv.stem}.metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    return {
        "candidate_label": "ideal_antiparallel_30deg_hexaplex",
        "source_path": str(args.source_pdb),
        "copied_path": str(args.copied_pdb),
        "filename": args.copied_pdb.name,
        "sha256": sha256(args.copied_pdb),
        "file_size_bytes": str(args.copied_pdb.stat().st_size),
        "total_pdb_atom_count": str(len(atoms)),
        "element_counts": counts_text(element_counts),
        "residue_counts": counts_text(residue_counts),
        "chain_ids": ";".join(chains),
        "chain_count": str(len(chains)),
        "coordinate_span_x_a": format_float(max_x - min_x),
        "coordinate_span_y_a": format_float(max_y - min_y),
        "coordinate_span_z_a": format_float(max_z - min_z),
        "hydrogens_present": "yes" if element_counts.get("H", 0) else "no",
        "duplicate_atom_coordinate_line_count": str(duplicates),
        "duplicate_atom_coordinate_lines_detected": "yes" if duplicates else "no",
        "deduped_atom_count": str(len(deduped_atoms)),
        "deduped_reference_path": str(args.deduped_reference),
        "deduped_reference_sha256": sha256(args.deduped_reference) if args.deduped_reference.exists() else "",
        "profile_path": str(args.profile_csv),
        "profile_source_path": str(args.existing_profile),
        "profile_method": "reused existing full-length baseline radial profile generated by the repo's established workflow",
        "provenance_note": "Hexaplex_AntiParallel_30deg_Ideal.pdb is treated as the ideal antiparallel 30-degree hexaplex model. Nick clarified that this is the file he meant by the 16-mer simulation benchmark shorthand.",
    }


def write_inventory(path: Path, row: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=INVENTORY_COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)


def write_report(path: Path, row: dict[str, str]) -> None:
    lines = [
        "# HXC590 S1 Ideal Antiparallel 30-Degree Hexaplex Inventory",
        "",
        "## Purpose",
        "",
        "This note inventories `Hexaplex_AntiParallel_30deg_Ideal.pdb` as the ideal antiparallel 30-degree hexaplex model for the corrected HXC590 S1 powder compatibility screen.",
        "",
        "Nick clarified that this is the file he meant when referring to the 16-mer simulation benchmark. For scoring and falsification outputs, the repo labels it by model provenance and geometry rather than by that informal shorthand. The file was copied unchanged into the scoring repo; no pNAB generation, structural model generation, or coordinate modification was performed.",
        "",
        "## Inventory",
        "",
        f"- Candidate label: `{row['candidate_label']}`",
        f"- Source path: `{row['source_path']}`",
        f"- Copied path: `{row['copied_path']}`",
        f"- SHA256: `{row['sha256']}`",
        f"- File size: {row['file_size_bytes']} bytes",
        f"- Total PDB atom count: {row['total_pdb_atom_count']}",
        f"- Element counts: {row['element_counts']}",
        f"- Residue counts: {row['residue_counts']}",
        f"- Chain IDs: {row['chain_ids']} ({row['chain_count']} chain entries)",
        f"- Coordinate spans: x {row['coordinate_span_x_a']} A, y {row['coordinate_span_y_a']} A, z {row['coordinate_span_z_a']} A",
        f"- Hydrogens present: {row['hydrogens_present']}",
        f"- Duplicate atom/coordinate records detected: {row['duplicate_atom_coordinate_lines_detected']} ({row['duplicate_atom_coordinate_line_count']})",
        f"- Deduped atom count by exact atom identity: {row['deduped_atom_count']}",
        "",
        "## Profile",
        "",
        f"- Scorer-compatible profile path: `{row['profile_path']}`",
        f"- Existing profile source: `{row['profile_source_path']}`",
        f"- Existing deduped reference: `{row['deduped_reference_path']}`",
        "",
        "The profile is reused from the repo's existing full-length baseline radial profile for this same original PDB family after the established deduped/full-length profile workflow. This is a falsification-style compatibility input, not a unique refined phase assignment.",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(args: argparse.Namespace) -> dict[str, str]:
    row = build_inventory_row(args)
    write_inventory(args.inventory_csv, row)
    write_report(args.report, row)
    return row


def main() -> None:
    row = run(parse_args())
    print(f"Inventoried {row['candidate_label']}: {row['total_pdb_atom_count']} PDB atoms")
    print(f"Element counts: {row['element_counts']}")
    print(f"Residue counts: {row['residue_counts']}")
    print(f"Duplicate atom/coordinate records: {row['duplicate_atom_coordinate_line_count']}")
    print(f"Wrote profile: {row['profile_path']}")


if __name__ == "__main__":
    main()
