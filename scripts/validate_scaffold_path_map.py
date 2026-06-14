#!/usr/bin/env python3
"""Validate a scaffold path map against scaffold PDB residue coverage."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter, OrderedDict
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from hexaplex_formation.geometry import group_atoms_by_residue, residue_centroid  # noqa: E402
from hexaplex_formation.pdb_utils import PDBAtom, load_pdb_atoms  # noqa: E402

try:  # noqa: SIM105
    from helical_order_fitted_axis import fitted_helical_rows, summary_row
except ImportError:  # pragma: no cover - exercised only if script layout changes
    fitted_helical_rows = None
    summary_row = None


REQUIRED_COLUMNS = [
    "map_name",
    "strand_id",
    "strand_label",
    "residue_index_in_pdb_order",
    "chain_id",
    "residue_name",
    "residue_number",
    "insertion_code",
    "residue_label",
    "source",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pdb",
        type=Path,
        default=Path("outputs/intermediates/normalized_structures/hexaplex_scaffold_only_complement_heavy_deduped.pdb"),
    )
    parser.add_argument("--map", type=Path, default=Path("inputs/metadata/scaffold_path_map_candidate.csv"))
    parser.add_argument("--out-json", type=Path, default=Path("outputs/reports/scaffold_path_map_validation.json"))
    parser.add_argument("--out-md", type=Path, default=Path("outputs/reports/scaffold_path_map_validation.md"))
    return parser.parse_args()


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def residue_label(chain_id: str, residue_name: str, residue_number: str, insertion_code: str) -> str:
    label = f"{residue_name}{residue_number}{insertion_code}"
    return f"{chain_id}:{label}" if chain_id else label


def map_residue_identity(row: dict[str, str]) -> tuple[str, str, str, str]:
    return (
        row.get("chain_id", "").strip(),
        row.get("residue_name", "").strip(),
        row.get("residue_number", "").strip(),
        row.get("insertion_code", "").strip(),
    )


def pdb_residue_identity(key: tuple[str, str, int | None, str]) -> tuple[str, str, str, str]:
    chain_id, residue_name, residue_number, insertion_code = key
    return (
        chain_id,
        residue_name,
        "" if residue_number is None else str(residue_number),
        insertion_code,
    )


def grouped_map_rows(rows: list[dict[str, str]]) -> "OrderedDict[str, list[dict[str, str]]]":
    grouped: "OrderedDict[str, list[dict[str, str]]]" = OrderedDict()
    for row in rows:
        strand_id = row.get("strand_id", "").strip()
        grouped.setdefault(strand_id, []).append(row)
    return grouped


def _atoms_for_map_rows(
    map_rows: list[dict[str, str]],
    atoms_by_identity: dict[tuple[str, str, str, str], list[PDBAtom]],
) -> list[PDBAtom]:
    atoms: list[PDBAtom] = []
    for row in map_rows:
        atoms.extend(atoms_by_identity.get(map_residue_identity(row), []))
    return atoms


def _z_span_for_atoms(atoms: list[PDBAtom]) -> str:
    centroids = []
    for residue_atoms in group_atoms_by_residue(atoms).values():
        center = residue_centroid(residue_atoms)
        if center is not None:
            centroids.append(center)
    if not centroids:
        return ""
    zs = [center[2] for center in centroids]
    return f"{max(zs) - min(zs):.6f}"


def _fitted_summary_for_atoms(atoms: list[PDBAtom]) -> dict[str, str]:
    if fitted_helical_rows is None or summary_row is None or not atoms:
        return {"fitted_angular_coverage_rad": "", "mean_radius_fitted": ""}
    rows, fit = fitted_helical_rows(atoms, use_heavy_only=False)
    summary = summary_row(Path("strand_subset.pdb"), rows, fit)
    return {
        "fitted_angular_coverage_rad": summary.get("angular_coverage_rad", ""),
        "mean_radius_fitted": summary.get("mean_radius_fitted", ""),
    }


def validate_map(pdb_path: Path, map_rows: list[dict[str, str]]) -> dict[str, object]:
    atoms = load_pdb_atoms(pdb_path)
    pdb_residue_items = list(group_atoms_by_residue(atoms).items())
    pdb_identities = [pdb_residue_identity(key) for key, _residue_atoms in pdb_residue_items]
    atoms_by_identity = {pdb_residue_identity(key): residue_atoms for key, residue_atoms in pdb_residue_items}
    map_identities = [map_residue_identity(row) for row in map_rows]

    errors: list[str] = []
    warnings: list[str] = []
    missing_columns = [column for column in REQUIRED_COLUMNS if map_rows and column not in map_rows[0]]
    if missing_columns:
        errors.append(f"Map is missing required columns: {', '.join(missing_columns)}")

    pdb_counts = Counter(pdb_identities)
    map_counts = Counter(map_identities)
    missing_from_map = [identity for identity in pdb_identities if map_counts[identity] == 0]
    mapped_missing_from_pdb = [identity for identity in map_identities if pdb_counts[identity] == 0]
    duplicates = [identity for identity, count in map_counts.items() if count > 1]
    blank_strand_ids = [row.get("residue_label", "") for row in map_rows if not row.get("strand_id", "").strip()]
    blank_strand_labels = [row.get("residue_label", "") for row in map_rows if not row.get("strand_label", "").strip()]

    if missing_from_map:
        errors.append(f"{len(missing_from_map)} PDB residue(s) are missing from the map.")
    if mapped_missing_from_pdb:
        errors.append(f"{len(mapped_missing_from_pdb)} mapped residue(s) are absent from the PDB.")
    if duplicates:
        errors.append(f"{len(duplicates)} duplicate mapped residue identity/identities were found.")
    if blank_strand_ids:
        errors.append(f"{len(blank_strand_ids)} map row(s) have blank strand_id.")
    if blank_strand_labels:
        errors.append(f"{len(blank_strand_labels)} map row(s) have blank strand_label.")

    strand_rows: list[dict[str, str]] = []
    for strand_id, rows in grouped_map_rows(map_rows).items():
        if not strand_id:
            continue
        strand_label = rows[0].get("strand_label", "")
        strand_atoms = _atoms_for_map_rows(rows, atoms_by_identity)
        fitted = _fitted_summary_for_atoms(strand_atoms)
        residue_names = [row.get("residue_name", "") for row in rows]
        strand_rows.append(
            {
                "strand_id": strand_id,
                "strand_label": strand_label,
                "residue_count": str(len(rows)),
                "first_residue_label": rows[0].get("residue_label", ""),
                "last_residue_label": rows[-1].get("residue_label", ""),
                "residue_name_pattern": "-".join(residue_names[:12]) + ("..." if len(residue_names) > 12 else ""),
                "z_span": _z_span_for_atoms(strand_atoms),
                "fitted_angular_coverage_rad": fitted["fitted_angular_coverage_rad"],
                "mean_radius_fitted": fitted["mean_radius_fitted"],
            }
        )

    if not strand_rows:
        warnings.append("No nonblank strand rows were available for per-strand summaries.")

    map_names = sorted({row.get("map_name", "") for row in map_rows if row.get("map_name", "")})
    sources = sorted({row.get("source", "") for row in map_rows if row.get("source", "")})
    return {
        "map_name": ", ".join(map_names),
        "source": ", ".join(sources),
        "validation_status": "pass" if not errors else "fail",
        "pdb_residue_count": len(pdb_identities),
        "mapped_residue_count": len(map_rows),
        "missing_from_map_count": len(missing_from_map),
        "mapped_missing_from_pdb_count": len(mapped_missing_from_pdb),
        "duplicate_residue_identity_count": len(duplicates),
        "blank_strand_id_count": len(blank_strand_ids),
        "blank_strand_label_count": len(blank_strand_labels),
        "strand_summaries": strand_rows,
        "errors": errors,
        "warnings": warnings,
        "caution": "Map validation checks internal consistency against the selected PDB; it does not establish biological truth or temporal assembly order.",
    }


def markdown_table(rows: list[dict[str, str]], columns: list[str]) -> str:
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(row.get(column, "") for column in columns) + " |")
    return "\n".join(lines)


def write_markdown_report(validation: dict[str, object], path: Path) -> None:
    strand_rows = validation["strand_summaries"]
    assert isinstance(strand_rows, list)
    columns = [
        "strand_id",
        "strand_label",
        "residue_count",
        "first_residue_label",
        "last_residue_label",
        "z_span",
        "fitted_angular_coverage_rad",
        "mean_radius_fitted",
    ]
    errors = validation["errors"]
    warnings = validation["warnings"]
    assert isinstance(errors, list)
    assert isinstance(warnings, list)
    lines = [
        "# Scaffold Path Map Validation",
        "",
        f"- Map name: {validation.get('map_name', '')}",
        f"- Source: {validation.get('source', '')}",
        f"- Validation status: {validation.get('validation_status', '')}",
        f"- PDB residues: {validation.get('pdb_residue_count', '')}",
        f"- Mapped residues: {validation.get('mapped_residue_count', '')}",
        "",
        "## Scientific caution",
        "",
        "Map validation checks internal consistency against the selected PDB. It does not establish biological truth, validate the PyMOL colored path interpretation by itself, or prove temporal assembly order.",
        "",
        "## Residue coverage summary",
        "",
        f"- Missing PDB residues in map: {validation.get('missing_from_map_count', '')}",
        f"- Mapped residues absent from PDB: {validation.get('mapped_missing_from_pdb_count', '')}",
        f"- Duplicate mapped residue identities: {validation.get('duplicate_residue_identity_count', '')}",
        f"- Blank strand_id rows: {validation.get('blank_strand_id_count', '')}",
        f"- Blank strand_label rows: {validation.get('blank_strand_label_count', '')}",
        "",
        "## Per-strand summary",
        "",
        markdown_table(strand_rows, columns),
        "",
        "## Warnings and errors",
        "",
    ]
    if errors:
        lines.extend(f"- ERROR: {error}" for error in errors)
    if warnings:
        lines.extend(f"- WARNING: {warning}" for warning in warnings)
    if not errors and not warnings:
        lines.append("- No warnings or errors.")
    lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    validation = validate_map(args.pdb, read_csv_rows(args.map))
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(validation, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_markdown_report(validation, args.out_md)
    print(f"Wrote {args.out_json}")
    print(f"Wrote {args.out_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
