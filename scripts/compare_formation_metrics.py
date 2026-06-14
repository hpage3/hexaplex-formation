#!/usr/bin/env python3
"""Aggregate formation metrics across normalized structures."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


METRIC_SUFFIXES = [
    "_heavy_deduped_contacts_4p5A",
    "_heavy_deduped_glu_motifs",
    "_heavy_deduped_helical_order_summary",
    "_allatom_deduped_hbond_candidates",
]

NORMALIZED_STRUCTURE_SUFFIXES = [
    "_allatom_deduped",
    "_heavy_deduped",
    "_heavy",
]

MOTIF_TYPES = [
    "GLU-GLU",
    "GLU-any",
    "backbone_O_to_GLU_sidechain_O",
    "GLU_sidechain_O_to_GLU_sidechain_O",
    "backbone_N_to_GLU_sidechain_O",
    "other_GLU_contact",
]

MOTIF_COLUMNS = {
    "GLU-GLU": "motif_GLU_GLU",
    "GLU-any": "motif_GLU_any",
    "backbone_O_to_GLU_sidechain_O": "motif_backbone_O_to_GLU_sidechain_O",
    "GLU_sidechain_O_to_GLU_sidechain_O": "motif_GLU_sidechain_O_to_GLU_sidechain_O",
    "backbone_N_to_GLU_sidechain_O": "motif_backbone_N_to_GLU_sidechain_O",
    "other_GLU_contact": "motif_other_GLU_contact",
}

NORMALIZATION_COLUMNS = [
    "original_atom_count",
    "allatom_deduped_atom_count",
    "heavy_atom_count",
    "heavy_deduped_atom_count",
    "hydrogens_detected",
]

HELICAL_COLUMNS = [
    "residue_count",
    "mean_radius_xy",
    "min_radius_xy",
    "max_radius_xy",
    "z_span",
    "angular_coverage_rad",
]

OUTPUT_COLUMNS = [
    "structure_base",
    *NORMALIZATION_COLUMNS,
    "residue_count",
    "contact_count_4p5A",
    "GLU_contact_count",
    "GLU_GLU_contact_count",
    "contact_count_per_residue",
    *MOTIF_COLUMNS.values(),
    "hbond_candidate_count",
    "hbond_candidates_per_residue",
    "GLU_involved_hbond_candidate_count",
    "mean_radius_xy",
    "min_radius_xy",
    "max_radius_xy",
    "z_span",
    "angular_coverage_rad",
]

REPORT_COLUMNS = [
    "structure_base",
    "residue_count",
    "heavy_deduped_atom_count",
    "contact_count_4p5A",
    "motif_GLU_GLU",
    "motif_GLU_any",
    "hbond_candidate_count",
    "GLU_involved_hbond_candidate_count",
    "mean_radius_xy",
    "z_span",
    "angular_coverage_rad",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--formation-dir", type=Path, default=Path("outputs/metrics/formation"))
    parser.add_argument(
        "--normalization-summary",
        type=Path,
        default=Path("outputs/metrics/structure_normalization_summary.csv"),
    )
    parser.add_argument("--out-csv", type=Path, default=Path("outputs/metrics/formation_comparison.csv"))
    parser.add_argument("--out-md", type=Path, default=Path("outputs/reports/formation_comparison.md"))
    return parser.parse_args()


def structure_base_from_name(path_or_name: str | Path) -> str:
    stem = Path(path_or_name).stem
    for suffix in METRIC_SUFFIXES:
        if stem.endswith(suffix):
            return stem[: -len(suffix)]
    for suffix in NORMALIZED_STRUCTURE_SUFFIXES:
        if stem.endswith(suffix):
            return stem[: -len(suffix)]
    return stem


def _empty_row(structure_base: str) -> dict[str, str]:
    row = {column: "" for column in OUTPUT_COLUMNS}
    row["structure_base"] = structure_base
    row["contact_count_4p5A"] = "0"
    row["GLU_contact_count"] = "0"
    row["GLU_GLU_contact_count"] = "0"
    row["hbond_candidate_count"] = "0"
    row["GLU_involved_hbond_candidate_count"] = "0"
    for column in MOTIF_COLUMNS.values():
        row[column] = "0"
    return row


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _float_ratio(numerator: int, denominator: str) -> str:
    if not denominator:
        return ""
    residue_count = int(denominator)
    if residue_count == 0:
        return ""
    return f"{numerator / residue_count:.6f}"


def read_motif_counts(path: Path) -> dict[str, int]:
    counts = {motif_type: 0 for motif_type in MOTIF_TYPES}
    for row in _rows(path):
        motif_type = row.get("motif_type", "")
        if motif_type in counts:
            counts[motif_type] = int(row.get("count", "0") or 0)
    return counts


def add_normalization_metrics(rows_by_base: dict[str, dict[str, str]], path: Path) -> None:
    if not path.exists():
        return
    for row in _rows(path):
        structure_base = structure_base_from_name(row["structure_file"])
        output_heavy = row.get("output_heavy_deduped_pdb", "")
        if output_heavy:
            structure_base = structure_base_from_name(output_heavy)
        summary = rows_by_base.setdefault(structure_base, _empty_row(structure_base))
        for column in NORMALIZATION_COLUMNS:
            summary[column] = row.get(column, "")


def add_contact_metrics(rows_by_base: dict[str, dict[str, str]], formation_dir: Path) -> None:
    for path in sorted(formation_dir.glob("*_contacts_4p5A.csv")):
        structure_base = structure_base_from_name(path)
        summary = rows_by_base.setdefault(structure_base, _empty_row(structure_base))
        contact_rows = _rows(path)
        glu_count = 0
        glu_glu_count = 0
        for row in contact_rows:
            res_i = row.get("residue_name_i", "").upper()
            res_j = row.get("residue_name_j", "").upper()
            if res_i == "GLU" or res_j == "GLU":
                glu_count += 1
            if res_i == "GLU" and res_j == "GLU":
                glu_glu_count += 1
        summary["contact_count_4p5A"] = str(len(contact_rows))
        summary["GLU_contact_count"] = str(glu_count)
        summary["GLU_GLU_contact_count"] = str(glu_glu_count)
        summary["contact_count_per_residue"] = _float_ratio(len(contact_rows), summary["residue_count"])


def add_motif_metrics(rows_by_base: dict[str, dict[str, str]], formation_dir: Path) -> None:
    for path in sorted(formation_dir.glob("*_glu_motifs.csv")):
        structure_base = structure_base_from_name(path)
        summary = rows_by_base.setdefault(structure_base, _empty_row(structure_base))
        counts = read_motif_counts(path)
        for motif_type, column in MOTIF_COLUMNS.items():
            summary[column] = str(counts[motif_type])


def add_helical_metrics(rows_by_base: dict[str, dict[str, str]], formation_dir: Path) -> None:
    for path in sorted(formation_dir.glob("*_helical_order_summary.csv")):
        structure_base = structure_base_from_name(path)
        summary = rows_by_base.setdefault(structure_base, _empty_row(structure_base))
        rows = _rows(path)
        if not rows:
            continue
        row = rows[0]
        for column in HELICAL_COLUMNS:
            summary[column] = row.get(column, "")
        summary["contact_count_per_residue"] = _float_ratio(
            int(summary["contact_count_4p5A"] or 0),
            summary["residue_count"],
        )
        summary["hbond_candidates_per_residue"] = _float_ratio(
            int(summary["hbond_candidate_count"] or 0),
            summary["residue_count"],
        )


def add_hbond_metrics(rows_by_base: dict[str, dict[str, str]], formation_dir: Path) -> None:
    for path in sorted(formation_dir.glob("*_hbond_candidates.csv")):
        structure_base = structure_base_from_name(path)
        summary = rows_by_base.setdefault(structure_base, _empty_row(structure_base))
        rows = _rows(path)
        glu_count = sum(
            1
            for row in rows
            if "GLU" in row.get("donor_residue", "").upper()
            or "GLU" in row.get("acceptor_residue", "").upper()
        )
        summary["hbond_candidate_count"] = str(len(rows))
        summary["GLU_involved_hbond_candidate_count"] = str(glu_count)
        summary["hbond_candidates_per_residue"] = _float_ratio(len(rows), summary["residue_count"])


def aggregate_metrics(formation_dir: Path, normalization_summary: Path) -> list[dict[str, str]]:
    rows_by_base: dict[str, dict[str, str]] = {}
    add_normalization_metrics(rows_by_base, normalization_summary)
    add_contact_metrics(rows_by_base, formation_dir)
    add_motif_metrics(rows_by_base, formation_dir)
    add_hbond_metrics(rows_by_base, formation_dir)
    add_helical_metrics(rows_by_base, formation_dir)
    return [rows_by_base[key] for key in sorted(rows_by_base)]


def write_csv(rows: list[dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def _as_int(row: dict[str, str], column: str) -> int:
    return int(row.get(column, "") or 0)


def _as_float(row: dict[str, str], column: str) -> float:
    return float(row.get(column, "") or 0.0)


def _markdown_table(rows: list[dict[str, str]], columns: list[str]) -> str:
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(row.get(column, "") for column in columns) + " |")
    return "\n".join(lines)


def report_notes(rows: list[dict[str, str]]) -> list[str]:
    if not rows:
        return ["No formation metric rows were available."]

    notes: list[str] = []
    largest_glu_glu = max(rows, key=lambda row: _as_int(row, "motif_GLU_GLU"))
    notes.append(
        f"The largest GLU-GLU motif count is in `{largest_glu_glu['structure_base']}` "
        f"({_as_int(largest_glu_glu, 'motif_GLU_GLU')}). This is consistent with a GLU-rich "
        "real-space contact pattern, but does not prove a specific stabilizing mechanism."
    )

    alanine_rows = [row for row in rows if row["structure_base"].startswith("alanine")]
    if alanine_rows and all(_as_int(row, "motif_GLU_any") == 0 for row in alanine_rows):
        notes.append("The alanine controls have zero GLU motif counts, as expected for GLU-free controls.")
    elif alanine_rows:
        notes.append("At least one alanine control has a nonzero GLU motif count, which requires follow-up.")
    else:
        notes.append("No alanine controls were present in the comparison table.")

    largest_coverage = max(rows, key=lambda row: _as_float(row, "angular_coverage_rad"))
    notes.append(
        f"The largest angular coverage is in `{largest_coverage['structure_base']}` "
        f"({_as_float(largest_coverage, 'angular_coverage_rad'):.6f} rad), supporting the working hypothesis "
        "only as a coarse z-axis ordering screen."
    )

    scaffold = next((row for row in rows if row["structure_base"] == "hexaplex_scaffold_only_complement"), None)
    if scaffold and alanine_rows:
        max_alanine_glu_any = max(_as_int(row, "motif_GLU_any") for row in alanine_rows)
        max_alanine_glu_glu = max(_as_int(row, "motif_GLU_GLU") for row in alanine_rows)
        scaffold_glu_any = _as_int(scaffold, "motif_GLU_any")
        scaffold_glu_glu = _as_int(scaffold, "motif_GLU_GLU")
        if scaffold_glu_any > max_alanine_glu_any and scaffold_glu_glu > max_alanine_glu_glu:
            notes.append(
                "The Hexaplex scaffold-only complement has higher GLU-any and GLU-GLU counts than the "
                "alanine controls. This supports the working hypothesis that GLU-rich scaffold contacts "
                "distinguish the organized scaffold, but requires follow-up and does not prove formation pathway."
            )
        else:
            notes.append(
                "The Hexaplex scaffold-only complement does not exceed alanine controls for both GLU-any and "
                "GLU-GLU counts in this table, which requires follow-up."
            )

    return notes


def write_markdown_report(rows: list[dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Formation Comparison Report",
        "",
        "## Scientific caution",
        "",
        "The d ~= 4.5 A diffraction feature is a reciprocal-space feature and should not be read as a literal 4.5 A atom contact. Contact maps and pair-distance summaries are real-space summaries. Hydrogen-bond outputs are rough geometric candidates, not definitive hydrogen bonds. Diffraction intensities from hexads-only and scaffold-only structures are comparative controls, not additive decompositions of the full intensity.",
        "",
        "## Comparison table",
        "",
        _markdown_table(rows, REPORT_COLUMNS),
        "",
        "## Automatically generated notes",
        "",
    ]
    lines.extend(f"- {note}" for note in report_notes(rows))
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    rows = aggregate_metrics(args.formation_dir, args.normalization_summary)
    write_csv(rows, args.out_csv)
    write_markdown_report(rows, args.out_md)
    print(f"Wrote {args.out_csv} with {len(rows)} structure row(s)")
    print(f"Wrote {args.out_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
