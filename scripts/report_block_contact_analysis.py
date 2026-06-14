#!/usr/bin/env python3
"""Generate a Markdown report for block contact decomposition outputs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


CATEGORIES = [
    "scaffold_within_block",
    "scaffold_between_blocks",
    "scaffold_hexad_or_other",
    "hexad_or_other_internal",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary-dir", type=Path, default=Path("outputs/metrics/block_contacts"))
    parser.add_argument("--out-md", type=Path, default=Path("outputs/reports/block_contact_analysis.md"))
    return parser.parse_args()


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def model_name_from_summary(path: Path) -> str:
    suffix = "_block_contact_summary"
    stem = path.stem
    return stem[: -len(suffix)] if stem.endswith(suffix) else stem


def load_category_summaries(summary_dir: Path) -> dict[str, dict[str, dict[str, str]]]:
    summaries: dict[str, dict[str, dict[str, str]]] = {}
    for path in sorted(summary_dir.glob("*_block_contact_summary.csv")):
        model_name = model_name_from_summary(path)
        summaries[model_name] = {row["contact_category"]: row for row in _rows(path)}
    return summaries


def _int(summary: dict[str, dict[str, str]], category: str, column: str) -> int:
    return int(summary.get(category, {}).get(column, "0") or 0)


def comparison_table_rows(summaries: dict[str, dict[str, dict[str, str]]]) -> list[dict[str, str]]:
    model_names = [
        "hexaplex_scaffold_only_complement_heavy_deduped",
        "full_hexaplex_anti_parallel_30deg_ideal_heavy_deduped",
    ]
    rows: list[dict[str, str]] = []
    for model_name in model_names:
        if model_name not in summaries:
            continue
        summary = summaries[model_name]
        row = {"model_name": model_name}
        for category in CATEGORIES:
            row[f"{category}_contacts"] = str(_int(summary, category, "contact_count"))
            row[f"{category}_GLU_involved"] = str(_int(summary, category, "GLU_involved_count"))
            row[f"{category}_GLU_GLU"] = str(_int(summary, category, "GLU_GLU_count"))
        rows.append(row)
    return rows


def _markdown_table(rows: list[dict[str, str]], columns: list[str]) -> str:
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(row.get(column, "") for column in columns) + " |")
    return "\n".join(lines)


def _block_count_from_model(model_name: str) -> int:
    prefix = "scaffold_blocks_"
    suffix = "_heavy_deduped"
    if not model_name.startswith(prefix) or not model_name.endswith(suffix):
        return 0
    block_text = model_name[len(prefix) : -len(suffix)]
    return len([part for part in block_text.split("_") if part])


def report_notes(summaries: dict[str, dict[str, dict[str, str]]]) -> list[str]:
    notes: list[str] = []
    scaffold = summaries.get("hexaplex_scaffold_only_complement_heavy_deduped")
    full = summaries.get("full_hexaplex_anti_parallel_30deg_ideal_heavy_deduped")

    if scaffold:
        within_glu_glu = _int(scaffold, "scaffold_within_block", "GLU_GLU_count")
        between_glu_glu = _int(scaffold, "scaffold_between_blocks", "GLU_GLU_count")
        if within_glu_glu >= between_glu_glu:
            notes.append(
                f"Full scaffold GLU-GLU contacts are mostly within-block ({within_glu_glu}) rather than "
                f"between-block ({between_glu_glu}) in the candidate map."
            )
        else:
            notes.append(
                f"Full scaffold GLU-GLU contacts are mostly between-block ({between_glu_glu}) rather than "
                f"within-block ({within_glu_glu}) in the candidate map."
            )

    if scaffold and full:
        scaffold_hexad_contacts = _int(full, "scaffold_hexad_or_other", "contact_count")
        scaffold_only_hexad_contacts = _int(scaffold, "scaffold_hexad_or_other", "contact_count")
        if scaffold_hexad_contacts > scaffold_only_hexad_contacts:
            notes.append(
                "Full Hexaplex adds scaffold-hexad/other contacts relative to scaffold-only, consistent with "
                "component coupling in the assembled model."
            )
        else:
            notes.append("Full Hexaplex does not add scaffold-hexad/other contacts relative to scaffold-only in this table.")

    ladder = [
        (model_name, summary)
        for model_name, summary in summaries.items()
        if model_name.startswith("scaffold_blocks_") and model_name.endswith("_heavy_deduped")
    ]
    ladder.sort(key=lambda item: _block_count_from_model(item[0]))
    between_counts = [
        (_block_count_from_model(model_name), _int(summary, "scaffold_between_blocks", "contact_count"))
        for model_name, summary in ladder
    ]
    increases = any(later > earlier for (_, earlier), (_, later) in zip(between_counts, between_counts[1:]))
    if increases:
        notes.append(
            "Inter-block contacts increase as scaffold blocks are added in the ladder, consistent with repeated "
            "contact-network stabilization during multi-block assembly."
        )
    else:
        notes.append("Inter-block contacts do not increase across the available scaffold ladder summaries.")

    notes.append(
        "Individual blocks appear geometrically helical, while multi-block assembly appears to add repeated "
        "contact-network stabilization; this does not prove temporal order."
    )
    return notes


def write_markdown_report(summaries: dict[str, dict[str, dict[str, str]]], path: Path) -> None:
    table_rows = comparison_table_rows(summaries)
    columns = ["model_name"]
    for category in CATEGORIES:
        columns.extend(
            [
                f"{category}_contacts",
                f"{category}_GLU_involved",
                f"{category}_GLU_GLU",
            ]
        )

    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Block Contact Analysis",
        "",
        "## Scientific cautions",
        "",
        "The d ~= 4.5 A diffraction feature is treated as a reciprocal-space scaffold signature in the working hypothesis. It does not imply literal 4.5 A real-space atom contacts. Contact maps are real-space proximity summaries only. Diffraction intensities from component structures are comparative controls, not additive decompositions.",
        "",
        "The block mapping is a candidate contiguous-residue mapping and must be validated against PyMOL colored strand paths. Current ladder output shows that one candidate block already spans nearly full angular coverage and axial length, so these blocks may represent long folded or twisted scaffold paths rather than angular wedge units.",
        "",
        "`hexad_or_other` assignments use base-like atom names such as N1/C2/N3/C4/N5/C6/OC2/OC4/OC6 when scaffold and hexad atoms share residue labels in full structures. This is an atom-level operational classification for comparison, not a final chemical annotation.",
        "",
        "## Scaffold And Full Hexaplex Summary",
        "",
        _markdown_table(table_rows, columns),
        "",
        "## Automatically generated notes",
        "",
    ]
    lines.extend(f"- {note}" for note in report_notes(summaries))
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    summaries = load_category_summaries(args.summary_dir)
    write_markdown_report(summaries, args.out_md)
    print(f"Wrote {args.out_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
