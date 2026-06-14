#!/usr/bin/env python3
"""Generate a PyMOL helper script for visual scaffold path map inspection."""

from __future__ import annotations

import argparse
import csv
from collections import OrderedDict
from pathlib import Path


COLORS = [
    "red",
    "orange",
    "yellow",
    "green",
    "cyan",
    "blue",
    "violet",
    "magenta",
    "salmon",
    "lime",
    "teal",
    "slate",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pdb",
        type=Path,
        default=Path("outputs/intermediates/normalized_structures/hexaplex_scaffold_only_complement_heavy_deduped.pdb"),
    )
    parser.add_argument("--map", type=Path, default=Path("inputs/metadata/scaffold_path_map_candidate.csv"))
    parser.add_argument("--out", type=Path, default=Path("outputs/reports/pymol_scaffold_path_map_helper.pml"))
    return parser.parse_args()


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def grouped_map_rows(rows: list[dict[str, str]]) -> "OrderedDict[tuple[str, str], list[dict[str, str]]]":
    grouped: "OrderedDict[tuple[str, str], list[dict[str, str]]]" = OrderedDict()
    for row in rows:
        strand_id = row.get("strand_id", "").strip()
        strand_label = row.get("strand_label", "").strip() or f"strand_{strand_id}"
        key = (strand_id, strand_label)
        grouped.setdefault(key, []).append(row)
    return grouped


def _selection_for_rows(rows: list[dict[str, str]]) -> str:
    by_chain: "OrderedDict[str, list[str]]" = OrderedDict()
    for row in rows:
        residue_number = row.get("residue_number", "").strip()
        insertion_code = row.get("insertion_code", "").strip()
        if not residue_number:
            continue
        resi = f"{residue_number}{insertion_code}"
        by_chain.setdefault(row.get("chain_id", "").strip(), []).append(resi)

    parts: list[str] = []
    for chain_id, residues in by_chain.items():
        unique_residues = list(dict.fromkeys(residues))
        resi_expr = "+".join(unique_residues)
        if chain_id:
            parts.append(f"(chain {chain_id} and resi {resi_expr})")
        else:
            parts.append(f"(resi {resi_expr})")
    return " or ".join(parts) if parts else "none"


def build_pml(pdb_path: Path, map_path: Path, rows: list[dict[str, str]]) -> str:
    lines = [
        "# PyMOL scaffold path map helper",
        "# Load this script in PyMOL, then compare these generated path colors",
        "# against the known colored strand/path representation.",
        "# This helper validates visual correspondence only; it does not prove biological truth.",
        f"load {pdb_path}, scaffold_path_map_model",
        "hide everything, scaffold_path_map_model",
        "show cartoon, scaffold_path_map_model",
        "show sticks, scaffold_path_map_model",
        "color gray70, scaffold_path_map_model",
        f"# Map CSV: {map_path}",
        "",
    ]
    for index, ((strand_id, strand_label), strand_rows) in enumerate(grouped_map_rows(rows).items(), start=1):
        if not strand_id:
            continue
        selection_name = f"path_{strand_label}".replace("-", "_").replace(" ", "_")
        color = COLORS[(index - 1) % len(COLORS)]
        selection = _selection_for_rows(strand_rows)
        lines.extend(
            [
                f"# strand_id={strand_id} strand_label={strand_label} residues={len(strand_rows)}",
                f"select {selection_name}, scaffold_path_map_model and ({selection})",
                f"color {color}, {selection_name}",
                f"show sticks, {selection_name}",
                f"label first {selection_name} and name CA, \"{strand_label}\"",
                f"print \"{strand_label}: {len(strand_rows)} residue(s), selection {selection_name}\"",
                "",
            ]
        )
    lines.extend(
        [
            "orient scaffold_path_map_model",
            "zoom scaffold_path_map_model",
            "# If generated colors do not match the reference colored paths, edit",
            "# inputs/metadata/scaffold_path_map_manual_template.csv into",
            "# inputs/metadata/scaffold_path_map_manual.csv and rerun the workflow.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    args = parse_args()
    text = build_pml(args.pdb, args.map, read_csv_rows(args.map))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(text, encoding="utf-8")
    print(f"Wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
