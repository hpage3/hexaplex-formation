#!/usr/bin/env python3
"""Summarize GLU-containing residue contacts from a contact-map CSV."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


MOTIF_TYPES = [
    "GLU-GLU",
    "GLU-any",
    "backbone_O_to_GLU_sidechain_O",
    "GLU_sidechain_O_to_GLU_sidechain_O",
    "backbone_N_to_GLU_sidechain_O",
    "other_GLU_contact",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contact-map", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    return parser.parse_args()


def _is_glu_sidechain_o(residue_name: str, atom_name: str) -> bool:
    return residue_name == "GLU" and atom_name in {"OE1", "OE2"}


def motif_counts_and_pairs(rows: list[dict[str, str]]) -> tuple[Counter[str], dict[str, tuple[int, float]]]:
    counts: Counter[str] = Counter({motif_type: 0 for motif_type in MOTIF_TYPES})
    residue_pairs: dict[str, tuple[int, float]] = {}

    for row in rows:
        res_i = row["residue_name_i"].upper()
        res_j = row["residue_name_j"].upper()
        atom_i = row["atom_i"].upper()
        atom_j = row["atom_j"].upper()
        has_glu = res_i == "GLU" or res_j == "GLU"
        if not has_glu:
            continue

        counted_specific = False
        if res_i == "GLU" and res_j == "GLU":
            counts["GLU-GLU"] += 1
        counts["GLU-any"] += 1

        i_side_o = _is_glu_sidechain_o(res_i, atom_i)
        j_side_o = _is_glu_sidechain_o(res_j, atom_j)
        if i_side_o and j_side_o:
            counts["GLU_sidechain_O_to_GLU_sidechain_O"] += 1
            counted_specific = True
        if (atom_i == "O" and j_side_o) or (atom_j == "O" and i_side_o):
            counts["backbone_O_to_GLU_sidechain_O"] += 1
            counted_specific = True
        if (atom_i == "N" and j_side_o) or (atom_j == "N" and i_side_o):
            counts["backbone_N_to_GLU_sidechain_O"] += 1
            counted_specific = True
        if not counted_specific:
            counts["other_GLU_contact"] += 1

        residue_pair = f"{row['residue_i']}--{row['residue_j']}"
        min_distance = float(row["min_distance_A"])
        pair_count, pair_min = residue_pairs.get(residue_pair, (0, min_distance))
        residue_pairs[residue_pair] = (pair_count + 1, min(pair_min, min_distance))

    return counts, residue_pairs


def main() -> int:
    args = parse_args()
    with args.contact_map.open("r", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    counts, residue_pairs = motif_counts_and_pairs(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["motif_type", "count"])
        writer.writeheader()
        for motif_type in MOTIF_TYPES:
            writer.writerow({"motif_type": motif_type, "count": counts[motif_type]})

    top_pairs_path = args.out.with_name(f"{args.out.stem}_top_residue_pairs.csv")
    sorted_pairs = sorted(residue_pairs.items(), key=lambda item: (-item[1][0], item[1][1], item[0]))
    with top_pairs_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["residue_pair", "count", "min_distance_A"])
        writer.writeheader()
        for residue_pair, (count, min_distance) in sorted_pairs:
            writer.writerow(
                {
                    "residue_pair": residue_pair,
                    "count": count,
                    "min_distance_A": f"{min_distance:.6f}",
                }
            )

    print(f"Wrote {args.out} and {top_pairs_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
