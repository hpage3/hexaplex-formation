#!/usr/bin/env python3
"""Attribute heavy-atom distance-window pairs to residue, chain, and atom classes."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from hexaplex_formation.pdb_utils import PDBAtom, chain_ids, dedupe_exact_atoms, heavy_atoms, load_pdb_atoms  # noqa: E402


BASE_RESIDUES = {"CYP", "MEP"}
BACKBONE_LIKE_ATOMS = {"N", "CA", "C", "O", "OXT"}
WINDOWS = {
    "3p0": (2.9, 3.1),
    "3p4": (3.3, 3.5),
    "4p5_5p0": (4.5, 5.0),
}

PAIR_COLUMNS = [
    "structure_id",
    "distance_window",
    "distance_A",
    "atom_i_serial",
    "atom_j_serial",
    "chain_i",
    "chain_j",
    "chain_pair_id",
    "residue_i",
    "residue_j",
    "residue_pair_class",
    "residue_sequence_relation",
    "chain_relation",
    "atom_i_name",
    "atom_j_name",
    "atom_pair_class",
    "unit_i",
    "unit_j",
]

SUMMARY_COLUMNS = [
    "structure_id",
    "source_pdb",
    "distance_window",
    "residue_pair_class",
    "chain_relation",
    "residue_sequence_relation",
    "atom_pair_class",
    "pair_count",
    "fraction_of_window",
    "mean_distance_A",
]

CHAIN_COLUMNS = [
    "structure_id",
    "distance_window",
    "chain_pair_id",
    "chain_relation",
    "pair_count",
    "fraction_of_window",
]

RESIDUE_CLASS_COLUMNS = [
    "structure_id",
    "distance_window",
    "residue_pair_class",
    "pair_count",
    "fraction_of_window",
]


@dataclass(frozen=True)
class StructureInput:
    structure_id: str
    path: Path


@dataclass(frozen=True)
class PreparedStructure:
    structure_id: str
    path: Path
    raw_atom_count: int
    deduped_atom_count: int
    heavy_atom_count: int
    atoms: tuple[PDBAtom, ...]
    unit_lookup: dict[tuple[str, str, int | None, str], int]
    warnings: tuple[str, ...]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--full-pdb",
        type=Path,
        default=Path("outputs/intermediates/ai_candidate_inputs/full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain.pdb"),
    )
    parser.add_argument(
        "--central6-pdb",
        type=Path,
        default=Path("outputs/mini_hexaplex/structures/mini_hexaplex_central6_units.pdb"),
    )
    parser.add_argument(
        "--central7-pdb",
        type=Path,
        default=Path("outputs/mini_hexaplex/structures/mini_hexaplex_central7_units.pdb"),
    )
    parser.add_argument("--summary-csv", type=Path, default=Path("outputs/metrics/distance_window_attribution_summary.csv"))
    parser.add_argument(
        "--pairs-sample-csv",
        type=Path,
        default=Path("outputs/metrics/distance_window_attribution_pairs_sample.csv"),
    )
    parser.add_argument(
        "--by-chain-pair-csv",
        type=Path,
        default=Path("outputs/metrics/distance_window_attribution_by_chain_pair.csv"),
    )
    parser.add_argument(
        "--by-residue-class-csv",
        type=Path,
        default=Path("outputs/metrics/distance_window_attribution_by_residue_class.csv"),
    )
    parser.add_argument("--report", type=Path, default=Path("outputs/reports/distance_window_attribution_report.md"))
    parser.add_argument("--plot-dir", type=Path, default=Path("outputs/plots/distance_window_attribution"))
    parser.add_argument("--pairs-sample-limit", type=int, default=2000)
    return parser.parse_args()


def residue_key(atom: PDBAtom) -> tuple[str, str, int | None, str]:
    return atom.chain_id, atom.residue_name, atom.residue_number, atom.insertion_code


def residue_order_by_chain(atoms: list[PDBAtom]) -> dict[str, list[tuple[str, int | None, str]]]:
    residues_by_chain: dict[str, list[tuple[str, int | None, str]]] = defaultdict(list)
    seen: set[tuple[str, str, int | None, str]] = set()
    for atom in atoms:
        key = residue_key(atom)
        if key in seen:
            continue
        seen.add(key)
        residues_by_chain[atom.chain_id].append((atom.residue_name, atom.residue_number, atom.insertion_code))
    return dict(residues_by_chain)


def infer_unit_lookup(atoms: list[PDBAtom]) -> tuple[dict[tuple[str, str, int | None, str], int], list[str]]:
    lookup: dict[tuple[str, str, int | None, str], int] = {}
    warnings: list[str] = []
    for chain_id, residues in residue_order_by_chain(atoms).items():
        if len(residues) % 2 != 0:
            warnings.append(f"chain {chain_id or '<blank>'} has odd residue count {len(residues)}; final unit assignment may be approximate")
        for offset, (residue_name, residue_number, insertion_code) in enumerate(residues):
            lookup[(chain_id, residue_name, residue_number, insertion_code)] = offset // 2 + 1
    return lookup, warnings


def prepare_structure(structure_id: str, path: Path) -> PreparedStructure:
    if not path.exists():
        raise FileNotFoundError(f"Missing PDB for {structure_id}: {path}")
    raw = load_pdb_atoms(path)
    deduped = dedupe_exact_atoms(raw)
    heavy = heavy_atoms(deduped)
    lookup, warnings = infer_unit_lookup(deduped)
    if len(raw) != len(deduped):
        warnings.append(f"removed {len(raw) - len(deduped)} exact duplicate atoms")
    chains = chain_ids(deduped)
    if len(chains) != 6:
        warnings.append(f"expected six chains; found {len(chains)} chain IDs: {','.join(chains)}")
    return PreparedStructure(
        structure_id=structure_id,
        path=path,
        raw_atom_count=len(raw),
        deduped_atom_count=len(deduped),
        heavy_atom_count=len(heavy),
        atoms=tuple(heavy),
        unit_lookup=lookup,
        warnings=tuple(warnings),
    )


def residue_group(atom: PDBAtom) -> str:
    if atom.residue_name in BASE_RESIDUES:
        return "CYP_MEP"
    if atom.residue_name == "GLU":
        return "GLU"
    return "other"


def residue_pair_class(atom_a: PDBAtom, atom_b: PDBAtom) -> str:
    left = residue_group(atom_a)
    right = residue_group(atom_b)
    if left == "CYP_MEP" and right == "CYP_MEP":
        return "CYP_MEP_vs_CYP_MEP"
    if {left, right} == {"CYP_MEP", "GLU"}:
        return "CYP_MEP_vs_GLU"
    if left == "GLU" and right == "GLU":
        return "GLU_vs_GLU"
    return "other"


def atom_role(atom: PDBAtom) -> str:
    if atom.atom_name.strip().upper() in BACKBONE_LIKE_ATOMS:
        return "backbone_like"
    if atom.residue_name in BASE_RESIDUES:
        return "base_like"
    if atom.residue_name == "GLU":
        return "scaffold_linker"
    return "other"


def atom_pair_class(atom_a: PDBAtom, atom_b: PDBAtom) -> str:
    roles = sorted([atom_role(atom_a), atom_role(atom_b)])
    return f"{roles[0]}_vs_{roles[1]}"


def chain_pair_id(atom_a: PDBAtom, atom_b: PDBAtom) -> str:
    if atom_a.chain_id == atom_b.chain_id:
        return f"{atom_a.chain_id}-{atom_b.chain_id}"
    left, right = sorted([atom_a.chain_id, atom_b.chain_id])
    return f"{left}-{right}"


def sequence_relation(unit_a: int | None, unit_b: int | None) -> str:
    if unit_a is None or unit_b is None:
        return "unknown"
    diff = abs(unit_a - unit_b)
    if diff == 0:
        return "same_unit"
    if diff == 1:
        return "adjacent_unit"
    if diff in {2, 3}:
        return "near_axial"
    return "longer_range"


def squared_distance(atom_a: PDBAtom, atom_b: PDBAtom) -> float:
    dx = atom_a.x - atom_b.x
    dy = atom_a.y - atom_b.y
    dz = atom_a.z - atom_b.z
    return dx * dx + dy * dy + dz * dz


def matching_window(distance_a: float) -> str | None:
    for window_id, (low, high) in WINDOWS.items():
        if low <= distance_a <= high:
            return window_id
    return None


def same_residue(atom_a: PDBAtom, atom_b: PDBAtom) -> bool:
    return residue_key(atom_a) == residue_key(atom_b)


def pair_rows_for_structure(structure: PreparedStructure) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    max_high = max(high for _, high in WINDOWS.values())
    max_high_sq = max_high * max_high
    atoms = structure.atoms
    for i, atom_a in enumerate(atoms):
        for j in range(i + 1, len(atoms)):
            atom_b = atoms[j]
            if same_residue(atom_a, atom_b):
                continue
            dist_sq = squared_distance(atom_a, atom_b)
            if dist_sq > max_high_sq:
                continue
            distance_a = math.sqrt(dist_sq)
            window_id = matching_window(distance_a)
            if window_id is None:
                continue
            unit_a = structure.unit_lookup.get(residue_key(atom_a))
            unit_b = structure.unit_lookup.get(residue_key(atom_b))
            rows.append(
                {
                    "structure_id": structure.structure_id,
                    "distance_window": window_id,
                    "distance_A": f"{distance_a:.4f}",
                    "atom_i_serial": str(atom_a.atom_serial or ""),
                    "atom_j_serial": str(atom_b.atom_serial or ""),
                    "chain_i": atom_a.chain_id,
                    "chain_j": atom_b.chain_id,
                    "chain_pair_id": chain_pair_id(atom_a, atom_b),
                    "residue_i": f"{atom_a.residue_name}{atom_a.residue_number or ''}",
                    "residue_j": f"{atom_b.residue_name}{atom_b.residue_number or ''}",
                    "residue_pair_class": residue_pair_class(atom_a, atom_b),
                    "residue_sequence_relation": sequence_relation(unit_a, unit_b),
                    "chain_relation": "same_chain" if atom_a.chain_id == atom_b.chain_id else "interchain",
                    "atom_i_name": atom_a.atom_name,
                    "atom_j_name": atom_b.atom_name,
                    "atom_pair_class": atom_pair_class(atom_a, atom_b),
                    "unit_i": str(unit_a or ""),
                    "unit_j": str(unit_b or ""),
                }
            )
    return rows


def summarize_pairs(rows: list[dict[str, str]], structures: dict[str, PreparedStructure]) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    totals_by_window = Counter((row["structure_id"], row["distance_window"]) for row in rows)
    grouped: dict[tuple[str, str, str, str, str, str], list[float]] = defaultdict(list)
    chain_counts = Counter()
    residue_counts = Counter()
    for row in rows:
        grouped[
            (
                row["structure_id"],
                row["distance_window"],
                row["residue_pair_class"],
                row["chain_relation"],
                row["residue_sequence_relation"],
                row["atom_pair_class"],
            )
        ].append(float(row["distance_A"]))
        chain_counts[(row["structure_id"], row["distance_window"], row["chain_pair_id"], row["chain_relation"])] += 1
        residue_counts[(row["structure_id"], row["distance_window"], row["residue_pair_class"])] += 1

    summary_rows: list[dict[str, str]] = []
    for key, distances in sorted(grouped.items()):
        structure_id, window_id, residue_class, chain_relation, relation, atom_class = key
        total = totals_by_window[(structure_id, window_id)]
        summary_rows.append(
            {
                "structure_id": structure_id,
                "source_pdb": str(structures[structure_id].path),
                "distance_window": window_id,
                "residue_pair_class": residue_class,
                "chain_relation": chain_relation,
                "residue_sequence_relation": relation,
                "atom_pair_class": atom_class,
                "pair_count": str(len(distances)),
                "fraction_of_window": f"{len(distances) / total:.6f}" if total else "",
                "mean_distance_A": f"{sum(distances) / len(distances):.4f}",
            }
        )

    chain_rows: list[dict[str, str]] = []
    for (structure_id, window_id, pair_id, relation), count in sorted(chain_counts.items()):
        total = totals_by_window[(structure_id, window_id)]
        chain_rows.append(
            {
                "structure_id": structure_id,
                "distance_window": window_id,
                "chain_pair_id": pair_id,
                "chain_relation": relation,
                "pair_count": str(count),
                "fraction_of_window": f"{count / total:.6f}" if total else "",
            }
        )

    residue_rows: list[dict[str, str]] = []
    for (structure_id, window_id, residue_class), count in sorted(residue_counts.items()):
        total = totals_by_window[(structure_id, window_id)]
        residue_rows.append(
            {
                "structure_id": structure_id,
                "distance_window": window_id,
                "residue_pair_class": residue_class,
                "pair_count": str(count),
                "fraction_of_window": f"{count / total:.6f}" if total else "",
            }
        )
    return summary_rows, chain_rows, residue_rows


def write_csv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def svg_bar_chart(title: str, groups: list[str], series: dict[str, list[int]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    width = 920
    height = 440
    left = 80
    bottom = 350
    max_value = max([value for values in series.values() for value in values] + [1])
    colors = ["#4c78a8", "#f58518", "#54a24b", "#b279a2", "#e45756"]
    bar_w = max(12, min(34, int(160 / max(1, len(series)))))
    group_gap = 50
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width / 2:.0f}" y="30" text-anchor="middle" font-family="Arial" font-size="18">{title}</text>',
        f'<line x1="{left}" y1="55" x2="{left}" y2="{bottom}" stroke="#333"/>',
        f'<line x1="{left}" y1="{bottom}" x2="{width - 30}" y2="{bottom}" stroke="#333"/>',
    ]
    series_names = list(series)
    for group_index, group in enumerate(groups):
        x0 = left + 30 + group_index * (len(series_names) * (bar_w + 4) + group_gap)
        for series_index, name in enumerate(series_names):
            value = series[name][group_index]
            h = value / max_value * 250 if max_value else 0
            x = x0 + series_index * (bar_w + 4)
            y = bottom - h
            parts.append(f'<rect x="{x}" y="{y:.2f}" width="{bar_w}" height="{h:.2f}" fill="{colors[series_index % len(colors)]}"/>')
            parts.append(f'<text x="{x + bar_w / 2:.0f}" y="{y - 5:.2f}" text-anchor="middle" font-family="Arial" font-size="9">{value}</text>')
        parts.append(f'<text x="{x0 + (len(series_names) * (bar_w + 4)) / 2:.0f}" y="{bottom + 18}" text-anchor="middle" font-family="Arial" font-size="12">{group}</text>')
    legend_x = 650
    for index, name in enumerate(series_names):
        y = 55 + index * 20
        parts.append(f'<rect x="{legend_x}" y="{y}" width="13" height="13" fill="{colors[index % len(colors)]}"/>')
        parts.append(f'<text x="{legend_x + 18}" y="{y + 11}" font-family="Arial" font-size="12">{name}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def write_plots(residue_rows: list[dict[str, str]], chain_rows: list[dict[str, str]], plot_dir: Path) -> list[Path]:
    written: list[Path] = []
    structures = ["central6", "central7", "full"]
    for window_id in WINDOWS:
        classes = sorted({row["residue_pair_class"] for row in residue_rows if row["distance_window"] == window_id})
        series = {
            residue_class: [
                sum(
                    int(row["pair_count"])
                    for row in residue_rows
                    if row["distance_window"] == window_id and row["structure_id"] == structure_id and row["residue_pair_class"] == residue_class
                )
                for structure_id in structures
            ]
            for residue_class in classes
        }
        out = plot_dir / f"{window_id}_residue_class_contributions.svg"
        svg_bar_chart(f"{window_id} residue-pair class contributions", structures, series, out)
        written.append(out)

    relation_series = {
        relation: [
            sum(
                int(row["pair_count"])
                for row in chain_rows
                if row["distance_window"] == "4p5_5p0" and row["structure_id"] == structure_id and row["chain_relation"] == relation
            )
            for structure_id in structures
        ]
        for relation in ["same_chain", "interchain"]
    }
    out = plot_dir / "4p5_5p0_chain_relation_contributions.svg"
    svg_bar_chart("4p5_5p0 same-chain vs interchain contributions", structures, relation_series, out)
    written.append(out)

    pairs = sorted({row["chain_pair_id"] for row in chain_rows if row["distance_window"] == "4p5_5p0" and row["chain_relation"] == "interchain"})
    top_pairs = pairs[:15]
    pair_series = {
        structure_id: [
            sum(
                int(row["pair_count"])
                for row in chain_rows
                if row["distance_window"] == "4p5_5p0" and row["structure_id"] == structure_id and row["chain_pair_id"] == pair_id
            )
            for pair_id in top_pairs
        ]
        for structure_id in structures
    }
    out = plot_dir / "4p5_5p0_chain_pair_contributions.svg"
    svg_bar_chart("4p5_5p0 interchain chain-pair contributions", top_pairs, pair_series, out)
    written.append(out)
    return written


def top_residue_class(residue_rows: list[dict[str, str]], structure_id: str, window_id: str) -> tuple[str, int] | None:
    matches = [row for row in residue_rows if row["structure_id"] == structure_id and row["distance_window"] == window_id]
    if not matches:
        return None
    best = max(matches, key=lambda row: int(row["pair_count"]))
    return best["residue_pair_class"], int(best["pair_count"])


def chain_relation_counts(chain_rows: list[dict[str, str]], structure_id: str, window_id: str) -> dict[str, int]:
    counts = Counter()
    for row in chain_rows:
        if row["structure_id"] == structure_id and row["distance_window"] == window_id:
            counts[row["chain_relation"]] += int(row["pair_count"])
    return dict(counts)


def residue_class_counts(residue_rows: list[dict[str, str]], structure_id: str, window_id: str) -> dict[str, int]:
    counts = Counter()
    for row in residue_rows:
        if row["structure_id"] == structure_id and row["distance_window"] == window_id:
            counts[row["residue_pair_class"]] += int(row["pair_count"])
    return dict(counts)


def sequence_relation_counts(summary_rows: list[dict[str, str]], structure_id: str, window_id: str) -> dict[str, int]:
    counts = Counter()
    for row in summary_rows:
        if row["structure_id"] == structure_id and row["distance_window"] == window_id:
            counts[row["residue_sequence_relation"]] += int(row["pair_count"])
    return dict(counts)


def report_text(
    structures: dict[str, PreparedStructure],
    summary_rows: list[dict[str, str]],
    residue_rows: list[dict[str, str]],
    chain_rows: list[dict[str, str]],
    plot_paths: list[Path],
    args: argparse.Namespace,
) -> str:
    lines = [
        "# Distance-Window Attribution Analysis",
        "",
        "This is a structural heavy-atom distance attribution analysis. It does not compute diffraction intensities directly. Pair counts identify enriched geometric contributors within X-ray-relevant distance windows in the molecular models.",
        "",
        "## Inputs",
        "",
    ]
    for structure_id, structure in structures.items():
        lines.append(
            f"- `{structure_id}`: `{structure.path}` ({structure.heavy_atom_count} heavy atoms; {structure.raw_atom_count} raw atom records)"
        )
    lines.extend(
        [
            "",
            "Selected full model: `outputs/intermediates/ai_candidate_inputs/full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain.pdb`, because it is the six-chain deduped full model used by prior workflows.",
            "",
            "Hydrogens and exact duplicate atom records are excluded. Same-residue atom pairs are excluded. No bond/topology information is used, so directly bonded or near-bonded pairs in different residues are not removed; this is geometry-only.",
            "",
            "Distance windows: 3p0 = 2.9-3.1 A, 3p4 = 3.3-3.5 A, 4p5_5p0 = 4.5-5.0 A.",
            "",
            "## Dominant Residue-Pair Classes",
            "",
            "| Window | Structure | Dominant class | Pair count |",
            "|---|---|---|---:|",
        ]
    )
    for window_id in WINDOWS:
        for structure_id in ["central6", "central7", "full"]:
            top = top_residue_class(residue_rows, structure_id, window_id)
            if top:
                lines.append(f"| {window_id} | {structure_id} | {top[0]} | {top[1]} |")

    lines.extend(
        [
            "",
            "## Growth Across 6, 7, And Full Models",
            "",
            "| Window | Structure | Total pairs | CYP/MEP-CYP/MEP | CYP/MEP-GLU | GLU-GLU | Dominant sequence relation |",
            "|---|---|---:|---:|---:|---:|---|",
        ]
    )
    for window_id in WINDOWS:
        for structure_id in ["central6", "central7", "full"]:
            residue_counts = residue_class_counts(residue_rows, structure_id, window_id)
            seq_counts = sequence_relation_counts(summary_rows, structure_id, window_id)
            total = sum(residue_counts.values())
            dominant_relation = "none"
            if seq_counts:
                dominant_relation = max(seq_counts.items(), key=lambda item: item[1])[0]
            lines.append(
                f"| {window_id} | {structure_id} | {total} | "
                f"{residue_counts.get('CYP_MEP_vs_CYP_MEP', 0)} | "
                f"{residue_counts.get('CYP_MEP_vs_GLU', 0)} | "
                f"{residue_counts.get('GLU_vs_GLU', 0)} | {dominant_relation} |"
            )

    lines.extend(
        [
            "",
            "## 4.5-5.0 A Chain Relation",
            "",
            "| Structure | Same-chain pairs | Interchain pairs | Interpretation |",
            "|---|---:|---:|---|",
        ]
    )
    for structure_id in ["central6", "central7", "full"]:
        counts = chain_relation_counts(chain_rows, structure_id, "4p5_5p0")
        same = counts.get("same_chain", 0)
        inter = counts.get("interchain", 0)
        if inter > same:
            interp = "mostly interchain"
        elif same > inter:
            interp = "mostly same-chain"
        else:
            interp = "mixed"
        lines.append(f"| {structure_id} | {same} | {inter} | {interp} |")

    lines.extend(
        [
            "",
            "## Interpretation For Nick",
            "",
            "- These pair counts are geometric contributors inside distance windows, not direct diffraction intensities and not proof that a class causes an experimental feature.",
            "- In these models, the 3.0 A window is dominated by CYP/MEP-GLU pairs, while the 3.4 A window is dominated by CYP/MEP-CYP/MEP pairs.",
            "- The 4.5-5.0 A window is dominated by CYP/MEP-CYP/MEP pairs in central6, central7, and full; CYP/MEP-GLU and GLU-GLU pairs are present but smaller contributors.",
            "- The 4.5-5.0 A window is mixed but slightly same-chain enriched in all three structures under this geometry-only counting scheme.",
            "- The 4.5-5.0 A count grows from central6 to central7 to full primarily because CYP/MEP-involving geometric opportunities grow with model length, with accompanying GLU/scaffold contributions.",
            "- Growth from 6 to 7 to full indicates increasing geometric opportunities for those pair classes in the model, not necessarily proportional intensity growth.",
            "",
            "## Outputs",
            "",
            f"- Summary CSV: `{args.summary_csv}`",
            f"- Pair sample CSV: `{args.pairs_sample_csv}`",
            f"- Chain-pair CSV: `{args.by_chain_pair_csv}`",
            f"- Residue-class CSV: `{args.by_residue_class_csv}`",
            f"- Plot directory: `{args.plot_dir}`",
        ]
    )
    for path in plot_paths:
        lines.append(f"- Plot: `{path}`")
    warnings = [warning for structure in structures.values() for warning in structure.warnings]
    if warnings:
        lines.extend(["", "## Warnings", ""])
        for warning in warnings:
            lines.append(f"- {warning}")
    return "\n".join(lines) + "\n"


def run(args: argparse.Namespace) -> dict[str, object]:
    inputs = [
        StructureInput("full", args.full_pdb),
        StructureInput("central6", args.central6_pdb),
        StructureInput("central7", args.central7_pdb),
    ]
    structures = {item.structure_id: prepare_structure(item.structure_id, item.path) for item in inputs}
    all_pair_rows: list[dict[str, str]] = []
    for structure in structures.values():
        all_pair_rows.extend(pair_rows_for_structure(structure))

    summary_rows, chain_rows, residue_rows = summarize_pairs(all_pair_rows, structures)
    write_csv(args.summary_csv, summary_rows, SUMMARY_COLUMNS)
    write_csv(args.pairs_sample_csv, all_pair_rows[: args.pairs_sample_limit], PAIR_COLUMNS)
    write_csv(args.by_chain_pair_csv, chain_rows, CHAIN_COLUMNS)
    write_csv(args.by_residue_class_csv, residue_rows, RESIDUE_CLASS_COLUMNS)
    plot_paths = write_plots(residue_rows, chain_rows, args.plot_dir)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(report_text(structures, summary_rows, residue_rows, chain_rows, plot_paths, args), encoding="utf-8")
    return {
        "pair_rows": len(all_pair_rows),
        "summary_rows": len(summary_rows),
        "chain_rows": len(chain_rows),
        "residue_rows": len(residue_rows),
        "plots": plot_paths,
        "structures": structures,
    }


def main() -> None:
    result = run(parse_args())
    print(f"Wrote {result['pair_rows']} attributed pair rows before sampling")
    print(f"Wrote {result['summary_rows']} summary rows")
    print(f"Wrote {result['chain_rows']} chain-pair rows")
    print(f"Wrote {result['residue_rows']} residue-class rows")
    print(f"Wrote {len(result['plots'])} plots")


if __name__ == "__main__":
    main()
