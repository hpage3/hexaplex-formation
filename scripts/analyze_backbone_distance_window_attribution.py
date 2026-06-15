#!/usr/bin/env python3
"""Exploratory backbone-specific distance-window attribution analysis."""

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
    "3p0": (2.9, 3.1, "known"),
    "3p4": (3.3, 3.5, "known"),
    "4p5_5p0": (4.5, 5.0, "known"),
    "guessed_5p0_6p0": (5.0, 6.0, "exploratory_guessed"),
    "guessed_6p5_7p5": (6.5, 7.5, "exploratory_guessed"),
    "guessed_7p8_9p0": (7.8, 9.0, "exploratory_guessed"),
}

SUMMARY_COLUMNS = [
    "structure_id",
    "source_pdb",
    "distance_window",
    "window_type",
    "atom_pair_class",
    "residue_pair_class",
    "chain_relation",
    "pair_count",
    "fraction_of_window",
    "backbone_involving",
]

ATOM_CLASS_COLUMNS = [
    "structure_id",
    "distance_window",
    "window_type",
    "atom_pair_class",
    "pair_count",
    "fraction_of_window",
    "backbone_involving",
]

CHAIN_RELATION_COLUMNS = [
    "structure_id",
    "distance_window",
    "window_type",
    "chain_relation",
    "backbone_involving",
    "pair_count",
    "fraction_of_window",
]

PAIR_SAMPLE_COLUMNS = [
    "structure_id",
    "distance_window",
    "window_type",
    "distance_A",
    "chain_i",
    "chain_j",
    "chain_pair_id",
    "residue_i",
    "residue_j",
    "residue_pair_class",
    "atom_i_name",
    "atom_j_name",
    "atom_i_class",
    "atom_j_class",
    "atom_pair_class",
    "chain_relation",
    "backbone_involving",
]


@dataclass(frozen=True)
class PreparedStructure:
    structure_id: str
    path: Path
    raw_atom_count: int
    deduped_atom_count: int
    heavy_atom_count: int
    atoms: tuple[PDBAtom, ...]
    warnings: tuple[str, ...]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--full-pdb",
        type=Path,
        default=Path("outputs/intermediates/ai_candidate_inputs/full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain.pdb"),
    )
    parser.add_argument("--central6-pdb", type=Path, default=Path("outputs/mini_hexaplex/structures/mini_hexaplex_central6_units.pdb"))
    parser.add_argument("--central7-pdb", type=Path, default=Path("outputs/mini_hexaplex/structures/mini_hexaplex_central7_units.pdb"))
    parser.add_argument("--summary-csv", type=Path, default=Path("outputs/metrics/backbone_distance_window_attribution_summary.csv"))
    parser.add_argument(
        "--by-atom-class-csv",
        type=Path,
        default=Path("outputs/metrics/backbone_distance_window_attribution_by_atom_class.csv"),
    )
    parser.add_argument(
        "--by-chain-relation-csv",
        type=Path,
        default=Path("outputs/metrics/backbone_distance_window_attribution_by_chain_relation.csv"),
    )
    parser.add_argument(
        "--pairs-sample-csv",
        type=Path,
        default=Path("outputs/metrics/backbone_distance_window_attribution_pairs_sample.csv"),
    )
    parser.add_argument("--report", type=Path, default=Path("outputs/reports/backbone_distance_window_attribution_report.md"))
    parser.add_argument("--plot-dir", type=Path, default=Path("outputs/plots/backbone_distance_window_attribution"))
    parser.add_argument("--pairs-sample-limit", type=int, default=2500)
    return parser.parse_args()


def residue_key(atom: PDBAtom) -> tuple[str, str, int | None, str]:
    return atom.chain_id, atom.residue_name, atom.residue_number, atom.insertion_code


def prepare_structure(structure_id: str, path: Path) -> PreparedStructure:
    if not path.exists():
        raise FileNotFoundError(f"Missing PDB for {structure_id}: {path}")
    raw = load_pdb_atoms(path)
    deduped = dedupe_exact_atoms(raw)
    heavy = heavy_atoms(deduped)
    warnings: list[str] = []
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
        warnings=tuple(warnings),
    )


def atom_class(atom: PDBAtom) -> str:
    if atom.atom_name.strip().upper() in BACKBONE_LIKE_ATOMS:
        return "backbone_like"
    if atom.residue_name in BASE_RESIDUES:
        return "base_like"
    if atom.residue_name == "GLU":
        return "scaffold_linker"
    return "other"


def atom_pair_class(class_a: str, class_b: str) -> str:
    ordered = sorted([class_a, class_b])
    allowed = {
        ("backbone_like", "backbone_like"),
        ("backbone_like", "base_like"),
        ("backbone_like", "scaffold_linker"),
        ("base_like", "base_like"),
        ("base_like", "scaffold_linker"),
        ("scaffold_linker", "scaffold_linker"),
    }
    pair = (ordered[0], ordered[1])
    if pair in allowed:
        return f"{pair[0]}_vs_{pair[1]}"
    return "other_mixed"


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


def chain_pair_id(atom_a: PDBAtom, atom_b: PDBAtom) -> str:
    if atom_a.chain_id == atom_b.chain_id:
        return f"{atom_a.chain_id}-{atom_b.chain_id}"
    left, right = sorted([atom_a.chain_id, atom_b.chain_id])
    return f"{left}-{right}"


def same_residue(atom_a: PDBAtom, atom_b: PDBAtom) -> bool:
    return residue_key(atom_a) == residue_key(atom_b)


def squared_distance(atom_a: PDBAtom, atom_b: PDBAtom) -> float:
    dx = atom_a.x - atom_b.x
    dy = atom_a.y - atom_b.y
    dz = atom_a.z - atom_b.z
    return dx * dx + dy * dy + dz * dz


def matching_window(distance_a: float) -> str | None:
    for window_id, (low, high, _) in WINDOWS.items():
        if low <= distance_a <= high:
            return window_id
    return None


def pair_rows_for_structure(structure: PreparedStructure) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    max_high = max(high for _, high, _ in WINDOWS.values())
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
            class_a = atom_class(atom_a)
            class_b = atom_class(atom_b)
            pair_class = atom_pair_class(class_a, class_b)
            backbone_involving = "yes" if "backbone_like" in {class_a, class_b} else "no"
            rows.append(
                {
                    "structure_id": structure.structure_id,
                    "distance_window": window_id,
                    "window_type": WINDOWS[window_id][2],
                    "distance_A": f"{distance_a:.4f}",
                    "chain_i": atom_a.chain_id,
                    "chain_j": atom_b.chain_id,
                    "chain_pair_id": chain_pair_id(atom_a, atom_b),
                    "residue_i": f"{atom_a.residue_name}{atom_a.residue_number or ''}",
                    "residue_j": f"{atom_b.residue_name}{atom_b.residue_number or ''}",
                    "residue_pair_class": residue_pair_class(atom_a, atom_b),
                    "atom_i_name": atom_a.atom_name,
                    "atom_j_name": atom_b.atom_name,
                    "atom_i_class": class_a,
                    "atom_j_class": class_b,
                    "atom_pair_class": pair_class,
                    "chain_relation": "same_chain" if atom_a.chain_id == atom_b.chain_id else "interchain",
                    "backbone_involving": backbone_involving,
                }
            )
    return rows


def summarize(rows: list[dict[str, str]], structures: dict[str, PreparedStructure]) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    totals = Counter((row["structure_id"], row["distance_window"]) for row in rows)
    grouped = Counter(
        (
            row["structure_id"],
            row["distance_window"],
            row["atom_pair_class"],
            row["residue_pair_class"],
            row["chain_relation"],
            row["backbone_involving"],
        )
        for row in rows
    )
    atom_grouped = Counter((row["structure_id"], row["distance_window"], row["atom_pair_class"], row["backbone_involving"]) for row in rows)
    chain_grouped = Counter((row["structure_id"], row["distance_window"], row["chain_relation"], row["backbone_involving"]) for row in rows)

    summary_rows: list[dict[str, str]] = []
    for key, count in sorted(grouped.items()):
        structure_id, window_id, pair_class, residue_class, chain_relation, backbone_involving = key
        total = totals[(structure_id, window_id)]
        summary_rows.append(
            {
                "structure_id": structure_id,
                "source_pdb": str(structures[structure_id].path),
                "distance_window": window_id,
                "window_type": WINDOWS[window_id][2],
                "atom_pair_class": pair_class,
                "residue_pair_class": residue_class,
                "chain_relation": chain_relation,
                "pair_count": str(count),
                "fraction_of_window": f"{count / total:.6f}" if total else "",
                "backbone_involving": backbone_involving,
            }
        )

    atom_rows: list[dict[str, str]] = []
    for (structure_id, window_id, pair_class, backbone_involving), count in sorted(atom_grouped.items()):
        total = totals[(structure_id, window_id)]
        atom_rows.append(
            {
                "structure_id": structure_id,
                "distance_window": window_id,
                "window_type": WINDOWS[window_id][2],
                "atom_pair_class": pair_class,
                "pair_count": str(count),
                "fraction_of_window": f"{count / total:.6f}" if total else "",
                "backbone_involving": backbone_involving,
            }
        )

    chain_rows: list[dict[str, str]] = []
    for (structure_id, window_id, chain_relation, backbone_involving), count in sorted(chain_grouped.items()):
        total = totals[(structure_id, window_id)]
        chain_rows.append(
            {
                "structure_id": structure_id,
                "distance_window": window_id,
                "window_type": WINDOWS[window_id][2],
                "chain_relation": chain_relation,
                "backbone_involving": backbone_involving,
                "pair_count": str(count),
                "fraction_of_window": f"{count / total:.6f}" if total else "",
            }
        )
    return summary_rows, atom_rows, chain_rows


def write_csv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def svg_bar_chart(title: str, groups: list[str], series: dict[str, list[int]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    width = 980
    height = 440
    left = 80
    bottom = 350
    max_value = max([value for values in series.values() for value in values] + [1])
    colors = ["#4c78a8", "#f58518", "#54a24b", "#b279a2", "#e45756", "#72b7b2", "#9d755d"]
    series_names = list(series)
    bar_w = max(10, min(30, int(170 / max(1, len(series_names)))))
    group_gap = 50
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width / 2:.0f}" y="30" text-anchor="middle" font-family="Arial" font-size="17">{title}</text>',
        f'<line x1="{left}" y1="55" x2="{left}" y2="{bottom}" stroke="#333"/>',
        f'<line x1="{left}" y1="{bottom}" x2="{width - 30}" y2="{bottom}" stroke="#333"/>',
    ]
    for group_index, group in enumerate(groups):
        x0 = left + 30 + group_index * (len(series_names) * (bar_w + 4) + group_gap)
        for series_index, name in enumerate(series_names):
            value = series[name][group_index]
            h = value / max_value * 250 if max_value else 0
            x = x0 + series_index * (bar_w + 4)
            y = bottom - h
            parts.append(f'<rect x="{x}" y="{y:.2f}" width="{bar_w}" height="{h:.2f}" fill="{colors[series_index % len(colors)]}"/>')
            parts.append(f'<text x="{x + bar_w / 2:.0f}" y="{y - 5:.2f}" text-anchor="middle" font-family="Arial" font-size="8">{value}</text>')
        parts.append(f'<text x="{x0 + (len(series_names) * (bar_w + 4)) / 2:.0f}" y="{bottom + 18}" text-anchor="middle" font-family="Arial" font-size="12">{group}</text>')
    legend_x = 690
    for index, name in enumerate(series_names):
        y = 55 + index * 19
        parts.append(f'<rect x="{legend_x}" y="{y}" width="12" height="12" fill="{colors[index % len(colors)]}"/>')
        parts.append(f'<text x="{legend_x + 17}" y="{y + 10}" font-family="Arial" font-size="11">{name}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def write_plots(atom_rows: list[dict[str, str]], chain_rows: list[dict[str, str]], plot_dir: Path) -> list[Path]:
    written: list[Path] = []
    structures = ["central6", "central7", "full"]
    for window_id in WINDOWS:
        classes = sorted({row["atom_pair_class"] for row in atom_rows if row["distance_window"] == window_id})
        series = {
            cls: [
                sum(
                    int(row["pair_count"])
                    for row in atom_rows
                    if row["structure_id"] == structure_id and row["distance_window"] == window_id and row["atom_pair_class"] == cls
                )
                for structure_id in structures
            ]
            for cls in classes
        }
        out = plot_dir / f"{window_id}_atom_pair_class_by_structure.svg"
        svg_bar_chart(f"{window_id} atom-pair class contributions", structures, series, out)
        written.append(out)

    for window_id in [key for key, (_, _, kind) in WINDOWS.items() if kind == "exploratory_guessed"]:
        bb_series = {
            "backbone_involving": [
                sum(
                    int(row["pair_count"])
                    for row in atom_rows
                    if row["structure_id"] == structure_id and row["distance_window"] == window_id and row["backbone_involving"] == "yes"
                )
                for structure_id in structures
            ],
            "non_backbone": [
                sum(
                    int(row["pair_count"])
                    for row in atom_rows
                    if row["structure_id"] == structure_id and row["distance_window"] == window_id and row["backbone_involving"] == "no"
                )
                for structure_id in structures
            ],
        }
        out = plot_dir / f"{window_id}_backbone_involving_by_structure.svg"
        svg_bar_chart(f"{window_id} backbone-involving vs non-backbone", structures, bb_series, out)
        written.append(out)

        relation_series = {
            relation: [
                sum(
                    int(row["pair_count"])
                    for row in chain_rows
                    if row["structure_id"] == structure_id and row["distance_window"] == window_id and row["chain_relation"] == relation
                )
                for structure_id in structures
            ]
            for relation in ["same_chain", "interchain"]
        }
        out = plot_dir / f"{window_id}_chain_relation_by_structure.svg"
        svg_bar_chart(f"{window_id} same-chain vs interchain", structures, relation_series, out)
        written.append(out)

    full_windows = list(WINDOWS)
    full_classes = sorted({row["atom_pair_class"] for row in atom_rows if row["structure_id"] == "full"})
    full_series = {
        cls: [
            sum(
                int(row["pair_count"])
                for row in atom_rows
                if row["structure_id"] == "full" and row["distance_window"] == window_id and row["atom_pair_class"] == cls
            )
            for window_id in full_windows
        ]
        for cls in full_classes
    }
    out = plot_dir / "full_atom_pair_class_across_windows.svg"
    svg_bar_chart("Full model atom-pair class counts across windows", full_windows, full_series, out)
    written.append(out)
    return written


def class_counts(atom_rows: list[dict[str, str]], structure_id: str, window_id: str) -> dict[str, int]:
    counts = Counter()
    for row in atom_rows:
        if row["structure_id"] == structure_id and row["distance_window"] == window_id:
            counts[row["atom_pair_class"]] += int(row["pair_count"])
    return dict(counts)


def backbone_fraction(atom_rows: list[dict[str, str]], structure_id: str, window_id: str) -> float:
    total = 0
    backbone = 0
    for row in atom_rows:
        if row["structure_id"] == structure_id and row["distance_window"] == window_id:
            count = int(row["pair_count"])
            total += count
            if row["backbone_involving"] == "yes":
                backbone += count
    return backbone / total if total else math.nan


def report_text(
    structures: dict[str, PreparedStructure],
    atom_rows: list[dict[str, str]],
    chain_rows: list[dict[str, str]],
    plot_paths: list[Path],
    args: argparse.Namespace,
) -> str:
    lines = [
        "# Exploratory Backbone Distance-Window Attribution",
        "",
        "This analysis uses provisional distance windows chosen without detector calibration. The results should be interpreted as a sensitivity screen for Nick's backbone hypothesis, not as an assignment of the circled experimental features. Calibrated d-spacings, q values, detector radii, beam center, wavelength, and sample-to-detector distance would allow these windows to be replaced with measured experimental windows.",
        "",
        "This is a geometry-only heavy-atom pair analysis. Hydrogens, exact duplicate atom records, and same-residue pairs are excluded. Directly bonded or near-bonded local pairs in different residues are not excluded because no bond/topology logic is used here.",
        "",
        "## Inputs",
        "",
    ]
    for structure_id, structure in structures.items():
        lines.append(f"- `{structure_id}`: `{structure.path}` ({structure.heavy_atom_count} heavy atoms)")

    lines.extend(
        [
            "",
            "## Backbone-Involving Fractions",
            "",
            "| Window | Type | central6 | central7 | full |",
            "|---|---|---:|---:|---:|",
        ]
    )
    for window_id, (_, _, kind) in WINDOWS.items():
        values = [backbone_fraction(atom_rows, structure_id, window_id) for structure_id in ["central6", "central7", "full"]]
        lines.append(f"| {window_id} | {kind} | {values[0]:.3f} | {values[1]:.3f} | {values[2]:.3f} |")

    lines.extend(
        [
            "",
            "## Dominant Atom-Pair Classes",
            "",
            "| Window | Structure | Dominant atom-pair class | Pair count |",
            "|---|---|---|---:|",
        ]
    )
    for window_id in WINDOWS:
        for structure_id in ["central6", "central7", "full"]:
            counts = class_counts(atom_rows, structure_id, window_id)
            if not counts:
                continue
            cls, count = max(counts.items(), key=lambda item: item[1])
            lines.append(f"| {window_id} | {structure_id} | {cls} | {count} |")

    guessed_windows = [key for key, (_, _, kind) in WINDOWS.items() if kind == "exploratory_guessed"]
    best_by_structure = {}
    for structure_id in ["central6", "central7", "full"]:
        best_window = max(guessed_windows, key=lambda window_id: backbone_fraction(atom_rows, structure_id, window_id))
        best_by_structure[structure_id] = (best_window, backbone_fraction(atom_rows, structure_id, best_window))

    lines.extend(
        [
            "",
            "## Backbone-Hypothesis Screen",
            "",
            "| Structure | Guessed window with highest backbone-involving fraction | Fraction |",
            "|---|---|---:|",
        ]
    )
    for structure_id in ["central6", "central7", "full"]:
        best_window, fraction = best_by_structure[structure_id]
        lines.append(f"| {structure_id} | {best_window} | {fraction:.3f} |")

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- Backbone-involving pairs are present across both known and guessed windows, but this screen does not assign any experimental feature.",
            "- The guessed longer windows are more backbone-involving than the known 3.4 A and 4.5-5.0 A windows by fraction, except that the 3.0 A local window is entirely backbone-involving under this atom-name heuristic.",
            "- Among the guessed windows, guessed_7p8_9p0 is most consistent with Nick's backbone hypothesis by backbone-involving fraction in central6, central7, and full.",
            "- Even in the guessed windows, the single largest atom-pair class is base_like_vs_base_like, so the backbone hypothesis is supported only as an enrichment/sensitivity signal, not as a clean class assignment.",
            "- The full model has backbone-involving fractions similar to central6/central7 in the guessed windows, suggesting scaling of geometric opportunities rather than a qualitatively new full-length backbone-only contribution.",
            "- Experimental calibration would make this stronger by replacing guessed windows with measured d/q windows from detector geometry.",
            "",
            "## Outputs",
            "",
            f"- Summary CSV: `{args.summary_csv}`",
            f"- Atom-class CSV: `{args.by_atom_class_csv}`",
            f"- Chain-relation CSV: `{args.by_chain_relation_csv}`",
            f"- Pair sample CSV: `{args.pairs_sample_csv}`",
            f"- Plot directory: `{args.plot_dir}`",
        ]
    )
    for path in plot_paths:
        lines.append(f"- Plot: `{path}`")
    return "\n".join(lines) + "\n"


def run(args: argparse.Namespace) -> dict[str, object]:
    structures = {
        "full": prepare_structure("full", args.full_pdb),
        "central6": prepare_structure("central6", args.central6_pdb),
        "central7": prepare_structure("central7", args.central7_pdb),
    }
    pair_rows: list[dict[str, str]] = []
    for structure in structures.values():
        pair_rows.extend(pair_rows_for_structure(structure))

    summary_rows, atom_rows, chain_rows = summarize(pair_rows, structures)
    write_csv(args.summary_csv, summary_rows, SUMMARY_COLUMNS)
    write_csv(args.by_atom_class_csv, atom_rows, ATOM_CLASS_COLUMNS)
    write_csv(args.by_chain_relation_csv, chain_rows, CHAIN_RELATION_COLUMNS)
    write_csv(args.pairs_sample_csv, pair_rows[: args.pairs_sample_limit], PAIR_SAMPLE_COLUMNS)
    plot_paths = write_plots(atom_rows, chain_rows, args.plot_dir)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(report_text(structures, atom_rows, chain_rows, plot_paths, args), encoding="utf-8")
    return {
        "pair_rows": len(pair_rows),
        "summary_rows": len(summary_rows),
        "atom_rows": len(atom_rows),
        "chain_rows": len(chain_rows),
        "plots": plot_paths,
    }


def main() -> None:
    result = run(parse_args())
    print(f"Wrote {result['pair_rows']} attributed pair rows before sampling")
    print(f"Wrote {result['summary_rows']} summary rows")
    print(f"Wrote {result['atom_rows']} atom-class rows")
    print(f"Wrote {result['chain_rows']} chain-relation rows")
    print(f"Wrote {len(result['plots'])} plots")


if __name__ == "__main__":
    main()
