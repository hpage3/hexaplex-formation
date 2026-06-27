"""Plot diagnostic/paper-draft figures for the final 20 Asem matched scan."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


GROUP_ORDER = [
    ("rise3p38", "sidechain"),
    ("rise3p38", "glu_ablated"),
    ("rise3p40", "sidechain"),
    ("rise3p40", "glu_ablated"),
]
GROUP_STYLE = {
    ("rise3p38", "sidechain"): {"label": "rise 3.38 sidechain", "color": "#1f77b4", "marker": "o"},
    ("rise3p38", "glu_ablated"): {"label": "rise 3.38 GLU-ablated", "color": "#ff7f0e", "marker": "s"},
    ("rise3p40", "sidechain"): {"label": "rise 3.40 sidechain", "color": "#2ca02c", "marker": "^"},
    ("rise3p40", "glu_ablated"): {"label": "rise 3.40 GLU-ablated", "color": "#d62728", "marker": "D"},
}
BANDS = ["base", "A", "B", "C", "D"]
TARGET_D = {"base": 3.38, "A": 3.8, "B": 4.4, "C": 5.65, "D": 7.3}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def by_group(rows: list[dict[str, str]], rise_label: str, variant: str) -> list[dict[str, str]]:
    selected = [
        row
        for row in rows
        if row["rise_label"] == rise_label and row["sidechain_variant"] == variant
    ]
    return sorted(selected, key=lambda row: float(row["twist_deg"]))


def plot_metric_by_twist(
    rows: list[dict[str, str]],
    metric: str,
    ylabel: str,
    title: str,
    output: Path,
    target_line: float | None = None,
    note: str | None = None,
) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for group in GROUP_ORDER:
        group_rows = by_group(rows, *group)
        style = GROUP_STYLE[group]
        ax.plot(
            [float(row["twist_deg"]) for row in group_rows],
            [float(row[metric]) for row in group_rows],
            color=style["color"],
            marker=style["marker"],
            linewidth=1.6,
            markersize=5,
            label=style["label"],
        )
    if target_line is not None:
        ax.axhline(target_line, color="#555555", linestyle="--", linewidth=1.1, label=f"target {target_line:g} A")
    ax.set_xlabel("Twist (degrees)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xticks([28, 29, 30, 31, 32])
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8, frameon=False)
    if note:
        fig.text(0.01, 0.01, note, fontsize=8, color="#444444")
        fig.tight_layout(rect=[0, 0.05, 1, 1])
    else:
        fig.tight_layout()
    fig.savefig(output, dpi=200)
    plt.close(fig)


def plot_band_positions(rows: list[dict[str, str]], output: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.4), sharex=True)
    axes_flat = axes.ravel()
    x_values = list(range(len(BANDS)))
    target_values = [TARGET_D[band] for band in BANDS]
    for ax, group in zip(axes_flat, GROUP_ORDER):
        group_rows = by_group(rows, *group)
        style = GROUP_STYLE[group]
        for row in group_rows:
            y_values = [float(row[f"observed_{band}_d_A"]) for band in BANDS]
            ax.plot(x_values, y_values, marker="o", linewidth=1.2, markersize=4, label=f"{float(row['twist_deg']):.0f} deg")
        ax.scatter(x_values, target_values, marker="x", s=42, color="#111111", label="target")
        ax.set_title(style["label"])
        ax.set_xticks(x_values)
        ax.set_xticklabels(BANDS)
        ax.set_ylabel("Observed d-spacing (A)")
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=7, frameon=False, ncol=2)
    fig.suptitle("Final 20 equalized set: observed band positions by twist", y=0.995)
    fig.tight_layout()
    fig.savefig(output, dpi=200)
    plt.close(fig)


def write_notes(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "# Figure Notes",
                "",
                "- The final 20 set is selected from 346 scored candidates.",
                "- Selection is one best-scoring candidate per rise x variant x twist.",
                "- Scoring uses the corrected HXC590 S1 target file via the rounded twist/rise scan target table.",
                "- Ranking is primarily driven by helical/backbone-associated bands, not base-stacking, because base_rmsd is nearly constant across candidates.",
                "- These are screening/model-family comparisons, not final structural assignments.",
                "",
            ]
        ),
        encoding="utf-8",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("outputs/asem_matched_twist_scan_28_32/final_20_equalized_best_candidates.csv"),
    )
    parser.add_argument(
        "--summary-by-twist",
        type=Path,
        default=Path("outputs/asem_matched_twist_scan_28_32/final_20_score_summary_by_twist.csv"),
    )
    parser.add_argument(
        "--summary-by-variant",
        type=Path,
        default=Path("outputs/asem_matched_twist_scan_28_32/final_20_score_summary_by_variant.csv"),
    )
    parser.add_argument(
        "--summary-by-rise-variant",
        type=Path,
        default=Path("outputs/asem_matched_twist_scan_28_32/final_20_score_summary_by_rise_variant.csv"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs/asem_matched_twist_scan_28_32/figures"),
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    rows = read_csv(args.input)
    if len(rows) != 20:
        raise ValueError(f"Expected 20 final rows, found {len(rows)}")
    # Read summary inputs to make the dependency explicit and fail fast if missing.
    for path in [args.summary_by_twist, args.summary_by_variant, args.summary_by_rise_variant]:
        if not read_csv(path):
            raise ValueError(f"No rows found in {path}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "combined": args.output_dir / "final_20_combined_rmsd_by_twist.png",
        "base": args.output_dir / "final_20_base_peak_by_twist.png",
        "helical": args.output_dir / "final_20_helical_rmsd_by_twist.png",
        "bands": args.output_dir / "final_20_band_positions_2x2.png",
        "notes": args.output_dir / "figure_notes.md",
    }

    plot_metric_by_twist(
        rows,
        metric="combined_rmsd",
        ylabel="combined RMSD (A)",
        title="Equalized final 20 set: combined RMSD by twist",
        output=outputs["combined"],
    )
    plot_metric_by_twist(
        rows,
        metric="observed_base_d_A",
        ylabel="observed base-stacking d-spacing (A)",
        title="Equalized final 20 set: base-stacking peak by twist",
        output=outputs["base"],
        target_line=3.38,
        note="Diagnostic note: base-stacking peak is weakly discriminating in this candidate set.",
    )
    plot_metric_by_twist(
        rows,
        metric="helical_rmsd",
        ylabel="helical/backbone-associated RMSD (A)",
        title="Equalized final 20 set: helical-band RMSD by twist",
        output=outputs["helical"],
    )
    plot_band_positions(rows, outputs["bands"])
    write_notes(outputs["notes"])

    for path in outputs.values():
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
