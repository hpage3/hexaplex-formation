#!/usr/bin/env python3
"""Compare Stage 2 with-GLU scores against the matched no-GLU ablation."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


COMPARISON_FIELDS = [
    "model_id", "twist_deg", "rise_A",
    "with_glu_base_rmsd", "no_glu_base_rmsd", "delta_base_rmsd",
    "with_glu_helical_rmsd", "no_glu_helical_rmsd", "delta_helical_rmsd",
    "with_glu_combined_rmsd", "no_glu_combined_rmsd", "delta_combined_rmsd",
    "with_glu_rank", "no_glu_rank", "rank_delta",
    "with_glu_observed_A_d_A", "no_glu_observed_A_d_A", "delta_observed_A_d_A",
    "with_glu_observed_B_d_A", "no_glu_observed_B_d_A", "delta_observed_B_d_A",
]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def number(row: dict[str, str], key: str) -> float:
    return float(row[key])


def difference(after: float, before: float) -> str:
    return f"{after - before:.9f}"


def build_comparison(with_rows: list[dict[str, str]], no_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    with_by_id = {row["model_id"]: row for row in with_rows}
    no_by_id = {row["model_id"]: row for row in no_rows}
    if set(with_by_id) != set(no_by_id):
        missing_no = sorted(set(with_by_id) - set(no_by_id))
        missing_with = sorted(set(no_by_id) - set(with_by_id))
        raise ValueError(f"Model sets differ; missing no-GLU={missing_no}, missing with-GLU={missing_with}")

    output = []
    for model_id in sorted(with_by_id, key=lambda key: (float(with_by_id[key]["twist_deg"]), key)):
        with_row = with_by_id[model_id]
        no_row = no_by_id[model_id]
        row = {"model_id": model_id, "twist_deg": with_row["twist_deg"], "rise_A": with_row["rise_A"]}
        for metric in ("base_rmsd", "helical_rmsd", "combined_rmsd"):
            before, after = number(with_row, metric), number(no_row, metric)
            row[f"with_glu_{metric}"] = f"{before:.9f}"
            row[f"no_glu_{metric}"] = f"{after:.9f}"
            row[f"delta_{metric}"] = difference(after, before)
        before_rank, after_rank = int(with_row["rank"]), int(no_row["rank"])
        row.update({
            "with_glu_rank": str(before_rank),
            "no_glu_rank": str(after_rank),
            "rank_delta": str(after_rank - before_rank),
        })
        for label in ("A", "B"):
            key = f"observed_{label}_d_A"
            before, after = number(with_row, key), number(no_row, key)
            row[f"with_glu_{key}"] = f"{before:.9f}"
            row[f"no_glu_{key}"] = f"{after:.9f}"
            row[f"delta_{key}"] = difference(after, before)
        output.append(row)
    return output


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=COMPARISON_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def means_by_twist(rows: list[dict[str, str]], key: str) -> list[tuple[float, float]]:
    grouped: dict[float, list[float]] = defaultdict(list)
    for row in rows:
        grouped[float(row["twist_deg"])].append(float(row[key]))
    return [(twist, sum(values) / len(values)) for twist, values in sorted(grouped.items())]


def paired_plot(path: Path, rows: list[dict[str, str]], metric: str, ylabel: str) -> None:
    with_points = means_by_twist(rows, f"with_glu_{metric}")
    no_points = means_by_twist(rows, f"no_glu_{metric}")
    fig, ax = plt.subplots(figsize=(7, 4.5), dpi=150)
    ax.plot(*zip(*with_points), marker="o", linewidth=1.5, label="with complete Glu")
    ax.plot(*zip(*no_points), marker="s", linewidth=1.5, label="without distal Glu branch")
    ax.set_xlabel("Twist (degrees)")
    ax.set_ylabel(ylabel)
    ax.set_title(f"3.40 A Stage 2 side-chain ablation: {ylabel}")
    ax.legend()
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def delta_plot(path: Path, rows: list[dict[str, str]]) -> None:
    points = means_by_twist(rows, "delta_combined_rmsd")
    fig, ax = plt.subplots(figsize=(7, 4.5), dpi=150)
    ax.axhline(0.0, color="black", linewidth=1)
    ax.plot(*zip(*points), marker="o", linewidth=1.5)
    ax.set_xlabel("Twist (degrees)")
    ax.set_ylabel("No-Glu minus with-Glu combined RMSD")
    ax.set_title("3.40 A Stage 2 side-chain ablation delta")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def load_profile(path: Path) -> tuple[list[float], list[float]]:
    rows = read_csv(path)
    points = [
        (float(row["d_A"]), float(row.get("intensity_mean") or row.get("mean_intensity") or row["intensity"]))
        for row in rows
        if 2.8 <= float(row["d_A"]) <= 10.0
    ]
    points.sort()
    maximum = max(intensity for _, intensity in points)
    return [d for d, _ in points], [intensity / maximum for _, intensity in points]


def representative_overlay(
    path: Path,
    rows: list[dict[str, str]],
    with_profile_dir: Path,
    no_profile_dir: Path,
) -> None:
    selected = {}
    for twist in (29.0, 30.0, 31.0):
        candidates = [row for row in rows if float(row["twist_deg"]) == twist]
        selected[twist] = min(candidates, key=lambda row: int(row["with_glu_rank"]))

    fig, axes = plt.subplots(3, 1, figsize=(8, 9), dpi=150, sharex=True)
    for ax, (twist, row) in zip(axes, selected.items()):
        model_id = row["model_id"]
        with_d, with_i = load_profile(with_profile_dir / f"{model_id}_radial.csv")
        no_d, no_i = load_profile(no_profile_dir / f"{model_id}_radial.csv")
        ax.plot(with_d, with_i, linewidth=1.5, label="with complete Glu")
        ax.plot(no_d, no_i, linewidth=1.5, label="without distal Glu branch")
        ax.set_ylabel("Normalized intensity")
        ax.set_title(f"{twist:.0f} degrees: {model_id}")
        ax.grid(alpha=0.2)
        ax.legend()
    axes[-1].set_xlabel("d-spacing (A)")
    axes[-1].invert_xaxis()
    fig.suptitle("3.40 A representative Stage 2 radial-profile ablations")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def write_summary(path: Path, rows: list[dict[str, str]]) -> None:
    mean = lambda key: sum(float(row[key]) for row in rows) / len(rows)
    improved = sum(float(row["delta_combined_rmsd"]) < 0 for row in rows)
    worsened = sum(float(row["delta_combined_rmsd"]) > 0 for row in rows)
    unchanged = len(rows) - improved - worsened
    lines = [
        "# Stage 2 Complete-Glu vs No-Glu Ablation",
        "",
        "This is a controlled 3.40 A sensitivity comparison using the same 30 Stage 2 candidates. "
        "It is not the future 3.38 A production comparison.",
        "",
        f"- Matched candidates: {len(rows)}",
        f"- Mean delta base RMSD (no-Glu minus with-Glu): {mean('delta_base_rmsd'):.6f} A",
        f"- Mean delta helical RMSD: {mean('delta_helical_rmsd'):.6f} A",
        f"- Mean delta combined RMSD: {mean('delta_combined_rmsd'):.6f} A",
        f"- Combined RMSD improved/worsened/unchanged after stripping: {improved}/{worsened}/{unchanged}",
        f"- Mean shift near 3.8 A (observed A): {mean('delta_observed_A_d_A'):.6f} A",
        f"- Mean shift near 4.5 A (observed B): {mean('delta_observed_B_d_A'):.6f} A",
        "",
        "Positive RMSD deltas mean the no-Glu ablation fits the target positions less well. "
        "These scores compare peak positions only; they do not establish a chemical mechanism or a final structural ranking.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--with-glu", type=Path, required=True)
    parser.add_argument("--no-glu", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--with-profile-dir", type=Path)
    parser.add_argument("--no-profile-dir", type=Path)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = build_comparison(read_csv(args.with_glu), read_csv(args.no_glu))
    write_csv(args.output_dir / "sidechain_ablation_with_vs_without_comparison.csv", rows)
    write_summary(args.output_dir / "sidechain_ablation_with_vs_without_summary.md", rows)
    paired_plot(
        args.output_dir / "sidechain_ablation_helical_rmsd_with_vs_without.png",
        rows, "helical_rmsd", "Helical RMSD (A)",
    )
    paired_plot(
        args.output_dir / "sidechain_ablation_combined_rmsd_with_vs_without.png",
        rows, "combined_rmsd", "Combined RMSD (A)",
    )
    delta_plot(args.output_dir / "sidechain_ablation_delta_combined_rmsd.png", rows)
    if args.with_profile_dir and args.no_profile_dir:
        representative_overlay(
            args.output_dir / "sidechain_ablation_representative_radial_overlays.png",
            rows,
            args.with_profile_dir,
            args.no_profile_dir,
        )
    print(f"Wrote paired comparison outputs for {len(rows)} candidates to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
