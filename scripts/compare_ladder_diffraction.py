#!/usr/bin/env python3
"""Compare ladder Debye radial-window scores."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


WINDOWS = ["d_8p4", "d_5p5_6p0", "d_4p5", "d_4p1", "d_3p4", "d_3p0"]
OUTPUT_COLUMNS = [
    "model_name",
    "included_blocks",
    "includes_hexads",
    "atom_mode",
    "atom_count",
    "residue_count",
]
for window in WINDOWS:
    OUTPUT_COLUMNS.extend([f"{window}_mean", f"{window}_fraction"])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--window-dir", type=Path, default=Path("outputs/metrics/ladder_diffraction/window_scores"))
    parser.add_argument("--ladder-summary", type=Path, default=Path("outputs/metrics/intermediate_ladder_summary.csv"))
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=Path("outputs/metrics/ladder_diffraction/ladder_diffraction_comparison.csv"),
    )
    parser.add_argument("--out-md", type=Path, default=Path("outputs/reports/ladder_diffraction_report.md"))
    parser.add_argument(
        "--ladder-metrics",
        type=Path,
        default=Path("outputs/metrics/intermediate_ladder/intermediate_ladder_comparison.csv"),
    )
    parser.add_argument("--profile-dir", type=Path, default=Path("outputs/metrics/ladder_diffraction/profiles"))
    return parser.parse_args()


def _rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def model_name_from_scores(path: Path) -> str:
    suffix = "_window_scores"
    stem = path.stem
    return stem[: -len(suffix)] if stem.endswith(suffix) else stem


def _score_by_window(path: Path) -> dict[str, dict[str, str]]:
    return {row["window_name"]: row for row in _rows(path)}


def aggregate_scores(window_dir: Path, ladder_summary: Path) -> list[dict[str, str]]:
    summary_by_model = {row["model_name"]: row for row in _rows(ladder_summary)}
    rows: list[dict[str, str]] = []
    for path in sorted(window_dir.glob("*_window_scores.csv")):
        model_name = model_name_from_scores(path)
        summary = summary_by_model.get(model_name, {})
        score_by_window = _score_by_window(path)
        row = {
            "model_name": model_name,
            "included_blocks": summary.get("included_blocks", ""),
            "includes_hexads": summary.get("includes_hexads", ""),
            "atom_mode": summary.get("atom_mode", ""),
            "atom_count": summary.get("atom_count", ""),
            "residue_count": summary.get("residue_count", ""),
        }
        for window in WINDOWS:
            score = score_by_window.get(window, {})
            row[f"{window}_mean"] = score.get("mean_intensity", "")
            row[f"{window}_fraction"] = score.get("integrated_intensity_fraction_of_total", "")
        rows.append(row)
    return rows


def write_csv(rows: list[dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _as_float(row: dict[str, str], column: str) -> float:
    return float(row.get(column, "") or 0.0)


def _block_count(row: dict[str, str]) -> int:
    included = row.get("included_blocks", "")
    return len([part for part in included.split(",") if part])


def _scaffold_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return sorted(
        [
            row
            for row in rows
            if row["model_name"].startswith("scaffold_blocks_")
            and row.get("includes_hexads") == "no"
            and row.get("atom_mode") == "heavy_deduped"
        ],
        key=_block_count,
    )


def _hexad_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return sorted(
        [
            row
            for row in rows
            if row["model_name"].startswith("hexads_plus_scaffold_blocks_")
            and row.get("includes_hexads") == "yes"
            and row.get("atom_mode") == "heavy_deduped"
        ],
        key=_block_count,
    )


def pearson(xs: list[float], ys: list[float]) -> float | None:
    if len(xs) < 3 or len(xs) != len(ys):
        return None
    x_mean = sum(xs) / len(xs)
    y_mean = sum(ys) / len(ys)
    numerator = sum((x - x_mean) * (y - y_mean) for x, y in zip(xs, ys))
    x_var = sum((x - x_mean) ** 2 for x in xs)
    y_var = sum((y - y_mean) ** 2 for y in ys)
    if x_var == 0 or y_var == 0:
        return None
    return numerator / math.sqrt(x_var * y_var)


def correlation_notes(rows: list[dict[str, str]], ladder_metrics_path: Path) -> list[str]:
    metric_rows = {row["model_name"]: row for row in _rows(ladder_metrics_path)}
    paired = [(row, metric_rows[row["model_name"]]) for row in rows if row["model_name"] in metric_rows]
    notes: list[str] = []
    for metric_column, label in [
        ("motif_GLU_GLU", "GLU-GLU motif count"),
        ("contact_count_4p5A", "real-space contact count"),
    ]:
        xs = [_as_float(row, "d_4p5_fraction") for row, _metric in paired]
        ys = [float(_metric.get(metric_column, "") or 0.0) for _row, _metric in paired]
        corr = pearson(xs, ys)
        if corr is not None:
            notes.append(f"Pearson correlation between d_4p5_fraction and {label}: {corr:.3f}.")
    xs = [_as_float(row, "d_3p4_fraction") for row, _metric in paired]
    ys = [1.0 if row.get("includes_hexads") == "yes" else 0.0 for row, _metric in paired]
    corr = pearson(xs, ys)
    if corr is not None:
        notes.append(f"Pearson correlation between d_3p4_fraction and hexad inclusion flag: {corr:.3f}.")
    return notes


def _metadata_notes(profile_dir: Path) -> list[str]:
    metadata_paths = sorted(profile_dir.glob("*_debye_profile.metadata.json"))
    metadata_items = []
    sampled = []
    first_n_truncated = []
    for path in metadata_paths:
        metadata = json.loads(path.read_text(encoding="utf-8"))
        metadata_items.append(metadata)
        if metadata.get("sampling_used") or metadata.get("sample_atoms") is not None:
            sampled.append(metadata)
        if metadata.get("first_n_truncated") or metadata.get("truncated"):
            first_n_truncated.append(metadata)
    if not metadata_items:
        return ["No profile metadata JSON files were found; runtime atom usage and method could not be verified."]

    methods = sorted({str(item.get("method", "unknown")) for item in metadata_items})
    bin_widths = sorted({str(item.get("distance_bin_width", "unknown")) for item in metadata_items})
    full_count = sum(
        1
        for item in metadata_items
        if not item.get("sampling_used")
        and item.get("sample_atoms") is None
        and not item.get("first_n_truncated")
        and not item.get("truncated")
    )
    notes = [
        f"Profile methods: {', '.join(methods)}.",
        f"Distance bin widths: {', '.join(bin_widths)}.",
        f"{full_count} of {len(metadata_items)} profile(s) used full post-filter atom sets.",
    ]
    if sampled:
        sample_values = sorted({str(item.get("sample_atoms")) for item in sampled})
        notes.append(
            "At least one profile used deterministic stratified sampling "
            f"(SAMPLE_ATOMS values: {', '.join(sample_values)}). Sampling is less biased than first-N truncation "
            "but remains a runtime-limited approximation."
        )
    if first_n_truncated:
        max_values = sorted({str(item.get("max_atoms")) for item in first_n_truncated})
        notes.append(
            "STRONG WARNING: at least one profile used legacy first-N atom truncation "
            f"(MAX_ATOMS values: {', '.join(max_values)}). For hexads-plus-scaffold structures, this can "
            "preferentially sample hexads before scaffold atoms and should not be used for scientific interpretation."
        )
    return notes


def report_notes(rows: list[dict[str, str]], ladder_metrics_path: Path) -> list[str]:
    notes: list[str] = []
    scaffold_rows = _scaffold_rows(rows)
    hexad_rows = _hexad_rows(rows)
    block1 = next((row for row in scaffold_rows if _block_count(row) == 1), None)
    if block1:
        notes.append(
            f"d_4p5_fraction is present in scaffold block 1 ({_as_float(block1, 'd_4p5_fraction'):.6f}), "
            "so this approximate profile does not require complete multi-block assembly to show a 4.5 A-window score."
        )
    if len(scaffold_rows) > 1:
        first = _as_float(scaffold_rows[0], "d_4p5_fraction")
        last = _as_float(scaffold_rows[-1], "d_4p5_fraction")
        direction = "increases" if last > first else "does not increase"
        notes.append(f"Across scaffold-only ladder rows, d_4p5_fraction {direction} from block 1 to all blocks.")
    paired_hexads = {row["included_blocks"]: row for row in hexad_rows}
    changes = []
    for scaffold in scaffold_rows:
        hexad = paired_hexads.get(scaffold["included_blocks"])
        if hexad:
            changes.append(_as_float(hexad, "d_3p4_fraction") - _as_float(scaffold, "d_3p4_fraction"))
    if changes:
        mean_change = sum(changes) / len(changes)
        notes.append(f"Adding hexads changes d_3p4_fraction by a mean paired difference of {mean_change:.6f}.")
    notes.extend(correlation_notes(rows, ladder_metrics_path))
    notes.append(
        "The d_4p5 window should be interpreted cautiously as a reciprocal-space-like comparative score, "
        "not as a literal 4.5 A contact count."
    )
    return notes


def _markdown_table(rows: list[dict[str, str]], columns: list[str]) -> str:
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(row.get(column, "") for column in columns) + " |")
    return "\n".join(lines)


def write_markdown_report(rows: list[dict[str, str]], path: Path, ladder_metrics_path: Path, profile_dir: Path) -> None:
    columns = ["model_name", "included_blocks", "includes_hexads"]
    for window in WINDOWS:
        columns.append(f"{window}_fraction")
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Ladder Diffraction Approximation Report",
        "",
        "## Scientific cautions",
        "",
        "This report uses a simplified Debye-style isotropic scattering approximation for comparative radial profiles. It is not a replacement for full fiber-diffraction simulations and does not model orientation, lattice effects, form factors, hydration, or experimental instrument response.",
        "",
        "The d ~= 4.5 A feature remains a reciprocal-space scaffold signature in the working hypothesis, not a literal 4.5 A atom-contact assignment. Component structure intensities are comparative controls, not additive decompositions of the full intensity.",
        "",
        "The candidate block map remains unvalidated against PyMOL colored strand paths, and ladder models do not prove temporal assembly order.",
        "",
        "## Runtime metadata",
        "",
    ]
    lines.extend(f"- {note}" for note in _metadata_notes(profile_dir))
    lines.extend(
        [
            "",
            "## Window fraction table",
            "",
            _markdown_table(rows, columns),
            "",
            "## Automatically generated notes",
            "",
        ]
    )
    lines.extend(f"- {note}" for note in report_notes(rows, ladder_metrics_path))
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    rows = aggregate_scores(args.window_dir, args.ladder_summary)
    write_csv(rows, args.out_csv)
    write_markdown_report(rows, args.out_md, args.ladder_metrics, args.profile_dir)
    print(f"Wrote {args.out_csv} with {len(rows)} model row(s)")
    print(f"Wrote {args.out_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
