#!/usr/bin/env python3
"""Analyze Asem-corrected compact twist diffraction radial profiles."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


TWISTS = [24, 26, 28, 30, 32, 34, 36]
FEATURE_WINDOWS = [
    ("d_3p4", "3.4 A", 3.3, 3.5),
    ("d_3p0", "3.0 A", 2.9, 3.1),
    ("d_4p5_5p0", "4.5-5.0 A", 4.5, 5.0),
    ("d_4p1", "4.1 A", 4.0, 4.2),
    ("d_5p5", "5.5 A", 5.4, 5.6),
    ("d_7p0", "7.0 A", 6.8, 7.2),
    ("d_8p4", "8.4 A", 8.2, 8.6),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--profile-dir",
        type=Path,
        default=Path("outputs/compact_twist_diffraction_asem_corrected/radial_profiles"),
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("outputs/compact_twist_diffraction_asem_corrected"),
    )
    parser.add_argument(
        "--summary-csv",
        type=Path,
        default=Path("outputs/metrics/compact_twist_feature_summary_asem_corrected.csv"),
    )
    parser.add_argument(
        "--ranking-csv",
        type=Path,
        default=Path("outputs/metrics/compact_twist_feature_rankings_asem_corrected.csv"),
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=Path("outputs/reports/compact_twist_diffraction_asem_corrected_report.md"),
    )
    parser.add_argument(
        "--old-summary",
        type=Path,
        default=Path(
            r"C:\Users\hpage3\OneDrive - Georgia Institute of Technology\Documents\GitHub"
            r"\fiber-diffraction-powder\outputs\metrics\compact_twist_feature_summary.csv"
        ),
    )
    parser.add_argument(
        "--old-ranking",
        type=Path,
        default=Path(
            r"C:\Users\hpage3\OneDrive - Georgia Institute of Technology\Documents\GitHub"
            r"\fiber-diffraction-powder\outputs\metrics\compact_twist_feature_rankings.csv"
        ),
    )
    parser.add_argument("--normalize-q-min", type=float, default=0.4)
    return parser.parse_args()


def safe_float(value: object) -> float | None:
    try:
        parsed = float(str(value))
    except (TypeError, ValueError):
        return None
    return parsed


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_profile(path: Path) -> list[dict[str, float]]:
    rows = []
    for row in read_csv(path):
        q_value = safe_float(row.get("q_Ainv") or row.get("q_center_inv_angstrom"))
        d_value = safe_float(row.get("d_A") or row.get("d_center_angstrom"))
        intensity = safe_float(row.get("mean_intensity") or row.get("intensity_mean") or row.get("intensity"))
        pixel_count = safe_float(row.get("pixel_count") or row.get("sample_count") or "1")
        if q_value is None or d_value is None or intensity is None or pixel_count is None:
            continue
        rows.append(
            {
                "q_Ainv": q_value,
                "d_A": d_value,
                "mean_intensity": intensity,
                "pixel_count": pixel_count,
            }
        )
    return rows


def normalize_profile(rows: list[dict[str, float]], normalize_q_min: float) -> list[dict[str, float]]:
    candidates = [
        row["mean_intensity"]
        for row in rows
        if row["q_Ainv"] >= normalize_q_min and row["pixel_count"] > 0
    ]
    max_intensity = max(candidates) if candidates else 0.0
    return [
        {
            **row,
            "normalized_intensity": row["mean_intensity"] / max_intensity if max_intensity else 0.0,
        }
        for row in rows
    ]


def profile_path(profile_dir: Path, twist: int) -> Path:
    return profile_dir / f"compact_hexaplex_twist_{twist}_radial.csv"


def load_profiles(profile_dir: Path, normalize_q_min: float) -> dict[int, list[dict[str, float]]]:
    profiles = {}
    for twist in TWISTS:
        path = profile_path(profile_dir, twist)
        if not path.exists():
            raise FileNotFoundError(f"missing radial profile: {path}")
        profiles[twist] = normalize_profile(read_profile(path), normalize_q_min)
    return profiles


def summarize_features(profiles: dict[int, list[dict[str, float]]]) -> list[dict[str, object]]:
    rows = []
    for twist in TWISTS:
        summary: dict[str, object] = {"twist": f"{twist:.1f}"}
        profile = profiles[twist]
        for key, _label, d_min, d_max in FEATURE_WINDOWS:
            points = [row for row in profile if d_min <= row["d_A"] <= d_max]
            if points:
                best = max(points, key=lambda row: row["normalized_intensity"])
                summary[f"{key}_max_norm_intensity"] = f"{best['normalized_intensity']:.12g}"
                summary[f"{key}_q_Ainv"] = f"{best['q_Ainv']:.12g}"
                summary[f"{key}_d_A"] = f"{best['d_A']:.12g}"
                summary[f"{key}_mean_intensity"] = f"{best['mean_intensity']:.12g}"
            else:
                summary[f"{key}_max_norm_intensity"] = ""
                summary[f"{key}_q_Ainv"] = ""
                summary[f"{key}_d_A"] = ""
                summary[f"{key}_mean_intensity"] = ""
        rows.append(summary)
    return rows


def summary_fieldnames() -> list[str]:
    fieldnames = ["twist"]
    for key, _label, _d_min, _d_max in FEATURE_WINDOWS:
        fieldnames.extend(
            [
                f"{key}_max_norm_intensity",
                f"{key}_q_Ainv",
                f"{key}_d_A",
                f"{key}_mean_intensity",
            ]
        )
    return fieldnames


def build_rankings(summary_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    rankings = []
    for key, label, _d_min, _d_max in FEATURE_WINDOWS:
        values = []
        for row in summary_rows:
            value = safe_float(row.get(f"{key}_max_norm_intensity"))
            if value is None:
                continue
            values.append(
                {
                    "feature_window": label,
                    "twist_deg": safe_float(row["twist"]),
                    "normalized_intensity": value,
                    "d_A_at_max": safe_float(row.get(f"{key}_d_A")),
                    "q_Ainv_at_max": safe_float(row.get(f"{key}_q_Ainv")),
                    "mean_intensity_at_max": safe_float(row.get(f"{key}_mean_intensity")),
                }
            )
        values.sort(key=lambda row: row["normalized_intensity"], reverse=True)
        for rank, row in enumerate(values, start=1):
            rankings.append(
                {
                    "feature_window": row["feature_window"],
                    "rank": rank,
                    "twist_deg": f"{row['twist_deg']:g}",
                    "normalized_intensity": f"{row['normalized_intensity']:.12g}",
                    "d_A_at_max": f"{row['d_A_at_max']:.12g}" if row["d_A_at_max"] is not None else "",
                    "q_Ainv_at_max": f"{row['q_Ainv_at_max']:.12g}" if row["q_Ainv_at_max"] is not None else "",
                    "mean_intensity_at_max": (
                        f"{row['mean_intensity_at_max']:.12g}" if row["mean_intensity_at_max"] is not None else ""
                    ),
                }
            )
    return rankings


def top_by_feature(rankings: list[dict[str, object]]) -> dict[str, dict[str, object]]:
    return {str(row["feature_window"]): row for row in rankings if str(row["rank"]) == "1" or row["rank"] == 1}


def old_top_by_feature(path: Path) -> dict[str, dict[str, str]]:
    rows = read_csv(path)
    return {row["feature_window"]: row for row in rows if row.get("rank") == "1"}


def plot_profiles_q(path: Path, profiles: dict[int, list[dict[str, float]]]) -> None:
    fig, ax = plt.subplots(figsize=(9, 5.5), dpi=160)
    for twist in TWISTS:
        rows = [row for row in profiles[twist] if row["q_Ainv"] >= 0.15]
        ax.plot(
            [row["q_Ainv"] for row in rows],
            [row["normalized_intensity"] for row in rows],
            linewidth=1.4,
            label=f"{twist} deg",
        )
    ax.set_xlabel("q (A^-1)")
    ax.set_ylabel("normalized mean intensity")
    ax.set_title("Asem-corrected compact twist radial profiles")
    ax.grid(True, alpha=0.25, linewidth=0.6)
    ax.legend(ncol=2)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_profiles_d(path: Path, profiles: dict[int, list[dict[str, float]]]) -> None:
    fig, ax = plt.subplots(figsize=(9, 5.5), dpi=160)
    for twist in TWISTS:
        rows = sorted(
            [row for row in profiles[twist] if 2.5 <= row["d_A"] <= 9.0],
            key=lambda row: row["d_A"],
        )
        ax.plot(
            [row["d_A"] for row in rows],
            [row["normalized_intensity"] for row in rows],
            linewidth=1.4,
            label=f"{twist} deg",
        )
    ax.invert_xaxis()
    ax.set_xlabel("d-spacing (A)")
    ax.set_ylabel("normalized mean intensity")
    ax.set_title("Asem-corrected compact twist radial profiles by d-spacing")
    ax.grid(True, alpha=0.25, linewidth=0.6)
    ax.legend(ncol=2)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_feature_response(path: Path, summary_rows: list[dict[str, object]]) -> None:
    fig, ax = plt.subplots(figsize=(9, 5.5), dpi=160)
    twists = [safe_float(row["twist"]) for row in summary_rows]
    for key, label, _d_min, _d_max in FEATURE_WINDOWS:
        values = [safe_float(row.get(f"{key}_max_norm_intensity")) for row in summary_rows]
        ax.plot(twists, values, marker="o", linewidth=1.4, label=label)
    ax.set_xlabel("target twist angle, degrees")
    ax.set_ylabel("window max normalized intensity")
    ax.set_title("Asem-corrected feature response versus compact twist")
    ax.grid(True, alpha=0.25, linewidth=0.6)
    ax.legend(ncol=2)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def markdown_table(rows: list[dict[str, object]], columns: list[str]) -> str:
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join("---" for _ in columns) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(str(row.get(column, "")) for column in columns) + " |")
    return "\n".join(lines)


def ranking_sentence(rankings: list[dict[str, object]], feature: str) -> str:
    rows = [row for row in rankings if row["feature_window"] == feature]
    return ", ".join(
        f"{row['twist_deg']} deg ({float(str(row['normalized_intensity'])):.6f})"
        for row in rows
    )


def write_report(
    path: Path,
    summary_rows: list[dict[str, object]],
    rankings: list[dict[str, object]],
    old_rankings: dict[str, dict[str, str]],
    output_root: Path,
    normalize_q_min: float,
) -> None:
    top = top_by_feature(rankings)
    comparison_rows = []
    changed = []
    for _key, label, _d_min, _d_max in FEATURE_WINDOWS:
        new_top = top.get(label, {})
        old_top = old_rankings.get(label, {})
        old_twist = old_top.get("twist_deg", "")
        new_twist = str(new_top.get("twist_deg", ""))
        changed_flag = "yes" if old_twist and new_twist and old_twist != new_twist else "no"
        if changed_flag == "yes":
            changed.append(label)
        comparison_rows.append(
            {
                "feature_window": label,
                "pre_Asem_fix_top_twist": old_twist,
                "corrected_top_twist": new_twist,
                "changed": changed_flag,
            }
        )
    primary_30 = all(str(top[label]["twist_deg"]) == "30" for label in ("3.4 A", "3.0 A", "4.5-5.0 A"))
    exception_8p4 = str(top["8.4 A"]["twist_deg"]) != "30"
    lines = [
        "# Compact Twist Diffraction Report: Asem-Corrected",
        "",
        "## Inputs",
        "",
        "Compact twist PDBs were copied into `inputs/compact_twist_variants` and the original source PDBs were not modified.",
        "",
        "## Corrected Workflow",
        "",
        "- Diffraction code path: `../fiber-diffraction/scripts.py` and `../fiber-diffraction/orientation_average.py`.",
        "- Correction: azimuthal rotations are applied to the same tilted coordinate frame rather than to already-rotated coordinates.",
        "- Vectorization: `numexpr` is installed in `.venv` and `scripts._HAVE_NUMEXPR` is true.",
        "- PDB conversion: heavy atoms only with exact duplicate removal.",
        "- Diffraction settings: grid size 41, grid limit 100 mm, theta-count 4, phi-count 8, psi-count 3, theta max 180 degrees, workers 1.",
        f"- Radial normalization: maximum mean intensity at q >= {normalize_q_min:g} A^-1 per twist.",
        "",
        "## Outputs",
        "",
        f"- Corrected output root: `{output_root}`",
        "- `outputs/metrics/compact_twist_feature_summary_asem_corrected.csv`",
        "- `outputs/metrics/compact_twist_feature_rankings_asem_corrected.csv`",
        "- `outputs/reports/compact_twist_diffraction_asem_corrected_report.md`",
        "",
        "## Corrected Feature Summary",
        "",
        markdown_table(
            summary_rows,
            [
                "twist",
                "d_3p4_max_norm_intensity",
                "d_3p0_max_norm_intensity",
                "d_4p5_5p0_max_norm_intensity",
                "d_4p1_max_norm_intensity",
                "d_5p5_max_norm_intensity",
                "d_7p0_max_norm_intensity",
                "d_8p4_max_norm_intensity",
            ],
        ),
        "",
        "## Corrected Top Feature Responses",
        "",
        markdown_table(
            [top[label] for _key, label, _d_min, _d_max in FEATURE_WINDOWS],
            ["feature_window", "rank", "twist_deg", "normalized_intensity", "d_A_at_max"],
        ),
        "",
        "## Pre/Post Asem-Fix Comparison",
        "",
        markdown_table(comparison_rows, ["feature_window", "pre_Asem_fix_top_twist", "corrected_top_twist", "changed"]),
        "",
        "## Interpretation",
        "",
        f"- Strongest twist per feature window changed after the correction: {'yes: ' + ', '.join(changed) if changed else 'no'}.",
        f"- 30 degrees remains strongest in the primary windows (3.4 A, 3.0 A, 4.5-5.0 A): {'yes' if primary_30 else 'no'}.",
        f"- The 8.4 A exception remains: {'yes' if exception_8p4 else 'no'}.",
        "- New corrected full rankings:",
    ]
    for _key, label, _d_min, _d_max in FEATURE_WINDOWS:
        lines.append(f"  - {label}: {ranking_sentence(rankings, label)}.")
    lines.extend(
        [
            "",
            "## Cautions",
            "",
            "- These are comparative powder-style radial summaries, not direct fits to oriented/fiber experimental images.",
            "- The corrected run uses the old pre-fix Track B report's 41 x 41 / 96-orientation / 160-bin settings for direct feature-window comparison.",
            "- The input structures are compact rigid twist variants, not pNAB-generated conformer ensembles.",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    profiles = load_profiles(args.profile_dir, args.normalize_q_min)
    summary_rows = summarize_features(profiles)
    rankings = build_rankings(summary_rows)
    plots_dir = args.output_root / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    plot_profiles_q(plots_dir / "compact_twist_radial_profiles_q_asem_corrected.png", profiles)
    plot_profiles_d(plots_dir / "compact_twist_radial_profiles_d_asem_corrected.png", profiles)
    plot_feature_response(plots_dir / "compact_twist_feature_response_asem_corrected.png", summary_rows)
    write_csv(args.summary_csv, summary_rows, summary_fieldnames())
    write_csv(
        args.ranking_csv,
        rankings,
        ["feature_window", "rank", "twist_deg", "normalized_intensity", "d_A_at_max", "q_Ainv_at_max", "mean_intensity_at_max"],
    )
    write_report(
        args.report,
        summary_rows,
        rankings,
        old_top_by_feature(args.old_ranking),
        args.output_root,
        args.normalize_q_min,
    )
    print(f"Wrote {args.summary_csv}")
    print(f"Wrote {args.ranking_csv}")
    print(f"Wrote {args.report}")
    print(f"Wrote {plots_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
