"""Score the matched Asem 28-32 twist/rise radial-profile screen.

This adapter joins the matched-coordinate, XYZ, and radial-profile manifests,
normalizes the radial CSVs into the existing twist/rise profile table shape,
then reuses the established peak-extraction and scoring helpers.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

from build_twist_rise_observed_peaks_table import (
    long_to_wide,
    output_fields as observed_output_fields,
    parse_target_order,
    validate_long_rows,
    write_csv as write_observed_csv,
)
from extract_twist_rise_profile_peaks import (
    extract_peaks,
    load_profiles,
    load_targets as load_extraction_targets,
    output_fields as extraction_output_fields,
    write_csv as write_extraction_csv,
)
from score_twist_rise_scan import (
    load_observed_peaks,
    load_targets as load_scoring_targets,
    merged_fields,
    read_csv,
    score_rows,
    write_csv,
)


EXPECTED_ROW_COUNT = 346
TARGET_ORDER = ["base", "A", "B", "C", "D"]


def resolve_path(path: Path, repo_root: Path) -> Path:
    if path.is_absolute():
        return path
    return repo_root / path


def require_unique(rows: list[dict[str, str]], key: str, path: Path) -> dict[str, dict[str, str]]:
    by_key: dict[str, dict[str, str]] = {}
    for index, row in enumerate(rows, start=2):
        value = row.get(key, "").strip()
        if not value:
            raise ValueError(f"{path} row {index} has empty {key}")
        if value in by_key:
            raise ValueError(f"Duplicate {key}={value!r} in {path}")
        by_key[value] = row
    return by_key


def radial_manifest_by_model(rows: list[dict[str, str]], path: Path) -> dict[str, dict[str, str]]:
    by_model: dict[str, dict[str, str]] = {}
    for index, row in enumerate(rows, start=2):
        model_id = row.get("model", "").strip() or row.get("model_id", "").strip()
        if not model_id:
            raise ValueError(f"{path} row {index} has empty model/model_id")
        if model_id in by_model:
            raise ValueError(f"Duplicate radial model {model_id!r} in {path}")
        by_model[model_id] = row
    return by_model


def join_manifests(
    primary_rows: list[dict[str, str]],
    xyz_by_model: dict[str, dict[str, str]],
    radial_by_model: dict[str, dict[str, str]],
    repo_root: Path,
) -> list[dict[str, str]]:
    joined: list[dict[str, str]] = []
    missing_xyz: list[str] = []
    missing_radial: list[str] = []
    missing_radial_files: list[str] = []

    for row in primary_rows:
        model_id = row["model_id"]
        xyz_row = xyz_by_model.get(model_id)
        radial_row = radial_by_model.get(model_id)

        if xyz_row is None:
            missing_xyz.append(model_id)
            continue
        if radial_row is None:
            missing_radial.append(model_id)
            continue

        radial_csv = radial_row.get("radial_csv", "").strip()
        if not radial_csv:
            missing_radial_files.append(f"{model_id}: blank radial_csv")
            continue
        radial_path = resolve_path(Path(radial_csv), repo_root)
        if not radial_path.exists():
            missing_radial_files.append(f"{model_id}: {radial_csv}")
            continue

        out = dict(row)
        out["xyz_file"] = xyz_row.get("xyz_file", "")
        out["xyz_status"] = xyz_row.get("xyz_status", "")
        out["radial_csv"] = radial_csv
        out["radial_status"] = radial_row.get("radial_status", radial_row.get("status", ""))
        out["diffraction_status"] = radial_row.get("status", out.get("diffraction_status", ""))
        out["diffraction_file"] = radial_row.get("npy_file", out.get("diffraction_file", ""))
        joined.append(out)

    problems = []
    if missing_xyz:
        problems.append(f"missing XYZ rows: {len(missing_xyz)}")
    if missing_radial:
        problems.append(f"missing radial rows: {len(missing_radial)}")
    if missing_radial_files:
        problems.append(f"missing radial CSV files: {len(missing_radial_files)}")
    if problems:
        detail = "; ".join(problems)
        examples = (missing_xyz + missing_radial + missing_radial_files)[:10]
        raise ValueError(f"Manifest join failed: {detail}. Examples: {examples}")

    if len(joined) != EXPECTED_ROW_COUNT:
        raise ValueError(f"Expected {EXPECTED_ROW_COUNT} joined rows, found {len(joined)}")

    return joined


def choose_column(fieldnames: list[str], candidates: list[str], path: Path) -> str:
    for name in candidates:
        if name in fieldnames:
            return name
    raise ValueError(f"{path} is missing all expected columns: {candidates}")


def adapt_radial_profile(path: Path, model_id: str) -> list[dict[str, str]]:
    rows = read_csv(path)
    if not rows:
        raise ValueError(f"No rows found in {path}")

    fieldnames = list(rows[0])
    d_column = choose_column(fieldnames, ["d_A", "d_center_angstrom", "d_angstrom"], path)
    intensity_column = choose_column(fieldnames, ["intensity_mean", "mean_intensity", "intensity"], path)
    q_column = next((name for name in ["q_Ainv", "q_center_inv_angstrom", "q_inv_angstrom"] if name in fieldnames), "")

    adapted: list[dict[str, str]] = []
    for index, row in enumerate(rows, start=2):
        d_text = row.get(d_column, "").strip()
        intensity_text = row.get(intensity_column, "").strip()
        if not d_text:
            raise ValueError(f"{path} row {index} has empty {d_column}")
        if not intensity_text:
            raise ValueError(f"{path} row {index} has empty {intensity_column}")

        out = {
            "model_id": model_id,
            "d_A": f"{float(d_text):.9f}",
            "intensity": f"{float(intensity_text):.9f}",
        }
        if q_column:
            q_text = row.get(q_column, "").strip()
            out["q"] = f"{float(q_text):.9f}" if q_text else ""
        else:
            out["q"] = ""
        adapted.append(out)

    return adapted


def combine_profiles(joined_rows: list[dict[str, str]], repo_root: Path) -> list[dict[str, str]]:
    combined: list[dict[str, str]] = []
    for row in joined_rows:
        radial_path = resolve_path(Path(row["radial_csv"]), repo_root)
        combined.extend(adapt_radial_profile(radial_path, row["model_id"]))
    return combined


def score_summary_rows(scored_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    fields = ["rise_label", "sidechain_variant", "twist_deg"]
    grouped: dict[tuple[str, str, str], list[dict[str, str]]] = {}
    for row in scored_rows:
        key = tuple(row[field] for field in fields)
        grouped.setdefault(key, []).append(row)

    summary: list[dict[str, str]] = []
    for key in sorted(grouped):
        rows = grouped[key]
        scored_count = sum(1 for row in rows if row.get("score_status") == "scored")
        missing_radial_count = sum(1 for row in rows if not row.get("radial_csv"))
        complete_count = sum(1 for row in rows if row.get("score_completeness") == "1.000000000")
        best = min(rows, key=lambda row: float(row["combined_rmsd"]) if row.get("combined_rmsd") else float("inf"))
        summary.append(
            {
                "rise_label": key[0],
                "sidechain_variant": key[1],
                "twist_deg": key[2],
                "candidate_count": str(len(rows)),
                "scored_count": str(scored_count),
                "complete_score_count": str(complete_count),
                "missing_radial_count": str(missing_radial_count),
                "best_model_id": best.get("model_id", ""),
                "best_candidate_id": best.get("candidate_id", ""),
                "best_combined_rmsd": best.get("combined_rmsd", ""),
                "best_base_rmsd": best.get("base_rmsd", ""),
                "best_helical_rmsd": best.get("helical_rmsd", ""),
            }
        )
    return summary


def summary_fields() -> list[str]:
    return [
        "rise_label",
        "sidechain_variant",
        "twist_deg",
        "candidate_count",
        "scored_count",
        "complete_score_count",
        "missing_radial_count",
        "best_model_id",
        "best_candidate_id",
        "best_combined_rmsd",
        "best_base_rmsd",
        "best_helical_rmsd",
    ]


def validate_scored_rows(scored_rows: list[dict[str, str]]) -> None:
    required_metadata = [
        "model_id",
        "rise_label",
        "rise_A",
        "twist_deg",
        "sidechain_variant",
        "candidate_id",
        "radial_csv",
    ]
    required_scores = [
        "observed_base_d_A",
        "observed_A_d_A",
        "observed_B_d_A",
        "observed_C_d_A",
        "observed_D_d_A",
        "base_rmsd",
        "helical_rmsd",
        "combined_rmsd",
        "observed_peak_count",
        "expected_peak_count",
        "missing_peak_count",
        "score_completeness",
    ]

    if len(scored_rows) != EXPECTED_ROW_COUNT:
        raise ValueError(f"Expected {EXPECTED_ROW_COUNT} scored rows, found {len(scored_rows)}")

    for field in required_metadata + required_scores:
        missing = [row.get("model_id", "<unknown>") for row in scored_rows if not row.get(field)]
        if missing:
            raise ValueError(f"Missing required field {field!r} for {len(missing)} rows; examples: {missing[:5]}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--xyz-manifest", type=Path, required=True)
    parser.add_argument("--radial-manifest", type=Path, required=True)
    parser.add_argument("--targets", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--window-half-width", type=float, default=0.20)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    repo_root = Path(__file__).resolve().parents[1]
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    primary_rows = read_csv(args.manifest)
    xyz_rows = read_csv(args.xyz_manifest)
    radial_rows = read_csv(args.radial_manifest)

    xyz_by_model = require_unique(xyz_rows, "model_id", args.xyz_manifest)
    radial_by_model = radial_manifest_by_model(radial_rows, args.radial_manifest)
    joined_rows = join_manifests(primary_rows, xyz_by_model, radial_by_model, repo_root)

    combined_profiles = combine_profiles(joined_rows, repo_root)
    combined_profiles_path = output_dir / "combined_profiles_for_scoring.csv"
    write_csv(combined_profiles_path, combined_profiles, ["model_id", "d_A", "q", "intensity"])

    targets = load_extraction_targets(args.targets)
    profiles_by_model = load_profiles(combined_profiles_path)
    observed_long = extract_peaks(
        profiles_by_model=profiles_by_model,
        targets=targets,
        window_half_width=args.window_half_width,
        include_missing=True,
    )
    observed_long_path = output_dir / "observed_peaks.csv"
    write_extraction_csv(observed_long_path, observed_long, extraction_output_fields())

    validate_long_rows(observed_long, observed_long_path)
    target_order = parse_target_order(",".join(TARGET_ORDER))
    observed_wide = long_to_wide(observed_long, target_order)
    observed_wide_path = output_dir / "observed_peaks_wide.csv"
    write_observed_csv(observed_wide_path, observed_wide, observed_output_fields(observed_wide, target_order))

    scoring_targets = load_scoring_targets(args.targets)
    observed_by_model = load_observed_peaks(observed_wide_path)
    scored_rows = score_rows(joined_rows, scoring_targets, observed_by_model)
    validate_scored_rows(scored_rows)

    scored_path = output_dir / "scored_radial_profiles.csv"
    write_csv(scored_path, scored_rows, merged_fields(scored_rows))

    summary = score_summary_rows(scored_rows)
    summary_path = output_dir / "score_summary.csv"
    write_csv(summary_path, summary, summary_fields())

    counts = Counter(
        (row["rise_label"], row["sidechain_variant"], row["twist_deg"])
        for row in scored_rows
    )

    print(f"Joined rows: {len(joined_rows)}")
    print(f"Combined profile rows: {len(combined_profiles)}")
    print(f"Observed peak rows: {len(observed_long)}")
    print(f"Wide observed rows: {len(observed_wide)}")
    print(f"Scored rows: {len(scored_rows)}")
    print(f"Score status counts: {dict(Counter(row.get('score_status', '') for row in scored_rows))}")
    print(f"Missing peak count distribution: {dict(Counter(row.get('missing_peak_count', '') for row in scored_rows))}")
    print("Counts by rise_label x sidechain_variant x twist_deg:")
    for key in sorted(counts):
        print(f"  {key[0]} {key[1]} twist {key[2]}: {counts[key]}")
    print("Outputs:")
    for path in [combined_profiles_path, observed_long_path, observed_wide_path, scored_path, summary_path]:
        print(f"  {path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
