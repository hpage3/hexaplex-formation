"""Run a twist/rise profile-extraction smoke workflow.

This is an end-to-end data-contract smoke test for the twist/rise scoring path.
It does not generate diffraction. It starts from an existing profile CSV with:

    model_id,d_A,intensity

Then it runs:

    profile peak extraction
      -> long-form observed peaks
      -> wide observed-peaks table
      -> scored twist/rise manifest

The purpose is to prove that profile-shaped data can flow through the
extraction and scoring pipeline before launching a large twist/rise scan.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.build_twist_rise_observed_peaks_table import (
    long_to_wide,
    output_fields as observed_output_fields,
    parse_target_order,
    validate_long_rows,
    write_csv as write_observed_csv,
)
from scripts.extract_twist_rise_profile_peaks import (
    extract_peaks,
    load_profiles,
    load_targets as load_extraction_targets,
    output_fields as extraction_output_fields,
    write_csv as write_extraction_csv,
)
from scripts.score_twist_rise_scan import (
    load_observed_peaks,
    load_targets as load_scoring_targets,
    merged_fields,
    read_csv,
    score_rows,
    write_csv as write_score_csv,
)


def count_status(rows: list[dict[str, str]], field: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        value = row.get(field, "") or "blank"
        counts[value] = counts.get(value, 0) + 1
    return counts


def write_manifest(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def write_demo_profiles(path: Path) -> None:
    """Write a tiny synthetic profile fixture.

    These values are only for pipeline validation. They are not real extracted
    diffraction results.
    """

    path.parent.mkdir(parents=True, exist_ok=True)
    rows = [
        ["model_id", "d_A", "intensity"],
        ["twist30p0_rise3p400", "3.30", "5"],
        ["twist30p0_rise3p400", "3.407", "100"],
        ["twist30p0_rise3p400", "3.853", "90"],
        ["twist30p0_rise3p400", "4.465", "80"],
        ["twist30p0_rise3p400", "5.730", "70"],
        ["twist30p0_rise3p400", "7.630", "60"],
        ["twist29p0_rise3p400", "3.396", "100"],
        ["twist29p0_rise3p400", "3.895", "90"],
        ["twist29p0_rise3p400", "4.460", "80"],
        ["twist29p0_rise3p400", "5.790", "70"],
        ["twist29p0_rise3p400", "7.465", "60"],
        ["twist31p0_rise3p400", "3.407", "100"],
        ["twist31p0_rise3p400", "3.770", "90"],
        ["twist31p0_rise3p400", "4.390", "80"],
        ["twist31p0_rise3p400", "5.550", "70"],
        ["twist31p0_rise3p400", "7.865", "60"],
    ]

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerows(rows)


def run_smoke(args: argparse.Namespace) -> dict[str, object]:
    out_dir: Path = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    profile_path: Path = args.profiles
    if args.write_demo_profiles:
        write_demo_profiles(profile_path)

    long_peaks_path = out_dir / "observed_peaks_long.csv"
    wide_peaks_path = out_dir / "observed_peaks_wide.csv"
    scored_manifest_path = out_dir / "twist_rise_scored_manifest.csv"
    smoke_manifest_path = out_dir / "smoke_manifest.json"

    extraction_targets = load_extraction_targets(args.targets)
    profiles_by_model = load_profiles(profile_path)

    extracted_rows = extract_peaks(
        profiles_by_model=profiles_by_model,
        targets=extraction_targets,
        window_half_width=args.window_half_width,
        include_missing=args.include_missing,
    )
    write_extraction_csv(long_peaks_path, extracted_rows, extraction_output_fields())

    long_rows = read_csv(long_peaks_path)
    validate_long_rows(long_rows, long_peaks_path)
    target_order = parse_target_order(args.target_order)
    wide_rows = long_to_wide(long_rows, target_order)
    write_observed_csv(wide_peaks_path, wide_rows, observed_output_fields(wide_rows, target_order))

    manifest_rows = read_csv(args.manifest)
    scoring_targets = load_scoring_targets(args.targets)
    observed_by_model = load_observed_peaks(wide_peaks_path)
    scored_rows = score_rows(manifest_rows, scoring_targets, observed_by_model)
    write_score_csv(scored_manifest_path, scored_rows, merged_fields(scored_rows))

    scored_subset = [
        row
        for row in scored_rows
        if row.get("score_status") == "scored"
    ]

    payload: dict[str, object] = {
        "status": "success",
        "profile_input": str(profile_path),
        "target_input": str(args.targets),
        "grid_manifest_input": str(args.manifest),
        "window_half_width_A": args.window_half_width,
        "include_missing": args.include_missing,
        "outputs": {
            "observed_peaks_long": str(long_peaks_path),
            "observed_peaks_wide": str(wide_peaks_path),
            "scored_manifest": str(scored_manifest_path),
            "smoke_manifest": str(smoke_manifest_path),
        },
        "counts": {
            "models_in_profile": len(profiles_by_model),
            "targets": len(extraction_targets),
            "extracted_peak_rows": len(extracted_rows),
            "wide_observed_rows": len(wide_rows),
            "scored_manifest_rows": len(scored_rows),
            "scored_rows": len(scored_subset),
            "peak_status": count_status(extracted_rows, "peak_status"),
        },
    }
    write_manifest(smoke_manifest_path, payload)
    return payload


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--profiles",
        type=Path,
        default=Path("outputs/twist_rise_profile_smoke/synthetic_profiles.csv"),
        help="Profile CSV with model_id,d_A,intensity. Can be auto-created with --write-demo-profiles.",
    )
    parser.add_argument(
        "--write-demo-profiles",
        action="store_true",
        help="Write a small synthetic profile CSV before running the smoke workflow.",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("outputs/twist_rise_scan/twist_rise_grid_manifest.csv"),
        help="Twist/rise grid manifest.",
    )
    parser.add_argument(
        "--targets",
        type=Path,
        default=Path("inputs/experimental_peak_windows/hxc590_twist_rise_scan_targets.csv"),
        help="Target peak CSV.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs/twist_rise_profile_smoke"),
        help="Directory for smoke workflow outputs.",
    )
    parser.add_argument(
        "--window-half-width",
        type=float,
        default=0.25,
        help="Half-width in Angstrom around each target d-spacing.",
    )
    parser.add_argument(
        "--include-missing",
        action="store_true",
        help="Write missing peak rows in the long-form output.",
    )
    parser.add_argument(
        "--target-order",
        default="base,A,B,C,D",
        help="Comma-separated target order for wide observed peak columns.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()

    if args.window_half_width <= 0:
        raise ValueError("--window-half-width must be positive")

    payload = run_smoke(args)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
