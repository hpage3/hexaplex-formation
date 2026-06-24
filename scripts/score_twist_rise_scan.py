"""Score a twist/rise scan manifest against HXC590 diffraction targets.

This scoring layer consumes:
  1. a twist/rise grid manifest,
  2. target peak definitions, and
  3. optionally, an observed peak table keyed by model_id.

The score is split into:
  - base_rmsd: base-stacking / rise-sensitive target(s)
  - helical_rmsd: backbone-associated / helical-geometry target(s)
  - combined_rmsd: all available target errors together
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Iterable


BASE_GROUP = "base_stacking"
HELICAL_GROUP = "backbone_associated"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def load_targets(path: Path) -> dict[str, dict[str, str]]:
    rows = read_csv(path)
    required = {"target_label", "target_d_A", "target_group"}
    if not rows:
        raise ValueError(f"No targets found in {path}")
    missing = required - set(rows[0])
    if missing:
        raise ValueError(f"Target file {path} is missing columns: {sorted(missing)}")

    targets: dict[str, dict[str, str]] = {}
    for row in rows:
        label = row["target_label"]
        if label in targets:
            raise ValueError(f"Duplicate target label {label!r} in {path}")
        float(row["target_d_A"])
        targets[label] = row
    return targets


def load_observed_peaks(path: Path | None) -> dict[str, dict[str, str]]:
    """Load observed peak assignments keyed by model_id.

    Expected columns include model_id plus any observed_<target_label>_d_A
    columns, for example:
      observed_base_d_A
      observed_A_d_A
      observed_B_d_A
    """

    if path is None:
        return {}

    rows = read_csv(path)
    if not rows:
        return {}

    if "model_id" not in rows[0]:
        raise ValueError(f"Observed peak file {path} is missing model_id column")

    observed: dict[str, dict[str, str]] = {}
    for row in rows:
        model_id = row["model_id"]
        if not model_id:
            raise ValueError(f"Observed peak file {path} has a row with empty model_id")
        if model_id in observed:
            raise ValueError(f"Duplicate model_id {model_id!r} in observed peak file {path}")
        observed[model_id] = row

    return observed


def merge_observed_peaks(
    manifest_rows: list[dict[str, str]],
    observed_by_model: dict[str, dict[str, str]],
) -> list[dict[str, str]]:
    """Merge observed peak columns into manifest rows by model_id."""

    if not observed_by_model:
        return [dict(row) for row in manifest_rows]

    merged: list[dict[str, str]] = []
    for row in manifest_rows:
        out = dict(row)
        model_id = out.get("model_id", "")
        observed = observed_by_model.get(model_id, {})
        for key, value in observed.items():
            if key == "model_id":
                continue
            if key.startswith("observed_"):
                out[key] = value
        merged.append(out)

    return merged


def rmsd(errors: Iterable[float]) -> str:
    values = list(errors)
    if not values:
        return ""
    return f"{math.sqrt(sum(v * v for v in values) / len(values)):.9f}"


def row_errors(row: dict[str, str], targets: dict[str, dict[str, str]], group: str) -> list[float]:
    """Return d-spacing errors for a target group."""

    errors: list[float] = []
    for label, target in targets.items():
        if target["target_group"] != group:
            continue
        observed_key = f"observed_{label}_d_A"
        observed_text = row.get(observed_key, "")
        if not observed_text:
            continue
        errors.append(float(observed_text) - float(target["target_d_A"]))
    return errors


def count_available_observations(row: dict[str, str], targets: dict[str, dict[str, str]]) -> int:
    count = 0
    for label in targets:
        if row.get(f"observed_{label}_d_A", ""):
            count += 1
    return count


def expected_observation_count(targets: dict[str, dict[str, str]]) -> int:
    return len(targets)


def completeness_fraction(observed_count: int, expected_count: int) -> str:
    if expected_count <= 0:
        return ""
    return f"{observed_count / expected_count:.9f}"


def score_rows(
    manifest_rows: list[dict[str, str]],
    targets: dict[str, dict[str, str]],
    observed_by_model: dict[str, dict[str, str]] | None = None,
) -> list[dict[str, str]]:
    rows_with_observed = merge_observed_peaks(manifest_rows, observed_by_model or {})
    scored: list[dict[str, str]] = []

    for row in rows_with_observed:
        out = dict(row)

        base_errors = row_errors(out, targets, BASE_GROUP)
        helical_errors = row_errors(out, targets, HELICAL_GROUP)

        base_score = rmsd(base_errors)
        helical_score = rmsd(helical_errors)
        combined_score = rmsd(base_errors + helical_errors)

        out["base_rmsd"] = base_score
        out["helical_rmsd"] = helical_score
        out["combined_rmsd"] = combined_score

        observed_count = count_available_observations(out, targets)
        expected_count = expected_observation_count(targets)
        missing_count = expected_count - observed_count

        out["observed_peak_count"] = str(observed_count)
        out["expected_peak_count"] = str(expected_count)
        out["missing_peak_count"] = str(missing_count)
        out["score_completeness"] = completeness_fraction(observed_count, expected_count)

        if combined_score:
            out["score_status"] = "scored"
        else:
            out["score_status"] = out.get("score_status") or "pending"

        scored.append(out)

    return scored


def merged_fields(rows: list[dict[str, str]]) -> list[str]:
    preferred = [
        "model_id",
        "twist_deg",
        "rise_A",
        "model_status",
        "coordinate_file",
        "diffraction_status",
        "diffraction_file",
        "score_status",
        "observed_peak_count",
        "expected_peak_count",
        "missing_peak_count",
        "score_completeness",
        "observed_base_d_A",
        "observed_A_d_A",
        "observed_B_d_A",
        "observed_C_d_A",
        "observed_D_d_A",
        "base_rmsd",
        "helical_rmsd",
        "combined_rmsd",
        "notes",
    ]

    seen = set()
    fields: list[str] = []
    for name in preferred:
        if any(name in row for row in rows):
            fields.append(name)
            seen.add(name)

    for row in rows:
        for name in row:
            if name not in seen:
                fields.append(name)
                seen.add(name)

    return fields


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("outputs/twist_rise_scan/twist_rise_grid_manifest.csv"),
    )
    parser.add_argument(
        "--targets",
        type=Path,
        default=Path("inputs/experimental_peak_windows/hxc590_twist_rise_scan_targets.csv"),
    )
    parser.add_argument(
        "--observed-peaks",
        type=Path,
        help="Optional observed peak table keyed by model_id.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("outputs/twist_rise_scan/twist_rise_scored_manifest.csv"),
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    manifest_rows = read_csv(args.manifest)
    targets = load_targets(args.targets)
    observed_by_model = load_observed_peaks(args.observed_peaks)
    scored = score_rows(manifest_rows, targets, observed_by_model)
    fields = merged_fields(scored)
    write_csv(args.output, scored, fields)
    print(f"Wrote {len(scored)} scored rows to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
