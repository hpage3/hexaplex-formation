"""Score a twist/rise scan manifest against HXC590 diffraction targets.

This is the scoring layer for the twist/rise model-family scan. It is designed
to consume a grid manifest plus later peak-extraction outputs and write base,
helical, and combined RMSD scores.

At this stage, rows without extracted peak assignments are preserved with
score_status=pending.
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
        writer = csv.DictWriter(handle, fieldnames=fields)
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


def rmsd(errors: Iterable[float]) -> str:
    values = list(errors)
    if not values:
        return ""
    return f"{math.sqrt(sum(v * v for v in values) / len(values)):.9f}"


def row_errors(row: dict[str, str], targets: dict[str, dict[str, str]], group: str) -> list[float]:
    """Return d-spacing errors for a target group.

    Later peak-extraction stages should add columns named:
      observed_<target_label>_d_A

    Example:
      observed_base_d_A
      observed_A_d_A
      observed_B_d_A
    """

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


def score_rows(
    manifest_rows: list[dict[str, str]],
    targets: dict[str, dict[str, str]],
) -> list[dict[str, str]]:
    scored: list[dict[str, str]] = []

    for row in manifest_rows:
        out = dict(row)

        base_errors = row_errors(out, targets, BASE_GROUP)
        helical_errors = row_errors(out, targets, HELICAL_GROUP)

        base_score = rmsd(base_errors)
        helical_score = rmsd(helical_errors)
        combined_score = rmsd(base_errors + helical_errors)

        out["base_rmsd"] = base_score
        out["helical_rmsd"] = helical_score
        out["combined_rmsd"] = combined_score

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
        "base_rmsd",
        "helical_rmsd",
        "combined_rmsd",
        "notes",
    ]

    seen = set()
    fields: list[str] = []
    for name in preferred:
        if rows and name in rows[0]:
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
        "--output",
        type=Path,
        default=Path("outputs/twist_rise_scan/twist_rise_scored_manifest.csv"),
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    manifest_rows = read_csv(args.manifest)
    targets = load_targets(args.targets)
    scored = score_rows(manifest_rows, targets)
    fields = merged_fields(scored)
    write_csv(args.output, scored, fields)
    print(f"Wrote {len(scored)} scored rows to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
