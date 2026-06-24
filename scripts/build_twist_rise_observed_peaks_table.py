"""Build the standardized observed-peaks table for twist/rise scoring.

This script converts a long-form peak assignment table:

    model_id,target_label,observed_d_A

into the wide table consumed by scripts/score_twist_rise_scan.py:

    model_id,observed_base_d_A,observed_A_d_A,...

The script does not perform peak picking. It only standardizes peak assignments
from whatever extraction method produced them.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


DEFAULT_TARGET_ORDER = ["base", "A", "B", "C", "D"]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_target_order(text: str | None) -> list[str]:
    if not text:
        return list(DEFAULT_TARGET_ORDER)
    targets = [item.strip() for item in text.split(",") if item.strip()]
    if not targets:
        raise ValueError("target order cannot be empty")
    if len(targets) != len(set(targets)):
        raise ValueError(f"target order contains duplicates: {targets}")
    return targets


def validate_long_rows(rows: list[dict[str, str]], source: Path) -> None:
    if not rows:
        raise ValueError(f"No rows found in {source}")

    required = {"model_id", "target_label", "observed_d_A"}
    missing = required - set(rows[0])
    if missing:
        raise ValueError(f"{source} is missing required columns: {sorted(missing)}")

    for index, row in enumerate(rows, start=2):
        if not row.get("model_id", "").strip():
            raise ValueError(f"{source} row {index} has empty model_id")
        if not row.get("target_label", "").strip():
            raise ValueError(f"{source} row {index} has empty target_label")
        if not row.get("observed_d_A", "").strip():
            raise ValueError(f"{source} row {index} has empty observed_d_A")
        float(row["observed_d_A"])


def long_to_wide(
    rows: list[dict[str, str]],
    target_order: list[str] | None = None,
) -> list[dict[str, str]]:
    """Convert long observed peaks to wide observed_<target>_d_A columns."""

    targets = target_order or list(DEFAULT_TARGET_ORDER)
    by_model: dict[str, dict[str, str]] = {}
    seen_pairs: set[tuple[str, str]] = set()

    for row in rows:
        model_id = row["model_id"].strip()
        target_label = row["target_label"].strip()
        observed_d = row["observed_d_A"].strip()

        pair = (model_id, target_label)
        if pair in seen_pairs:
            raise ValueError(f"Duplicate observed peak for model_id={model_id!r}, target_label={target_label!r}")
        seen_pairs.add(pair)

        if model_id not in by_model:
            by_model[model_id] = {"model_id": model_id}

        by_model[model_id][f"observed_{target_label}_d_A"] = observed_d

    # Preserve deterministic output order by model_id. This keeps diffs stable.
    wide_rows: list[dict[str, str]] = []
    for model_id in sorted(by_model):
        row = {"model_id": model_id}
        for target in targets:
            row[f"observed_{target}_d_A"] = by_model[model_id].get(f"observed_{target}_d_A", "")
        # Preserve non-standard targets, if any, after the default target columns.
        for key in sorted(by_model[model_id]):
            if key != "model_id" and key not in row:
                row[key] = by_model[model_id][key]
        wide_rows.append(row)

    return wide_rows


def output_fields(rows: list[dict[str, str]], target_order: list[str]) -> list[str]:
    fields = ["model_id"] + [f"observed_{target}_d_A" for target in target_order]
    seen = set(fields)

    for row in rows:
        for key in row:
            if key not in seen:
                fields.append(key)
                seen.add(key)

    return fields


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Long-form observed peaks CSV with model_id,target_label,observed_d_A.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("outputs/twist_rise_scan/observed_peaks.csv"),
        help="Wide observed-peaks CSV for score_twist_rise_scan.py.",
    )
    parser.add_argument(
        "--target-order",
        default="base,A,B,C,D",
        help="Comma-separated target order for wide output columns.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    target_order = parse_target_order(args.target_order)
    rows = read_csv(args.input)
    validate_long_rows(rows, args.input)
    wide_rows = long_to_wide(rows, target_order)
    fields = output_fields(wide_rows, target_order)
    write_csv(args.output, wide_rows, fields)
    print(f"Wrote {len(wide_rows)} observed-peak rows to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
