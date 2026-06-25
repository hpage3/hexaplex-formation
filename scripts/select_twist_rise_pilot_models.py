"""Select a small twist/rise pilot manifest from the full scan manifest.

Default pilot grid:
  twists: 28,29,30,31,32
  rises: 3.35,3.40,3.45,3.50

The output keeps the same manifest columns as the input and adds:
  pilot_selection = true
  pilot_notes
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


DEFAULT_TWISTS = "28,29,30,31,32"
DEFAULT_RISES = "3.35,3.40,3.45,3.50"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_float_set(text: str) -> set[float]:
    values: set[float] = set()
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        values.add(round(float(item), 6))

    if not values:
        raise ValueError("Selection list cannot be empty")

    return values


def require_columns(rows: list[dict[str, str]], path: Path) -> None:
    if not rows:
        raise ValueError(f"No rows found in {path}")

    required = {"model_id", "twist_deg", "rise_A"}
    missing = required - set(rows[0])
    if missing:
        raise ValueError(f"{path} is missing required columns: {sorted(missing)}")


def select_pilot_rows(
    rows: list[dict[str, str]],
    twists: set[float],
    rises: set[float],
) -> list[dict[str, str]]:
    selected: list[dict[str, str]] = []

    for row in rows:
        twist = round(float(row["twist_deg"]), 6)
        rise = round(float(row["rise_A"]), 6)

        if twist in twists and rise in rises:
            out = dict(row)
            out["pilot_selection"] = "true"
            out["pilot_notes"] = "twist/rise pilot grid"
            selected.append(out)

    return selected


def output_fields(input_rows: list[dict[str, str]], selected_rows: list[dict[str, str]]) -> list[str]:
    fields: list[str] = []
    seen: set[str] = set()

    for row in input_rows + selected_rows:
        for field in row:
            if field not in seen:
                fields.append(field)
                seen.add(field)

    for field in ["pilot_selection", "pilot_notes"]:
        if field not in seen:
            fields.append(field)
            seen.add(field)

    return fields


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("outputs/twist_rise_scan/twist_rise_grid_manifest.csv"),
        help="Full twist/rise scan manifest.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("outputs/twist_rise_scan/twist_rise_pilot_manifest.csv"),
        help="Pilot manifest output path.",
    )
    parser.add_argument(
        "--twists",
        default=DEFAULT_TWISTS,
        help="Comma-separated twist values to select.",
    )
    parser.add_argument(
        "--rises",
        default=DEFAULT_RISES,
        help="Comma-separated rise values to select.",
    )
    parser.add_argument(
        "--expected-count",
        type=int,
        default=20,
        help="Fail unless this many rows are selected. Use 0 to disable.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()

    rows = read_csv(args.manifest)
    require_columns(rows, args.manifest)

    twists = parse_float_set(args.twists)
    rises = parse_float_set(args.rises)

    selected = select_pilot_rows(rows, twists, rises)

    if args.expected_count and len(selected) != args.expected_count:
        raise ValueError(
            f"Expected {args.expected_count} pilot rows, selected {len(selected)}. "
            f"Twists={sorted(twists)}, rises={sorted(rises)}"
        )

    fields = output_fields(rows, selected)
    write_csv(args.output, selected, fields)

    print(f"Wrote {len(selected)} pilot rows to {args.output}")
    print(f"Twists: {', '.join(str(v) for v in sorted(twists))}")
    print(f"Rises: {', '.join(str(v) for v in sorted(rises))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
