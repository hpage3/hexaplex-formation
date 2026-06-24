"""Build a twist/rise scan manifest for Hexaplex model-family searches.

This script does not generate coordinates or diffraction outputs. It only creates
a deterministic CSV manifest describing the twist/rise grid to be evaluated by
later pipeline stages.
"""

from __future__ import annotations

import argparse
import csv
from decimal import Decimal
from pathlib import Path


FIELDS = [
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


def decimal_range(start: str, stop: str, step: str) -> list[Decimal]:
    start_d = Decimal(str(start))
    stop_d = Decimal(str(stop))
    step_d = Decimal(str(step))

    if step_d <= 0:
        raise ValueError("step must be positive")
    if stop_d < start_d:
        raise ValueError("stop must be greater than or equal to start")

    values: list[Decimal] = []
    current = start_d
    epsilon = step_d / Decimal("1000000")

    while current <= stop_d + epsilon:
        values.append(current)
        current += step_d

    return values


def format_decimal(value: Decimal, places: int = 3) -> str:
    quant = Decimal("1").scaleb(-places)
    return str(value.quantize(quant))


def model_id(twist: Decimal, rise: Decimal) -> str:
    twist_text = format_decimal(twist, 1).replace("-", "m").replace(".", "p")
    rise_text = format_decimal(rise, 3).replace("-", "m").replace(".", "p")
    return f"twist{twist_text}_rise{rise_text}"


def build_rows(
    twist_min: str,
    twist_max: str,
    twist_step: str,
    rise_min: str,
    rise_max: str,
    rise_step: str,
) -> list[dict[str, str]]:
    twists = decimal_range(twist_min, twist_max, twist_step)
    rises = decimal_range(rise_min, rise_max, rise_step)

    rows: list[dict[str, str]] = []
    for rise in rises:
        for twist in twists:
            rows.append(
                {
                    "model_id": model_id(twist, rise),
                    "twist_deg": format_decimal(twist, 3),
                    "rise_A": format_decimal(rise, 3),
                    "model_status": "pending",
                    "coordinate_file": "",
                    "diffraction_status": "pending",
                    "diffraction_file": "",
                    "score_status": "pending",
                    "base_rmsd": "",
                    "helical_rmsd": "",
                    "combined_rmsd": "",
                    "notes": "",
                }
            )
    return rows


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--twist-min", default="20.0")
    parser.add_argument("--twist-max", default="40.0")
    parser.add_argument("--twist-step", default="0.5")
    parser.add_argument("--rise-min", default="3.0")
    parser.add_argument("--rise-max", default="4.0")
    parser.add_argument("--rise-step", default="0.05")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("outputs/twist_rise_scan/twist_rise_grid_manifest.csv"),
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    rows = build_rows(
        twist_min=args.twist_min,
        twist_max=args.twist_max,
        twist_step=args.twist_step,
        rise_min=args.rise_min,
        rise_max=args.rise_max,
        rise_step=args.rise_step,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
