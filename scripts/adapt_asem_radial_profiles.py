"""Adapt Asem-corrected radial profile CSVs to the twist/rise profile format.

Asem radial profile files already contain calibrated d-spacing and intensity
columns, for example:

    d_A,intensity_mean

This adapter writes the standard profile table consumed by
extract_twist_rise_profile_peaks.py:

    model_id,d_A,intensity

By default, model_id is derived from the input filename by removing
"_radial" and the extension:

    compact_hexaplex_twist_30_radial.csv -> compact_hexaplex_twist_30
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


DEFAULT_D_COLUMN = "d_A"
DEFAULT_INTENSITY_COLUMN = "intensity_mean"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def derive_model_id(path: Path) -> str:
    stem = path.stem
    if stem.endswith("_radial"):
        stem = stem[: -len("_radial")]
    return stem


def adapt_profile_rows(
    rows: list[dict[str, str]],
    model_id: str,
    source: Path,
    d_column: str = DEFAULT_D_COLUMN,
    intensity_column: str = DEFAULT_INTENSITY_COLUMN,
) -> list[dict[str, str]]:
    if not rows:
        raise ValueError(f"No rows found in {source}")

    missing = {d_column, intensity_column} - set(rows[0])
    if missing:
        raise ValueError(f"{source} is missing required columns: {sorted(missing)}")

    adapted: list[dict[str, str]] = []

    for index, row in enumerate(rows, start=2):
        d_text = row.get(d_column, "").strip()
        intensity_text = row.get(intensity_column, "").strip()

        if not d_text:
            raise ValueError(f"{source} row {index} has empty {d_column}")
        if not intensity_text:
            raise ValueError(f"{source} row {index} has empty {intensity_column}")

        d_A = float(d_text)
        intensity = float(intensity_text)

        adapted.append(
            {
                "model_id": model_id,
                "d_A": f"{d_A:.9f}",
                "intensity": f"{intensity:.9f}",
            }
        )

    return adapted


def adapt_files(
    input_files: list[Path],
    d_column: str = DEFAULT_D_COLUMN,
    intensity_column: str = DEFAULT_INTENSITY_COLUMN,
) -> list[dict[str, str]]:
    all_rows: list[dict[str, str]] = []

    for path in sorted(input_files):
        model_id = derive_model_id(path)
        rows = read_csv(path)
        all_rows.extend(
            adapt_profile_rows(
                rows=rows,
                model_id=model_id,
                source=path,
                d_column=d_column,
                intensity_column=intensity_column,
            )
        )

    return all_rows


def collect_input_files(input_path: Path, pattern: str) -> list[Path]:
    if input_path.is_file():
        return [input_path]

    if input_path.is_dir():
        files = sorted(input_path.glob(pattern))
        if not files:
            raise ValueError(f"No files matching {pattern!r} found in {input_path}")
        return files

    raise ValueError(f"Input path does not exist: {input_path}")


def output_fields() -> list[str]:
    return ["model_id", "d_A", "intensity"]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Asem radial CSV file or directory of radial CSV files.",
    )
    parser.add_argument(
        "--pattern",
        default="*_radial.csv",
        help="Glob pattern used when --input is a directory.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("outputs/twist_rise_profile_smoke/asem_profiles.csv"),
        help="Standard profile CSV with model_id,d_A,intensity.",
    )
    parser.add_argument(
        "--d-column",
        default=DEFAULT_D_COLUMN,
        help="D-spacing column in the Asem radial CSV.",
    )
    parser.add_argument(
        "--intensity-column",
        default=DEFAULT_INTENSITY_COLUMN,
        help="Intensity column in the Asem radial CSV.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    input_files = collect_input_files(args.input, args.pattern)
    rows = adapt_files(
        input_files=input_files,
        d_column=args.d_column,
        intensity_column=args.intensity_column,
    )
    write_csv(args.output, rows, output_fields())
    print(f"Wrote {len(rows)} adapted profile rows from {len(input_files)} files to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
