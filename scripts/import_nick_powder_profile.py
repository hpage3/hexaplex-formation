"""Import Nick/Asem experimental powder profile from the preserved RTF file."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


RTF_PREFIX = re.compile(r"^(?:[{}]\s*|\\[a-zA-Z*]+-?\d* ?)+")
NUMERIC_ROW = re.compile(r"([0-9]+(?:\.[0-9]+)?)\s+([0-9]+(?:\.[0-9]+)?)\\?")


def parse_profile(path: Path) -> list[tuple[float, float]]:
    text = path.read_text(encoding="utf-8", errors="replace")
    rows = []
    for raw_line in text.splitlines():
        line = RTF_PREFIX.sub("", raw_line.strip()).strip()
        match = NUMERIC_ROW.fullmatch(line)
        if match:
            d_spacing, intensity = match.groups()
            rows.append((float(d_spacing), float(intensity)))
    if not rows:
        raise ValueError(f"No numeric powder profile rows found in {path}")
    return rows


def write_csv(rows: list[tuple[float, float]], path: Path) -> None:
    max_intensity = max(intensity for _, intensity in rows)
    if max_intensity <= 0:
        raise ValueError("Maximum intensity must be positive for normalization")

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["d_A", "intensity_raw", "intensity_normalized"])
        for d_spacing, intensity in rows:
            writer.writerow(
                [
                    f"{d_spacing:.15g}",
                    f"{intensity:.15g}",
                    f"{intensity / max_intensity:.15g}",
                ]
            )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("reference/nick_asem_powder_archive/Experimental_Powder_Diffracxtion_Pattern.rtf"),
        help="Preserved RTF experimental powder profile.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("inputs/experimental/nick_powder_profile.csv"),
        help="Output CSV with d_A, intensity_raw, and intensity_normalized.",
    )
    args = parser.parse_args()

    rows = parse_profile(args.input)
    write_csv(rows, args.output)
    d_values = [row[0] for row in rows]
    max_norm = max(row[1] for row in rows) / max(row[1] for row in rows)
    print(
        f"Wrote {len(rows)} rows to {args.output}; "
        f"d_A range {min(d_values):.6g}-{max(d_values):.6g}; "
        f"max normalized intensity {max_norm:.6g}"
    )


if __name__ == "__main__":
    main()
