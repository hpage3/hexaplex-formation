#!/usr/bin/env python3
"""Score predefined d-spacing windows in a radial profile."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from hexaplex_formation.scattering import integrate_window, q_from_d  # noqa: E402


WINDOWS = [
    ("d_8p4", 8.0, 8.8),
    ("d_5p5_6p0", 5.4, 6.1),
    ("d_4p5", 4.35, 4.65),
    ("d_4p1", 3.95, 4.25),
    ("d_3p4", 3.25, 3.50),
    ("d_3p0", 2.90, 3.10),
]

FIELDNAMES = [
    "window_name",
    "d_min_A",
    "d_max_A",
    "q_min_Ainv",
    "q_max_Ainv",
    "point_count",
    "mean_intensity",
    "max_intensity",
    "integrated_intensity",
    "mean_intensity_fraction_of_total",
    "integrated_intensity_fraction_of_total",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    return parser.parse_args()


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def score_windows(profile_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    total_integrated = sum(float(row["intensity"]) for row in profile_rows)
    total_mean = total_integrated / len(profile_rows) if profile_rows else 0.0
    rows: list[dict[str, str]] = []
    for name, d_min, d_max in WINDOWS:
        summary = integrate_window(profile_rows, d_min=d_min, d_max=d_max)
        mean_intensity = float(summary["mean_intensity"] or 0.0)
        integrated = float(summary["integrated_intensity"] or 0.0)
        rows.append(
            {
                "window_name": name,
                "d_min_A": f"{d_min:.2f}",
                "d_max_A": f"{d_max:.2f}",
                "q_min_Ainv": f"{q_from_d(d_max):.6f}",
                "q_max_Ainv": f"{q_from_d(d_min):.6f}",
                **summary,
                "mean_intensity_fraction_of_total": f"{mean_intensity / total_mean:.6f}" if total_mean else "",
                "integrated_intensity_fraction_of_total": f"{integrated / total_integrated:.6f}" if total_integrated else "",
            }
        )
    return rows


def main() -> int:
    args = parse_args()
    rows = score_windows(_rows(args.profile))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
