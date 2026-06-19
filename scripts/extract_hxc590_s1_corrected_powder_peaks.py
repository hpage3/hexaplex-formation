#!/usr/bin/env python3
"""Extract conservative HXC590 S1 powder targets from the corrected trace."""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path


SAMPLE_ID = "HXC590_S1_TM_TC_without_salt_powder_corrected"
DEFAULT_TRACE = Path("inputs/experimental/hxc590_s1_powder_corrected_trace.csv")
DEFAULT_TARGETS = Path("outputs/metrics/hxc590_s1_powder_corrected_peak_targets.csv")
DEFAULT_SHIFT = Path("outputs/metrics/hxc590_s1_powder_peak_correction_shift.csv")

TARGET_COLUMNS = [
    "sample_id",
    "peak_label",
    "target_id",
    "distance_a",
    "d_angstrom",
    "q_inv_angstrom",
    "intensity",
    "region",
    "diagnostic_role",
    "confidence",
    "d_window_half_width_angstrom",
    "notes",
]

SHIFT_COLUMNS = [
    "feature_label",
    "previous_distance_a",
    "corrected_distance_a",
    "shift_a",
    "interpretation_note",
]


@dataclass(frozen=True)
class TracePoint:
    row_index: int
    distance_a: float
    intensity: float


@dataclass(frozen=True)
class TargetRegion:
    peak_label: str
    lower_a: float
    upper_a: float
    region: str
    diagnostic_role: str
    confidence: str
    window_half_width_a: float
    notes: str


REGIONS = [
    TargetRegion(
        "corrected_broad_7p3",
        7.05,
        7.55,
        "broad ~7.3 A feature",
        "supporting",
        "lower",
        0.25,
        "Corrected broad high-d feature; useful as a supporting spacing window, not a strong diagnostic by itself.",
    ),
    TargetRegion(
        "corrected_5p56",
        5.55,
        5.60,
        "~5.55-5.60 A feature",
        "supporting",
        "lower",
        0.15,
        "Corrected mid-d feature. It is retained as a supporting window because the local background is broad.",
    ),
    TargetRegion(
        "corrected_4p40",
        4.35,
        4.45,
        "strong ~4.40 A feature",
        "diagnostic",
        "medium_high",
        0.10,
        "Strong corrected feature in the 4.4 A diagnostic region.",
    ),
    TargetRegion(
        "corrected_3p80",
        3.75,
        3.85,
        "~3.75-3.85 A feature",
        "diagnostic",
        "medium_high",
        0.08,
        "Corrected feature in the 3.7-3.9 A diagnostic region.",
    ),
    TargetRegion(
        "corrected_3p39_stacking",
        3.38,
        3.40,
        "~3.38-3.40 A base-stacking feature",
        "diagnostic",
        "medium_high",
        0.06,
        "Corrected base-stacking feature close to the expected ~3.4 A stacked-heterocycle separation.",
    ),
]


SHIFT_ROWS = [
    (
        "broad_7p3_feature",
        7.30,
        "corrected_broad_7p3",
        "The corrected high-d feature shifts only modestly relative to the previous broad 7.30 A target.",
    ),
    (
        "5p5_feature",
        5.50,
        "corrected_5p56",
        "The corrected mid-d feature remains close to the previous 5.50 A supporting target.",
    ),
    (
        "4p3_to_4p4_feature",
        4.33,
        "corrected_4p40",
        "The strong corrected feature is modestly higher in d-spacing than the previous 4.33 A target.",
    ),
    (
        "3p9_side_of_3p7_3p9_region",
        3.90,
        "corrected_3p80",
        "The corrected 3.75-3.85 A feature sits between the previous 3.90 A and 3.71 A targets.",
    ),
    (
        "3p71_side_of_3p7_3p9_region",
        3.71,
        "corrected_3p80",
        "The corrected 3.75-3.85 A feature sits between the previous 3.90 A and 3.71 A targets.",
    ),
    (
        "base_stacking_feature",
        3.35,
        "corrected_3p39_stacking",
        "The corrected base-stacking feature moves closer to the expected ~3.4 A region; this strengthens compatibility but does not uniquely determine phase, exact twist, or exact rise.",
    ),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trace-csv", type=Path, default=DEFAULT_TRACE)
    parser.add_argument("--targets-csv", type=Path, default=DEFAULT_TARGETS)
    parser.add_argument("--shift-csv", type=Path, default=DEFAULT_SHIFT)
    return parser.parse_args()


def q_from_d(distance_a: float) -> float:
    if distance_a <= 0:
        raise ValueError("distance must be positive")
    return 2.0 * math.pi / distance_a


def format_float(value: float, digits: int = 6) -> str:
    return f"{value:.{digits}f}"


def read_trace(path: Path) -> list[TracePoint]:
    if not path.exists():
        raise FileNotFoundError(f"Corrected trace file not found: {path}")
    rows: list[TracePoint] = []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != ["distance_a", "intensity"]:
            raise ValueError(f"Expected columns distance_a,intensity in {path}")
        for index, row in enumerate(reader, start=1):
            rows.append(TracePoint(index, float(row["distance_a"]), float(row["intensity"])))
    return rows


def sorted_for_peak_detection(points: list[TracePoint]) -> list[TracePoint]:
    return sorted(points, key=lambda point: (point.distance_a, point.row_index))


def local_maxima(points: list[TracePoint]) -> list[TracePoint]:
    sorted_points = sorted_for_peak_detection(points)
    maxima: list[TracePoint] = []
    for index in range(1, len(sorted_points) - 1):
        previous_i = sorted_points[index - 1].intensity
        current_i = sorted_points[index].intensity
        next_i = sorted_points[index + 1].intensity
        if current_i >= previous_i and current_i >= next_i and (current_i > previous_i or current_i > next_i):
            maxima.append(sorted_points[index])
    return maxima


def choose_region_peak(points: list[TracePoint], region: TargetRegion) -> TracePoint:
    candidates = [point for point in points if region.lower_a <= point.distance_a <= region.upper_a]
    if not candidates:
        raise ValueError(f"No corrected trace points found for {region.peak_label}")
    maxima_in_region = [point for point in local_maxima(points) if region.lower_a <= point.distance_a <= region.upper_a]
    pool = maxima_in_region or candidates
    return max(pool, key=lambda point: (point.intensity, -point.row_index))


def target_id_for_label(label: str) -> str:
    return label.replace("corrected_", "d_")


def extract_targets(points: list[TracePoint]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for region in REGIONS:
        peak = choose_region_peak(points, region)
        rows.append(
            {
                "sample_id": SAMPLE_ID,
                "peak_label": region.peak_label,
                "target_id": target_id_for_label(region.peak_label),
                "distance_a": format_float(peak.distance_a, 3),
                "d_angstrom": format_float(peak.distance_a, 3),
                "q_inv_angstrom": format_float(q_from_d(peak.distance_a)),
                "intensity": format_float(peak.intensity, 3),
                "region": region.region,
                "diagnostic_role": region.diagnostic_role,
                "confidence": region.confidence,
                "d_window_half_width_angstrom": format_float(region.window_half_width_a, 2),
                "notes": region.notes + " Extracted from the raw corrected trace after internal sorting for peak detection; raw trace rows were not deduplicated or rewritten.",
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_shift_table(path: Path, target_rows: list[dict[str, str]]) -> None:
    by_label = {row["peak_label"]: row for row in target_rows}
    rows: list[dict[str, str]] = []
    for feature_label, previous, corrected_label, note in SHIFT_ROWS:
        corrected = float(by_label[corrected_label]["distance_a"])
        rows.append(
            {
                "feature_label": feature_label,
                "previous_distance_a": format_float(previous, 3),
                "corrected_distance_a": format_float(corrected, 3),
                "shift_a": format_float(corrected - previous, 3),
                "interpretation_note": note,
            }
        )
    write_csv(path, rows, SHIFT_COLUMNS)


def run(args: argparse.Namespace) -> dict[str, object]:
    points = read_trace(args.trace_csv)
    targets = extract_targets(points)
    write_csv(args.targets_csv, targets, TARGET_COLUMNS)
    write_shift_table(args.shift_csv, targets)
    return {
        "trace_rows": len(points),
        "target_rows": len(targets),
        "stacking_distance": next(row["distance_a"] for row in targets if row["peak_label"] == "corrected_3p39_stacking"),
    }


def main() -> None:
    result = run(parse_args())
    print(f"Read {result['trace_rows']} corrected trace rows")
    print(f"Wrote {result['target_rows']} corrected peak targets")
    print(f"Corrected stacking target: {result['stacking_distance']} A")


if __name__ == "__main__":
    main()
