"""Score theoretical peak positions against experimental d-spacing targets."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


SUMMARY_FIELDS = [
    "rank",
    "model_id",
    "twist_deg",
    "target_count",
    "assigned_peak_count",
    "missing_peak_count",
    "root_sum_relative_error",
    "rms_relative_error",
    "weighted_rms_relative_error",
    "max_abs_relative_error",
    "mean_abs_relative_error",
    "missing_peak_penalty",
]

PER_PEAK_FIELDS = [
    "model_id",
    "twist_deg",
    "target_label",
    "experimental_d_A",
    "theoretical_d_A",
    "relative_error",
    "abs_relative_error",
    "weight",
    "is_missing",
    "assignment_method",
    "matched_peak_d_A",
    "peak_intensity_optional",
    "notes",
]

ASSIGNMENT_FIELDS = [
    "model_id",
    "twist_deg",
    "target_label",
    "theoretical_d_A",
    "assignment_method",
    "notes",
]


@dataclass(frozen=True)
class TargetPeak:
    target_label: str
    experimental_d_A: float
    default_weight: float
    structural_note: str = ""


@dataclass(frozen=True)
class PeakAssignment:
    model_id: str
    twist_deg: str
    target_label: str
    theoretical_d_A: float
    assignment_method: str = ""
    notes: str = ""
    matched_peak_d_A: str = ""
    peak_intensity_optional: str = ""


def _as_float(value: str, field_name: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid {field_name}: {value!r}") from exc


def relative_error(experimental_d_A: float, theoretical_d_A: float) -> float:
    """Signed relative error, using (experimental - theoretical) / experimental."""

    if experimental_d_A == 0:
        raise ValueError("experimental_d_A must be non-zero")
    return (experimental_d_A - theoretical_d_A) / experimental_d_A


def root_sum_relative_error(errors: Iterable[float]) -> float:
    values = list(errors)
    return math.sqrt(sum(error * error for error in values))


def rms_relative_error(errors: Iterable[float]) -> float:
    values = list(errors)
    if not values:
        return 0.0
    return math.sqrt(sum(error * error for error in values) / len(values))


def weighted_rms_relative_error(errors: Iterable[float], weights: Iterable[float]) -> float:
    error_values = list(errors)
    weight_values = list(weights)
    if len(error_values) != len(weight_values):
        raise ValueError("errors and weights must have the same length")
    weight_sum = sum(weight_values)
    if weight_sum <= 0:
        return 0.0
    return math.sqrt(
        sum(weight * error * error for error, weight in zip(error_values, weight_values))
        / weight_sum
    )


def read_targets(path: Path) -> list[TargetPeak]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"target_label", "experimental_d_A", "default_weight"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing required columns: {sorted(missing)}")
        return [
            TargetPeak(
                target_label=row["target_label"],
                experimental_d_A=_as_float(row["experimental_d_A"], "experimental_d_A"),
                default_weight=_as_float(row["default_weight"], "default_weight"),
                structural_note=row.get("structural_note", ""),
            )
            for row in reader
        ]


def read_assignments(path: Path) -> list[PeakAssignment]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"model_id", "twist_deg", "target_label", "theoretical_d_A"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing required columns: {sorted(missing)}")
        assignments = []
        for row in reader:
            if not any((value or "").strip() for value in row.values()):
                continue
            assignments.append(
                PeakAssignment(
                    model_id=row["model_id"],
                    twist_deg=row["twist_deg"],
                    target_label=row["target_label"],
                    theoretical_d_A=_as_float(row["theoretical_d_A"], "theoretical_d_A"),
                    assignment_method=row.get("assignment_method", ""),
                    notes=row.get("notes", ""),
                )
            )
        return assignments


def match_peak_list(
    targets: list[TargetPeak],
    peak_rows: Iterable[dict[str, str]],
    tolerance_A: float,
) -> list[PeakAssignment]:
    """Match each target to the nearest listed peak for each model within tolerance."""

    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in peak_rows:
        grouped[(row["model_id"], row["twist_deg"])].append(row)

    assignments: list[PeakAssignment] = []
    for (model_id, twist_deg), rows in grouped.items():
        available = list(rows)
        for target in targets:
            best_index = None
            best_delta = None
            best_peak = None
            for index, peak in enumerate(available):
                peak_d_A = _as_float(peak["peak_d_A"], "peak_d_A")
                delta = abs(peak_d_A - target.experimental_d_A)
                if best_delta is None or delta < best_delta:
                    best_index = index
                    best_delta = delta
                    best_peak = peak
            if best_peak is None or best_delta is None or best_delta > tolerance_A:
                continue
            available.pop(best_index)
            peak_d_A = _as_float(best_peak["peak_d_A"], "peak_d_A")
            assignments.append(
                PeakAssignment(
                    model_id=model_id,
                    twist_deg=twist_deg,
                    target_label=target.target_label,
                    theoretical_d_A=peak_d_A,
                    assignment_method=f"nearest_within_{tolerance_A:g}A",
                    notes=best_peak.get("notes", ""),
                    matched_peak_d_A=f"{peak_d_A:.6g}",
                    peak_intensity_optional=best_peak.get("peak_intensity_optional", ""),
                )
            )
    return assignments


def read_peak_list(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"model_id", "twist_deg", "peak_d_A"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing required columns: {sorted(missing)}")
        return list(reader)


def score_models(
    targets: list[TargetPeak],
    assignments: list[PeakAssignment],
    missing_peak_penalty: float = 0.25,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    grouped: dict[tuple[str, str], dict[str, PeakAssignment]] = defaultdict(dict)
    for assignment in assignments:
        grouped[(assignment.model_id, assignment.twist_deg)][assignment.target_label] = assignment

    summary_rows: list[dict[str, object]] = []
    per_peak_rows: list[dict[str, object]] = []
    for (model_id, twist_deg), by_label in grouped.items():
        errors: list[float] = []
        weights: list[float] = []
        missing_count = 0
        assigned_count = 0
        for target in targets:
            assignment = by_label.get(target.target_label)
            weights.append(target.default_weight)
            if assignment is None:
                error = missing_peak_penalty
                missing_count += 1
                per_peak_rows.append(
                    {
                        "model_id": model_id,
                        "twist_deg": twist_deg,
                        "target_label": target.target_label,
                        "experimental_d_A": f"{target.experimental_d_A:.8g}",
                        "theoretical_d_A": "",
                        "relative_error": f"{error:.12g}",
                        "abs_relative_error": f"{abs(error):.12g}",
                        "weight": f"{target.default_weight:.8g}",
                        "is_missing": "true",
                        "assignment_method": "",
                        "matched_peak_d_A": "",
                        "peak_intensity_optional": "",
                        "notes": "missing target; penalty applied",
                    }
                )
            else:
                error = relative_error(target.experimental_d_A, assignment.theoretical_d_A)
                assigned_count += 1
                per_peak_rows.append(
                    {
                        "model_id": model_id,
                        "twist_deg": twist_deg,
                        "target_label": target.target_label,
                        "experimental_d_A": f"{target.experimental_d_A:.8g}",
                        "theoretical_d_A": f"{assignment.theoretical_d_A:.8g}",
                        "relative_error": f"{error:.12g}",
                        "abs_relative_error": f"{abs(error):.12g}",
                        "weight": f"{target.default_weight:.8g}",
                        "is_missing": "false",
                        "assignment_method": assignment.assignment_method,
                        "matched_peak_d_A": assignment.matched_peak_d_A,
                        "peak_intensity_optional": assignment.peak_intensity_optional,
                        "notes": assignment.notes,
                    }
                )
            errors.append(error)

        abs_errors = [abs(error) for error in errors]
        summary_rows.append(
            {
                "rank": "",
                "model_id": model_id,
                "twist_deg": twist_deg,
                "target_count": len(targets),
                "assigned_peak_count": assigned_count,
                "missing_peak_count": missing_count,
                "root_sum_relative_error": f"{root_sum_relative_error(errors):.12g}",
                "rms_relative_error": f"{rms_relative_error(errors):.12g}",
                "weighted_rms_relative_error": f"{weighted_rms_relative_error(errors, weights):.12g}",
                "max_abs_relative_error": f"{(max(abs_errors) if abs_errors else 0.0):.12g}",
                "mean_abs_relative_error": f"{(sum(abs_errors) / len(abs_errors) if abs_errors else 0.0):.12g}",
                "missing_peak_penalty": f"{missing_peak_penalty:.12g}",
            }
        )

    summary_rows.sort(
        key=lambda row: (
            float(row["weighted_rms_relative_error"]),
            int(row["missing_peak_count"]),
            str(row["model_id"]),
        )
    )
    for index, row in enumerate(summary_rows, start=1):
        row["rank"] = index
    return summary_rows, per_peak_rows


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Score theoretical or matched peak positions against experimental d-spacing targets."
    )
    parser.add_argument(
        "--targets",
        default="inputs/experimental_peak_windows/hxc590_backbone_peak_targets.csv",
        type=Path,
        help="CSV of experimental target peaks.",
    )
    parser.add_argument(
        "--assignments",
        default="inputs/theoretical_peak_assignments/nick_example_29_30_31_peak_assignments.csv",
        type=Path,
        help="CSV of assigned theoretical peaks.",
    )
    parser.add_argument(
        "--peak-list",
        type=Path,
        help="Optional peak-list CSV for automatic nearest matching.",
    )
    parser.add_argument(
        "--auto-match",
        action="store_true",
        help="Build assignments by matching peak-list rows to target peaks.",
    )
    parser.add_argument(
        "--match-tolerance",
        default=0.35,
        type=float,
        help="Maximum d-spacing difference in Angstrom for automatic matching.",
    )
    parser.add_argument(
        "--missing-peak-penalty",
        default=0.25,
        type=float,
        help="Relative-error penalty used for missing target assignments.",
    )
    parser.add_argument(
        "--summary-output",
        default="outputs/metrics/hxc590_peak_position_fit_summary.csv",
        type=Path,
        help="Per-model summary output CSV.",
    )
    parser.add_argument(
        "--per-peak-output",
        default="outputs/metrics/hxc590_peak_position_fit_per_peak_errors.csv",
        type=Path,
        help="Per-target error output CSV.",
    )
    parser.add_argument(
        "--matched-assignments-output",
        type=Path,
        help="Optional CSV output for auto-matched assignments.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    targets = read_targets(args.targets)

    if args.auto_match:
        if args.peak_list is None:
            raise SystemExit("--auto-match requires --peak-list")
        assignments = match_peak_list(targets, read_peak_list(args.peak_list), args.match_tolerance)
        if args.matched_assignments_output:
            write_csv(
                args.matched_assignments_output,
                ASSIGNMENT_FIELDS,
                [
                    {
                        "model_id": assignment.model_id,
                        "twist_deg": assignment.twist_deg,
                        "target_label": assignment.target_label,
                        "theoretical_d_A": f"{assignment.theoretical_d_A:.8g}",
                        "assignment_method": assignment.assignment_method,
                        "notes": assignment.notes,
                    }
                    for assignment in assignments
                ],
            )
    else:
        assignments = read_assignments(args.assignments)

    summary_rows, per_peak_rows = score_models(
        targets,
        assignments,
        missing_peak_penalty=args.missing_peak_penalty,
    )
    write_csv(args.summary_output, SUMMARY_FIELDS, summary_rows)
    write_csv(args.per_peak_output, PER_PEAK_FIELDS, per_peak_rows)

    print(f"Read {len(targets)} targets and {len(assignments)} assigned peaks.")
    print(f"Wrote {len(summary_rows)} model rows to {args.summary_output}.")
    print(f"Wrote {len(per_peak_rows)} peak-error rows to {args.per_peak_output}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
