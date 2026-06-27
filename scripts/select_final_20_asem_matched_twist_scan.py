"""Select the final equalized 20-row Asem matched twist/rise comparison set.

Selection is one best candidate per rise_label x sidechain_variant x twist_deg.
Lower combined_rmsd is primary; deterministic tie-breakers are lower
helical_rmsd, lower base_rmsd, lower missing_peak_count, then lower candidate
numeric index.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path


GROUP_FIELDS = ["rise_label", "sidechain_variant", "twist_deg"]
REQUIRED_OUTPUT_FIELDS = [
    "model_id",
    "rise_label",
    "rise_A",
    "twist_deg",
    "sidechain_variant",
    "candidate_id",
    "radial_csv",
    "observed_base_d_A",
    "observed_A_d_A",
    "observed_B_d_A",
    "observed_C_d_A",
    "observed_D_d_A",
    "base_rmsd",
    "helical_rmsd",
    "combined_rmsd",
    "observed_peak_count",
    "missing_peak_count",
    "score_completeness",
]
PER_BAND_ERROR_PREFIXES = ("delta_", "error_", "abs_error_")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def numeric(value: str, field: str, model_id: str) -> float:
    text = (value or "").strip()
    if not text:
        raise ValueError(f"{model_id} has empty {field}")
    parsed = float(text)
    if not math.isfinite(parsed):
        raise ValueError(f"{model_id} has non-finite {field}: {value!r}")
    return parsed


def candidate_index(candidate_id: str) -> int:
    match = re.search(r"(\d+)$", candidate_id or "")
    if not match:
        return 10**9
    return int(match.group(1))


def sort_key(row: dict[str, str]) -> tuple[float, float, float, int, int]:
    model_id = row.get("model_id", "<unknown>")
    return (
        numeric(row.get("combined_rmsd", ""), "combined_rmsd", model_id),
        numeric(row.get("helical_rmsd", ""), "helical_rmsd", model_id),
        numeric(row.get("base_rmsd", ""), "base_rmsd", model_id),
        int(numeric(row.get("missing_peak_count", ""), "missing_peak_count", model_id)),
        candidate_index(row.get("candidate_id", "")),
    )


def group_key(row: dict[str, str]) -> tuple[str, str, str]:
    return tuple(row[field] for field in GROUP_FIELDS)


def select_best(rows: list[dict[str, str]]) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    grouped: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[group_key(row)].append(row)

    selected: list[dict[str, str]] = []
    tie_rows: list[dict[str, str]] = []

    for key in sorted(grouped):
        candidates = grouped[key]
        ranked = sorted(candidates, key=sort_key)
        winner = ranked[0]
        selected.append(winner)

        best_combined = numeric(winner["combined_rmsd"], "combined_rmsd", winner["model_id"])
        tied = [
            row
            for row in candidates
            if math.isclose(
                numeric(row["combined_rmsd"], "combined_rmsd", row["model_id"]),
                best_combined,
                rel_tol=0.0,
                abs_tol=1e-12,
            )
        ]
        if len(tied) > 1:
            tie_rows.append(
                {
                    "rise_label": key[0],
                    "sidechain_variant": key[1],
                    "twist_deg": key[2],
                    "tie_count": str(len(tied)),
                    "selected_model_id": winner["model_id"],
                    "selected_candidate_id": winner["candidate_id"],
                    "tie_breaker": tie_breaker_used(winner, tied),
                    "tied_model_ids": ";".join(row["model_id"] for row in sorted(tied, key=sort_key)),
                }
            )

    return selected, tie_rows


def tie_breaker_used(winner: dict[str, str], tied: list[dict[str, str]]) -> str:
    checks = [
        ("helical_rmsd", lambda row: numeric(row["helical_rmsd"], "helical_rmsd", row["model_id"])),
        ("base_rmsd", lambda row: numeric(row["base_rmsd"], "base_rmsd", row["model_id"])),
        (
            "missing_peak_count",
            lambda row: int(numeric(row["missing_peak_count"], "missing_peak_count", row["model_id"])),
        ),
        ("candidate_numeric_index", lambda row: candidate_index(row["candidate_id"])),
    ]
    remaining = list(tied)
    for name, getter in checks:
        best = min(getter(row) for row in remaining)
        narrowed = [row for row in remaining if getter(row) == best]
        if len(narrowed) < len(remaining):
            return name
        remaining = narrowed
    return "model_id_order_after_equal_tie_breakers"


def output_fields(input_fields: list[str]) -> list[str]:
    fields = list(REQUIRED_OUTPUT_FIELDS)
    seen = set(fields)
    for name in input_fields:
        if name in seen:
            continue
        if name.startswith(PER_BAND_ERROR_PREFIXES):
            fields.append(name)
            seen.add(name)
    for name in input_fields:
        if name not in seen:
            fields.append(name)
            seen.add(name)
    return fields


def validate_selection(selected: list[dict[str, str]], all_rows: list[dict[str, str]]) -> None:
    if len(selected) != 20:
        raise ValueError(f"Expected 20 selected rows, found {len(selected)}")
    keys = [group_key(row) for row in selected]
    if len(keys) != len(set(keys)):
        raise ValueError("Selected rows contain duplicate group keys")

    by_group: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in all_rows:
        by_group[group_key(row)].append(row)

    for selected_row in selected:
        selected_key = group_key(selected_row)
        selected_score = numeric(selected_row["combined_rmsd"], "combined_rmsd", selected_row["model_id"])
        lower = [
            row["model_id"]
            for row in by_group[selected_key]
            if numeric(row["combined_rmsd"], "combined_rmsd", row["model_id"]) < selected_score
        ]
        if lower:
            raise ValueError(f"{selected_row['model_id']} was selected despite lower combined_rmsd rows: {lower}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("outputs/asem_matched_twist_scan_28_32/scored_radial_profiles.csv"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("outputs/asem_matched_twist_scan_28_32/final_20_equalized_best_candidates.csv"),
    )
    parser.add_argument(
        "--tie-report",
        type=Path,
        default=Path("outputs/asem_matched_twist_scan_28_32/final_20_tie_break_report.csv"),
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    rows = read_csv(args.input)
    if not rows:
        raise ValueError(f"No scored rows found in {args.input}")

    selected, tie_rows = select_best(rows)
    validate_selection(selected, rows)

    selected = sorted(selected, key=lambda row: (row["rise_label"], row["sidechain_variant"], float(row["twist_deg"])))
    write_csv(args.output, selected, output_fields(list(rows[0])))

    tie_fields = [
        "rise_label",
        "sidechain_variant",
        "twist_deg",
        "tie_count",
        "selected_model_id",
        "selected_candidate_id",
        "tie_breaker",
        "tied_model_ids",
    ]
    write_csv(args.tie_report, tie_rows, tie_fields)

    print(f"Input scored rows: {len(rows)}")
    print(f"Selected rows: {len(selected)}")
    print(f"Tie groups: {len(tie_rows)}")
    print(f"Output: {args.output}")
    print(f"Tie report: {args.tie_report}")
    if tie_rows:
        print("Tie groups encountered:")
        for row in tie_rows:
            print(
                f"  {row['rise_label']} {row['sidechain_variant']} twist {row['twist_deg']}: "
                f"{row['tie_count']} tied, selected {row['selected_model_id']} by {row['tie_breaker']}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
