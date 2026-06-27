"""Summarize the final 20 Asem matched twist/rise screening candidates."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


SCORE_FIELDS = ["base_rmsd", "helical_rmsd", "combined_rmsd"]
OBSERVED_FIELDS = [
    "observed_base_d_A",
    "observed_A_d_A",
    "observed_B_d_A",
    "observed_C_d_A",
    "observed_D_d_A",
]
TARGET_D = {
    "base": 3.38,
    "A": 3.8,
    "B": 4.4,
    "C": 5.65,
    "D": 7.3,
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def f(value: str) -> float:
    return float(value)


def fmt(value: float) -> str:
    if not math.isfinite(value):
        return ""
    return f"{value:.9f}"


def mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else float("nan")


def best_row(rows: list[dict[str, str]]) -> dict[str, str]:
    return min(rows, key=lambda row: f(row["combined_rmsd"]))


def tie_lookup(rows: list[dict[str, str]]) -> dict[tuple[str, str, str], dict[str, str]]:
    return {
        (row["rise_label"], row["sidechain_variant"], row["twist_deg"]): row
        for row in rows
    }


def add_band_errors(row: dict[str, str]) -> dict[str, str]:
    out = dict(row)
    for label in ["base", "A", "B", "C", "D"]:
        observed_key = f"observed_{label}_d_A"
        error_key = f"{label}_error_d_A"
        out[error_key] = fmt(f(out[observed_key]) - TARGET_D[label])
    return out


def model_selection_summary(rows: list[dict[str, str]], ties: dict[tuple[str, str, str], dict[str, str]]) -> list[dict[str, str]]:
    out_rows: list[dict[str, str]] = []
    for row in sorted(rows, key=lambda r: (r["rise_label"], r["sidechain_variant"], f(r["twist_deg"]))):
        out = add_band_errors(row)
        tie = ties.get((row["rise_label"], row["sidechain_variant"], row["twist_deg"]), {})
        out["tie_count"] = tie.get("tie_count", "1")
        out["tie_breaker"] = tie.get("tie_breaker", "")
        out_rows.append(out)
    return out_rows


def summarize_groups(rows: list[dict[str, str]], group_fields: list[str]) -> list[dict[str, str]]:
    grouped: dict[tuple[str, ...], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[tuple(row[field] for field in group_fields)].append(row)

    out_rows: list[dict[str, str]] = []
    for key in sorted(grouped, key=lambda item: tuple(float(part) if part.replace(".", "", 1).isdigit() else part for part in item)):
        group_rows = grouped[key]
        best = best_row(group_rows)
        out = {field: value for field, value in zip(group_fields, key)}
        out["row_count"] = str(len(group_rows))
        out["mean_base_rmsd"] = fmt(mean([f(row["base_rmsd"]) for row in group_rows]))
        out["mean_helical_rmsd"] = fmt(mean([f(row["helical_rmsd"]) for row in group_rows]))
        out["mean_combined_rmsd"] = fmt(mean([f(row["combined_rmsd"]) for row in group_rows]))
        out["best_model_id"] = best["model_id"]
        out["best_candidate_id"] = best["candidate_id"]
        out["best_rise_label"] = best["rise_label"]
        out["best_sidechain_variant"] = best["sidechain_variant"]
        out["best_twist_deg"] = best["twist_deg"]
        out["best_base_rmsd"] = best["base_rmsd"]
        out["best_helical_rmsd"] = best["helical_rmsd"]
        out["best_combined_rmsd"] = best["combined_rmsd"]
        out["min_score_completeness"] = fmt(min(f(row["score_completeness"]) for row in group_rows))
        out["max_missing_peak_count"] = str(max(int(f(row["missing_peak_count"])) for row in group_rows))
        for label in ["base", "A", "B", "C", "D"]:
            observed_key = f"observed_{label}_d_A"
            out[f"mean_observed_{label}_d_A"] = fmt(mean([f(row[observed_key]) for row in group_rows]))
            out[f"mean_{label}_error_d_A"] = fmt(mean([f(row[observed_key]) - TARGET_D[label] for row in group_rows]))
        out_rows.append(out)
    return out_rows


def summary_fields(group_fields: list[str]) -> list[str]:
    fields = list(group_fields)
    fields.extend(
        [
            "row_count",
            "mean_base_rmsd",
            "mean_helical_rmsd",
            "mean_combined_rmsd",
            "best_model_id",
            "best_candidate_id",
            "best_rise_label",
            "best_sidechain_variant",
            "best_twist_deg",
            "best_base_rmsd",
            "best_helical_rmsd",
            "best_combined_rmsd",
            "min_score_completeness",
            "max_missing_peak_count",
        ]
    )
    for label in ["base", "A", "B", "C", "D"]:
        fields.append(f"mean_observed_{label}_d_A")
        fields.append(f"mean_{label}_error_d_A")
    return fields


def model_fields(input_fields: list[str]) -> list[str]:
    preferred = [
        "model_id",
        "rise_label",
        "rise_A",
        "twist_deg",
        "sidechain_variant",
        "candidate_id",
        "radial_csv",
        "observed_base_d_A",
        "base_error_d_A",
        "observed_A_d_A",
        "A_error_d_A",
        "observed_B_d_A",
        "B_error_d_A",
        "observed_C_d_A",
        "C_error_d_A",
        "observed_D_d_A",
        "D_error_d_A",
        "base_rmsd",
        "helical_rmsd",
        "combined_rmsd",
        "observed_peak_count",
        "missing_peak_count",
        "score_completeness",
        "tie_count",
        "tie_breaker",
    ]
    fields = list(preferred)
    seen = set(fields)
    for name in input_fields:
        if name not in seen:
            fields.append(name)
            seen.add(name)
    return fields


def validate_rows(rows: list[dict[str, str]]) -> None:
    if len(rows) != 20:
        raise ValueError(f"Expected 20 final rows, found {len(rows)}")
    required = [
        "model_id",
        "rise_label",
        "rise_A",
        "twist_deg",
        "sidechain_variant",
        "candidate_id",
        "base_rmsd",
        "helical_rmsd",
        "combined_rmsd",
        "score_completeness",
        "missing_peak_count",
        *OBSERVED_FIELDS,
    ]
    for field in required:
        missing = [row.get("model_id", "<unknown>") for row in rows if not row.get(field)]
        if missing:
            raise ValueError(f"Missing {field} for {len(missing)} selected rows: {missing[:5]}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("outputs/asem_matched_twist_scan_28_32/final_20_equalized_best_candidates.csv"),
    )
    parser.add_argument(
        "--tie-report",
        type=Path,
        default=Path("outputs/asem_matched_twist_scan_28_32/final_20_tie_break_report.csv"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs/asem_matched_twist_scan_28_32"),
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    rows = read_csv(args.input)
    validate_rows(rows)
    ties = tie_lookup(read_csv(args.tie_report)) if args.tie_report.exists() else {}

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    by_twist = summarize_groups(rows, ["twist_deg"])
    by_variant = summarize_groups(rows, ["sidechain_variant"])
    by_rise_variant = summarize_groups(rows, ["rise_label", "sidechain_variant"])
    model_summary = model_selection_summary(rows, ties)

    paths = {
        "by_twist": output_dir / "final_20_score_summary_by_twist.csv",
        "by_variant": output_dir / "final_20_score_summary_by_variant.csv",
        "by_rise_variant": output_dir / "final_20_score_summary_by_rise_variant.csv",
        "model_selection": output_dir / "final_20_model_selection_summary.csv",
    }

    write_csv(paths["by_twist"], by_twist, summary_fields(["twist_deg"]))
    write_csv(paths["by_variant"], by_variant, summary_fields(["sidechain_variant"]))
    write_csv(paths["by_rise_variant"], by_rise_variant, summary_fields(["rise_label", "sidechain_variant"]))
    write_csv(paths["model_selection"], model_summary, model_fields(list(rows[0])))

    print(f"Input final rows: {len(rows)}")
    for name, path in paths.items():
        print(f"{name}: {path}")
    print(f"Twists summarized: {', '.join(row['twist_deg'] for row in by_twist)}")
    print(f"Variants summarized: {', '.join(row['sidechain_variant'] for row in by_variant)}")
    print(
        "Rise x variant groups summarized: "
        + ", ".join(f"{row['rise_label']} {row['sidechain_variant']}" for row in by_rise_variant)
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
