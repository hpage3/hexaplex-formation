"""Summarize and rank scored twist/rise manifest rows.

Ranking order:
  1. score_completeness descending
  2. combined_rmsd ascending
  3. helical_rmsd ascending

This keeps incomplete rows from looking artificially better merely because
difficult peaks were not observed.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


REQUIRED_COLUMNS = {
    "model_id",
    "twist_deg",
    "rise_A",
    "score_status",
    "score_completeness",
    "combined_rmsd",
    "helical_rmsd",
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


def require_columns(rows: list[dict[str, str]], path: Path) -> None:
    if not rows:
        raise ValueError(f"No rows found in {path}")
    missing = REQUIRED_COLUMNS - set(rows[0])
    if missing:
        raise ValueError(f"{path} is missing required columns: {sorted(missing)}")


def parse_float(value: str, default: float = float("inf")) -> float:
    text = (value or "").strip()
    if not text:
        return default
    return float(text)


def scored_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return [row for row in rows if row.get("score_status") == "scored"]


def rank_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    ranked = sorted(
        rows,
        key=lambda row: (
            -parse_float(row.get("score_completeness", ""), default=0.0),
            parse_float(row.get("combined_rmsd", "")),
            parse_float(row.get("helical_rmsd", "")),
        ),
    )

    output: list[dict[str, str]] = []
    for index, row in enumerate(ranked, start=1):
        out = dict(row)
        out["rank"] = str(index)
        output.append(out)

    return output


def selected_fields(rows: list[dict[str, str]]) -> list[str]:
    preferred = [
        "rank",
        "model_id",
        "twist_deg",
        "rise_A",
        "score_status",
        "observed_peak_count",
        "expected_peak_count",
        "missing_peak_count",
        "score_completeness",
        "base_rmsd",
        "helical_rmsd",
        "combined_rmsd",
        "observed_base_d_A",
        "observed_A_d_A",
        "observed_B_d_A",
        "observed_C_d_A",
        "observed_D_d_A",
        "notes",
    ]

    fields: list[str] = []
    seen: set[str] = set()

    for field in preferred:
        if any(field in row for row in rows):
            fields.append(field)
            seen.add(field)

    for row in rows:
        for field in row:
            if field not in seen:
                fields.append(field)
                seen.add(field)

    return fields


def format_value(row: dict[str, str], key: str) -> str:
    value = row.get(key, "")
    if value == "":
        return ""
    try:
        return f"{float(value):.6f}"
    except ValueError:
        return value


def write_markdown_summary(path: Path, ranked: list[dict[str, str]], top_n: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    lines: list[str] = []
    lines.append("# Twist/Rise Scored Manifest Summary")
    lines.append("")
    lines.append("Rows are ranked by score completeness descending, then combined RMSD ascending, then helical RMSD ascending.")
    lines.append("")
    lines.append(f"Scored rows: {len(ranked)}")
    lines.append("")

    if not ranked:
        lines.append("No scored rows were found.")
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return

    lines.append(f"## Top {min(top_n, len(ranked))} ranked rows")
    lines.append("")
    lines.append("| Rank | Model ID | Twist | Rise | Completeness | Missing | Base RMSD | Helical RMSD | Combined RMSD |")
    lines.append("|---:|---|---:|---:|---:|---:|---:|---:|---:|")

    for row in ranked[:top_n]:
        lines.append(
            "| "
            + " | ".join(
                [
                    row.get("rank", ""),
                    row.get("model_id", ""),
                    format_value(row, "twist_deg"),
                    format_value(row, "rise_A"),
                    format_value(row, "score_completeness"),
                    row.get("missing_peak_count", ""),
                    format_value(row, "base_rmsd"),
                    format_value(row, "helical_rmsd"),
                    format_value(row, "combined_rmsd"),
                ]
            )
            + " |"
        )

    lines.append("")
    lines.append("## Interpretation guardrail")
    lines.append("")
    lines.append(
        "A lower RMSD should be interpreted together with score completeness. "
        "Rows with missing peaks may look artificially favorable because difficult target windows were not scored."
    )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def summarize_manifest(input_path: Path, ranked_csv: Path, summary_md: Path, top_n: int) -> list[dict[str, str]]:
    rows = read_csv(input_path)
    require_columns(rows, input_path)
    ranked = rank_rows(scored_rows(rows))
    fields = selected_fields(ranked or rows)
    if "rank" not in fields:
        fields = ["rank"] + fields
    write_csv(ranked_csv, ranked, fields)
    write_markdown_summary(summary_md, ranked, top_n)
    return ranked


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Scored twist/rise manifest CSV.",
    )
    parser.add_argument(
        "--ranked-csv",
        type=Path,
        default=Path("outputs/twist_rise_scan_reports/scored_manifest_ranked.csv"),
    )
    parser.add_argument(
        "--summary-md",
        type=Path,
        default=Path("outputs/twist_rise_scan_reports/scored_manifest_summary.md"),
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=20,
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.top_n <= 0:
        raise ValueError("--top-n must be positive")

    ranked = summarize_manifest(
        input_path=args.input,
        ranked_csv=args.ranked_csv,
        summary_md=args.summary_md,
        top_n=args.top_n,
    )
    print(f"Wrote {len(ranked)} ranked scored rows to {args.ranked_csv}")
    print(f"Wrote summary to {args.summary_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
