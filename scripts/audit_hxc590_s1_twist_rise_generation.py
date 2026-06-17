#!/usr/bin/env python3
"""Audit whether HXC590 S1 twist/rise sensitivity candidates can be generated."""

from __future__ import annotations

import argparse
import csv
import importlib.util
from pathlib import Path


TWISTS_DEG = [24, 26, 28, 30, 32, 34, 36]
RISES_A = [3.20, 3.30, 3.35, 3.40, 3.50, 3.60]

DEFAULT_SENSITIVITY_CSV = Path("outputs/metrics/hxc590_s1_twist_rise_sensitivity.csv")
DEFAULT_GAP_CSV = Path("outputs/metrics/hxc590_s1_twist_rise_generation_gap.csv")
DEFAULT_REPORT = Path("outputs/reports/hxc590_s1_twist_rise_generation_gap_report.md")

CSV_COLUMNS = [
    "candidate_id",
    "candidate_family",
    "parameter_value",
    "model_path",
    "profile_path",
    "status",
    "reason",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sensitivity-csv", type=Path, default=DEFAULT_SENSITIVITY_CSV)
    parser.add_argument("--gap-csv", type=Path, default=DEFAULT_GAP_CSV)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    return parser.parse_args()


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def pnab_importable() -> bool:
    return importlib.util.find_spec("pnab") is not None


def current_model_yaml_files(root: Path) -> list[Path]:
    yaml_files = list(root.glob("*.yaml")) + list(root.glob("*.yml"))
    candidates = []
    for path in yaml_files:
        try:
            text = path.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        if "HelicalParameters" in text and "h_twist" in text and "CYP" in text:
            candidates.append(path)
    return candidates


def status_for_pair(model_path: Path, profile_path: Path, unavailable_reason: str) -> tuple[str, str]:
    model_exists = (repo_root() / model_path).is_file()
    profile_exists = (repo_root() / profile_path).is_file()
    if model_exists and profile_exists:
        return "located", "Coordinate and radial-profile files are present."
    if model_exists:
        return "located_needs_profile", "Coordinate file is present but no matching radial profile was found."
    return "unavailable", unavailable_reason


def twist_rows() -> list[dict[str, str]]:
    yaml_candidates = current_model_yaml_files(repo_root())
    pnab_note = "pNAB is importable." if pnab_importable() else "python cannot import pnab."
    current_yaml_note = (
        "current-model baseline YAML candidate found."
        if yaml_candidates
        else "no current-model baseline YAML with HelicalParameters.h_twist and CYP content was found."
    )
    unavailable_reason = (
        "No coordinate/profile files were found for this twist. The repo contains a pNAB twist helper, "
        f"but {current_yaml_note} Also, {pnab_note}"
    )
    rows = []
    for twist in TWISTS_DEG:
        if twist == 30:
            model_path = Path("outputs/intermediates/ai_candidate_inputs/full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain.pdb")
            profile_path = Path("outputs/mini_hexaplex/radial_profiles/full_length_baseline_radial.csv")
            status, reason = status_for_pair(
                model_path,
                profile_path,
                "The full-length 30-degree baseline was expected but not found.",
            )
        else:
            model_path = Path(f"outputs/length_twist_diffraction/structures/full_length_twist_{twist}.pdb")
            profile_path = Path(f"outputs/mini_hexaplex/radial_profiles/full_length_twist_{twist}_radial.csv")
            status, reason = status_for_pair(model_path, profile_path, unavailable_reason)
        rows.append(
            {
                "candidate_id": f"full_length_twist_{twist}",
                "candidate_family": "full_length_twist_variant",
                "parameter_value": f"{twist} deg",
                "model_path": str(model_path),
                "profile_path": str(profile_path),
                "status": status,
                "reason": reason,
            }
        )
    return rows


def rise_token(rise: float) -> str:
    return f"{rise:.2f}".replace(".", "p")


def rise_rows() -> list[dict[str, str]]:
    unavailable_reason = (
        "No coordinate/profile files were found for this rise. No safe existing rise-generation workflow "
        "or audited stack-axis transform was found in this repo."
    )
    rows = []
    for rise in RISES_A:
        token = rise_token(rise)
        model_path = Path(f"outputs/rise_variants/structures/hexaflex_rise_{token}.pdb")
        profile_path = Path(f"outputs/rise_variants/radial_profiles/hexaflex_rise_{token}_radial.csv")
        status, reason = status_for_pair(model_path, profile_path, unavailable_reason)
        rows.append(
            {
                "candidate_id": f"rise_{token}_synthetic_control",
                "candidate_family": "rise_variant",
                "parameter_value": f"{rise:.2f} A",
                "model_path": str(model_path),
                "profile_path": str(profile_path),
                "status": status,
                "reason": reason,
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def markdown_table(rows: list[dict[str, str]], columns: list[str]) -> list[str]:
    lines = ["| " + " | ".join(columns) + " |", "|" + "|".join("---" for _ in columns) + "|"]
    for row in rows:
        lines.append("| " + " | ".join(row.get(column, "") for column in columns) + " |")
    return lines


def write_report(path: Path, rows: list[dict[str, str]]) -> None:
    twist = [row for row in rows if row["candidate_family"] == "full_length_twist_variant"]
    rise = [row for row in rows if row["candidate_family"] == "rise_variant"]
    current_yaml = current_model_yaml_files(repo_root())
    lines = [
        "# HXC590 S1 Twist/Rise Generation Gap Audit",
        "",
        "## Purpose",
        "",
        "This audit checks whether the repository can safely generate the missing HXC590 S1 twist and rise sensitivity candidates before rerunning the powder falsification screen.",
        "",
        "This remains a falsification-style screen, not a definitive phase assignment.",
        "",
        "Synthetic twist/rise variants are controls for diffraction sensitivity, not chemically optimized structures.",
        "",
        "## Generator inspection",
        "",
        "- Existing pNAB twist helper: `scripts/generate_pnab_twist_variants.py`.",
        f"- pNAB Python import status: {'available' if pnab_importable() else 'not available'}.",
        f"- Current-model pNAB baseline YAML candidates found: {len(current_yaml)}.",
        "- Existing length/twist sensitivity workflow keeps non-30-degree full-length twist rows as pending placeholders when the official builder inputs are absent.",
        "- No audited rise-generation workflow or stack-axis rise-transform workflow was found for the requested rise scan.",
        "",
        "## Twist candidate status",
        "",
    ]
    lines.extend(markdown_table(twist, ["candidate_id", "parameter_value", "status", "reason"]))
    lines.extend(["", "## Rise candidate status", ""])
    lines.extend(markdown_table(rise, ["candidate_id", "parameter_value", "status", "reason"]))
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The full-length 30-degree baseline is available, but the missing non-30-degree twist variants remain unavailable in this checkout because the current-model pNAB baseline YAML and matching runtime inputs are absent.",
            "",
            "Rise variants remain unavailable because the repository does not contain a safe, audited rise-generation path for the current candidate model.",
            "",
            "If nearby twist or rise variants survive under current tolerances, the current powder peak list supports the conformation family but does not uniquely determine those parameters.",
            "",
            "Because nearby twist and rise variants are not currently generated or profiled, they cannot be used to refine or falsify the twist/rise parameters in this pass.",
            "",
            "## Outputs",
            "",
            f"- `{DEFAULT_SENSITIVITY_CSV}`",
            f"- `{DEFAULT_GAP_CSV}`",
            f"- `{DEFAULT_REPORT}`",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(args: argparse.Namespace) -> dict[str, int]:
    rows = twist_rows() + rise_rows()
    write_csv(args.sensitivity_csv, rows)
    write_csv(args.gap_csv, [row for row in rows if row["status"] != "located"])
    write_report(args.report, rows)
    return {
        "rows": len(rows),
        "located": sum(1 for row in rows if row["status"] == "located"),
        "unavailable": sum(1 for row in rows if row["status"] == "unavailable"),
    }


def main() -> None:
    result = run(parse_args())
    print(f"Wrote {result['rows']} twist/rise candidate audit rows")
    print(f"Located candidates with profiles: {result['located']}")
    print(f"Unavailable candidates: {result['unavailable']}")


if __name__ == "__main__":
    main()
