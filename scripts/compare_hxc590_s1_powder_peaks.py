#!/usr/bin/env python3
"""Compare HXC590 S1 powder peak windows with existing model radial profiles."""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path


SAMPLE_ID = "HXC590_S1_TM_TC_without_salt_powder"

DEFAULT_EXPERIMENTAL_PEAKS = Path("inputs/experimental/hxc590_s1_powder_peaks.csv")
DEFAULT_TARGETS_CSV = Path("outputs/metrics/hxc590_s1_powder_peak_targets.csv")
DEFAULT_SCORES_CSV = Path("outputs/metrics/hxc590_s1_powder_peak_match_scores.csv")
DEFAULT_REPORT = Path("outputs/reports/hxc590_s1_powder_peak_comparison_report.md")
DEFAULT_PLOT_DIR = Path("outputs/plots/hxc590_s1_powder_peak_comparison")

WINDOW_HALF_WIDTHS_A = {
    3.35: 0.06,
    3.71: 0.08,
    3.90: 0.08,
    4.33: 0.10,
    5.50: 0.15,
    6.50: 0.20,
    7.30: 0.25,
    7.90: 0.25,
}

TARGET_COLUMNS = [
    "sample_id",
    "target_id",
    "d_angstrom",
    "q_inv_angstrom",
    "d_window_half_width_angstrom",
    "d_min_angstrom",
    "d_max_angstrom",
    "q_min_inv_angstrom",
    "q_max_inv_angstrom",
    "confidence",
    "diagnostic_window",
    "note",
]

SCORE_COLUMNS = [
    "candidate_id",
    "candidate_label",
    "model_path",
    "profile_path",
    "target_id",
    "target_d_angstrom",
    "target_q_inv_angstrom",
    "confidence",
    "diagnostic_window",
    "window_half_width_angstrom",
    "matched",
    "matched_d_angstrom",
    "matched_q_inv_angstrom",
    "abs_d_error_angstrom",
    "fractional_d_error",
    "window_point_count",
    "window_max_intensity",
    "candidate_match_count",
    "candidate_diagnostic_match_count",
    "candidate_mean_abs_d_error_angstrom",
    "candidate_mean_fractional_d_error",
    "candidate_rank_score",
]

SUMMARY_COLUMNS = [
    "candidate_id",
    "candidate_label",
    "model_path",
    "profile_path",
    "match_count",
    "diagnostic_match_count",
    "mean_abs_d_error_angstrom",
    "mean_fractional_d_error",
    "rank_score",
]


@dataclass(frozen=True)
class PeakTarget:
    sample_id: str
    target_id: str
    d_angstrom: float
    q_inv_angstrom: float
    window_half_width: float
    confidence: str
    note: str

    @property
    def d_min(self) -> float:
        return self.d_angstrom - self.window_half_width

    @property
    def d_max(self) -> float:
        return self.d_angstrom + self.window_half_width

    @property
    def q_min(self) -> float:
        return q_from_d(self.d_max)

    @property
    def q_max(self) -> float:
        return q_from_d(self.d_min)

    @property
    def diagnostic_window(self) -> bool:
        return self.confidence != "lower"


@dataclass(frozen=True)
class CandidateModel:
    candidate_id: str
    candidate_label: str
    model_path: Path
    profile_path: Path
    category: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--experimental-peaks", type=Path, default=DEFAULT_EXPERIMENTAL_PEAKS)
    parser.add_argument("--targets-csv", type=Path, default=DEFAULT_TARGETS_CSV)
    parser.add_argument("--scores-csv", type=Path, default=DEFAULT_SCORES_CSV)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--plot-dir", type=Path, default=DEFAULT_PLOT_DIR)
    parser.add_argument("--include-nick-8hexad", action="store_true", help="Include Nick's provided Hexaplex_8Hexads.xyz 8-hexad candidate.")
    parser.add_argument("--include-ideal-hexaplex", action="store_true", help="Include the ideal antiparallel 30-degree Hexaplex_AntiParallel_30deg_Ideal.pdb model.")
    return parser.parse_args()


def q_from_d(d_angstrom: float) -> float:
    if d_angstrom <= 0:
        raise ValueError("d-spacing must be greater than zero")
    return 2.0 * math.pi / d_angstrom


def format_float(value: float | None, digits: int = 6) -> str:
    if value is None or not math.isfinite(value):
        return ""
    return f"{value:.{digits}f}"


def target_id_for_d(d_angstrom: float) -> str:
    return f"d_{str(f'{d_angstrom:.2f}').replace('.', 'p')}"


def truthy_text(value: str) -> bool:
    return value.strip().lower() in {"yes", "true", "1", "diagnostic", "primary"}


def read_experimental_targets(path: Path) -> list[PeakTarget]:
    if not path.exists():
        raise FileNotFoundError(f"Experimental peak file not found: {path}")
    targets: list[PeakTarget] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            d_value = float(row.get("d_angstrom") or row.get("distance_a") or row["distance_angstrom"])
            q_text = row.get("q_inv_angstrom") or row.get("q_a_inv") or row.get("q_Ainv")
            q_value = float(q_text) if q_text else q_from_d(d_value)
            expected_q = q_from_d(d_value)
            if abs(q_value - expected_q) > 5e-6:
                raise ValueError(f"q conversion mismatch for d={d_value}: got {q_value}, expected {expected_q:.6f}")
            width_text = row.get("d_window_half_width_angstrom") or row.get("window_half_width_a")
            rounded_d = round(d_value, 2)
            if width_text:
                window_half_width = float(width_text)
            elif rounded_d in WINDOW_HALF_WIDTHS_A:
                window_half_width = WINDOW_HALF_WIDTHS_A[rounded_d]
            else:
                raise ValueError(f"No d-window tolerance configured for {d_value}")
            confidence = row.get("confidence") or ("medium_high" if truthy_text(row.get("diagnostic_role", "")) else "lower")
            targets.append(
                PeakTarget(
                    sample_id=row.get("sample_id") or SAMPLE_ID,
                    target_id=row.get("target_id") or row.get("peak_label") or target_id_for_d(d_value),
                    d_angstrom=d_value,
                    q_inv_angstrom=q_value,
                    window_half_width=window_half_width,
                    confidence=confidence,
                    note=row.get("note") or row.get("notes", ""),
                )
            )
    return sorted(targets, key=lambda target: target.d_angstrom, reverse=True)


def write_targets_csv(path: Path, targets: list[PeakTarget]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=TARGET_COLUMNS, lineterminator="\n")
        writer.writeheader()
        for target in targets:
            writer.writerow(
                {
                    "sample_id": target.sample_id,
                    "target_id": target.target_id,
                    "d_angstrom": format_float(target.d_angstrom, 2),
                    "q_inv_angstrom": format_float(target.q_inv_angstrom),
                    "d_window_half_width_angstrom": format_float(target.window_half_width, 2),
                    "d_min_angstrom": format_float(target.d_min, 2),
                    "d_max_angstrom": format_float(target.d_max, 2),
                    "q_min_inv_angstrom": format_float(target.q_min),
                    "q_max_inv_angstrom": format_float(target.q_max),
                    "confidence": target.confidence,
                    "diagnostic_window": "yes" if target.diagnostic_window else "no",
                    "note": target.note,
                }
            )


def default_candidate_models(include_nick_8hexad: bool = False, include_ideal_hexaplex: bool = False) -> list[CandidateModel]:
    candidates = [
        CandidateModel(
            "central6_units_30deg",
            "6-unit formed endpoint, central 30-degree model",
            Path("outputs/mini_hexaplex/structures/mini_hexaplex_central6_units.pdb"),
            Path("outputs/mini_hexaplex/radial_profiles/central6_units_radial.csv"),
            "central_formed_endpoint",
        ),
        CandidateModel(
            "central7_units_30deg",
            "7-unit formed endpoint, central 30-degree model",
            Path("outputs/mini_hexaplex/structures/mini_hexaplex_central7_units.pdb"),
            Path("outputs/mini_hexaplex/radial_profiles/central7_units_radial.csv"),
            "central_formed_endpoint",
        ),
        CandidateModel(
            "central8_units_30deg",
            "8-unit formed endpoint, central 30-degree model",
            Path("outputs/mini_hexaplex/structures/mini_hexaplex_central8_units.pdb"),
            Path("outputs/mini_hexaplex/radial_profiles/central8_units_radial.csv"),
            "central_formed_endpoint",
        ),
        CandidateModel(
            "central12_units_30deg",
            "12-unit central 30-degree model",
            Path("outputs/mini_hexaplex/structures/mini_hexaplex_central12_units.pdb"),
            Path("outputs/mini_hexaplex/radial_profiles/central12_units_radial.csv"),
            "central_formed_endpoint",
        ),
        CandidateModel(
            "full_length_twist_30",
            "Full-length anti-parallel 30-degree model",
            Path("outputs/intermediates/ai_candidate_inputs/full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain.pdb"),
            Path("outputs/mini_hexaplex/radial_profiles/full_length_baseline_radial.csv"),
            "full_length_30deg",
        ),
    ]
    for scale in ["0p85", "0p90", "0p95", "1p00", "1p05", "1p10", "1p15"]:
        candidates.append(
            CandidateModel(
                f"base_length_scale_{scale}",
                f"Base-length variant scale {scale.replace('p', '.')}",
                Path(f"outputs/base_length_sweep/structures/hexaplex_base_length_scale_{scale}.pdb"),
                Path(f"outputs/base_length_sweep/radial_profiles/hexaplex_base_length_scale_{scale}_radial.csv"),
                "base_length_variant",
            )
        )
    if include_nick_8hexad:
        candidates.append(
            CandidateModel(
                "nick_hexaplex_8hexads",
                "Nick-provided Hexaplex_8Hexads.xyz 8-hexad candidate",
                Path("inputs/candidates/nick_hexaplex_8hexads.xyz"),
                Path("outputs/metrics/hxc590_s1_nick_hexaplex_8hexads_profile.csv"),
                "nick_provided_8hexad",
            )
        )
    if include_ideal_hexaplex:
        candidates.append(
            CandidateModel(
                "ideal_antiparallel_30deg_hexaplex",
                "Ideal antiparallel 30-degree hexaplex model (Hexaplex_AntiParallel_30deg_Ideal.pdb)",
                Path("inputs/candidates/ideal_antiparallel_30deg_hexaplex.pdb"),
                Path("outputs/metrics/hxc590_s1_ideal_antiparallel_30deg_hexaplex_profile.csv"),
                "ideal_antiparallel_30deg_hexaplex",
            )
        )
    return [candidate for candidate in candidates if candidate.profile_path.exists()]


def unavailable_candidate_notes() -> list[str]:
    return [
        "Length/twist manifest rows for 24, 26, 28, 32, 34, and 36 degree full-length twists are marked pending locally and were not scored because their coordinate/profile files are not present.",
        "No local 15 or 26.7 degree coordinate/profile outputs were found in this branch.",
    ]


def read_profile(path: Path) -> list[dict[str, float]]:
    if not path.exists():
        raise FileNotFoundError(f"Candidate profile not found: {path}")
    rows: list[dict[str, float]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            d_text = row.get("d_A") or row.get("d_center_angstrom")
            q_text = row.get("q_Ainv") or row.get("q_center_inv_angstrom")
            intensity_text = row.get("intensity") or row.get("mean_intensity") or row.get("intensity_mean")
            if not d_text or not q_text or not intensity_text:
                continue
            d_value = float(d_text)
            q_value = float(q_text)
            intensity = float(intensity_text)
            if d_value > 0 and q_value > 0 and math.isfinite(intensity):
                rows.append({"d": d_value, "q": q_value, "intensity": intensity})
    return sorted(rows, key=lambda row: row["d"])


def local_feature_rows(profile_rows: list[dict[str, float]]) -> list[dict[str, float]]:
    """Return local maxima so matching asks for a nearby profile feature."""
    if len(profile_rows) < 3:
        return profile_rows
    features: list[dict[str, float]] = []
    for index in range(1, len(profile_rows) - 1):
        previous_i = profile_rows[index - 1]["intensity"]
        current_i = profile_rows[index]["intensity"]
        next_i = profile_rows[index + 1]["intensity"]
        if current_i >= previous_i and current_i >= next_i and (current_i > previous_i or current_i > next_i):
            features.append(profile_rows[index])
    return features


def match_target_to_profile(target: PeakTarget, profile_rows: list[dict[str, float]]) -> dict[str, object]:
    feature_rows = local_feature_rows(profile_rows)
    window_rows = [row for row in feature_rows if target.d_min <= row["d"] <= target.d_max]
    nearest = min(profile_rows, key=lambda row: abs(row["d"] - target.d_angstrom)) if profile_rows else None
    if not window_rows:
        return {
            "matched": False,
            "matched_row": nearest,
            "point_count": 0,
            "max_intensity": None,
            "abs_error": abs(nearest["d"] - target.d_angstrom) if nearest else None,
            "fractional_error": abs(nearest["d"] - target.d_angstrom) / target.d_angstrom if nearest else None,
        }
    # Choose the closest local feature as the spacing match. Intensity is reported but not used for ranking.
    matched = min(window_rows, key=lambda row: abs(row["d"] - target.d_angstrom))
    max_intensity = max(row["intensity"] for row in window_rows)
    abs_error = abs(matched["d"] - target.d_angstrom)
    return {
        "matched": True,
        "matched_row": matched,
        "point_count": len(window_rows),
        "max_intensity": max_intensity,
        "abs_error": abs_error,
        "fractional_error": abs_error / target.d_angstrom,
    }


def summarize_candidate(matches: list[dict[str, object]], targets: list[PeakTarget]) -> dict[str, float | int | None]:
    matched = [match for match in matches if match["matched"]]
    diagnostic = [match for match, target in zip(matches, targets) if match["matched"] and target.diagnostic_window]
    abs_errors = [float(match["abs_error"]) for match in matched if match["abs_error"] is not None]
    frac_errors = [float(match["fractional_error"]) for match in matched if match["fractional_error"] is not None]
    mean_abs = sum(abs_errors) / len(abs_errors) if abs_errors else None
    mean_frac = sum(frac_errors) / len(frac_errors) if frac_errors else None
    rank_score = len(matched) + 0.5 * len(diagnostic) - (mean_frac or 1.0)
    return {
        "match_count": len(matched),
        "diagnostic_match_count": len(diagnostic),
        "mean_abs_error": mean_abs,
        "mean_fractional_error": mean_frac,
        "rank_score": rank_score,
    }


def score_candidates(targets: list[PeakTarget], candidates: list[CandidateModel]) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    score_rows: list[dict[str, str]] = []
    summary_rows: list[dict[str, str]] = []
    for candidate in candidates:
        profile_rows = read_profile(candidate.profile_path)
        matches = [match_target_to_profile(target, profile_rows) for target in targets]
        summary = summarize_candidate(matches, targets)
        summary_row = {
            "candidate_id": candidate.candidate_id,
            "candidate_label": candidate.candidate_label,
            "model_path": str(candidate.model_path),
            "profile_path": str(candidate.profile_path),
            "match_count": str(summary["match_count"]),
            "diagnostic_match_count": str(summary["diagnostic_match_count"]),
            "mean_abs_d_error_angstrom": format_float(summary["mean_abs_error"]),
            "mean_fractional_d_error": format_float(summary["mean_fractional_error"]),
            "rank_score": format_float(summary["rank_score"]),
        }
        summary_rows.append(summary_row)
        for target, match in zip(targets, matches):
            matched_row = match["matched_row"]
            assert matched_row is None or isinstance(matched_row, dict)
            score_rows.append(
                {
                    "candidate_id": candidate.candidate_id,
                    "candidate_label": candidate.candidate_label,
                    "model_path": str(candidate.model_path),
                    "profile_path": str(candidate.profile_path),
                    "target_id": target.target_id,
                    "target_d_angstrom": format_float(target.d_angstrom, 2),
                    "target_q_inv_angstrom": format_float(target.q_inv_angstrom),
                    "confidence": target.confidence,
                    "diagnostic_window": "yes" if target.diagnostic_window else "no",
                    "window_half_width_angstrom": format_float(target.window_half_width, 2),
                    "matched": "yes" if match["matched"] else "no",
                    "matched_d_angstrom": format_float(matched_row["d"] if matched_row else None),
                    "matched_q_inv_angstrom": format_float(matched_row["q"] if matched_row else None),
                    "abs_d_error_angstrom": format_float(match["abs_error"]),
                    "fractional_d_error": format_float(match["fractional_error"]),
                    "window_point_count": str(match["point_count"]),
                    "window_max_intensity": format_float(match["max_intensity"]),
                    "candidate_match_count": summary_row["match_count"],
                    "candidate_diagnostic_match_count": summary_row["diagnostic_match_count"],
                    "candidate_mean_abs_d_error_angstrom": summary_row["mean_abs_d_error_angstrom"],
                    "candidate_mean_fractional_d_error": summary_row["mean_fractional_d_error"],
                    "candidate_rank_score": summary_row["rank_score"],
                }
            )
    summary_rows.sort(
        key=lambda row: (
            -float(row["rank_score"] or "0"),
            -int(row["diagnostic_match_count"]),
            -int(row["match_count"]),
            float(row["mean_abs_d_error_angstrom"] or "999"),
            row["candidate_id"],
        )
    )
    return score_rows, summary_rows


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def markdown_table(rows: list[dict[str, str]], columns: list[str], limit: int | None = None) -> list[str]:
    selected = rows[:limit] if limit is not None else rows
    lines = ["| " + " | ".join(columns) + " |", "|" + "|".join("---" for _ in columns) + "|"]
    for row in selected:
        lines.append("| " + " | ".join(row.get(column, "") for column in columns) + " |")
    return lines


def write_coverage_plot(path: Path, summary_rows: list[dict[str, str]]) -> None:
    rows = summary_rows[:10]
    width = 980
    height = 420
    left = 70
    bottom = 310
    bar_w = 28
    gap = 18
    max_matches = 8
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="490" y="28" text-anchor="middle" font-family="Arial" font-size="17">HXC590 S1 powder peak-window coverage</text>',
        f'<line x1="{left}" y1="70" x2="{left}" y2="{bottom}" stroke="#333"/>',
        f'<line x1="{left}" y1="{bottom}" x2="{width - 20}" y2="{bottom}" stroke="#333"/>',
    ]
    for index, row in enumerate(rows):
        x = left + 24 + index * (bar_w + gap)
        match_count = int(row["match_count"])
        diagnostic_count = int(row["diagnostic_match_count"])
        h = match_count / max_matches * 220
        diagnostic_h = diagnostic_count / 4 * 220
        parts.append(f'<rect x="{x}" y="{bottom - h:.2f}" width="{bar_w}" height="{h:.2f}" fill="#8fb3ff"/>')
        parts.append(f'<rect x="{x + 8}" y="{bottom - diagnostic_h:.2f}" width="12" height="{diagnostic_h:.2f}" fill="#f28e2b"/>')
        label = row["candidate_id"].replace("_", " ")
        parts.append(f'<text x="{x + 12}" y="{bottom + 14}" text-anchor="end" font-family="Arial" font-size="8" transform="rotate(-42 {x + 12},{bottom + 14})">{label}</text>')
    parts.append('<rect x="740" y="52" width="14" height="14" fill="#8fb3ff"/>')
    parts.append('<text x="760" y="64" font-family="Arial" font-size="12">all peak windows, out of 8</text>')
    parts.append('<rect x="740" y="74" width="14" height="14" fill="#f28e2b"/>')
    parts.append('<text x="760" y="86" font-family="Arial" font-size="12">diagnostic windows, out of 4</text>')
    parts.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def write_target_match_plot(path: Path, score_rows: list[dict[str, str]], best_candidate_id: str) -> None:
    rows = [row for row in score_rows if row["candidate_id"] == best_candidate_id]
    width = 760
    height = 360
    left = 80
    bottom = 280
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width / 2:.0f}" y="28" text-anchor="middle" font-family="Arial" font-size="16">Best candidate d-spacing errors</text>',
        f'<line x1="{left}" y1="60" x2="{left}" y2="{bottom}" stroke="#333"/>',
        f'<line x1="{left}" y1="{bottom}" x2="{width - 30}" y2="{bottom}" stroke="#333"/>',
        f'<line x1="{left}" y1="{bottom}" x2="{width - 30}" y2="{bottom}" stroke="#777" stroke-dasharray="4 3"/>',
    ]
    max_error = max([float(row["abs_d_error_angstrom"] or "0") for row in rows] + [0.25])
    for index, row in enumerate(rows):
        x = left + 35 + index * 75
        error = float(row["abs_d_error_angstrom"] or "0")
        h = error / max_error * 180
        color = "#59a14f" if row["matched"] == "yes" else "#e15759"
        parts.append(f'<rect x="{x}" y="{bottom - h:.2f}" width="30" height="{h:.2f}" fill="{color}"/>')
        parts.append(f'<text x="{x + 15}" y="{bottom + 16}" text-anchor="middle" font-family="Arial" font-size="10">{row["target_d_angstrom"]}</text>')
    parts.append('<text x="22" y="168" font-family="Arial" font-size="11" transform="rotate(-90 22,168)">absolute d error (A)</text>')
    parts.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def write_report(
    path: Path,
    targets: list[PeakTarget],
    candidates: list[CandidateModel],
    score_rows: list[dict[str, str]],
    summary_rows: list[dict[str, str]],
    plot_dir: Path,
) -> None:
    corrected = any("corrected" in target.note.lower() or "corrected" in target.sample_id.lower() for target in targets)
    diagnostic_rows = [
        row
        for row in score_rows
        if row["candidate_id"] == summary_rows[0]["candidate_id"] and row["diagnostic_window"] == "yes"
    ]
    diagnostic_summary = [
        f"{row['target_d_angstrom']} A = {row['matched']}"
        for row in diagnostic_rows
    ]
    purpose_text = (
        "This report compares John Bacsa / Nick's corrected HXC590 S1 powder distance scale against existing simulated radial profiles for hexaflex / stacked-hexad candidate coordinate models. It is a peak-window diagnostic-spacing comparison, not a full Rietveld refinement or phase-refinement workflow."
        if corrected
        else "This report compares John Bacsa's approximate HXC590 S1 powder peak list against existing simulated radial profiles for hexaflex / stacked-hexad candidate coordinate models. It is a peak-window diagnostic-spacing comparison, not a full Rietveld refinement or phase-refinement workflow."
    )
    experimental_context = (
        "John Bacsa provided a corrected distance scale for the HXC590 S1 powder data. The correction shifts the target peak positions modestly and places the base-stacking feature close to the expected 3.4 A region. Relative intensities are preserved in the target table but are treated as approximate rather than reliable scoring constraints."
        if corrected
        else "John reported that merged raw frames show uniform powder rings rather than oriented fiber arcs, with large broad background under peaks. Peak positions are approximate and relative intensities are treated as approximate rather than reliable constraints."
    )
    interpretation_text = (
        "The corrected HXC590 S1 powder targets strengthen compatibility with a related stacked-hexad/hexaflex-like model family if the corrected diagnostic windows, including the near-3.4 A stacking feature, are reproduced by the simulated model profile."
        if corrected
        else "The HXC590 S1 powder peak list is consistent with a related stacked-hexad/hexaflex phase if the diagnostic 3.35 A, 4.33 A, and 3.7-3.9 A windows are reproduced by the simulated model profile."
    )
    lines = [
        "# HXC590 S1 Powder Peak Comparison",
        "",
        "## Purpose",
        "",
        purpose_text,
        "",
        "## Experimental Context",
        "",
        "Sample: `TM_TC without salt`, `HXC590 S1 powder`.",
        "",
        experimental_context,
        "",
        "Uniform rings are expected for an unoriented powder sample, whereas oriented fibers give arcs. Therefore, the ring-versus-arc difference is not by itself evidence for a different molecular phase.",
        "",
        "## Experimental Peaks and q Conversion",
        "",
    ]
    target_rows = [
        {
            "d_A": format_float(target.d_angstrom, 2),
            "q_Ainv": format_float(target.q_inv_angstrom),
            "window_A": f"+/-{target.window_half_width:.2f}",
            "confidence": target.confidence,
        }
        for target in targets
    ]
    lines.extend(markdown_table(target_rows, ["d_A", "q_Ainv", "window_A", "confidence"]))
    lines.extend(["", "The q conversion used `q = 2*pi/d`. Window matching was performed in d-space, with q-space bounds recorded in `outputs/metrics/hxc590_s1_powder_peak_targets.csv`.", ""])
    lines.extend(["## Candidate Models", ""])
    candidate_rows = [
        {
            "candidate_id": candidate.candidate_id,
            "category": candidate.category,
            "profile_path": str(candidate.profile_path),
        }
        for candidate in candidates
    ]
    lines.extend(markdown_table(candidate_rows, ["candidate_id", "category", "profile_path"]))
    if any(candidate.candidate_id == "nick_hexaplex_8hexads" for candidate in candidates):
        nick_row = next((row for row in summary_rows if row["candidate_id"] == "nick_hexaplex_8hexads"), None)
        central12_row = next((row for row in summary_rows if row["candidate_id"] == "central12_units_30deg"), None)
        central8_row = next((row for row in summary_rows if row["candidate_id"] == "central8_units_30deg"), None)
        lines.extend(
            [
                "",
                "## Nick 8-Hexad Candidate",
                "",
                "The `nick_hexaplex_8hexads` row is Nick's included `Hexaplex_8Hexads.xyz` 8-hexad candidate.",
            ]
        )
        if nick_row:
            lines.append(
                f"It scores with {nick_row['match_count']} of {len(targets)} corrected windows and {nick_row['diagnostic_match_count']} diagnostic windows matched."
            )
        if central12_row and nick_row:
            relation = "above" if float(nick_row["rank_score"]) > float(central12_row["rank_score"]) else "below"
            lines.append(f"By the existing rank score, Nick's 8-hexad candidate is {relation} `central12_units_30deg` in this corrected screen.")
        if central8_row and nick_row:
            relation = "above" if float(nick_row["rank_score"]) > float(central8_row["rank_score"]) else "below"
            lines.append(f"By the existing rank score, Nick's 8-hexad candidate is {relation} `central8_units_30deg` in this corrected screen.")
    if any(candidate.candidate_id == "ideal_antiparallel_30deg_hexaplex" for candidate in candidates):
        ideal_row = next((row for row in summary_rows if row["candidate_id"] == "ideal_antiparallel_30deg_hexaplex"), None)
        central12_row = next((row for row in summary_rows if row["candidate_id"] == "central12_units_30deg"), None)
        central8_row = next((row for row in summary_rows if row["candidate_id"] == "central8_units_30deg"), None)
        nick8_row = next((row for row in summary_rows if row["candidate_id"] == "nick_hexaplex_8hexads"), None)
        lines.extend(
            [
                "",
                "## Ideal Antiparallel 30-Degree Hexaplex Model",
                "",
                "The `ideal_antiparallel_30deg_hexaplex` row is `Hexaplex_AntiParallel_30deg_Ideal.pdb`, labeled by ideal antiparallel 30-degree model geometry for this screen.",
                "Nick clarified that this is the file he meant when referring to the 16-mer simulation benchmark; the scoring outputs use model-provenance naming instead of that informal shorthand.",
            ]
        )
        if ideal_row:
            lines.append(
                f"It scores with {ideal_row['match_count']} of {len(targets)} corrected windows and {ideal_row['diagnostic_match_count']} diagnostic windows matched."
            )
        for label, row in [
            ("central12_units_30deg", central12_row),
            ("central8_units_30deg", central8_row),
            ("nick_hexaplex_8hexads", nick8_row),
        ]:
            if row and ideal_row:
                relation = "above" if float(ideal_row["rank_score"]) > float(row["rank_score"]) else "below"
                lines.append(f"By the existing rank score, the ideal antiparallel 30-degree hexaplex model is {relation} `{label}` in this corrected screen.")
    lines.extend(["", "Unavailable twist variants:", ""])
    lines.extend(f"- {note}" for note in unavailable_candidate_notes())
    lines.extend(
        [
            "",
            "## Scoring Method",
            "",
            "Each experimental d-spacing was assigned a conservative d-window tolerance from the prompt. A candidate was counted as matching that peak when its radial profile contained at least one local profile maximum inside the d-window. The reported matched d-spacing is the closest local feature in that window.",
            "",
            "Ranking uses peak-window coverage, diagnostic-window coverage, and spacing error. Relative intensity is not used in the ranking because the experimental background is broad and relative intensities were not treated as reliable constraints.",
            "",
            "## Best-Matching Candidates",
            "",
        ]
    )
    lines.extend(markdown_table(summary_rows, SUMMARY_COLUMNS, limit=8))
    if corrected:
        central8_row = next((row for row in summary_rows if row["candidate_id"] == "central8_units_30deg"), None)
        central8_note = (
            f" with {central8_row['match_count']} of 5 total windows and {central8_row['diagnostic_match_count']} diagnostic windows matched."
            if central8_row
            else "."
        )
        lines.extend(
            [
                "",
                "Compared with the earlier approximate-target run, the corrected five-target comparison changes the top-ranked available candidate from `central8_units_30deg` to "
                f"`{summary_rows[0]['candidate_id']}` under the existing scoring logic. This is a rank change within the available screened profiles, not a unique refined phase assignment.",
                "",
                f"The corrected near-3.4 A stacking target is included as a diagnostic window. `central8_units_30deg` remains scored in the corrected comparison{central8_note}",
            ]
        )
    lines.extend(["", "## Diagnostic Peak-Window Match Table", ""])
    lines.extend(
        markdown_table(
            diagnostic_rows,
            [
                "candidate_id",
                "target_d_angstrom",
                "matched",
                "matched_d_angstrom",
                "abs_d_error_angstrom",
                "window_point_count",
            ],
        )
    )
    best = summary_rows[0]
    best_rows = [row for row in score_rows if row["candidate_id"] == best["candidate_id"]]
    matched_by_target = {row["target_d_angstrom"]: row["matched"] == "yes" for row in best_rows}
    lines.extend(
        [
            "",
            "## Conservative Same-Phase Interpretation",
            "",
            interpretation_text,
            "",
            f"For the current top-ranked candidate, `{best['candidate_id']}`, the diagnostic windows are reproduced as follows: {', '.join(diagnostic_summary)}.",
            "",
            "A cautious answer to John's same-phase question is that the available powder peak list is compatible with a related stacked-hexad/hexaflex-like model family when these diagnostic spacing regions are matched, but the current comparison is not sufficient by itself to assign the phase.",
            "",
            "## Limitations",
            "",
            "This comparison does not establish a definitive phase assignment. The experimental background is broad, peak positions are approximate, and relative intensities were not treated as reliable constraints.",
            "",
            "The model profiles are simplified radial-profile diagnostics. They are useful for comparing spacing windows, but they are not a full powder diffraction refinement, not a full fiber diffraction simulator, and not a calibrated experimental fit.",
            "",
            "Improved comparison would benefit from calibrated q values, detector radii or pixel positions, beam center, wavelength, sample-to-detector distance, background subtraction details, and uncertainty estimates for each peak.",
            "",
            "## Outputs",
            "",
            "- `outputs/metrics/hxc590_s1_powder_peak_targets.csv`",
            "- `outputs/metrics/hxc590_s1_powder_peak_match_scores.csv`",
            f"- `{plot_dir / 'peak_window_coverage.svg'}`",
            f"- `{plot_dir / 'best_candidate_d_errors.svg'}`",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(args: argparse.Namespace) -> dict[str, object]:
    targets = read_experimental_targets(args.experimental_peaks)
    candidates = default_candidate_models(
        include_nick_8hexad=args.include_nick_8hexad,
        include_ideal_hexaplex=args.include_ideal_hexaplex,
    )
    if not candidates:
        raise ValueError("No candidate profiles found")
    write_targets_csv(args.targets_csv, targets)
    score_rows, summary_rows = score_candidates(targets, candidates)
    write_csv(args.scores_csv, score_rows, SCORE_COLUMNS)
    args.plot_dir.mkdir(parents=True, exist_ok=True)
    write_coverage_plot(args.plot_dir / "peak_window_coverage.svg", summary_rows)
    write_target_match_plot(args.plot_dir / "best_candidate_d_errors.svg", score_rows, summary_rows[0]["candidate_id"])
    write_report(args.report, targets, candidates, score_rows, summary_rows, args.plot_dir)
    return {
        "targets": len(targets),
        "candidates": len(candidates),
        "score_rows": len(score_rows),
        "best_candidate": summary_rows[0]["candidate_id"],
        "best_match_count": summary_rows[0]["match_count"],
        "best_diagnostic_match_count": summary_rows[0]["diagnostic_match_count"],
    }


def main() -> None:
    result = run(parse_args())
    print(f"Wrote {result['targets']} powder peak targets")
    print(f"Compared {result['candidates']} candidate profiles")
    print(f"Wrote {result['score_rows']} peak-match score rows")
    print(
        f"Best candidate: {result['best_candidate']} "
        f"({result['best_match_count']} total matches; {result['best_diagnostic_match_count']} diagnostic matches)"
    )


if __name__ == "__main__":
    main()
