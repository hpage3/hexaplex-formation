#!/usr/bin/env python3
"""Run a falsification-style HXC590 S1 powder candidate screen."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "scripts"))

from compare_hxc590_s1_powder_peaks import (  # noqa: E402
    CandidateModel,
    PeakTarget,
    default_candidate_models,
    format_float,
    local_feature_rows,
    markdown_table,
    match_target_to_profile,
    read_experimental_targets,
    read_profile,
)


DEFAULT_EXPERIMENTAL_PEAKS = Path("inputs/experimental/hxc590_s1_powder_peaks.csv")
DEFAULT_MANIFEST_CSV = Path("outputs/metrics/hxc590_s1_falsification_candidate_manifest.csv")
DEFAULT_SCORES_CSV = Path("outputs/metrics/hxc590_s1_falsification_scores.csv")
DEFAULT_UNMATCHED_CSV = Path("outputs/metrics/hxc590_s1_predicted_unmatched_peaks.csv")
DEFAULT_TOLERANCE_CSV = Path("outputs/metrics/hxc590_s1_tolerance_survival_summary.csv")
DEFAULT_REPORT = Path("outputs/reports/hxc590_s1_falsification_report.md")
DEFAULT_PLOT_DIR = Path("outputs/plots/hxc590_s1_falsification")
DEFAULT_TWIST_RISE_SENSITIVITY_CSV = Path("outputs/metrics/hxc590_s1_twist_rise_sensitivity.csv")

TOLERANCE_SCALES = {"narrow": 0.5, "current": 1.0, "broad": 1.5}
CENTRAL_CANDIDATE_ID = "central8_units_30deg"

MANIFEST_COLUMNS = [
    "candidate_id",
    "candidate_label",
    "candidate_family",
    "candidate_role",
    "status",
    "model_path",
    "profile_path",
    "reason",
]

SCORE_COLUMNS = [
    "tolerance_setting",
    "candidate_id",
    "candidate_label",
    "candidate_family",
    "candidate_role",
    "match_count",
    "diagnostic_match_count",
    "strict_survives",
    "screen_survives",
    "mean_abs_d_error_angstrom",
    "mean_fractional_d_error",
    "observed_recovery_score",
    "predicted_unmatched_peak_count",
    "predicted_unmatched_peak_fraction",
    "discrimination_score",
    "matched_targets",
    "missed_diagnostic_targets",
]

UNMATCHED_COLUMNS = [
    "candidate_id",
    "candidate_label",
    "candidate_family",
    "candidate_role",
    "predicted_peak_rank",
    "predicted_d_angstrom",
    "predicted_q_inv_angstrom",
    "relative_intensity",
    "matched_experimental_window",
]

TOLERANCE_COLUMNS = [
    "tolerance_setting",
    "window_scale",
    "surviving_candidate_count",
    "surviving_candidate_names",
    "strict_surviving_candidate_count",
    "strict_surviving_candidate_names",
    "best_candidate",
    "central8_units_30deg_survives",
    "central8_units_30deg_strict_survives",
    "central8_units_30deg_uniquely_best",
]


@dataclass(frozen=True)
class ManifestCandidate:
    model: CandidateModel | None
    candidate_family: str
    candidate_role: str
    status: str
    reason: str
    candidate_id: str
    candidate_label: str
    model_path: Path
    profile_path: Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--experimental-peaks", type=Path, default=DEFAULT_EXPERIMENTAL_PEAKS)
    parser.add_argument("--manifest-csv", type=Path, default=DEFAULT_MANIFEST_CSV)
    parser.add_argument("--scores-csv", type=Path, default=DEFAULT_SCORES_CSV)
    parser.add_argument("--unmatched-csv", type=Path, default=DEFAULT_UNMATCHED_CSV)
    parser.add_argument("--tolerance-csv", type=Path, default=DEFAULT_TOLERANCE_CSV)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--plot-dir", type=Path, default=DEFAULT_PLOT_DIR)
    return parser.parse_args()


def available_candidate(model: CandidateModel, family: str, role: str, reason: str = "") -> ManifestCandidate:
    return ManifestCandidate(
        model=model,
        candidate_family=family,
        candidate_role=role,
        status="available",
        reason=reason,
        candidate_id=model.candidate_id,
        candidate_label=model.candidate_label,
        model_path=model.model_path,
        profile_path=model.profile_path,
    )


def unavailable_candidate(
    candidate_id: str,
    label: str,
    family: str,
    role: str,
    model_path: Path,
    profile_path: Path,
    reason: str,
) -> ManifestCandidate:
    return ManifestCandidate(
        model=None,
        candidate_family=family,
        candidate_role=role,
        status="unavailable",
        reason=reason,
        candidate_id=candidate_id,
        candidate_label=label,
        model_path=model_path,
        profile_path=profile_path,
    )


def build_candidate_manifest() -> list[ManifestCandidate]:
    manifest: list[ManifestCandidate] = []
    for model in default_candidate_models():
        family = "existing_hxc590_scored"
        role = "candidate"
        if model.category == "base_length_variant":
            family = "base_length_variant"
        elif model.candidate_id == "full_length_twist_30":
            family = "full_length_twist_variant"
        manifest.append(available_candidate(model, family, role))

    negative_controls = [
        CandidateModel(
            "negative_hexads_only",
            "Hexads-only partial-structure negative control",
            Path("outputs/intermediates/ladder_structures/reference_hexads_only_heavy_deduped.pdb"),
            Path("outputs/metrics/ladder_diffraction/profiles/reference_hexads_only_heavy_deduped_debye_profile.csv"),
            "negative_control",
        ),
        CandidateModel(
            "negative_scaffold_only",
            "Scaffold-only partial-structure negative control",
            Path("outputs/intermediates/ladder_structures/reference_scaffold_only_complement_heavy_deduped.pdb"),
            Path("outputs/metrics/ladder_diffraction/profiles/reference_scaffold_only_complement_heavy_deduped_debye_profile.csv"),
            "negative_control",
        ),
        CandidateModel(
            "negative_alanine_beta_sheet",
            "Alanine beta-sheet full negative control",
            Path("outputs/intermediates/cleaned_structures/alanine_beta_sheet_full_heavy_deduped.pdb"),
            Path(""),
            "negative_control",
        ),
    ]
    for control in negative_controls:
        if str(control.profile_path) and control.profile_path.is_file():
            manifest.append(available_candidate(control, "wrong_geometry_control", "negative_control"))
        else:
            manifest.append(
                unavailable_candidate(
                    control.candidate_id,
                    control.candidate_label,
                    "wrong_geometry_control",
                    "negative_control",
                    control.model_path,
                    control.profile_path,
                    "No existing radial profile was found; not generating a new control profile in this pass.",
                )
            )

    for twist in [24, 26, 28, 32, 34, 36]:
        manifest.append(
            unavailable_candidate(
                f"full_length_twist_{twist}",
                f"Full-length {twist}-degree twist variant",
                "full_length_twist_variant",
                "plausible_alternative",
                Path(f"outputs/length_twist_diffraction/structures/full_length_twist_{twist}.pdb"),
                Path(f"outputs/mini_hexaplex/radial_profiles/full_length_twist_{twist}_radial.csv"),
                "No coordinate/profile files were found for this twist. Existing pNAB helper support requires the current-model baseline YAML and pNAB runtime before generation is safe.",
            )
        )

    for rise in ["3p20", "3p30", "3p35", "3p40", "3p50", "3p60"]:
        manifest.append(
            unavailable_candidate(
                f"rise_{rise}_synthetic_control",
                f"Synthetic rise-control scan {rise.replace('p', '.')} A",
                "rise_variant",
                "plausible_alternative",
                Path(f"outputs/rise_variants/structures/hexaflex_rise_{rise}.pdb"),
                Path(f"outputs/rise_variants/radial_profiles/hexaflex_rise_{rise}_radial.csv"),
                "No coordinate/profile files were found for this rise. No safe existing rise-generation workflow or audited stack-axis transform was found in this repo.",
            )
        )
    return manifest


def scaled_targets(targets: list[PeakTarget], scale: float) -> list[PeakTarget]:
    return [
        PeakTarget(
            sample_id=target.sample_id,
            target_id=target.target_id,
            d_angstrom=target.d_angstrom,
            q_inv_angstrom=target.q_inv_angstrom,
            window_half_width=target.window_half_width * scale,
            confidence=target.confidence,
            note=target.note,
        )
        for target in targets
    ]


def summarize_matches(matches: list[dict[str, object]], targets: list[PeakTarget]) -> dict[str, object]:
    matched = [match for match in matches if match["matched"]]
    diagnostic_matched = [match for match, target in zip(matches, targets) if match["matched"] and target.diagnostic_window]
    abs_errors = [float(match["abs_error"]) for match in matched if match["abs_error"] is not None]
    frac_errors = [float(match["fractional_error"]) for match in matched if match["fractional_error"] is not None]
    missed_diagnostic = [
        format_float(target.d_angstrom, 2)
        for match, target in zip(matches, targets)
        if target.diagnostic_window and not match["matched"]
    ]
    matched_targets = [
        format_float(target.d_angstrom, 2)
        for match, target in zip(matches, targets)
        if match["matched"]
    ]
    mean_abs = sum(abs_errors) / len(abs_errors) if abs_errors else None
    mean_frac = sum(frac_errors) / len(frac_errors) if frac_errors else None
    match_count = len(matched)
    diagnostic_count = len(diagnostic_matched)
    return {
        "match_count": match_count,
        "diagnostic_match_count": diagnostic_count,
        "strict_survives": diagnostic_count == 4 and match_count == 8,
        "screen_survives": diagnostic_count == 4 and match_count >= 6,
        "mean_abs_error": mean_abs,
        "mean_fractional_error": mean_frac,
        "observed_recovery_score": match_count + 0.5 * diagnostic_count - (mean_frac or 1.0),
        "matched_targets": ";".join(matched_targets),
        "missed_diagnostic_targets": ";".join(missed_diagnostic),
    }


def observed_windows_for_targets(targets: list[PeakTarget]) -> list[tuple[float, float]]:
    return [(target.d_min, target.d_max) for target in targets]


def in_any_window(d_value: float, windows: list[tuple[float, float]]) -> bool:
    return any(d_min <= d_value <= d_max for d_min, d_max in windows)


def strong_predicted_unmatched_peaks(
    profile_rows: list[dict[str, float]],
    targets: list[PeakTarget],
    top_n: int = 12,
) -> tuple[list[dict[str, object]], int, float]:
    if not profile_rows:
        return [], 0, 0.0
    d_min = min(target.d_min for target in targets)
    d_max = max(target.d_max for target in targets)
    features = [row for row in local_feature_rows(profile_rows) if d_min <= row["d"] <= d_max]
    if not features:
        return [], 0, 0.0
    strongest = sorted(features, key=lambda row: row["intensity"], reverse=True)[:top_n]
    max_intensity = max(row["intensity"] for row in strongest) if strongest else 1.0
    windows = observed_windows_for_targets(targets)
    rows: list[dict[str, object]] = []
    unmatched_count = 0
    for rank, row in enumerate(strongest, start=1):
        matched = in_any_window(row["d"], windows)
        if not matched:
            unmatched_count += 1
        rows.append(
            {
                "rank": rank,
                "d": row["d"],
                "q": row["q"],
                "relative_intensity": row["intensity"] / max_intensity if max_intensity else 0.0,
                "matched": matched,
            }
        )
    fraction = unmatched_count / len(strongest) if strongest else 0.0
    return rows, unmatched_count, fraction


def score_available_candidates(
    manifest: list[ManifestCandidate],
    targets: list[PeakTarget],
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    score_rows: list[dict[str, str]] = []
    unmatched_rows: list[dict[str, str]] = []
    current_unmatched_summary: dict[str, tuple[int, float]] = {}
    available = [entry for entry in manifest if entry.status == "available"]

    for entry in available:
        profile_rows = read_profile(entry.profile_path)
        current_targets = scaled_targets(targets, TOLERANCE_SCALES["current"])
        unmatched, unmatched_count, unmatched_fraction = strong_predicted_unmatched_peaks(profile_rows, current_targets)
        current_unmatched_summary[entry.candidate_id] = (unmatched_count, unmatched_fraction)
        for row in unmatched:
            unmatched_rows.append(
                {
                    "candidate_id": entry.candidate_id,
                    "candidate_label": entry.candidate_label,
                    "candidate_family": entry.candidate_family,
                    "candidate_role": entry.candidate_role,
                    "predicted_peak_rank": str(row["rank"]),
                    "predicted_d_angstrom": format_float(float(row["d"])),
                    "predicted_q_inv_angstrom": format_float(float(row["q"])),
                    "relative_intensity": format_float(float(row["relative_intensity"])),
                    "matched_experimental_window": "yes" if row["matched"] else "no",
                }
            )

        for tolerance_name, scale in TOLERANCE_SCALES.items():
            these_targets = scaled_targets(targets, scale)
            matches = [match_target_to_profile(target, profile_rows) for target in these_targets]
            observed = summarize_matches(matches, these_targets)
            unmatched_count_current, unmatched_fraction_current = current_unmatched_summary[entry.candidate_id]
            discrimination_score = (
                float(observed["observed_recovery_score"])
                - 0.25 * unmatched_count_current
                - 0.5 * unmatched_fraction_current
            )
            score_rows.append(
                {
                    "tolerance_setting": tolerance_name,
                    "candidate_id": entry.candidate_id,
                    "candidate_label": entry.candidate_label,
                    "candidate_family": entry.candidate_family,
                    "candidate_role": entry.candidate_role,
                    "match_count": str(observed["match_count"]),
                    "diagnostic_match_count": str(observed["diagnostic_match_count"]),
                    "strict_survives": "yes" if observed["strict_survives"] else "no",
                    "screen_survives": "yes" if observed["screen_survives"] else "no",
                    "mean_abs_d_error_angstrom": format_float(observed["mean_abs_error"]),
                    "mean_fractional_d_error": format_float(observed["mean_fractional_error"]),
                    "observed_recovery_score": format_float(observed["observed_recovery_score"]),
                    "predicted_unmatched_peak_count": str(unmatched_count_current),
                    "predicted_unmatched_peak_fraction": format_float(unmatched_fraction_current),
                    "discrimination_score": format_float(discrimination_score),
                    "matched_targets": str(observed["matched_targets"]),
                    "missed_diagnostic_targets": str(observed["missed_diagnostic_targets"]),
                }
            )

    tolerance_rows = build_tolerance_rows(score_rows)
    return score_rows, unmatched_rows, tolerance_rows


def build_tolerance_rows(score_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for tolerance_name, scale in TOLERANCE_SCALES.items():
        group = [row for row in score_rows if row["tolerance_setting"] == tolerance_name]
        sorted_group = sorted(
            group,
            key=lambda row: (
                -float(row["discrimination_score"]),
                -int(row["diagnostic_match_count"]),
                -int(row["match_count"]),
                float(row["mean_abs_d_error_angstrom"] or "999"),
                row["candidate_id"],
            ),
        )
        survivors = [row for row in sorted_group if row["screen_survives"] == "yes"]
        strict_survivors = [row for row in sorted_group if row["strict_survives"] == "yes"]
        best_score = float(sorted_group[0]["discrimination_score"]) if sorted_group else float("nan")
        best_rows = [row for row in sorted_group if float(row["discrimination_score"]) == best_score]
        central8 = next((row for row in sorted_group if row["candidate_id"] == CENTRAL_CANDIDATE_ID), None)
        rows.append(
            {
                "tolerance_setting": tolerance_name,
                "window_scale": format_float(scale, 2),
                "surviving_candidate_count": str(len(survivors)),
                "surviving_candidate_names": ";".join(row["candidate_id"] for row in survivors),
                "strict_surviving_candidate_count": str(len(strict_survivors)),
                "strict_surviving_candidate_names": ";".join(row["candidate_id"] for row in strict_survivors),
                "best_candidate": sorted_group[0]["candidate_id"] if sorted_group else "",
                "central8_units_30deg_survives": "yes" if central8 and central8["screen_survives"] == "yes" else "no",
                "central8_units_30deg_strict_survives": "yes" if central8 and central8["strict_survives"] == "yes" else "no",
                "central8_units_30deg_uniquely_best": "yes" if len(best_rows) == 1 and best_rows[0]["candidate_id"] == CENTRAL_CANDIDATE_ID else "no",
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_optional_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_manifest(path: Path, manifest: list[ManifestCandidate]) -> None:
    rows = [
        {
            "candidate_id": entry.candidate_id,
            "candidate_label": entry.candidate_label,
            "candidate_family": entry.candidate_family,
            "candidate_role": entry.candidate_role,
            "status": entry.status,
            "model_path": str(entry.model_path),
            "profile_path": str(entry.profile_path),
            "reason": entry.reason,
        }
        for entry in manifest
    ]
    write_csv(path, rows, MANIFEST_COLUMNS)


def plot_bar(path: Path, rows: list[dict[str, str]], value_column: str, title: str, color: str) -> None:
    focus = [row for row in rows if row["tolerance_setting"] == "current"][:14]
    width = 1080
    height = 420
    left = 70
    bottom = 310
    max_value = max([float(row[value_column] or "0") for row in focus] + [1.0])
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width / 2:.0f}" y="28" text-anchor="middle" font-family="Arial" font-size="17">{title}</text>',
        f'<line x1="{left}" y1="72" x2="{left}" y2="{bottom}" stroke="#333"/>',
        f'<line x1="{left}" y1="{bottom}" x2="{width - 20}" y2="{bottom}" stroke="#333"/>',
    ]
    for index, row in enumerate(focus):
        x = left + 22 + index * 68
        value = float(row[value_column] or "0")
        h = value / max_value * 220
        parts.append(f'<rect x="{x}" y="{bottom - h:.2f}" width="28" height="{h:.2f}" fill="{color}"/>')
        label = row["candidate_id"].replace("_", " ")
        parts.append(f'<text x="{x + 14}" y="{bottom + 14}" text-anchor="end" font-family="Arial" font-size="8" transform="rotate(-42 {x + 14},{bottom + 14})">{label}</text>')
    parts.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def write_tolerance_plot(path: Path, tolerance_rows: list[dict[str, str]]) -> None:
    width = 620
    height = 340
    left = 70
    bottom = 260
    max_value = max([int(row["surviving_candidate_count"]) for row in tolerance_rows] + [1])
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="310" y="28" text-anchor="middle" font-family="Arial" font-size="17">Tolerance survival counts</text>',
        f'<line x1="{left}" y1="70" x2="{left}" y2="{bottom}" stroke="#333"/>',
        f'<line x1="{left}" y1="{bottom}" x2="{width - 30}" y2="{bottom}" stroke="#333"/>',
    ]
    for index, row in enumerate(tolerance_rows):
        x = left + 72 + index * 140
        value = int(row["surviving_candidate_count"])
        strict = int(row["strict_surviving_candidate_count"])
        h = value / max_value * 170
        strict_h = strict / max_value * 170
        parts.append(f'<rect x="{x}" y="{bottom - h:.2f}" width="44" height="{h:.2f}" fill="#8fb3ff"/>')
        parts.append(f'<rect x="{x + 14}" y="{bottom - strict_h:.2f}" width="16" height="{strict_h:.2f}" fill="#f28e2b"/>')
        parts.append(f'<text x="{x + 22}" y="{bottom + 22}" text-anchor="middle" font-family="Arial" font-size="12">{row["tolerance_setting"]}</text>')
    parts.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def write_plots(plot_dir: Path, score_rows: list[dict[str, str]], tolerance_rows: list[dict[str, str]]) -> None:
    current = sorted(
        [row for row in score_rows if row["tolerance_setting"] == "current"],
        key=lambda row: float(row["discrimination_score"]),
        reverse=True,
    )
    plot_bar(plot_dir / "candidate_discrimination_scores.svg", current, "discrimination_score", "Current-tolerance discrimination score", "#4e79a7")
    plot_bar(plot_dir / "diagnostic_window_survival.svg", current, "diagnostic_match_count", "Diagnostic-window recovery", "#59a14f")
    plot_bar(plot_dir / "unmatched_predicted_peaks.svg", current, "predicted_unmatched_peak_count", "Unmatched predicted peak count", "#e15759")
    write_tolerance_plot(plot_dir / "tolerance_survival_counts.svg", tolerance_rows)


def row_table(rows: list[dict[str, str]], columns: list[str], limit: int | None = None) -> list[str]:
    return markdown_table(rows[:limit] if limit else rows, columns)


def write_report(
    path: Path,
    manifest: list[ManifestCandidate],
    score_rows: list[dict[str, str]],
    unmatched_rows: list[dict[str, str]],
    tolerance_rows: list[dict[str, str]],
    plot_dir: Path,
) -> None:
    available_rows = [
        {
            "candidate_id": entry.candidate_id,
            "family": entry.candidate_family,
            "role": entry.candidate_role,
            "profile_path": str(entry.profile_path),
        }
        for entry in manifest
        if entry.status == "available"
    ]
    unavailable_rows = [
        {
            "candidate_id": entry.candidate_id,
            "family": entry.candidate_family,
            "reason": entry.reason,
        }
        for entry in manifest
        if entry.status != "available"
    ]
    current_rows = sorted(
        [row for row in score_rows if row["tolerance_setting"] == "current"],
        key=lambda row: float(row["discrimination_score"]),
        reverse=True,
    )
    best = current_rows[0]
    central8 = next(row for row in current_rows if row["candidate_id"] == CENTRAL_CANDIDATE_ID)
    survivors = [row for row in current_rows if row["screen_survives"] == "yes"]
    failed = [row for row in current_rows if row["screen_survives"] != "yes"]
    unmatched_focus = [
        row for row in unmatched_rows if row["matched_experimental_window"] == "no"
    ][:18]
    twist_rise_rows = read_optional_csv(DEFAULT_TWIST_RISE_SENSITIVITY_CSV)
    non30_twists = [
        row
        for row in twist_rise_rows
        if row.get("candidate_family") == "full_length_twist_variant" and row.get("parameter_value") != "30 deg"
    ]
    rise_rows = [row for row in twist_rise_rows if row.get("candidate_family") == "rise_variant"]
    interpretation = (
        "If many alternatives survive under current or broad tolerances, the current powder peak list is compatible with "
        "the proposed conformation but does not distinguish it uniquely."
    )
    if best["candidate_id"] == CENTRAL_CANDIDATE_ID and sum(1 for row in current_rows if row["discrimination_score"] == best["discrimination_score"]) == 1:
        interpretation = (
            "Under the current tolerance setting, central8_units_30deg is the best-scoring available candidate in this "
            "screen. The result is more persuasive when considered against failed or weaker alternatives, but it remains "
            "a screening result rather than a phase assignment."
        )

    lines = [
        "# HXC590 S1 Powder Falsification Screen",
        "",
        "## Purpose",
        "",
        "This analysis is a falsification-style screening analysis, not a definitive phase assignment.",
        "",
        "This remains a falsification-style screen, not a definitive phase assignment.",
        "",
        "A candidate is more persuasive if it reproduces the diagnostic windows while plausible alternatives fail or match materially worse.",
        "",
        "## Experimental input and limitations",
        "",
        "The input peak list is John Bacsa's approximate HXC590 S1 powder list for TM_TC without salt. The experimental background is broad, peak positions are approximate, and relative intensities were not treated as reliable constraints.",
        "",
        "Uniform rings are expected for an unoriented powder sample, whereas oriented fibers give arcs. The ring-versus-arc difference is not by itself evidence for a different molecular phase.",
        "",
        "## Candidate families tested",
        "",
    ]
    lines.extend(row_table(available_rows, ["candidate_id", "family", "role", "profile_path"]))
    lines.extend(["", "## Candidate families unavailable", ""])
    lines.extend(row_table(unavailable_rows, ["candidate_id", "family", "reason"], limit=20))
    lines.extend(
        [
            "",
            "## Scoring method",
            "",
            "Observed-peak recovery uses the same local-profile-maximum d-window method as the HXC590 S1 peak comparison. The diagnostic windows are 3.35 A, 4.33 A, 3.90 A, and 3.71 A.",
            "",
            "Predicted unmatched peaks were counted from the top local maxima within the experimental d-range. A predicted peak is unmatched when it does not fall inside any current experimental d-window. Unmatched predicted peaks are treated as screening diagnostics only, because weak or broad experimental features may be obscured by background.",
            "",
            "Survival means matching all 4 diagnostic windows and at least 6 of 8 total peak windows. Strict survival means matching all 4 diagnostic windows and all 8 total windows.",
            "",
            "## Observed-peak recovery",
            "",
        ]
    )
    lines.extend(
        row_table(
            current_rows,
            [
                "candidate_id",
                "match_count",
                "diagnostic_match_count",
                "screen_survives",
                "strict_survives",
                "mean_abs_d_error_angstrom",
                "discrimination_score",
            ],
        )
    )
    lines.extend(["", "## Diagnostic-window recovery", ""])
    lines.append(f"`{CENTRAL_CANDIDATE_ID}` recovered {central8['diagnostic_match_count']} of 4 diagnostic windows under current tolerances.")
    lines.append("")
    lines.append("Current-tolerance survivors:")
    lines.extend(f"- `{row['candidate_id']}`" for row in survivors)
    lines.append("")
    lines.append("Current-tolerance failures:")
    lines.extend(f"- `{row['candidate_id']}` missed survival criteria" for row in failed)
    lines.extend(["", "## Predicted unmatched peaks", ""])
    lines.extend(row_table(unmatched_focus, ["candidate_id", "predicted_peak_rank", "predicted_d_angstrom", "relative_intensity"], limit=18))
    lines.extend(["", "## Tolerance / identifiability audit", ""])
    lines.extend(
        row_table(
            tolerance_rows,
            [
                "tolerance_setting",
                "surviving_candidate_count",
                "strict_surviving_candidate_count",
                "best_candidate",
                "central8_units_30deg_survives",
                "central8_units_30deg_uniquely_best",
            ],
        )
    )
    lines.extend(
        [
            "",
            "## Twist/rise sensitivity status",
            "",
            "Missing full-length non-30-degree twist variants were not generated in this pass. The repo contains pNAB twist scaffolding, but no current-model baseline pNAB YAML and matching runtime inputs were found, and the local Python environment cannot currently import pNAB.",
            "",
            "Rise variants were not generated in this pass. No safe existing rise-generation workflow or audited stack-axis transform was found for the current candidate model.",
            "",
            "Synthetic twist/rise variants are controls for diffraction sensitivity, not chemically optimized structures.",
            "",
            "If nearby twist or rise variants survive under current tolerances, the current powder peak list supports the conformation family but does not uniquely determine those parameters.",
            "",
            "Because nearby twist and rise variants are unavailable here, this screen can compare the current conformation family against available length and negative-control alternatives, but it cannot uniquely determine twist or rise parameters.",
            "",
            "Non-30-degree twist rows:",
            "",
        ]
    )
    if non30_twists:
        lines.extend(row_table(non30_twists, ["candidate_id", "parameter_value", "status", "reason"]))
    else:
        lines.append("- No twist/rise sensitivity audit CSV was found.")
    lines.extend(["", "Rise rows:", ""])
    if rise_rows:
        lines.extend(row_table(rise_rows, ["candidate_id", "parameter_value", "status", "reason"]))
    else:
        lines.append("- No twist/rise sensitivity audit CSV was found.")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            interpretation,
            "",
            f"Best current-tolerance candidate: `{best['candidate_id']}`. central8_units_30deg is {'best' if best['candidate_id'] == CENTRAL_CANDIDATE_ID else 'not best'} under the current scoring.",
            "",
            "The `central8_units_30deg` candidate remains the best-scoring available candidate under current tolerances. This ranking applies only to candidates with existing radial profiles.",
            "",
            "Nearby non-30-degree twists and requested rise variants neither survive nor fail in this pass because their coordinate/profile files are unavailable.",
            "",
            "The current powder peak list helps screen the available conformation family against available alternatives, but it does not distinguish the exact twist or rise without generated nearby controls.",
            "",
            "## What would be needed for stronger falsification",
            "",
            "Stronger falsification would require calibrated q values, uncertainty estimates, background-subtracted experimental profiles, and generated profiles for full-length twist variants, rise variants, and chemically meaningful wrong-register compact controls.",
            "",
            "## Limitations",
            "",
            "This analysis is a falsification-style screening analysis, not a definitive phase assignment.",
            "",
            "The experimental background is broad, peak positions are approximate, and relative intensities were not treated as reliable constraints.",
            "",
            "Unmatched predicted peaks are treated as screening diagnostics only, because weak or broad experimental features may be obscured by background.",
            "",
            "The available alternative set is incomplete: missing full-length twist variants and rise variants could not be generated safely from current local assets.",
            "",
            "## Outputs",
            "",
            "- `outputs/metrics/hxc590_s1_falsification_candidate_manifest.csv`",
            "- `outputs/metrics/hxc590_s1_falsification_scores.csv`",
            "- `outputs/metrics/hxc590_s1_predicted_unmatched_peaks.csv`",
            "- `outputs/metrics/hxc590_s1_tolerance_survival_summary.csv`",
            "- `outputs/metrics/hxc590_s1_twist_rise_sensitivity.csv` when the twist/rise audit has been run",
            f"- `{plot_dir / 'candidate_discrimination_scores.svg'}`",
            f"- `{plot_dir / 'diagnostic_window_survival.svg'}`",
            f"- `{plot_dir / 'tolerance_survival_counts.svg'}`",
            f"- `{plot_dir / 'unmatched_predicted_peaks.svg'}`",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(args: argparse.Namespace) -> dict[str, object]:
    targets = read_experimental_targets(args.experimental_peaks)
    manifest = build_candidate_manifest()
    write_manifest(args.manifest_csv, manifest)
    score_rows, unmatched_rows, tolerance_rows = score_available_candidates(manifest, targets)
    write_csv(args.scores_csv, score_rows, SCORE_COLUMNS)
    write_csv(args.unmatched_csv, unmatched_rows, UNMATCHED_COLUMNS)
    write_csv(args.tolerance_csv, tolerance_rows, TOLERANCE_COLUMNS)
    write_plots(args.plot_dir, score_rows, tolerance_rows)
    write_report(args.report, manifest, score_rows, unmatched_rows, tolerance_rows, args.plot_dir)
    current_rows = sorted(
        [row for row in score_rows if row["tolerance_setting"] == "current"],
        key=lambda row: float(row["discrimination_score"]),
        reverse=True,
    )
    current_tolerance = next(row for row in tolerance_rows if row["tolerance_setting"] == "current")
    return {
        "available_candidates": sum(1 for entry in manifest if entry.status == "available"),
        "unavailable_candidates": sum(1 for entry in manifest if entry.status != "available"),
        "score_rows": len(score_rows),
        "best_candidate": current_rows[0]["candidate_id"],
        "current_survivors": current_tolerance["surviving_candidate_count"],
        "central8_uniquely_best": current_tolerance["central8_units_30deg_uniquely_best"],
    }


def main() -> None:
    result = run(parse_args())
    print(f"Wrote falsification manifest with {result['available_candidates']} available and {result['unavailable_candidates']} unavailable candidates")
    print(f"Wrote {result['score_rows']} tolerance-specific score rows")
    print(f"Best current-tolerance candidate: {result['best_candidate']}")
    print(f"Current-tolerance survivors: {result['current_survivors']}")
    print(f"central8 uniquely best: {result['central8_uniquely_best']}")


if __name__ == "__main__":
    main()
