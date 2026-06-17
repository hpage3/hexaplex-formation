import importlib.util
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
FORBIDDEN_TERMS = [
    "proves",
    "proof",
    "confirmed same phase",
    "actual phase",
    "true structure",
    "nailed",
    "beyond doubt",
]


def _load_script_module(name: str, relative_path: str):
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


falsify = _load_script_module("falsify_hxc590_s1_powder_candidates", "scripts/falsify_hxc590_s1_powder_candidates.py")


def target(d_value: float, half_width: float, confidence: str = "medium_high") -> falsify.PeakTarget:
    return falsify.PeakTarget("sample", f"d_{d_value}", d_value, 2.0, half_width, confidence, "")


def test_candidate_manifest_marks_available_and_unavailable_candidates():
    manifest = falsify.build_candidate_manifest()

    statuses = {entry.status for entry in manifest}
    assert "available" in statuses
    assert "unavailable" in statuses
    assert any(entry.candidate_family == "rise_variant" and entry.status == "unavailable" for entry in manifest)
    assert any(entry.candidate_id == "central8_units_30deg" and entry.status == "available" for entry in manifest)


def test_observed_window_recovery_scoring():
    targets = [target(3.35, 0.06), target(4.33, 0.10, "lower")]
    matches = [
        {"matched": True, "abs_error": 0.01, "fractional_error": 0.01 / 3.35},
        {"matched": False, "abs_error": 0.20, "fractional_error": 0.20 / 4.33},
    ]

    summary = falsify.summarize_matches(matches, targets)

    assert summary["match_count"] == 1
    assert summary["diagnostic_match_count"] == 1
    assert summary["strict_survives"] is False


def test_diagnostic_survival_logic():
    targets = [target(3.35, 0.06), target(4.33, 0.10), target(3.90, 0.08), target(3.71, 0.08)] + [
        target(7.90, 0.25, "lower"),
        target(7.30, 0.25, "lower"),
    ]
    matches = [{"matched": True, "abs_error": 0.01, "fractional_error": 0.001} for _ in targets]

    summary = falsify.summarize_matches(matches, targets)

    assert summary["screen_survives"] is True
    assert summary["strict_survives"] is False


def test_predicted_unmatched_peak_detection():
    targets = [target(3.35, 0.06), target(4.33, 0.10)]
    profile = [
        {"d": 3.00, "q": 2.09, "intensity": 1.0},
        {"d": 3.35, "q": 1.88, "intensity": 10.0},
        {"d": 3.80, "q": 1.65, "intensity": 1.0},
        {"d": 4.00, "q": 1.57, "intensity": 8.0},
        {"d": 4.33, "q": 1.45, "intensity": 2.0},
        {"d": 4.60, "q": 1.37, "intensity": 1.0},
    ]

    rows, count, fraction = falsify.strong_predicted_unmatched_peaks(profile, targets, top_n=2)

    assert count == 1
    assert fraction == 0.5
    assert any(row["matched"] is False for row in rows)


def test_tolerance_scaling_changes_survival_counts():
    score_rows = [
        {
            "tolerance_setting": "narrow",
            "candidate_id": "a",
            "screen_survives": "no",
            "strict_survives": "no",
            "discrimination_score": "1",
            "diagnostic_match_count": "3",
            "match_count": "5",
            "mean_abs_d_error_angstrom": "0.1",
        },
        {
            "tolerance_setting": "current",
            "candidate_id": "a",
            "screen_survives": "yes",
            "strict_survives": "no",
            "discrimination_score": "2",
            "diagnostic_match_count": "4",
            "match_count": "6",
            "mean_abs_d_error_angstrom": "0.1",
        },
    ]

    rows = falsify.build_tolerance_rows(score_rows)

    by_setting = {row["tolerance_setting"]: row for row in rows}
    assert by_setting["narrow"]["surviving_candidate_count"] == "0"
    assert by_setting["current"]["surviving_candidate_count"] == "1"


def test_relative_experimental_intensities_are_not_required():
    targets = [target(3.35, 0.06)]
    profile = [
        {"d": 3.20, "q": 1.96, "intensity": 1.0},
        {"d": 3.35, "q": 1.88, "intensity": 5.0},
        {"d": 3.50, "q": 1.79, "intensity": 1.0},
    ]

    match = falsify.match_target_to_profile(targets[0], profile)
    summary = falsify.summarize_matches([match], targets)

    assert summary["match_count"] == 1
    assert "observed_recovery_score" in summary


def test_report_includes_conservative_language_and_avoids_forbidden_terms(tmp_path):
    manifest = [
        falsify.unavailable_candidate(
            "missing",
            "Missing candidate",
            "rise_variant",
            "plausible_alternative",
            Path("missing.pdb"),
            Path("missing.csv"),
            "missing by design",
        )
    ]
    score_rows = [
        {
            "tolerance_setting": "current",
            "candidate_id": "central8_units_30deg",
            "candidate_label": "central8",
            "candidate_family": "existing",
            "candidate_role": "candidate",
            "match_count": "8",
            "diagnostic_match_count": "4",
            "screen_survives": "yes",
            "strict_survives": "yes",
            "mean_abs_d_error_angstrom": "0.010000",
            "mean_fractional_d_error": "0.001000",
            "observed_recovery_score": "9.999000",
            "predicted_unmatched_peak_count": "0",
            "predicted_unmatched_peak_fraction": "0.000000",
            "discrimination_score": "9.999000",
            "matched_targets": "3.35",
            "missed_diagnostic_targets": "",
        }
    ]
    unmatched_rows = []
    tolerance_rows = [
        {
            "tolerance_setting": "current",
            "surviving_candidate_count": "1",
            "strict_surviving_candidate_count": "1",
            "best_candidate": "central8_units_30deg",
            "central8_units_30deg_survives": "yes",
            "central8_units_30deg_uniquely_best": "yes",
        }
    ]

    report = tmp_path / "report.md"
    falsify.write_report(report, manifest, score_rows, unmatched_rows, tolerance_rows, tmp_path / "plots")

    text = report.read_text(encoding="utf-8")
    assert "falsification-style screening analysis, not a definitive phase assignment" in text
    assert "Unmatched predicted peaks are treated as screening diagnostics only" in text
    lowered = text.lower()
    assert not any(term in lowered for term in FORBIDDEN_TERMS)
