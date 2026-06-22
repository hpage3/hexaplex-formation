import importlib.util
import math
import sys
from pathlib import Path


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "score_peak_position_fit.py"
SPEC = importlib.util.spec_from_file_location("score_peak_position_fit", SCRIPT_PATH)
score_peak_position_fit = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = score_peak_position_fit
SPEC.loader.exec_module(score_peak_position_fit)


TargetPeak = score_peak_position_fit.TargetPeak
PeakAssignment = score_peak_position_fit.PeakAssignment


def test_relative_error_and_unweighted_metrics():
    errors = [
        score_peak_position_fit.relative_error(7.3, 7.1),
        score_peak_position_fit.relative_error(5.65, 5.75),
    ]

    expected_root = math.sqrt(sum(error * error for error in errors))
    expected_rms = math.sqrt(sum(error * error for error in errors) / 2)

    assert errors[0] == (7.3 - 7.1) / 7.3
    assert score_peak_position_fit.root_sum_relative_error(errors) == expected_root
    assert score_peak_position_fit.rms_relative_error(errors) == expected_rms


def test_weighted_rms_relative_error():
    errors = [0.10, 0.20, 0.30]
    weights = [1.0, 2.0, 3.0]
    expected = math.sqrt((1.0 * 0.01 + 2.0 * 0.04 + 3.0 * 0.09) / 6.0)

    assert score_peak_position_fit.weighted_rms_relative_error(errors, weights) == expected


def test_missing_peak_penalty_and_ranking():
    targets = [
        TargetPeak("D", 7.3, 1.0),
        TargetPeak("C", 5.65, 2.0),
    ]
    assignments = [
        PeakAssignment("good", "29", "D", 7.3),
        PeakAssignment("good", "29", "C", 5.65),
        PeakAssignment("missing", "31", "D", 7.3),
    ]

    summary_rows, per_peak_rows = score_peak_position_fit.score_models(
        targets,
        assignments,
        missing_peak_penalty=0.25,
    )

    assert [row["model_id"] for row in summary_rows] == ["good", "missing"]
    missing_summary = summary_rows[1]
    assert missing_summary["missing_peak_count"] == 1
    assert math.isclose(float(missing_summary["weighted_rms_relative_error"]), math.sqrt(
        (1.0 * 0.0 + 2.0 * 0.25 * 0.25) / 3.0
    ))
    missing_peak_rows = [
        row for row in per_peak_rows if row["model_id"] == "missing" and row["target_label"] == "C"
    ]
    assert missing_peak_rows[0]["is_missing"] == "true"
    assert float(missing_peak_rows[0]["relative_error"]) == 0.25


def test_nearest_peak_matching_respects_tolerance():
    targets = [
        TargetPeak("D", 7.3, 1.0),
        TargetPeak("C", 5.65, 1.0),
        TargetPeak("B", 4.4, 1.0),
    ]
    peak_rows = [
        {"model_id": "model-a", "twist_deg": "30", "peak_d_A": "7.28", "peak_intensity_optional": "10"},
        {"model_id": "model-a", "twist_deg": "30", "peak_d_A": "5.90", "peak_intensity_optional": "8"},
        {"model_id": "model-a", "twist_deg": "30", "peak_d_A": "4.90", "peak_intensity_optional": "7"},
    ]

    assignments = score_peak_position_fit.match_peak_list(targets, peak_rows, tolerance_A=0.35)

    assert [(item.target_label, item.theoretical_d_A) for item in assignments] == [
        ("D", 7.28),
        ("C", 5.90),
    ]
    assert assignments[0].assignment_method == "nearest_within_0.35A"
