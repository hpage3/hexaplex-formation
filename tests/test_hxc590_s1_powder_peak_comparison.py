import csv
import importlib.util
import math
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_script_module(name: str, relative_path: str):
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


powder = _load_script_module("compare_hxc590_s1_powder_peaks", "scripts/compare_hxc590_s1_powder_peaks.py")


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def test_d_to_q_conversion():
    assert powder.q_from_d(3.35) == pytest.approx(2.0 * math.pi / 3.35)


def test_d_window_matching_uses_configured_tolerance():
    target = powder.PeakTarget("sample", "d_3p35", 3.35, powder.q_from_d(3.35), 0.06, "medium_high", "")
    profile = [
        {"d": 3.20, "q": powder.q_from_d(3.20), "intensity": 1000.0},
        {"d": 3.36, "q": powder.q_from_d(3.36), "intensity": 10.0},
    ]

    match = powder.match_target_to_profile(target, profile)

    assert match["matched"] is True
    assert match["matched_row"]["d"] == pytest.approx(3.36)
    assert match["abs_error"] == pytest.approx(0.01)


def test_broad_background_confidence_marks_long_spacing_as_non_diagnostic(tmp_path):
    peaks = tmp_path / "peaks.csv"
    write_csv(
        peaks,
        [
            {
                "sample_id": powder.SAMPLE_ID,
                "d_angstrom": "7.90",
                "q_inv_angstrom": f"{powder.q_from_d(7.90):.6f}",
                "confidence": "lower",
                "note": "broad background; relative intensities approximate",
            },
            {
                "sample_id": powder.SAMPLE_ID,
                "d_angstrom": "3.35",
                "q_inv_angstrom": f"{powder.q_from_d(3.35):.6f}",
                "confidence": "medium_high",
                "note": "diagnostic spacing",
            },
        ],
    )

    targets = powder.read_experimental_targets(peaks)

    by_d = {round(target.d_angstrom, 2): target for target in targets}
    assert by_d[7.90].diagnostic_window is False
    assert by_d[3.35].diagnostic_window is True


def test_ranking_does_not_depend_on_relative_intensity():
    targets = [powder.PeakTarget("sample", "d_4p33", 4.33, powder.q_from_d(4.33), 0.10, "medium_high", "")]
    low_intensity_match = powder.match_target_to_profile(
        targets[0],
        [{"d": 4.32, "q": powder.q_from_d(4.32), "intensity": 1.0}],
    )
    high_intensity_match = powder.match_target_to_profile(
        targets[0],
        [{"d": 4.32, "q": powder.q_from_d(4.32), "intensity": 999999.0}],
    )

    low_summary = powder.summarize_candidate([low_intensity_match], targets)
    high_summary = powder.summarize_candidate([high_intensity_match], targets)

    assert low_summary["rank_score"] == pytest.approx(high_summary["rank_score"])


def test_report_includes_conservative_limitation_language(tmp_path):
    targets = [
        powder.PeakTarget("sample", "d_4p33", 4.33, powder.q_from_d(4.33), 0.10, "medium_high", ""),
        powder.PeakTarget("sample", "d_3p90", 3.90, powder.q_from_d(3.90), 0.08, "medium_high", ""),
        powder.PeakTarget("sample", "d_3p71", 3.71, powder.q_from_d(3.71), 0.08, "medium_high", ""),
        powder.PeakTarget("sample", "d_3p35", 3.35, powder.q_from_d(3.35), 0.06, "medium_high", ""),
    ]
    candidate = powder.CandidateModel("fixture", "Fixture candidate", Path("fixture.pdb"), Path("fixture.csv"), "fixture")
    score_rows = []
    for target in targets:
        score_rows.append(
            {
                "candidate_id": "fixture",
                "target_d_angstrom": f"{target.d_angstrom:.2f}",
                "diagnostic_window": "yes",
                "matched": "yes",
                "matched_d_angstrom": f"{target.d_angstrom:.6f}",
                "abs_d_error_angstrom": "0.000000",
                "window_point_count": "1",
            }
        )
    summary_rows = [
        {
            "candidate_id": "fixture",
            "candidate_label": "Fixture candidate",
            "model_path": "fixture.pdb",
            "profile_path": "fixture.csv",
            "match_count": "4",
            "diagnostic_match_count": "4",
            "mean_abs_d_error_angstrom": "0.000000",
            "mean_fractional_d_error": "0.000000",
            "rank_score": "6.000000",
        }
    ]

    report = tmp_path / "report.md"
    powder.write_report(report, targets, [candidate], score_rows, summary_rows, tmp_path / "plots")

    text = report.read_text(encoding="utf-8")
    assert "does not establish a definitive phase assignment" in text
    assert "relative intensities were not treated as reliable constraints" in text
    assert "ring-versus-arc difference is not by itself evidence" in text
