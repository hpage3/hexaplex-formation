import csv
import importlib.util
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
TRACE_CSV = REPO_ROOT / "inputs/experimental/hxc590_s1_powder_corrected_trace.csv"
TARGETS_CSV = REPO_ROOT / "outputs/metrics/hxc590_s1_powder_corrected_peak_targets.csv"
CORRECTED_REPORTS = [
    REPO_ROOT / "outputs/reports/hxc590_s1_powder_corrected_peak_comparison_report.md",
    REPO_ROOT / "outputs/reports/hxc590_s1_corrected_falsification_report.md",
]
FORBIDDEN_REPORT_TERMS = [
    "proof",
    "proves",
    "confirmed",
    "actual phase",
    "true structure",
    "exact structure",
]


def _load_script_module(name: str, relative_path: str):
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


extract = _load_script_module(
    "extract_hxc590_s1_corrected_powder_peaks",
    "scripts/extract_hxc590_s1_corrected_powder_peaks.py",
)


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def test_corrected_trace_csv_exists_with_expected_columns_and_raw_values():
    rows = read_rows(TRACE_CSV)

    assert rows
    assert set(rows[0]) == {"distance_a", "intensity"}
    assert len(rows) > 300
    assert rows[0] == {"distance_a": "9.828", "intensity": "1.518"}
    assert rows[-1] == {"distance_a": "2.222", "intensity": "0.030"}


def test_extraction_preserves_raw_trace_and_finds_stacking_feature_near_3p4():
    points = extract.read_trace(TRACE_CSV)
    targets = extract.extract_targets(points)

    by_label = {row["peak_label"]: row for row in targets}
    stacking = by_label["corrected_3p39_stacking"]

    assert len(points) > 300
    assert 3.38 <= float(stacking["distance_a"]) <= 3.40
    assert stacking["diagnostic_role"] == "diagnostic"
    assert "raw trace rows were not deduplicated or rewritten" in stacking["notes"]


def test_corrected_target_outputs_are_created():
    assert TARGETS_CSV.exists()
    rows = read_rows(TARGETS_CSV)
    labels = {row["peak_label"] for row in rows}

    assert {
        "corrected_broad_7p3",
        "corrected_5p56",
        "corrected_4p40",
        "corrected_3p80",
        "corrected_3p39_stacking",
    } <= labels


def test_corrected_reports_avoid_overclaiming_terms():
    for report in CORRECTED_REPORTS:
        assert report.exists()
        text = report.read_text(encoding="utf-8").lower()
        assert not any(term in text for term in FORBIDDEN_REPORT_TERMS)
