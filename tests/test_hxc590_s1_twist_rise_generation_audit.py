import importlib.util
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_script_module(name: str, relative_path: str):
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


audit = _load_script_module("audit_hxc590_s1_twist_rise_generation", "scripts/audit_hxc590_s1_twist_rise_generation.py")


def test_rise_token_is_deterministic():
    assert audit.rise_token(3.2) == "3p20"
    assert audit.rise_token(3.35) == "3p35"


def test_audit_candidate_counts_and_labels():
    rows = audit.twist_rows() + audit.rise_rows()

    assert len(rows) == 13
    assert {row["candidate_id"] for row in rows} >= {
        "full_length_twist_24",
        "full_length_twist_30",
        "full_length_twist_36",
        "rise_3p20_synthetic_control",
        "rise_3p60_synthetic_control",
    }


def test_non30_twists_and_rises_are_marked_unavailable_without_profiles():
    rows = audit.twist_rows() + audit.rise_rows()
    by_id = {row["candidate_id"]: row for row in rows}

    assert by_id["full_length_twist_24"]["status"] == "unavailable"
    assert by_id["rise_3p40_synthetic_control"]["status"] == "unavailable"
    assert "baseline YAML" in by_id["full_length_twist_24"]["reason"]
    assert "rise-generation workflow" in by_id["rise_3p40_synthetic_control"]["reason"]


def test_report_includes_conservative_language(tmp_path):
    rows = audit.twist_rows() + audit.rise_rows()
    report = tmp_path / "gap.md"

    audit.write_report(report, rows)

    text = report.read_text(encoding="utf-8")
    assert "falsification-style screen, not a definitive phase assignment" in text
    assert "controls for diffraction sensitivity, not chemically optimized structures" in text
    assert "does not uniquely determine those parameters" in text
