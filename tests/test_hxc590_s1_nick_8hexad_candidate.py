import csv
import importlib.util
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
INVENTORY_CSV = REPO_ROOT / "outputs/metrics/hxc590_s1_nick_8hexad_candidate_inventory.csv"
PROFILE_CSV = REPO_ROOT / "outputs/metrics/hxc590_s1_nick_hexaplex_8hexads_profile.csv"
COMPARISON_SCORES = REPO_ROOT / "outputs/metrics/hxc590_s1_powder_corrected_with_nick_8hexad_peak_match_scores.csv"
FALSIFICATION_SCORES = REPO_ROOT / "outputs/metrics/hxc590_s1_corrected_with_nick_8hexad_falsification_scores.csv"
WITH_NICK_REPORTS = [
    REPO_ROOT / "outputs/reports/hxc590_s1_nick_8hexad_candidate_inventory_report.md",
    REPO_ROOT / "outputs/reports/hxc590_s1_powder_corrected_with_nick_8hexad_peak_comparison_report.md",
    REPO_ROOT / "outputs/reports/hxc590_s1_corrected_with_nick_8hexad_falsification_report.md",
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


inventory = _load_script_module(
    "inventory_hxc590_s1_nick_8hexad_candidate",
    "scripts/inventory_hxc590_s1_nick_8hexad_candidate.py",
)
compare = _load_script_module("compare_hxc590_s1_powder_peaks", "scripts/compare_hxc590_s1_powder_peaks.py")
falsify = _load_script_module("falsify_hxc590_s1_powder_candidates", "scripts/falsify_hxc590_s1_powder_candidates.py")


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def test_inventory_output_records_checksum_atom_count_and_provenance():
    row = read_rows(INVENTORY_CSV)[0]

    assert row["candidate_label"] == "nick_hexaplex_8hexads"
    assert row["sha256"] == "1cacfc20ed6e6660f9832d9d7c5c68a01386d468651e693d71559c818a632702"
    assert row["xyz_atom_count_observed"] == "2064"
    assert row["element_counts"] == "C:1032;N:556;O:476"
    assert "Nick-provided Hexaplex_8Hexads.xyz" in row["provenance_note"]
    assert "not a 16-mer" in row["not_16mer_note"]


def test_xyz_parser_detects_duplicate_lines_without_modifying_coordinates():
    _header, _comment, atoms = inventory.read_xyz(REPO_ROOT / "inputs/candidates/nick_hexaplex_8hexads.xyz")
    duplicate_count = len(atoms) - len({atom.raw_line for atom in atoms})

    assert len(atoms) == 2064
    assert duplicate_count == 1032


def test_profile_and_corrected_scores_include_nick_8hexad_candidate():
    profile_rows = read_rows(PROFILE_CSV)
    comparison_rows = read_rows(COMPARISON_SCORES)
    nick_rows = [row for row in comparison_rows if row["candidate_id"] == "nick_hexaplex_8hexads"]

    assert profile_rows[0] == {"q_Ainv": "0.200000", "d_A": "31.415927", "intensity": profile_rows[0]["intensity"]}
    assert len(nick_rows) == 5
    stacking = next(row for row in nick_rows if row["target_d_angstrom"] == "3.38")
    assert stacking["matched"] == "yes"
    assert float(stacking["matched_d_angstrom"]) == 3.378057


def test_optional_candidate_manifests_include_nick_candidate():
    candidates = compare.default_candidate_models(include_nick_8hexad=True)
    manifest = falsify.build_candidate_manifest(include_nick_8hexad=True)

    assert any(candidate.candidate_id == "nick_hexaplex_8hexads" for candidate in candidates)
    nick_entry = next(entry for entry in manifest if entry.candidate_id == "nick_hexaplex_8hexads")
    assert nick_entry.candidate_family == "nick_provided_8hexad"
    assert "8-hexad candidate" in nick_entry.candidate_label


def test_corrected_with_nick_falsification_rows_show_no_current_survival():
    rows = read_rows(FALSIFICATION_SCORES)
    current = next(row for row in rows if row["candidate_id"] == "nick_hexaplex_8hexads" and row["tolerance_setting"] == "current")

    assert current["match_count"] == "3"
    assert current["diagnostic_match_count"] == "2"
    assert current["screen_survives"] == "no"
    assert current["strict_survives"] == "no"


def test_with_nick_reports_avoid_overclaiming_and_do_not_label_candidate_as_16mer():
    for report in WITH_NICK_REPORTS:
        text = report.read_text(encoding="utf-8").lower()
        assert not any(term in text for term in FORBIDDEN_REPORT_TERMS)
        if "16-mer" in text:
            assert "not a 16-mer" in text
