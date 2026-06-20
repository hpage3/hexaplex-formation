import csv
import importlib.util
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
INVENTORY_CSV = REPO_ROOT / "outputs/metrics/hxc590_s1_ideal_antiparallel_30deg_hexaplex_inventory.csv"
PROFILE_CSV = REPO_ROOT / "outputs/metrics/hxc590_s1_ideal_antiparallel_30deg_hexaplex_profile.csv"
COMPARISON_SCORES = (
    REPO_ROOT / "outputs/metrics/hxc590_s1_powder_corrected_with_ideal_hexaplex_peak_match_scores.csv"
)
FALSIFICATION_SCORES = (
    REPO_ROOT / "outputs/metrics/hxc590_s1_corrected_with_ideal_hexaplex_falsification_scores.csv"
)
TOLERANCE_SUMMARY = (
    REPO_ROOT / "outputs/metrics/hxc590_s1_corrected_with_ideal_hexaplex_tolerance_survival_summary.csv"
)
WITH_IDEAL_HEXAPLEX_REPORTS = [
    REPO_ROOT / "outputs/reports/hxc590_s1_ideal_antiparallel_30deg_hexaplex_inventory_report.md",
    REPO_ROOT / "outputs/reports/hxc590_s1_powder_corrected_with_ideal_hexaplex_peak_comparison_report.md",
    REPO_ROOT / "outputs/reports/hxc590_s1_corrected_with_ideal_hexaplex_falsification_report.md",
]
FORBIDDEN_OVERCLAIMING_PHRASES = [
    "proves",
    "proof",
    "confirmed structure",
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
    "inventory_hxc590_s1_ideal_antiparallel_30deg_hexaplex",
    "scripts/inventory_hxc590_s1_ideal_antiparallel_30deg_hexaplex.py",
)
compare = _load_script_module("compare_hxc590_s1_powder_peaks", "scripts/compare_hxc590_s1_powder_peaks.py")
falsify = _load_script_module("falsify_hxc590_s1_powder_candidates", "scripts/falsify_hxc590_s1_powder_candidates.py")


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def test_inventory_records_ideal_hexaplex_model_and_duplicate_status():
    row = read_rows(INVENTORY_CSV)[0]

    assert row["candidate_label"] == "ideal_antiparallel_30deg_hexaplex"
    assert row["sha256"] == "9a880544a551b3d16f9024e7897a23af0286820ed69df9e3cb805c94113b6aca"
    assert row["total_pdb_atom_count"] == "7146"
    assert row["deduped_atom_count"] == "3573"
    assert row["element_counts"] == "C:2160;H:2814;N:1170;O:1002"
    assert row["residue_counts"] == "CYP:2082;GLU:2712;MEP:2352"
    assert row["hydrogens_present"] == "yes"
    assert row["duplicate_atom_coordinate_line_count"] == "3573"
    assert row["duplicate_atom_coordinate_lines_detected"] == "yes"
    assert "ideal antiparallel 30-degree hexaplex model" in row["provenance_note"]
    assert "16-mer simulation benchmark shorthand" in row["provenance_note"]


def test_pdb_parser_preserves_raw_copied_coordinates_and_detects_duplicates():
    atoms = inventory.load_pdb_atoms(REPO_ROOT / "inputs/candidates/ideal_antiparallel_30deg_hexaplex.pdb")
    unique_identities = {inventory.atom_identity_key(atom) for atom in atoms}

    assert len(atoms) == 7146
    assert len(unique_identities) == 3573
    assert atoms[0].record_type == "ATOM"
    assert atoms[0].residue_name == "CYP"
    assert atoms[0].element == "N"


def test_reused_profile_and_corrected_scores_include_ideal_hexaplex_candidate():
    profile_rows = read_rows(PROFILE_CSV)
    comparison_rows = read_rows(COMPARISON_SCORES)
    ideal_rows = [row for row in comparison_rows if row["candidate_id"] == "ideal_antiparallel_30deg_hexaplex"]

    assert {"q_Ainv", "d_A", "mean_intensity"}.issubset(profile_rows[0])
    assert len(profile_rows) > 100
    assert len(ideal_rows) == 5
    assert {row["matched"] for row in ideal_rows} == {"yes", "no"}
    stacking = next(row for row in ideal_rows if row["target_id"] == "d_3p39_stacking")
    assert stacking["matched"] == "yes"
    assert float(stacking["matched_d_angstrom"]) == 3.328208
    assert stacking["candidate_match_count"] == "3"
    assert stacking["candidate_diagnostic_match_count"] == "2"


def test_optional_candidate_lists_include_nick_8hexad_and_ideal_hexaplex_candidates():
    candidates = compare.default_candidate_models(include_nick_8hexad=True, include_ideal_hexaplex=True)
    manifest = falsify.build_candidate_manifest(include_nick_8hexad=True, include_ideal_hexaplex=True)

    assert any(candidate.candidate_id == "nick_hexaplex_8hexads" for candidate in candidates)
    assert any(candidate.candidate_id == "ideal_antiparallel_30deg_hexaplex" for candidate in candidates)
    ideal_entry = next(entry for entry in manifest if entry.candidate_id == "ideal_antiparallel_30deg_hexaplex")
    assert ideal_entry.candidate_family == "ideal_antiparallel_30deg_hexaplex"
    assert "Ideal antiparallel 30-degree hexaplex model" in ideal_entry.candidate_label


def test_corrected_with_ideal_hexaplex_falsification_rows_do_not_survive_current_screen():
    rows = read_rows(FALSIFICATION_SCORES)
    current = next(
        row
        for row in rows
        if row["candidate_id"] == "ideal_antiparallel_30deg_hexaplex"
        and row["tolerance_setting"] == "current"
    )

    assert current["match_count"] == "3"
    assert current["diagnostic_match_count"] == "2"
    assert current["matched_targets"] == "5.58;4.40;3.38"
    assert current["missed_diagnostic_targets"] == "3.79"
    assert current["screen_survives"] == "no"
    assert current["strict_survives"] == "no"


def test_survivor_counts_remain_unchanged_with_ideal_hexaplex_candidate_added():
    rows = read_rows(TOLERANCE_SUMMARY)
    by_setting = {row["tolerance_setting"]: row for row in rows}

    assert by_setting["narrow"]["surviving_candidate_count"] == "0"
    assert by_setting["current"]["surviving_candidate_count"] == "1"
    assert by_setting["current"]["surviving_candidate_names"] == "central12_units_30deg"
    assert by_setting["broad"]["surviving_candidate_count"] == "2"
    assert by_setting["broad"]["surviving_candidate_names"] == "central12_units_30deg;central8_units_30deg"


def test_with_ideal_hexaplex_reports_avoid_overclaiming_and_constrain_16mer_language():
    for report in WITH_IDEAL_HEXAPLEX_REPORTS:
        text = report.read_text(encoding="utf-8")
        lowered = text.lower()
        assert not any(term in lowered for term in FORBIDDEN_OVERCLAIMING_PHRASES)
        assert "16-mer" not in lowered or "16-mer simulation benchmark" in lowered
