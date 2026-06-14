import csv
import importlib.util
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_script_module(name: str, relative_path: str):
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


combined = _load_script_module("combined_evidence_table", "scripts/combined_evidence_table.py")


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def test_model_name_normalization():
    assert combined.normalize_model_name("outputs/foo/model_a_heavy_deduped.pdb") == "model_a_heavy_deduped"
    assert combined.normalize_model_name("model_a_block_contact_summary") == "model_a"
    assert combined.normalize_model_name("model_a_debye_profile") == "model_a"


def test_block_count_extraction_from_included_blocks():
    assert combined.block_count_from_included_blocks("1,2,3") == "3"
    assert combined.block_count_from_included_blocks("") == ""
    assert combined.block_count_from_included_blocks(" 1, 3 ") == "2"


def test_safe_numeric_parsing_with_blanks():
    assert combined.safe_float("") is None
    assert combined.safe_float(None) is None
    assert combined.safe_float("1.25") == pytest.approx(1.25)
    assert combined.safe_int("") == 0


def test_pearson_helper():
    assert combined.pearson([1.0, 2.0, 3.0], [2.0, 4.0, 6.0]) == pytest.approx(1.0)
    assert combined.pearson([1.0, 2.0], [1.0, 2.0]) is None


def test_joining_small_fixture_csvs_and_derived_fields(tmp_path):
    block_dir = tmp_path / "block_contacts"
    _write_csv(
        block_dir / "scaffold_blocks_1_2_heavy_deduped_block_contact_summary.csv",
        ["contact_category", "contact_count", "GLU_involved_count", "GLU_GLU_count"],
        [
            {
                "contact_category": "scaffold_within_block",
                "contact_count": "10",
                "GLU_involved_count": "6",
                "GLU_GLU_count": "2",
            },
            {
                "contact_category": "scaffold_between_blocks",
                "contact_count": "5",
                "GLU_involved_count": "4",
                "GLU_GLU_count": "3",
            },
        ],
    )
    _write_csv(
        block_dir / "scaffold_blocks_1_2_heavy_deduped_block_pair_summary.csv",
        ["block_pair", "contact_count", "GLU_involved_count", "GLU_GLU_count", "min_distance_A"],
        [
            {
                "block_pair": "1--2",
                "contact_count": "5",
                "GLU_involved_count": "4",
                "GLU_GLU_count": "3",
                "min_distance_A": "3.2",
            }
        ],
    )
    ladder_rows = [
        {
            "model_name": "scaffold_blocks_1_2_heavy_deduped",
            "atom_mode": "heavy_deduped",
            "included_blocks": "1,2",
            "includes_hexads": "no",
            "atom_count": "100",
            "residue_count": "20",
            "contact_count_4p5A": "20",
            "motif_GLU_GLU": "4",
            "motif_GLU_any": "8",
            "angular_coverage_rad": "6.0",
        }
    ]
    fitted_rows = [
        {
            "model": "scaffold_blocks_1_2_heavy_deduped",
            "angular_coverage_rad": "6.1",
            "approximate_turns": "1.2",
            "approximate_pitch_per_turn": "10.0",
            "axial_span": "12.0",
            "z_axis_angle_degrees": "5.0",
        }
    ]
    diffraction_rows = [
        {
            "model_name": "scaffold_blocks_1_2_heavy_deduped",
            "d_4p5_fraction": "0.02",
            "d_3p4_fraction": "0.01",
        }
    ]

    rows = combined.build_combined_rows(ladder_rows, fitted_rows, block_dir, diffraction_rows)

    assert len(rows) == 1
    row = rows[0]
    assert row["block_count"] == "2"
    assert row["angular_coverage_rad_fitted"] == "6.1"
    assert row["scaffold_between_block_contacts"] == "5"
    assert row["top_block_pair_by_contacts"] == "1--2"
    assert row["scaffold_between_fraction"] == "0.333333"
    assert row["GLU_GLU_between_fraction"] == "0.600000"
    assert row["d4p5_per_GLU_GLU"] == "0.005000"
    assert row["d4p5_per_contact"] == "0.001000"


def test_markdown_report_includes_cautions_and_interpretation(tmp_path):
    out = tmp_path / "combined.md"
    rows = [
        {
            column: ""
            for column in combined.OUTPUT_COLUMNS
        }
    ]
    rows[0].update(
        {
            "model_name": "scaffold_blocks_1_heavy_deduped",
            "included_blocks": "1",
            "block_count": "1",
            "includes_hexads": "no",
            "angular_coverage_rad_fitted": "6.2",
            "approximate_turns": "1.0",
            "contact_count_4p5A": "57",
            "motif_GLU_GLU": "14",
            "scaffold_between_block_contacts": "0",
            "scaffold_hexad_or_other_contacts": "0",
            "d_4p5_fraction": "0.014",
            "d_3p4_fraction": "0.015",
        }
    )

    combined.write_markdown_report(rows, out, diffraction_available=True)

    text = out.read_text(encoding="utf-8")
    assert "Combined Evidence Report: Hexaplex Formation Metrics" in text
    assert "reciprocal-space feature, not a literal atom-distance" in text
    assert "Debye scores are simplified comparative approximations" in text
    assert "## Mechanistic interpretation" in text
    assert "individual scaffold paths are already folded/twisted" in text
