import csv
import importlib.util
import sys
from pathlib import Path

import pytest

from hexaplex_formation.pdb_utils import PDBAtom, write_pdb_atoms


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_ROOT = REPO_ROOT / "scripts"
if str(SCRIPTS_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_ROOT))


def _load_script_module(name: str, relative_path: str):
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


convert_map = _load_script_module("convert_strand_map_to_scaffold_path_map", "scripts/convert_strand_map_to_scaffold_path_map.py")
validate_map = _load_script_module("validate_scaffold_path_map", "scripts/validate_scaffold_path_map.py")
pymol_helper = _load_script_module("generate_pymol_strand_map_helper", "scripts/generate_pymol_strand_map_helper.py")
compare_maps = _load_script_module("compare_scaffold_path_maps", "scripts/compare_scaffold_path_maps.py")


def atom(serial: int, residue_name: str, residue_number: int, x: float, y: float, z: float) -> PDBAtom:
    return PDBAtom(
        record_type="ATOM",
        atom_serial=serial,
        atom_name="CA",
        alt_loc="",
        residue_name=residue_name,
        chain_id="",
        residue_number=residue_number,
        insertion_code="",
        x=x,
        y=y,
        z=z,
        occupancy=1.0,
        temp_factor=0.0,
        element="C",
    )


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _fixture_pdb(path: Path) -> None:
    write_pdb_atoms(
        [
            atom(1, "GLU", 1, 0.0, 0.0, 0.0),
            atom(2, "CYP", 2, 1.0, 0.0, 1.0),
            atom(3, "GLU", 3, 0.0, 1.0, 2.0),
        ],
        path,
    )


def _canonical_rows() -> list[dict[str, str]]:
    return [
        {
            "map_name": "fixture",
            "strand_id": "1",
            "strand_label": "path_1",
            "residue_index_in_pdb_order": "1",
            "chain_id": "",
            "residue_name": "GLU",
            "residue_number": "1",
            "insertion_code": "",
            "residue_label": "GLU1",
            "source": "test",
        },
        {
            "map_name": "fixture",
            "strand_id": "1",
            "strand_label": "path_1",
            "residue_index_in_pdb_order": "2",
            "chain_id": "",
            "residue_name": "CYP",
            "residue_number": "2",
            "insertion_code": "",
            "residue_label": "CYP2",
            "source": "test",
        },
        {
            "map_name": "fixture",
            "strand_id": "2",
            "strand_label": "path_2",
            "residue_index_in_pdb_order": "3",
            "chain_id": "",
            "residue_name": "GLU",
            "residue_number": "3",
            "insertion_code": "",
            "residue_label": "GLU3",
            "source": "test",
        },
    ]


def test_conversion_from_old_block_map_to_canonical_map():
    rows = [
        {
            "block_id": "2",
            "residue_index_in_pdb_order": "31",
            "chain_id": "",
            "residue_name": "GLU",
            "residue_number": "31",
            "insertion_code": "",
            "residue_label": "GLU31",
        }
    ]

    converted = convert_map.convert_rows(rows, "candidate", "legacy")

    assert converted[0]["map_name"] == "candidate"
    assert converted[0]["strand_id"] == "2"
    assert converted[0]["strand_label"] == "block_2"
    assert converted[0]["source"] == "legacy"


def test_map_validation_passes_complete_fixture(tmp_path):
    pdb_path = tmp_path / "fixture.pdb"
    _fixture_pdb(pdb_path)

    result = validate_map.validate_map(pdb_path, _canonical_rows())

    assert result["validation_status"] == "pass"
    assert result["pdb_residue_count"] == 3
    assert result["mapped_residue_count"] == 3
    assert len(result["strand_summaries"]) == 2


def test_map_validation_catches_missing_residue(tmp_path):
    pdb_path = tmp_path / "fixture.pdb"
    _fixture_pdb(pdb_path)

    result = validate_map.validate_map(pdb_path, _canonical_rows()[:2])

    assert result["validation_status"] == "fail"
    assert result["missing_from_map_count"] == 1
    assert any("missing from the map" in error for error in result["errors"])


def test_map_validation_catches_duplicate_residue(tmp_path):
    pdb_path = tmp_path / "fixture.pdb"
    _fixture_pdb(pdb_path)
    rows = _canonical_rows()
    duplicate = dict(rows[1])

    result = validate_map.validate_map(pdb_path, rows + [duplicate])

    assert result["validation_status"] == "fail"
    assert result["duplicate_residue_identity_count"] == 1


def test_pymol_helper_contains_selections_for_each_strand(tmp_path):
    text = pymol_helper.build_pml(tmp_path / "fixture.pdb", tmp_path / "map.csv", _canonical_rows())

    assert "select path_path_1" in text
    assert "select path_path_2" in text
    assert "color red" in text
    assert "compare these generated path colors" in text


def test_comparison_handles_missing_manual_map(tmp_path):
    report = compare_maps.comparison_report(tmp_path / "candidate.csv", tmp_path / "missing_manual.csv")

    assert "manual/PyMOL map is not available yet" in report
    assert "manual template" in report


def test_comparison_handles_changed_residue_assignment(tmp_path):
    map_a = tmp_path / "a.csv"
    map_b = tmp_path / "b.csv"
    rows_a = _canonical_rows()
    rows_b = [dict(row) for row in rows_a]
    rows_b[1]["strand_id"] = "2"
    rows_b[1]["strand_label"] = "path_2"
    _write_csv(map_a, convert_map.FIELDNAMES, rows_a)
    _write_csv(map_b, convert_map.FIELDNAMES, rows_b)

    report = compare_maps.comparison_report(map_a, map_b)

    assert "Different or missing assignments: 1" in report
    assert "CYP2" in report
    assert "| 1 | 2 | 1 |" in report
