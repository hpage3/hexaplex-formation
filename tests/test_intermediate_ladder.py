import csv
import importlib.util
from pathlib import Path

from hexaplex_formation.pdb_utils import PDBAtom


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_script_module(name: str, relative_path: str):
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


strand_map = _load_script_module("build_strand_map_candidate", "scripts/build_strand_map_candidate.py")
ladder = _load_script_module("build_intermediate_ladder", "scripts/build_intermediate_ladder.py")
compare_ladder = _load_script_module("compare_intermediate_ladder", "scripts/compare_intermediate_ladder.py")


def atom(serial: int, residue_name: str, residue_number: int, chain_id: str = "A") -> PDBAtom:
    return PDBAtom(
        record_type="ATOM",
        atom_serial=serial,
        atom_name="CA",
        alt_loc="",
        residue_name=residue_name,
        chain_id=chain_id,
        residue_number=residue_number,
        insertion_code="",
        x=float(serial),
        y=0.0,
        z=0.0,
        occupancy=1.0,
        temp_factor=0.0,
        element="C",
    )


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def test_block_assignment_from_twelve_residues_into_three_blocks():
    residue_keys = [("A", "GLU", residue_number, "") for residue_number in range(1, 13)]

    rows = strand_map.build_rows_from_residue_keys(residue_keys, blocks=3)

    assert [row["block_id"] for row in rows] == ["1"] * 4 + ["2"] * 4 + ["3"] * 4
    assert rows[0]["residue_label"] == "A:GLU1"
    assert rows[-1]["residue_label"] == "A:GLU12"
    assert strand_map.residues_per_block(rows) == {"1": 4, "2": 4, "3": 4}


def test_parse_strand_map_candidate_csv(tmp_path):
    path = tmp_path / "strand_map_candidate.csv"
    _write_csv(
        path,
        strand_map.FIELDNAMES,
        [
            {
                "block_id": "1",
                "residue_index_in_pdb_order": "1",
                "chain_id": "A",
                "residue_name": "GLU",
                "residue_number": "10",
                "insertion_code": "",
                "residue_label": "A:GLU10",
            },
            {
                "block_id": "2",
                "residue_index_in_pdb_order": "2",
                "chain_id": "A",
                "residue_name": "ALA",
                "residue_number": "11",
                "insertion_code": "",
                "residue_label": "A:ALA11",
            },
        ],
    )

    parsed = ladder.read_strand_map(path)

    assert parsed == {
        1: {("A", "GLU", 10, "")},
        2: {("A", "ALA", 11, "")},
    }


def test_select_residues_by_block_preserves_atom_order():
    atoms = [
        atom(1, "GLU", 1),
        atom(2, "ALA", 2),
        atom(3, "GLU", 3),
        atom(4, "ALA", 4),
    ]
    blocks_by_id = {
        1: {("A", "GLU", 1, ""), ("A", "ALA", 2, "")},
        2: {("A", "GLU", 3, ""), ("A", "ALA", 4, "")},
    }

    selected = ladder.select_atoms_by_blocks(atoms, blocks_by_id, [2])

    assert [selected_atom.atom_serial for selected_atom in selected] == [3, 4]


def test_ladder_model_naming():
    assert (
        ladder.model_blocks_name("hexads_plus_scaffold_blocks", [1, 2, 3], "heavy_deduped")
        == "hexads_plus_scaffold_blocks_1_2_3_heavy_deduped"
    )


def test_ladder_report_generation_includes_caution_language(tmp_path):
    out = tmp_path / "intermediate_ladder_report.md"
    rows = [
        {
            "model_name": "scaffold_blocks_1_heavy_deduped",
            "atom_mode": "heavy_deduped",
            "included_blocks": "1",
            "includes_hexads": "no",
            "atom_count": "10",
            "residue_count": "3",
            "contact_count_4p5A": "1",
            "motif_GLU_GLU": "0",
            "motif_GLU_any": "2",
            "hbond_candidate_count": "0",
            "mean_radius_xy": "1.0",
            "z_span": "2.0",
            "angular_coverage_rad": "2.5",
        },
        {
            "model_name": "scaffold_blocks_1_2_heavy_deduped",
            "atom_mode": "heavy_deduped",
            "included_blocks": "1,2",
            "includes_hexads": "no",
            "atom_count": "20",
            "residue_count": "6",
            "contact_count_4p5A": "3",
            "motif_GLU_GLU": "1",
            "motif_GLU_any": "4",
            "hbond_candidate_count": "0",
            "mean_radius_xy": "1.0",
            "z_span": "5.0",
            "angular_coverage_rad": "5.6",
        },
        {
            "model_name": "hexads_plus_scaffold_blocks_1_2_heavy_deduped",
            "atom_mode": "heavy_deduped",
            "included_blocks": "1,2",
            "includes_hexads": "yes",
            "atom_count": "30",
            "residue_count": "10",
            "contact_count_4p5A": "7",
            "motif_GLU_GLU": "1",
            "motif_GLU_any": "4",
            "hbond_candidate_count": "0",
            "mean_radius_xy": "2.0",
            "z_span": "6.0",
            "angular_coverage_rad": "5.8",
        },
    ]

    compare_ladder.write_markdown_report(rows, out)

    text = out.read_text(encoding="utf-8")
    assert "candidate intermediate" in text
    assert "consistent with" in text
    assert "does not prove temporal assembly order" in text
    assert "requires validation against PyMOL strand mapping and/or simulation" in text
