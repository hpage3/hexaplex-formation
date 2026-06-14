import csv
import importlib.util
import sys
from pathlib import Path

from hexaplex_formation.pdb_utils import PDBAtom, write_pdb_atoms


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_script_module(name: str, relative_path: str):
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


attrib = _load_script_module("analyze_distance_window_attribution", "scripts/analyze_distance_window_attribution.py")


def atom(serial: int, name: str, residue: str, residue_number: int, chain: str, x: float, y: float, z: float, element: str = "C") -> PDBAtom:
    return PDBAtom(
        record_type="ATOM",
        atom_serial=serial,
        atom_name=name,
        alt_loc="",
        residue_name=residue,
        chain_id=chain,
        residue_number=residue_number,
        insertion_code="",
        x=x,
        y=y,
        z=z,
        occupancy=1.0,
        temp_factor=0.0,
        element=element,
    )


def synthetic_structure(path: Path) -> None:
    atoms = [
        atom(1, "CB", "CYP", 1, "A", 0.0, 0.0, 0.0),
        atom(2, "CG", "MEP", 1, "B", 3.0, 0.0, 0.0),
        atom(3, "OE1", "GLU", 2, "B", 0.0, 3.4, 0.0, "O"),
        atom(4, "N", "GLU", 2, "A", 4.7, 0.0, 0.0, "N"),
        atom(5, "H", "CYP", 1, "A", 3.0, 0.0, 0.0, "H"),
    ]
    write_pdb_atoms(atoms, path)


def test_residue_and_atom_classes():
    cyp = atom(1, "CB", "CYP", 1, "A", 0, 0, 0)
    mep = atom(2, "CG", "MEP", 1, "B", 0, 0, 0)
    glu = atom(3, "OE1", "GLU", 2, "B", 0, 0, 0, "O")
    backbone = atom(4, "CA", "CYP", 1, "A", 0, 0, 0)

    assert attrib.residue_pair_class(cyp, mep) == "CYP_MEP_vs_CYP_MEP"
    assert attrib.residue_pair_class(cyp, glu) == "CYP_MEP_vs_GLU"
    assert attrib.residue_pair_class(glu, glu) == "GLU_vs_GLU"
    assert attrib.atom_role(cyp) == "base_like"
    assert attrib.atom_role(glu) == "scaffold_linker"
    assert attrib.atom_role(backbone) == "backbone_like"


def test_sequence_relation():
    assert attrib.sequence_relation(1, 1) == "same_unit"
    assert attrib.sequence_relation(1, 2) == "adjacent_unit"
    assert attrib.sequence_relation(1, 3) == "near_axial"
    assert attrib.sequence_relation(1, 5) == "longer_range"
    assert attrib.sequence_relation(None, 1) == "unknown"


def test_pair_rows_exclude_hydrogen_and_same_residue(tmp_path):
    path = tmp_path / "structure.pdb"
    synthetic_structure(path)
    prepared = attrib.prepare_structure("test", path)

    rows = attrib.pair_rows_for_structure(prepared)

    assert {row["distance_window"] for row in rows} >= {"3p0", "3p4", "4p5_5p0"}
    assert all(row["atom_i_name"] != "H" and row["atom_j_name"] != "H" for row in rows)
    assert any(row["residue_pair_class"] == "CYP_MEP_vs_CYP_MEP" for row in rows)
    assert any(row["residue_pair_class"] == "CYP_MEP_vs_GLU" for row in rows)


def test_run_writes_expected_outputs(tmp_path):
    full = tmp_path / "full.pdb"
    c6 = tmp_path / "c6.pdb"
    c7 = tmp_path / "c7.pdb"
    for path in [full, c6, c7]:
        synthetic_structure(path)

    args = type(
        "Args",
        (),
        {
            "full_pdb": full,
            "central6_pdb": c6,
            "central7_pdb": c7,
            "summary_csv": tmp_path / "summary.csv",
            "pairs_sample_csv": tmp_path / "pairs.csv",
            "by_chain_pair_csv": tmp_path / "chain.csv",
            "by_residue_class_csv": tmp_path / "residue.csv",
            "report": tmp_path / "report.md",
            "plot_dir": tmp_path / "plots",
            "pairs_sample_limit": 10,
        },
    )()

    result = attrib.run(args)

    assert result["pair_rows"] > 0
    assert args.summary_csv.exists()
    assert args.pairs_sample_csv.exists()
    assert args.by_chain_pair_csv.exists()
    assert args.by_residue_class_csv.exists()
    assert args.report.exists()
    assert result["plots"]
    with args.summary_csv.open(newline="", encoding="utf-8") as handle:
        assert csv.DictReader(handle).fieldnames == attrib.SUMMARY_COLUMNS
