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


attrib = _load_script_module("analyze_backbone_distance_window_attribution", "scripts/analyze_backbone_distance_window_attribution.py")


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
        atom(1, "CA", "CYP", 1, "A", 0.0, 0.0, 0.0),
        atom(2, "CB", "MEP", 1, "B", 3.0, 0.0, 0.0),
        atom(3, "OE1", "GLU", 2, "B", 0.0, 5.5, 0.0, "O"),
        atom(4, "N", "GLU", 2, "A", 7.0, 0.0, 0.0, "N"),
        atom(5, "CG", "CYP", 3, "A", 8.4, 0.0, 0.0),
        atom(6, "H", "CYP", 3, "A", 8.4, 0.0, 0.1, "H"),
    ]
    write_pdb_atoms(atoms, path)


def test_atom_and_pair_classification():
    backbone = atom(1, "CA", "CYP", 1, "A", 0, 0, 0)
    base = atom(2, "CB", "MEP", 1, "B", 0, 0, 0)
    scaffold = atom(3, "OE1", "GLU", 2, "B", 0, 0, 0, "O")

    assert attrib.atom_class(backbone) == "backbone_like"
    assert attrib.atom_class(base) == "base_like"
    assert attrib.atom_class(scaffold) == "scaffold_linker"
    assert attrib.atom_pair_class("backbone_like", "base_like") == "backbone_like_vs_base_like"
    assert attrib.atom_pair_class("base_like", "scaffold_linker") == "base_like_vs_scaffold_linker"


def test_pair_rows_cover_known_and_guessed_windows(tmp_path):
    path = tmp_path / "structure.pdb"
    synthetic_structure(path)
    prepared = attrib.prepare_structure("test", path)
    rows = attrib.pair_rows_for_structure(prepared)

    windows = {row["distance_window"] for row in rows}
    assert "3p0" in windows
    assert "guessed_5p0_6p0" in windows
    assert "guessed_6p5_7p5" in windows
    assert "guessed_7p8_9p0" in windows
    assert all(row["atom_i_name"] != "H" and row["atom_j_name"] != "H" for row in rows)
    assert any(row["backbone_involving"] == "yes" for row in rows)


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
            "by_atom_class_csv": tmp_path / "atom.csv",
            "by_chain_relation_csv": tmp_path / "chain.csv",
            "pairs_sample_csv": tmp_path / "pairs.csv",
            "report": tmp_path / "report.md",
            "plot_dir": tmp_path / "plots",
            "pairs_sample_limit": 20,
        },
    )()

    result = attrib.run(args)

    assert result["pair_rows"] > 0
    assert args.summary_csv.exists()
    assert args.by_atom_class_csv.exists()
    assert args.by_chain_relation_csv.exists()
    assert args.pairs_sample_csv.exists()
    assert args.report.exists()
    assert result["plots"]
    with args.summary_csv.open(newline="", encoding="utf-8") as handle:
        assert csv.DictReader(handle).fieldnames == attrib.SUMMARY_COLUMNS
