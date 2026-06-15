import csv
import importlib.util
import math
import sys
from pathlib import Path

import pytest

from hexaplex_formation.pdb_utils import PDBAtom, write_pdb_atoms


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_script_module(name: str, relative_path: str):
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


aleph = _load_script_module("analyze_aleph_fingerprint", "scripts/analyze_aleph_fingerprint.py")


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


def synthetic_hexaplex(path: Path, unit_count: int = 4) -> None:
    atoms = []
    serial = 1
    for chain_index, chain in enumerate("ABCDEF"):
        chain_phase = chain_index * math.tau / 6.0
        for unit_index in range(unit_count):
            phase = chain_phase + unit_index * math.radians(30.0)
            z = unit_index * 3.4
            x = 5.0 * math.cos(phase)
            y = 5.0 * math.sin(phase)
            residue_number = unit_index * 2 + 1
            residue = "CYP" if chain_index % 2 == 0 else "MEP"
            atoms.append(atom(serial, "CA", residue, residue_number, chain, x, y, z))
            serial += 1
            atoms.append(atom(serial, "CB", residue, residue_number, chain, x + 0.2, y, z))
            serial += 1
            atoms.append(atom(serial, "N", "GLU", residue_number + 1, chain, x * 1.1, y * 1.1, z + 0.2, "N"))
            serial += 1
            atoms.append(atom(serial, "OE1", "GLU", residue_number + 1, chain, x * 1.1 + 0.1, y * 1.1, z + 0.2, "O"))
            serial += 1
    write_pdb_atoms(atoms, path)


def test_atom_classification():
    assert aleph.atom_class(atom(1, "CA", "CYP", 1, "A", 0, 0, 0)) == "backbone_like"
    assert aleph.atom_class(atom(2, "CB", "MEP", 1, "A", 0, 0, 0)) == "base_like"
    assert aleph.atom_class(atom(3, "OE1", "GLU", 2, "A", 0, 0, 0, "O")) == "scaffold_linker"


def test_axis_fitting_returns_unit_vector():
    _, axis = aleph.infer_axis([(0, 0, 0), (0, 0, 1), (0, 0, 2), (0.1, 0, 3)])
    assert math.sqrt(sum(value * value for value in axis)) == pytest.approx(1.0)


def test_local_twist_and_rise_on_simple_geometry(tmp_path):
    path = tmp_path / "synthetic.pdb"
    synthetic_hexaplex(path, 5)
    result = aleph.analyze_structure("synthetic", path)
    twists = [float(row["local_twist_deg"]) for row in result.per_unit_rows if row["local_twist_deg"]]
    rises = [float(row["local_rise_A"]) for row in result.per_unit_rows if row["local_rise_A"]]

    assert len(twists) == 4
    assert sum(twists) / len(twists) == pytest.approx(30.0, abs=5.0)
    assert sum(rises) / len(rises) == pytest.approx(3.4, abs=0.5)


def test_fft_summary_handles_short_signal():
    row = aleph.fft_summary("short", "local_twist_deg", [1.0, 2.0, 3.0])

    assert row["too_short_for_reliable_interpretation"] == "true"
    assert "too short" in row["warnings"]


def test_output_schemas_are_written(tmp_path):
    full = tmp_path / "full.pdb"
    c6 = tmp_path / "c6.pdb"
    c7 = tmp_path / "c7.pdb"
    for path in [full, c6, c7]:
        synthetic_hexaplex(path, 4)
    args = type(
        "Args",
        (),
        {
            "per_unit_csv": tmp_path / "per_unit.csv",
            "summary_csv": tmp_path / "summary.csv",
            "fft_csv": tmp_path / "fft.csv",
            "report": tmp_path / "report.md",
            "plot_dir": tmp_path / "plots",
            "include_optional_twists": False,
        },
    )()
    original = aleph.DEFAULT_STRUCTURES
    aleph.DEFAULT_STRUCTURES = [("full", full), ("central6", c6), ("central7", c7)]
    try:
        result = aleph.run(args)
    finally:
        aleph.DEFAULT_STRUCTURES = original

    assert result["per_unit_rows"] == 12
    assert args.per_unit_csv.exists()
    assert args.summary_csv.exists()
    assert args.fft_csv.exists()
    assert args.report.exists()
    with args.per_unit_csv.open(newline="", encoding="utf-8") as handle:
        assert csv.DictReader(handle).fieldnames == aleph.PER_UNIT_COLUMNS
