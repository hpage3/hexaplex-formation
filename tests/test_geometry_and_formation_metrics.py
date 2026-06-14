import importlib.util
from pathlib import Path

import pytest

from hexaplex_formation.geometry import (
    distance,
    group_atoms_by_residue,
    min_inter_residue_distance,
)
from hexaplex_formation.pdb_utils import PDBAtom


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_script_module(name: str, relative_path: str):
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


contact_map = _load_script_module("contact_map", "scripts/contact_map.py")
hbond_candidates = _load_script_module("hbond_candidates", "scripts/hbond_candidates.py")


def atom(
    serial: int,
    atom_name: str,
    residue_name: str,
    residue_number: int,
    x: float,
    y: float,
    z: float,
    element: str,
    chain_id: str = "A",
) -> PDBAtom:
    return PDBAtom(
        record_type="ATOM",
        atom_serial=serial,
        atom_name=atom_name,
        alt_loc="",
        residue_name=residue_name,
        chain_id=chain_id,
        residue_number=residue_number,
        insertion_code="",
        x=x,
        y=y,
        z=z,
        occupancy=1.0,
        temp_factor=0.0,
        element=element,
    )


def test_geometry_distance():
    assert distance((0.0, 0.0, 0.0), (1.0, 2.0, 2.0)) == pytest.approx(3.0)


def test_group_atoms_by_residue():
    atoms = [
        atom(1, "N", "GLU", 1, 0.0, 0.0, 0.0, "N"),
        atom(2, "CA", "GLU", 1, 1.0, 0.0, 0.0, "C"),
        atom(3, "CA", "ALA", 2, 5.0, 0.0, 0.0, "C"),
    ]
    grouped = group_atoms_by_residue(atoms)
    assert list(grouped) == [("A", "GLU", 1, ""), ("A", "ALA", 2, "")]
    assert [len(values) for values in grouped.values()] == [2, 1]


def test_min_inter_residue_distance():
    residue_a = [
        atom(1, "N", "GLU", 1, 0.0, 0.0, 0.0, "N"),
        atom(2, "CA", "GLU", 1, 1.0, 0.0, 0.0, "C"),
    ]
    residue_b = [
        atom(3, "CA", "ALA", 2, 5.0, 0.0, 0.0, "C"),
        atom(4, "O", "ALA", 2, 2.5, 0.0, 0.0, "O"),
    ]
    min_distance, atom_a, atom_b = min_inter_residue_distance(residue_a, residue_b)
    assert min_distance == pytest.approx(1.5)
    assert atom_a.atom_name == "CA"
    assert atom_b.atom_name == "O"


def test_contact_map_logic_on_tiny_fixture():
    atoms = [
        atom(1, "OE1", "GLU", 1, 0.0, 0.0, 0.0, "O"),
        atom(2, "CA", "ALA", 2, 3.0, 0.0, 0.0, "C"),
        atom(3, "CA", "ALA", 3, 10.0, 0.0, 0.0, "C"),
    ]
    rows = contact_map.contact_rows(atoms, cutoff=4.5, use_heavy_only=True)
    assert len(rows) == 1
    assert rows[0]["residue_i"] == "A:GLU1"
    assert rows[0]["residue_j"] == "A:ALA2"
    assert rows[0]["atom_i"] == "OE1"
    assert rows[0]["atom_j"] == "CA"


def test_hydrogen_parent_inference_detects_simple_candidate():
    atoms = [
        atom(1, "N", "GLU", 1, 0.0, 0.0, 0.0, "N"),
        atom(2, "H", "GLU", 1, 1.0, 0.0, 0.0, "H"),
        atom(3, "O", "ALA", 2, 2.6, 0.0, 0.0, "O"),
    ]
    donor = hbond_candidates.infer_parent_donor(atoms[1], hbond_candidates.donor_heavy_atoms(atoms))
    assert donor == atoms[0]

    rows = hbond_candidates.hbond_candidate_rows(atoms)
    assert len(rows) == 1
    assert rows[0]["donor_residue"] == "A:GLU1"
    assert rows[0]["donor_atom"] == "N"
    assert rows[0]["hydrogen_atom"] == "H"
    assert rows[0]["acceptor_residue"] == "A:ALA2"
