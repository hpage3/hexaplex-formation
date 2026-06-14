import importlib.util
import sys
from pathlib import Path

import pytest

from hexaplex_formation.pdb_utils import PDBAtom


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


inspect_arm = _load_script_module("inspect_hexad_arm_geometry", "scripts/inspect_hexad_arm_geometry.py")


def atom(serial: int, name: str, residue_name: str, residue_number: int, x: float, y: float, z: float) -> PDBAtom:
    return PDBAtom(
        record_type="ATOM",
        atom_serial=serial,
        atom_name=name,
        alt_loc="",
        residue_name=residue_name,
        chain_id="A",
        residue_number=residue_number,
        insertion_code="",
        x=x,
        y=y,
        z=z,
        occupancy=1.0,
        temp_factor=0.0,
        element=name[0],
    )


def test_atom_classification_distinguishes_backbone_glu_and_arm_atoms():
    assert inspect_arm.atom_classification(atom(1, "CA", "CYP", 1, 0.0, 0.0, 0.0)) == "backbone_like"
    assert inspect_arm.atom_classification(atom(2, "CB", "CYP", 1, 0.0, 0.0, 0.0)) == "candidate_sidechain_or_arm"
    assert inspect_arm.atom_classification(atom(3, "OE1", "GLU", 2, 0.0, 0.0, 0.0)) == "glu_scaffold"


def test_inventory_counts_atom_names_by_residue_type():
    atoms = [
        atom(1, "N", "CYP", 1, 0.0, 0.0, 0.0),
        atom(2, "CB", "CYP", 1, 3.0, 0.0, 0.0),
        atom(3, "N", "CYP", 2, 0.0, 0.0, 1.0),
        atom(4, "CB", "CYP", 2, 3.0, 0.0, 1.0),
    ]

    rows = {(row["residue_name"], row["atom_name"]): row for row in inspect_arm.build_inventory_rows(atoms)}

    assert rows[("CYP", "N")]["atom_count"] == "2"
    assert rows[("CYP", "N")]["residue_occurrence_count"] == "2"
    assert rows[("CYP", "CB")]["classification"] == "candidate_sidechain_or_arm"


def test_geometry_rows_flag_consistently_outward_arm_atoms():
    atoms = [
        atom(1, "N", "CYP", 1, 1.0, 0.0, 0.0),
        atom(2, "CA", "CYP", 1, 1.0, 0.0, 0.5),
        atom(3, "C", "CYP", 1, 1.0, 0.0, 1.0),
        atom(4, "CB", "CYP", 1, 4.0, 0.0, 0.5),
        atom(5, "N", "CYP", 2, 1.0, 0.0, 3.0),
        atom(6, "CA", "CYP", 2, 1.0, 0.0, 3.5),
        atom(7, "C", "CYP", 2, 1.0, 0.0, 4.0),
        atom(8, "CB", "CYP", 2, 4.0, 0.0, 3.5),
    ]
    origin = (0.0, 0.0, 0.0)
    axis = (0.0, 0.0, 1.0)
    radial = inspect_arm.radial_distances(atoms, origin, axis)

    rows = {(row["residue_name"], row["atom_name"]): row for row in inspect_arm.build_geometry_rows(atoms, radial)}

    assert rows[("CYP", "CB")]["outward_candidate"] == "yes"
    assert rows[("CYP", "CB")]["axis_facing_candidate"] == "no"
    assert rows[("CYP", "CB")]["selection_recommendation"] == "transformable_candidate"
    assert rows[("CYP", "N")]["selection_recommendation"] == "fixed_candidate"
    assert float(rows[("CYP", "CB")]["mean_delta_from_residue_backbone_A"]) == pytest.approx(3.0)


def test_geometry_rows_flag_consistently_axis_facing_arm_atoms():
    atoms = [
        atom(1, "N", "MEP", 1, 4.0, 0.0, 0.0),
        atom(2, "CA", "MEP", 1, 4.0, 0.0, 0.5),
        atom(3, "C", "MEP", 1, 4.0, 0.0, 1.0),
        atom(4, "NX1", "MEP", 1, 1.0, 0.0, 0.5),
        atom(5, "N", "MEP", 2, 4.0, 0.0, 3.0),
        atom(6, "CA", "MEP", 2, 4.0, 0.0, 3.5),
        atom(7, "C", "MEP", 2, 4.0, 0.0, 4.0),
        atom(8, "NX1", "MEP", 2, 1.0, 0.0, 3.5),
    ]
    radial = inspect_arm.radial_distances(atoms, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0))

    rows = {(row["residue_name"], row["atom_name"]): row for row in inspect_arm.build_geometry_rows(atoms, radial)}

    assert rows[("MEP", "NX1")]["outward_candidate"] == "no"
    assert rows[("MEP", "NX1")]["axis_facing_candidate"] == "yes"
    assert rows[("MEP", "NX1")]["selection_recommendation"] == "transformable_candidate"


def test_infer_axis_uses_backbone_like_atoms():
    atoms = [
        atom(1, "N", "GLU", 1, 0.0, 0.0, 0.0),
        atom(2, "CA", "GLU", 1, 0.0, 0.0, 1.0),
        atom(3, "C", "GLU", 1, 0.0, 0.0, 2.0),
        atom(4, "CB", "CYP", 2, 5.0, 0.0, 1.0),
    ]

    _origin, axis, source = inspect_arm.infer_axis(atoms)

    assert source == "backbone_like_heavy_atoms"
    assert abs(axis[2]) > 0.9
