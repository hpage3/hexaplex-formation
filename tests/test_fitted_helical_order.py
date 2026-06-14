import importlib.util
import math
from pathlib import Path

import pytest

from hexaplex_formation.geometry import (
    covariance_matrix_3d,
    cross,
    dot,
    norm,
    power_iteration_principal_axis,
)
from hexaplex_formation.pdb_utils import PDBAtom


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_script_module(name: str, relative_path: str):
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


fitted = _load_script_module("helical_order_fitted_axis", "scripts/helical_order_fitted_axis.py")


def atom(serial: int, residue_number: int, x: float, y: float, z: float) -> PDBAtom:
    return PDBAtom(
        record_type="ATOM",
        atom_serial=serial,
        atom_name="CA",
        alt_loc="",
        residue_name="GLU",
        chain_id="A",
        residue_number=residue_number,
        insertion_code="",
        x=x,
        y=y,
        z=z,
        occupancy=1.0,
        temp_factor=0.0,
        element="C",
    )


def test_covariance_matrix_on_simple_points():
    cov = covariance_matrix_3d([(-1.0, 0.0, 0.0), (1.0, 0.0, 0.0)])

    assert cov[0][0] == pytest.approx(1.0)
    assert cov[1][1] == pytest.approx(0.0)
    assert cov[2][2] == pytest.approx(0.0)


def test_principal_axis_for_points_along_z():
    cov = covariance_matrix_3d([(0.0, 0.0, -2.0), (0.0, 0.0, 0.0), (0.0, 0.0, 2.0)])

    axis = power_iteration_principal_axis(cov)

    assert axis[0] == pytest.approx(0.0, abs=1e-6)
    assert axis[1] == pytest.approx(0.0, abs=1e-6)
    assert axis[2] == pytest.approx(1.0, abs=1e-6)


def test_principal_axis_for_points_along_x():
    cov = covariance_matrix_3d([(-2.0, 0.0, 0.0), (0.0, 0.0, 0.0), (2.0, 0.0, 0.0)])

    axis = power_iteration_principal_axis(cov)

    assert axis[0] == pytest.approx(1.0, abs=1e-6)
    assert axis[1] == pytest.approx(0.0, abs=1e-6)
    assert axis[2] == pytest.approx(0.0, abs=1e-6)


def test_perpendicular_basis_orthogonality():
    basis_u, basis_v = fitted.build_perpendicular_basis((0.0, 0.0, 1.0))

    assert norm(basis_u) == pytest.approx(1.0)
    assert norm(basis_v) == pytest.approx(1.0)
    assert dot(basis_u, (0.0, 0.0, 1.0)) == pytest.approx(0.0)
    assert dot(basis_v, (0.0, 0.0, 1.0)) == pytest.approx(0.0)
    assert dot(basis_u, basis_v) == pytest.approx(0.0)
    assert norm(cross(basis_u, basis_v)) == pytest.approx(1.0)


def test_angle_unwrapping_behavior():
    angles = [3.0, -3.0, -2.5]

    unwrapped = fitted.unwrap_angles(angles)

    assert unwrapped[1] > unwrapped[0]
    assert unwrapped[1] == pytest.approx(-3.0 + 2.0 * math.pi)
    assert unwrapped[2] == pytest.approx(-2.5 + 2.0 * math.pi)


def test_fitted_axis_summary_on_synthetic_helix_like_fixture():
    atoms = []
    for index in range(12):
        theta = index * math.pi / 3.0
        atoms.append(atom(index + 1, index + 1, math.cos(theta), math.sin(theta), index * 0.5))

    rows, fit = fitted.fitted_helical_rows(atoms, use_heavy_only=True)
    summary = fitted.summary_row(Path("synthetic.pdb"), rows, fit)

    assert len(rows) == 12
    assert summary["residue_count"] == "12"
    assert float(summary["mean_radius_fitted"]) > 0.5
    assert float(summary["axial_span"]) > 4.0
    assert float(summary["approximate_turns"]) > 1.0
    assert summary["approximate_pitch_per_turn"]
