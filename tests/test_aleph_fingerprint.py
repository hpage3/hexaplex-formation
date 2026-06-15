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


def test_angle_wrapping_and_unwrapping():
    assert aleph.wrap_degrees(190.0) == pytest.approx(-170.0)
    assert aleph.wrap_degrees(-190.0) == pytest.approx(170.0)
    raw = [math.radians(170.0), math.radians(-170.0), math.radians(-150.0)]
    unwrapped = aleph.unwrap_radians(raw)

    assert aleph.phase_unwrap_ok(raw, unwrapped)
    assert math.degrees(unwrapped[1] - unwrapped[0]) == pytest.approx(20.0)
    assert math.degrees(unwrapped[2] - unwrapped[1]) == pytest.approx(20.0)


def test_axis_flip_logic_makes_unit_order_positive():
    axis, flipped = aleph.orient_axis_by_unit_order((0.0, 0.0, -1.0), (0.0, 0.0, 0.0), [(0.0, 0.0, 0.0), (0.0, 0.0, 3.4)])

    assert flipped is True
    assert axis == pytest.approx((0.0, 0.0, 1.0))


def test_circular_dispersion_is_bounded_and_resultant_is_interpretable():
    angles = [0.0, math.pi / 2, math.pi, 3 * math.pi / 2]
    dispersion = aleph.circular_dispersion_deg(angles)
    resultant = aleph.circular_mean_resultant_length(angles)

    assert dispersion is not None
    assert 0.0 <= dispersion <= 180.0
    assert resultant == pytest.approx(0.0, abs=1e-12)


def test_local_twist_and_rise_on_simple_geometry(tmp_path):
    path = tmp_path / "synthetic.pdb"
    synthetic_hexaplex(path, 5)
    result = aleph.analyze_structure("synthetic", path)
    twists = [float(row["local_twist_deg"]) for row in result.per_unit_rows if row["local_twist_deg"]]
    rises = [float(row["local_rise_A"]) for row in result.per_unit_rows if row["local_rise_A"]]

    assert len(twists) == 4
    assert sum(twists) / len(twists) == pytest.approx(30.0, abs=5.0)
    assert sum(rises) / len(rises) == pytest.approx(3.4, abs=0.5)
    assert all(rise > 0.0 for rise in rises)


def test_fft_summary_handles_short_signal():
    row = aleph.fft_summary("short", "local_twist_deg", [1.0, 2.0, 3.0])

    assert row["too_short_for_reliable_interpretation"] == "true"
    assert "too short" in row["warnings"]


def test_feature_normalization_handles_missing_values():
    assert aleph.normalize_feature_value("", 0.0, 1.0) is None
    assert aleph.normalize_feature_value("bad", 0.0, 1.0) is None
    assert aleph.normalize_feature_value("5", 0.0, 10.0) == pytest.approx(0.5)
    assert aleph.normalize_feature_value("370", None, None, cyclic=True) == pytest.approx(10.0 / 360.0)


def test_qc_warning_row_is_created():
    row = {"local_twist_warning": "", "local_rise_warning": "negative rise", "plane_fit_warning": "", "warnings": "", "missing_chain_count": "0"}

    assert aleph.row_has_qc_warning(row)


def test_fingerprint_strip_plot_writes_svg(tmp_path):
    rows = [
        {
            "structure_id": "central7",
            "unit_index": "1",
            "local_abs_twist_deg": "30",
            "local_twist_deg": "-30",
            "local_rise_A": "3.4",
            "base_radial_spread_A": "1.0",
            "aleph_phase_deg": "45",
            "aleph_base_plane_bend_deg": "10",
            "aleph_scaffold_plane_bend_deg": "12",
            "chain_mean_resultant_length": "0.7",
            "local_twist_warning": "",
            "local_rise_warning": "",
            "plane_fit_warning": "",
            "warnings": "",
            "missing_chain_count": "0",
        }
    ]
    path = tmp_path / "strip.svg"

    aleph.svg_fingerprint_strip("central7", rows, path)

    assert path.exists()
    assert path.read_text(encoding="utf-8").startswith("<svg")


def test_fingerprint_comparison_plot_writes_svg(tmp_path):
    rows = []
    for structure_id in ["full", "central6", "central7"]:
        rows.append(
            {
                "structure_id": structure_id,
                "unit_index": "1",
                "local_abs_twist_deg": "30",
                "local_twist_deg": "-30",
                "local_rise_A": "3.4",
                "base_radial_spread_A": "1.0",
                "aleph_phase_deg": "45",
                "aleph_base_plane_bend_deg": "10",
                "aleph_scaffold_plane_bend_deg": "12",
                "chain_mean_resultant_length": "0.7",
                "local_twist_warning": "",
                "local_rise_warning": "",
                "plane_fit_warning": "",
                "warnings": "",
                "missing_chain_count": "0",
            }
        )
    path = tmp_path / "comparison.svg"

    aleph.svg_fingerprint_comparison(rows, path)

    assert path.exists()
    assert "Aleph fingerprint comparison" in path.read_text(encoding="utf-8")


def test_series_points_use_signed_local_twist():
    rows = [
        {
            "structure_id": "central7",
            "unit_index": "1",
            "local_twist_deg": "-30",
            "local_abs_twist_deg": "30",
            "local_rise_A": "3.4",
            "local_twist_warning": "",
            "local_rise_warning": "",
            "plane_fit_warning": "",
            "warnings": "",
            "missing_chain_count": "0",
        }
    ]

    points = aleph.series_points(rows, "central7", "local_twist_deg")

    assert points == [(1.0, -30.0, False)]


def test_series_fingerprint_plot_writes_svg(tmp_path):
    rows = [
        {
            "structure_id": "central7",
            "unit_index": "1",
            "local_twist_deg": "-30",
            "local_abs_twist_deg": "30",
            "local_rise_A": "3.4",
            "aleph_base_plane_bend_deg": "10",
            "aleph_scaffold_plane_bend_deg": "12",
            "local_twist_warning": "",
            "local_rise_warning": "",
            "plane_fit_warning": "",
            "warnings": "",
            "missing_chain_count": "0",
        }
    ]
    path = tmp_path / "series.svg"

    aleph.svg_series_fingerprint("central7", rows, path)

    assert path.exists()
    assert "Aleph series fingerprint" in path.read_text(encoding="utf-8")


def test_companion_traces_handle_missing_short_data(tmp_path):
    rows = [
        {
            "structure_id": "full",
            "unit_index": "1",
            "local_twist_deg": "-12",
            "local_abs_twist_deg": "12",
            "local_rise_A": "",
            "aleph_base_plane_bend_deg": "20",
            "aleph_scaffold_plane_bend_deg": "",
            "local_twist_warning": "",
            "local_rise_warning": "",
            "plane_fit_warning": "",
            "warnings": "",
            "missing_chain_count": "0",
        }
    ]
    path = tmp_path / "companions.svg"

    aleph.svg_companion_traces("full", rows, path)

    assert path.exists()
    assert "companion traces" in path.read_text(encoding="utf-8")


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
            "qc_csv": tmp_path / "qc.csv",
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
    assert args.qc_csv.exists()
    assert args.report.exists()
    with args.per_unit_csv.open(newline="", encoding="utf-8") as handle:
        assert csv.DictReader(handle).fieldnames == aleph.PER_UNIT_COLUMNS
    with args.qc_csv.open(newline="", encoding="utf-8") as handle:
        assert csv.DictReader(handle).fieldnames == aleph.QC_COLUMNS
