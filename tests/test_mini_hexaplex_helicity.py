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


helicity = _load_script_module("analyze_mini_hexaplex_helicity", "scripts/analyze_mini_hexaplex_helicity.py")


def atom(
    serial: int,
    name: str,
    residue_name: str,
    residue_number: int,
    chain_id: str,
    x: float,
    y: float,
    z: float,
) -> PDBAtom:
    return PDBAtom(
        record_type="ATOM",
        atom_serial=serial,
        atom_name=name,
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
        element=name[0],
    )


def synthetic_hexaplex(unit_count: int = 6, helical: bool = True) -> list[PDBAtom]:
    atoms: list[PDBAtom] = []
    serial = 1
    for chain_index, chain_id in enumerate("ABCDEF"):
        base_name = "CYP" if chain_index % 2 == 0 else "MEP"
        phase = chain_index * (2.0 * math.pi / 6.0)
        for unit_index in range(unit_count):
            z = unit_index * 3.2
            if helical:
                theta = phase + unit_index * math.radians(28.0)
            else:
                theta = phase + math.sin(unit_index * 1.7) * 1.1
            x = 10.0 * math.cos(theta)
            y = 10.0 * math.sin(theta)
            residue_number = unit_index * 2 + 1
            atoms.append(atom(serial, "CA", base_name, residue_number, chain_id, x, y, z))
            serial += 1
            atoms.append(atom(serial, "C", base_name, residue_number, chain_id, x + 0.2, y, z))
            serial += 1
            atoms.append(atom(serial, "N", "GLU", residue_number + 1, chain_id, x, y + 0.2, z + 0.2))
            serial += 1
    return atoms


def test_cylindrical_coordinate_transform_on_synthetic_helix_point():
    origin = (0.0, 0.0, 0.0)
    axis = (0.0, 0.0, 1.0)
    basis_u, basis_v = helicity.build_perpendicular_basis(axis)
    point = (2.0 * basis_u[0] + 5.0 * axis[0], 2.0 * basis_u[1] + 5.0 * axis[1], 2.0 * basis_u[2] + 5.0 * axis[2])

    z, theta, radius = helicity.cylindrical_coordinates(point, origin, axis, basis_u, basis_v)

    assert z == pytest.approx(5.0)
    assert radius == pytest.approx(2.0)
    assert theta == pytest.approx(0.0, abs=1e-6)


def test_unwrap_angles_across_two_pi_boundary():
    angles = [math.radians(170), math.radians(179), math.radians(-175), math.radians(-166)]

    unwrapped = helicity.unwrap_angles(angles)

    assert [math.degrees(value) for value in unwrapped] == pytest.approx([170, 179, 185, 194])


def test_linear_theta_z_fit_is_high_for_synthetic_helix_and_lower_for_nonhelix():
    z_values = [float(index) for index in range(8)]
    theta_values = [0.4 * z + 0.2 for z in z_values]
    noisy_values = [math.sin(index * 1.7) for index in range(8)]

    assert helicity.coefficient_of_determination(z_values, theta_values) > 0.999
    assert helicity.coefficient_of_determination(z_values, noisy_values) < 0.4


def test_script_writes_helicity_summary_csv(tmp_path):
    structures_dir = tmp_path / "structures"
    structures_dir.mkdir()
    baseline_pdb = tmp_path / "baseline.pdb"
    variant_pdb = structures_dir / "mini_hexaplex_central6_units.pdb"
    manifest = tmp_path / "manifest.csv"
    geometry = tmp_path / "geometry.csv"
    feature_summary = tmp_path / "features.csv"
    out_csv = tmp_path / "helicity.csv"
    plot_dir = tmp_path / "plots"
    report = tmp_path / "report.md"
    write_pdb_atoms(synthetic_hexaplex(unit_count=8), baseline_pdb)
    write_pdb_atoms(synthetic_hexaplex(unit_count=6), variant_pdb)
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "variant_id",
                "truncation_rule",
                "units_per_chain",
                "residues_per_chain",
                "total_residue_count",
                "total_atom_count",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "variant_id": "central6_units",
                "truncation_rule": "synthetic central six units",
                "units_per_chain": "6",
                "residues_per_chain": "A:12;B:12;C:12;D:12;E:12;F:12",
                "total_residue_count": "72",
                "total_atom_count": "108",
            }
        )
    with geometry.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["variant_id", "units_per_chain", "axial_extent_A", "structural_coherence_flag"])
        writer.writeheader()
        writer.writerow(
            {
                "variant_id": "central6_units",
                "units_per_chain": "6",
                "axial_extent_A": "19.200000",
                "structural_coherence_flag": "borderline",
            }
        )
    with feature_summary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["variant_id", "ratio_to_full_length_4p5_5A"])
        writer.writeheader()
        writer.writerow({"variant_id": "central6_units", "ratio_to_full_length_4p5_5A": "0.5"})

    axis_origin, axis = helicity.infer_axis(helicity.load_pdb_atoms(baseline_pdb))
    rows = [
        helicity.analyze_structure(
            "full_length_baseline",
            baseline_pdb,
            axis_origin,
            axis,
            {},
            {},
        ),
        helicity.analyze_structure(
            "central6_units",
            variant_pdb,
            axis_origin,
            axis,
            {"central6_units": {"units_per_chain": "6", "truncation_rule": "synthetic central six units"}},
            {"central6_units": {"axial_extent_A": "19.200000", "structural_coherence_flag": "borderline"}},
        ),
    ]
    helicity.add_normalized_extent_metrics(rows)
    helicity.write_csv(out_csv, rows, helicity.SUMMARY_COLUMNS)

    assert out_csv.exists()
    with out_csv.open("r", newline="", encoding="utf-8") as handle:
        written = list(csv.DictReader(handle))
    central_row = next(row for row in written if row["variant_id"] == "central6_units")
    assert central_row["units_per_chain"] == "6"
    assert float(central_row["helical_coherence_score"]) > 0.95
    assert central_row["axial_extent_A"] == "19.200000"
    assert float(central_row["coherent_helical_turns"]) > 0.0
    assert 0.0 < float(central_row["normalized_axial_extent_vs_full"]) < 1.0

    plot_paths = helicity.write_plots(rows, list(csv.DictReader(feature_summary.open())), plot_dir)
    helicity.write_report(rows, list(csv.DictReader(feature_summary.open())), plot_paths, report)
    assert (plot_dir / "mini_hexaplex_units_vs_helical_coherence_score.png").exists()
    assert (plot_dir / "mini_hexaplex_units_vs_coherent_helical_turns.png").exists()
    assert (plot_dir / "mini_hexaplex_coherent_turns_and_4p5_5p0_response.png").exists()
    assert report.exists()
    assert "coherent_helical_turns = axial_extent_A / mean_pitch_A" in report.read_text(encoding="utf-8")
