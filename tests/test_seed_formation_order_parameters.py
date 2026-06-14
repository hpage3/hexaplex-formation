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


seed = _load_script_module("analyze_seed_formation_order_parameters", "scripts/analyze_seed_formation_order_parameters.py")


def atom(
    serial: int,
    name: str,
    residue_name: str,
    residue_number: int,
    chain_id: str,
    x: float,
    y: float,
    z: float,
    element: str | None = None,
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
        element=element or name[0],
    )


def synthetic_seed(unit_count: int = 4, radius: float = 6.0) -> list[PDBAtom]:
    atoms: list[PDBAtom] = []
    serial = 1
    for chain_index, chain_id in enumerate("ABCDEF"):
        phase = chain_index * (2.0 * math.pi / 6.0)
        base_name = "CYP" if chain_index % 2 == 0 else "MEP"
        for unit_index in range(unit_count):
            theta = phase + unit_index * math.radians(24.0)
            z = unit_index * 3.2
            x = radius * math.cos(theta)
            y = radius * math.sin(theta)
            residue_number = unit_index * 2 + 1
            atoms.append(atom(serial, "CA", base_name, residue_number, chain_id, x, y, z, "C"))
            serial += 1
            atoms.append(atom(serial, "C", base_name, residue_number, chain_id, x + 0.3, y, z, "C"))
            serial += 1
            atoms.append(atom(serial, "N", "GLU", residue_number + 1, chain_id, x, y + 0.3, z + 0.2, "N"))
            serial += 1
            atoms.append(atom(serial, "O", "GLU", residue_number + 1, chain_id, x - 0.2, y, z + 0.2, "O"))
            serial += 1
    return atoms


def small_contact_atoms() -> list[PDBAtom]:
    return [
        atom(1, "CA", "CYP", 1, "A", 0.0, 0.0, 0.0, "C"),
        atom(2, "N", "GLU", 2, "A", 8.0, 0.0, 0.0, "N"),
        atom(3, "CA", "MEP", 1, "B", 3.0, 0.0, 0.0, "C"),
        atom(4, "O", "GLU", 2, "B", 20.0, 0.0, 0.0, "O"),
    ]


def test_six_chains_and_units_are_preserved(tmp_path):
    path = tmp_path / "mini_hexaplex_central4_units.pdb"
    write_pdb_atoms(synthetic_seed(4), path)

    reference, warnings = seed.load_seed_reference(path, 4, 4.5)

    assert list(reference.chain_atoms) == list("ABCDEF")
    assert all(len(seed.residue_keys_for_chain(atoms)) == 8 for atoms in reference.chain_atoms.values())
    assert warnings == []


def test_rigid_body_transform_preserves_intrachain_distances():
    chain_atoms = tuple(synthetic_seed(2)[:8])
    before = seed.distance(chain_atoms[0], chain_atoms[-1])
    rotation = seed.rotation_matrix_from_axis_angle((1.0, 2.0, 3.0), math.radians(47.0))

    transformed = seed.transform_chain(chain_atoms, rotation, (4.0, -2.0, 1.5))

    assert seed.distance(transformed[0], transformed[-1]) == pytest.approx(before)


def test_contact_fraction_works_on_synthetic_example():
    target_contacts, cyp_mep_contacts, backbone_contacts = seed.find_interchain_contacts(small_contact_atoms(), 4.5)
    sample_atoms = small_contact_atoms()
    moved_apart = [
        sample_atoms[0],
        sample_atoms[1],
        atom(3, "CA", "MEP", 1, "B", 10.0, 0.0, 0.0, "C"),
        sample_atoms[3],
    ]
    sample_contacts, sample_cyp_mep, sample_backbone = seed.find_interchain_contacts(moved_apart, 4.5)

    fraction, present = seed.contact_fraction(frozenset(target_contacts), sample_contacts)
    cyp_fraction, _ = seed.contact_fraction(frozenset(cyp_mep_contacts), sample_cyp_mep)
    backbone_fraction, _ = seed.contact_fraction(frozenset(backbone_contacts), sample_backbone)

    assert len(target_contacts) == 1
    assert present == 0
    assert fraction == 0.0
    assert cyp_fraction == 0.0
    assert backbone_fraction == 0.0


def test_rmsd_to_identical_structure_is_near_zero(tmp_path):
    path = tmp_path / "mini_hexaplex_central4_units.pdb"
    atoms = synthetic_seed(4)
    write_pdb_atoms(atoms, path)
    reference, _ = seed.load_seed_reference(path, 4, 4.5)

    assert seed.rmsd_to_reference(list(reference.atoms), reference.atoms) == pytest.approx(0.0, abs=1e-8)


def test_loose_ensemble_has_larger_radius_of_gyration_than_formed_reference(tmp_path):
    path = tmp_path / "mini_hexaplex_central4_units.pdb"
    write_pdb_atoms(synthetic_seed(4), path)
    reference, _ = seed.load_seed_reference(path, 4, 4.5)
    rng = seed.random.Random(12)

    loose = seed.generate_loose_initial(reference, 0, rng)

    assert seed.radius_of_gyration(loose) > reference.radius_of_gyration


def test_angular_randomized_loose_preserves_intrachain_distances_and_reduces_phase_score(tmp_path):
    path = tmp_path / "mini_hexaplex_central4_units.pdb"
    write_pdb_atoms(synthetic_seed(4), path)
    reference, _ = seed.load_seed_reference(path, 4, 4.5)
    rng = seed.random.Random(14)
    chain_id = "A"
    before_atoms = reference.chain_atoms[chain_id]
    before_distance = seed.distance(before_atoms[0], before_atoms[-1])

    randomized = seed.generate_angular_randomized_loose_initial(reference, 0, rng)
    randomized_chain = tuple(atom for atom in randomized if atom.chain_id == chain_id)
    _, _, centroids = seed.chain_centroid_metrics(seed.group_by_chain(randomized))
    phase_score = seed.angular_phase_order_score(reference, centroids)

    assert seed.distance(randomized_chain[0], randomized_chain[-1]) == pytest.approx(before_distance)
    assert phase_score is not None
    assert phase_score < 0.9


def test_output_csv_is_written_with_expected_columns(tmp_path):
    structures_dir = tmp_path / "structures"
    structures_dir.mkdir()
    write_pdb_atoms(synthetic_seed(4), structures_dir / "mini_hexaplex_central4_units.pdb")
    out_csv = tmp_path / "metrics" / "seed.csv"
    args = type(
        "Args",
        (),
        {
            "structures_dir": structures_dir,
            "unit_counts": "4",
            "samples_per_ensemble": 1,
            "loose_mode": "angular_randomized",
            "random_seed": 1,
            "contact_cutoff": 4.5,
            "ensemble_dir": tmp_path / "ensembles",
            "plot_dir": tmp_path / "plots",
            "out_csv": out_csv,
            "endpoint_means_csv": tmp_path / "means.csv",
            "out_report": tmp_path / "report.md",
            "save_examples": 0,
        },
    )()

    rows = seed.run(args)

    assert len(rows) == 2
    with out_csv.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        assert reader.fieldnames == seed.SUMMARY_COLUMNS
        written = list(reader)
    assert {row["ensemble_type"] for row in written} == {"formed_perturbed", "angular_randomized_loose_initial"}
    assert (tmp_path / "report.md").exists()
    assert (tmp_path / "means.csv").exists()


def test_invalid_or_missing_pdb_path_fails_clearly(tmp_path):
    missing = tmp_path / "missing.pdb"

    with pytest.raises(FileNotFoundError, match="Mini-hexaplex PDB not found"):
        seed.load_seed_reference(missing, 4, 4.5)

    invalid = tmp_path / "invalid.pdb"
    invalid.write_text("END\n", encoding="utf-8")
    with pytest.raises(ValueError, match="No ATOM/HETATM records"):
        seed.load_seed_reference(invalid, 4, 4.5)
