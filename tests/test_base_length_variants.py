import csv
import importlib.util
import sys
from pathlib import Path

import pytest

from hexaplex_formation.pdb_utils import PDBAtom, load_pdb_atoms, write_pdb_atoms


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


variants = _load_script_module("generate_base_length_variants", "scripts/generate_base_length_variants.py")


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


def fixture_atoms() -> list[PDBAtom]:
    return [
        atom(1, "N", "CYP", 1, "A", 0.0, 0.0, 0.0),
        atom(2, "CA", "CYP", 1, "A", 1.0, 0.0, 0.0),
        atom(3, "C", "CYP", 1, "A", 2.0, 0.0, 0.0),
        atom(4, "O", "CYP", 1, "A", 3.0, 0.0, 0.0),
        atom(5, "CB", "CYP", 1, "A", 1.0, 2.0, 0.0),
        atom(6, "N", "GLU", 2, "A", 0.0, 5.0, 0.0),
        atom(7, "CA", "GLU", 2, "A", 1.0, 5.0, 0.0),
        atom(8, "CB", "GLU", 2, "A", 1.0, 7.0, 0.0),
        atom(9, "N", "MEP", 3, "B", 0.0, 0.0, 5.0),
        atom(10, "CA", "MEP", 3, "B", 1.0, 0.0, 5.0),
        atom(11, "C", "MEP", 3, "B", 2.0, 0.0, 5.0),
        atom(12, "NX1", "MEP", 3, "B", 1.0, 0.0, 8.0),
    ]


def transformable_names() -> dict[str, set[str]]:
    return {"CYP": {"CB"}, "MEP": {"NX1"}}


def test_variant_generator_preserves_atom_count_and_chain_count():
    atoms = fixture_atoms()
    variant_atoms, transformed, warnings = variants.generate_variant_atoms(atoms, transformable_names(), 1.20)

    assert len(variant_atoms) == len(atoms)
    assert len(transformed) == 2
    assert warnings == []
    assert variants.chain_count(variant_atoms) == 2


def test_scale_1p00_leaves_coordinates_unchanged():
    atoms = fixture_atoms()
    variant_atoms, _transformed, _warnings = variants.generate_variant_atoms(atoms, transformable_names(), 1.00)

    assert variants.coordinates_match(atoms, variant_atoms)


def test_fixed_atoms_and_glu_atoms_remain_unchanged_for_all_scales():
    atoms = fixture_atoms()
    variant_atoms, transformed, _warnings = variants.generate_variant_atoms(atoms, transformable_names(), 0.85)

    for index, (old_atom, new_atom) in enumerate(zip(atoms, variant_atoms)):
        if index not in transformed:
            assert variants.distance(old_atom, new_atom) == pytest.approx(0.0)
        if old_atom.residue_name == "GLU":
            assert variants.distance(old_atom, new_atom) == pytest.approx(0.0)


def test_transformed_atoms_follow_anchor_scaling_rule():
    atoms = fixture_atoms()
    variant_atoms, transformed, _warnings = variants.generate_variant_atoms(atoms, transformable_names(), 1.20)

    cyp_cb = variant_atoms[4]
    mep_nx1 = variant_atoms[11]
    assert 4 in transformed
    assert 11 in transformed
    assert (cyp_cb.x, cyp_cb.y, cyp_cb.z) == pytest.approx((1.0, 2.4, 0.0))
    assert (mep_nx1.x, mep_nx1.y, mep_nx1.z) == pytest.approx((1.0, 0.0, 8.6))


def test_manifest_and_geometry_csvs_are_created(tmp_path):
    atoms = fixture_atoms()
    input_pdb = tmp_path / "fixture.pdb"
    out_dir = tmp_path / "structures"
    manifest = out_dir / "base_length_variant_manifest.csv"
    geometry = tmp_path / "geometry.csv"
    report = tmp_path / "report.md"
    write_pdb_atoms(atoms, input_pdb)

    manifest_rows, geometry_rows, _summaries = variants.build_variant_rows(
        input_pdb,
        out_dir,
        [1.00, 1.20],
        load_pdb_atoms(input_pdb),
        transformable_names(),
        0.5,
    )
    variants.write_csv(manifest_rows, manifest, variants.MANIFEST_FIELDNAMES)
    variants.write_csv(geometry_rows, geometry, variants.GEOMETRY_FIELDNAMES)
    variants.write_report(report, input_pdb, manifest_rows, geometry_rows, transformable_names())

    assert manifest.exists()
    assert geometry.exists()
    assert report.exists()
    assert (out_dir / "hexaplex_base_length_scale_1p00.pdb").exists()
    assert (out_dir / "hexaplex_base_length_scale_1p20.pdb").exists()

    with manifest.open("r", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["transformed_atom_count"] == "2"
    assert rows[0]["notes"] == "scale 1.00 baseline coordinate match"

    with geometry.open("r", newline="", encoding="utf-8") as handle:
        geometry_csv_rows = list(csv.DictReader(handle))
    assert geometry_csv_rows[0]["chain_count"] == "2"
    assert geometry_csv_rows[0]["residue_count"] == "3"
    assert "geometry sensitivity study" in report.read_text(encoding="utf-8")
