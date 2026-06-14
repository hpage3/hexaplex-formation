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
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


mini = _load_script_module("generate_mini_hexaplex_variants", "scripts/generate_mini_hexaplex_variants.py")


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


def six_chain_fixture(unit_count: int = 10) -> list[PDBAtom]:
    atoms: list[PDBAtom] = []
    serial = 1
    for chain_index, chain_id in enumerate("ABCDEF"):
        base_name = "CYP" if chain_index % 2 == 0 else "MEP"
        x = float(chain_index * 4)
        increasing = chain_index % 2 == 0
        for unit_index in range(1, unit_count + 1):
            axial = float(unit_index if increasing else unit_count + 1 - unit_index)
            residue_number = (unit_index - 1) * 2 + 1
            atoms.append(atom(serial, "C", base_name, residue_number, chain_id, x, 0.0, axial))
            serial += 1
            atoms.append(atom(serial, "N", "GLU", residue_number + 1, chain_id, x, 1.0, axial + 0.05))
            serial += 1
    return atoms


def _units_by_chain(atoms: list[PDBAtom]):
    origin = (10.0, 0.5, 5.5)
    axis = (0.0, 0.0, 1.0)
    residues = mini.build_residues(atoms, origin, axis)
    units_by_chain, warnings = mini.group_repeat_units(residues)
    return residues, units_by_chain, warnings


def test_six_chains_and_base_glu_units_are_confirmed():
    atoms = six_chain_fixture()
    residues, units_by_chain, warnings = _units_by_chain(atoms)

    assert mini.confirm_six_chains(atoms) == list("ABCDEF")
    assert warnings == []
    assert set(units_by_chain) == set("ABCDEF")
    assert all(len(units) == 10 for units in units_by_chain.values())
    assert all(len(residues[chain_id]) == 20 for chain_id in "ABCDEF")
    for units in units_by_chain.values():
        for unit in units:
            assert unit.base_residue.residue_name in {"CYP", "MEP"}
            assert unit.glu_residue.residue_name == "GLU"


def test_literal_first8_uses_residue_order():
    atoms = six_chain_fixture()
    _residues, units_by_chain, _warnings = _units_by_chain(atoms)

    selection = mini.select_literal_first(units_by_chain, 8)

    assert {chain_id: [unit.unit_index for unit in units] for chain_id, units in selection.items()} == {
        chain_id: list(range(1, 9)) for chain_id in "ABCDEF"
    }


def test_same_physical_lower_end_uses_consistent_axial_end():
    atoms = six_chain_fixture()
    _residues, units_by_chain, _warnings = _units_by_chain(atoms)

    selection = mini.select_physical_end(units_by_chain, 8, "lower")

    assert [unit.unit_index for unit in selection["A"]] == list(range(1, 9))
    assert [unit.unit_index for unit in selection["B"]] == list(range(3, 11))
    for units in selection.values():
        assert len(units) == 8
        assert max(unit.axial_t for unit in units) < 9.2


def test_central8_uses_axial_coordinates():
    atoms = six_chain_fixture()
    _residues, units_by_chain, _warnings = _units_by_chain(atoms)

    selection = mini.select_central(units_by_chain, 8)

    assert all([unit.unit_index for unit in units] == list(range(2, 10)) for units in selection.values())


def test_12_unit_variants_preserve_six_chains_and_use_axial_coordinates():
    atoms = six_chain_fixture(unit_count=15)
    _residues, units_by_chain, _warnings = _units_by_chain(atoms)

    literal = mini.select_literal_first(units_by_chain, 12)
    lower = mini.select_physical_end(units_by_chain, 12, "lower")
    central = mini.select_central(units_by_chain, 12)

    assert set(literal) == set("ABCDEF")
    assert all(len(units) == 12 for units in literal.values())
    assert all(len(units) == 12 for units in lower.values())
    assert all(len(units) == 12 for units in central.values())
    assert [unit.unit_index for unit in literal["A"]] == list(range(1, 13))
    assert [unit.unit_index for unit in lower["A"]] == list(range(1, 13))
    assert [unit.unit_index for unit in lower["B"]] == list(range(4, 16))
    assert [unit.unit_index for unit in central["A"]] == list(range(2, 14))
    for selection in (literal, lower, central):
        for units in selection.values():
            for unit in units:
                assert unit.base_residue.residue_name in {"CYP", "MEP"}
                assert unit.glu_residue.residue_name == "GLU"


def test_short_unit_counts_4_to_7_are_accepted_and_preserve_geometry():
    atoms = six_chain_fixture(unit_count=15)
    _residues, units_by_chain, _warnings = _units_by_chain(atoms)

    for count in (4, 5, 6, 7):
        lower = mini.select_physical_end(units_by_chain, count, "lower")
        central = mini.select_central(units_by_chain, count)

        assert all(len(units) == count for units in lower.values())
        assert all(len(units) == count for units in central.values())
        assert len({unit.unit_index for unit in lower["A"]}) == count
        assert len({unit.unit_index for unit in central["A"]}) == count
        for selection in (lower, central):
            for units in selection.values():
                for unit in units:
                    assert unit.base_residue.residue_name in {"CYP", "MEP"}
                    assert unit.glu_residue.residue_name == "GLU"


def test_variant_outputs_preserve_chains_coordinates_and_write_csvs(tmp_path):
    atoms = six_chain_fixture(unit_count=15)
    input_pdb = tmp_path / "fixture.pdb"
    out_dir = tmp_path / "structures"
    manifest = tmp_path / "manifest.csv"
    geometry = tmp_path / "geometry.csv"
    write_pdb_atoms(atoms, input_pdb)
    loaded_atoms = load_pdb_atoms(input_pdb)
    origin = (10.0, 0.5, 5.5)
    axis = (0.0, 0.0, 1.0)
    residues = mini.build_residues(loaded_atoms, origin, axis)
    units_by_chain, _warnings = mini.group_repeat_units(residues)

    manifest_rows, geometry_rows, _selections = mini.build_variant_rows(
        input_pdb,
        out_dir,
        loaded_atoms,
        units_by_chain,
        origin,
        axis,
        8,
        0.01,
    )
    mini.write_csv(manifest_rows, manifest, mini.MANIFEST_FIELDNAMES)
    mini.write_csv(geometry_rows, geometry, mini.GEOMETRY_FIELDNAMES)

    assert manifest.exists()
    assert geometry.exists()
    assert (out_dir / "mini_hexaplex_literal_first8_units.pdb").exists()
    with manifest.open("r", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert {row["variant_id"] for row in rows} == {
        "literal_first8_units",
        "lower_end_first8_units",
        "central8_units",
    }
    assert all(row["chains_included"] == "A,B,C,D,E,F" for row in rows)
    assert all(row["units_per_chain"] == "8" for row in rows)
    assert all(row["total_residue_count"] == "96" for row in rows)

    selected_atoms = load_pdb_atoms(out_dir / "mini_hexaplex_literal_first8_units.pdb")
    source_by_key = {
        (a.chain_id, a.residue_name, a.residue_number, a.atom_name): (a.x, a.y, a.z) for a in loaded_atoms
    }
    assert len({atom.chain_id for atom in selected_atoms}) == 6
    assert mini.duplicate_atom_count(selected_atoms) == 0
    for selected_atom in selected_atoms:
        key = (
            selected_atom.chain_id,
            selected_atom.residue_name,
            selected_atom.residue_number,
            selected_atom.atom_name,
        )
        assert (selected_atom.x, selected_atom.y, selected_atom.z) == pytest.approx(source_by_key[key])


def test_12_unit_variant_outputs_and_invalid_16_unit_request(tmp_path, monkeypatch):
    atoms = six_chain_fixture(unit_count=15)
    input_pdb = tmp_path / "fixture.pdb"
    out_dir = tmp_path / "structures"
    manifest = tmp_path / "manifest.csv"
    geometry = tmp_path / "geometry.csv"
    report = tmp_path / "diagnostics.md"
    write_pdb_atoms(atoms, input_pdb)

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "generate_mini_hexaplex_variants.py",
            "--pdb",
            str(input_pdb),
            "--units",
            "12",
            "--out-dir",
            str(out_dir),
            "--manifest",
            str(manifest),
            "--geometry-out",
            str(geometry),
            "--diagnostic-report",
            str(report),
        ],
    )
    assert mini.main() == 0
    assert (out_dir / "mini_hexaplex_literal_first12_units.pdb").exists()
    assert (out_dir / "mini_hexaplex_lower_end_first12_units.pdb").exists()
    assert (out_dir / "mini_hexaplex_central12_units.pdb").exists()

    with manifest.open("r", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert {row["variant_id"] for row in rows} == {
        "literal_first12_units",
        "lower_end_first12_units",
        "central12_units",
    }
    assert all(row["units_per_chain"] == "12" for row in rows)
    assert all(row["total_residue_count"] == "144" for row in rows)

    selected_atoms = load_pdb_atoms(out_dir / "mini_hexaplex_central12_units.pdb")
    assert len({atom.chain_id for atom in selected_atoms}) == 6
    assert mini.duplicate_atom_count(selected_atoms) == 0

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "generate_mini_hexaplex_variants.py",
            "--pdb",
            str(input_pdb),
            "--units",
            "16",
            "--out-dir",
            str(out_dir),
        ],
    )
    with pytest.raises(SystemExit, match="only has 15 units; cannot select 16"):
        mini.main()
