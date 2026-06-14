from dataclasses import replace
from pathlib import Path

import pytest

from hexaplex_formation.pdb_utils import (
    atom_count,
    atom_identity_key,
    bounding_box,
    centroid,
    chain_ids,
    dedupe_exact_atoms,
    heavy_atoms,
    is_hydrogen,
    load_pdb_atoms,
    residue_count,
    residue_names,
    write_pdb_atoms,
)


PDB_FIXTURE = """\
ATOM      1  N   GLU A   1       0.000   1.000   2.000  1.00 10.00           N
ATOM      2  CA  GLU A   1       1.000   2.000   3.000  1.00 11.00
ATOM      3  H1  GLU A   1       2.000   3.000   4.000  1.00 12.00
ATOM      4  C   ALA B   2      -1.000   0.000   5.000  0.50 13.00           C
HETATM    5  O   HOH B   3       3.000  -1.000   1.000  1.00 20.00           O
END
"""


@pytest.fixture
def atoms(tmp_path: Path):
    pdb_path = tmp_path / "tiny.pdb"
    pdb_path.write_text(PDB_FIXTURE, encoding="utf-8")
    return load_pdb_atoms(pdb_path)


def test_counts_and_identifiers(atoms):
    assert atom_count(atoms) == 5
    assert residue_count(atoms) == 3
    assert chain_ids(atoms) == ["A", "B"]
    assert residue_names(atoms) == ["ALA", "GLU", "HOH"]


def test_centroid_and_bounding_box(atoms):
    assert centroid(atoms) == pytest.approx((1.0, 1.0, 3.0))
    assert bounding_box(atoms) == pytest.approx((-1.0, 3.0, -1.0, 3.0, 1.0, 5.0))


def test_element_fallback_behavior(atoms):
    assert atoms[1].atom_name == "CA"
    assert atoms[1].element == "C"
    assert atoms[2].atom_name == "H1"
    assert atoms[2].element == "H"


def test_hydrogen_detection(atoms):
    assert is_hydrogen(atoms[2])
    assert is_hydrogen(replace(atoms[1], atom_name=" H2 ", element=""))
    assert not is_hydrogen(atoms[1])


def test_heavy_atom_filtering(atoms):
    filtered = heavy_atoms(atoms)
    assert atom_count(filtered) == 4
    assert [atom.atom_name for atom in filtered] == ["N", "CA", "C", "O"]


def test_all_atom_dedupe_preserves_hydrogens_and_first_atom_order(atoms):
    duplicate_heavy = replace(atoms[1], atom_serial=99, occupancy=0.25, temp_factor=99.0)
    duplicate_hydrogen = replace(atoms[2], atom_serial=98, occupancy=0.25, temp_factor=99.0)
    shifted = replace(atoms[1], atom_serial=100, x=1.004)
    deduped = dedupe_exact_atoms([atoms[0], atoms[1], duplicate_heavy, atoms[2], duplicate_hydrogen, atoms[3], shifted])

    assert [atom.atom_serial for atom in deduped] == [1, 2, 3, 4, 100]
    assert atom_identity_key(atoms[1]) == atom_identity_key(duplicate_heavy)
    assert atom_identity_key(atoms[2]) == atom_identity_key(duplicate_hydrogen)
    assert atom_identity_key(atoms[1]) != atom_identity_key(shifted)


def test_heavy_only_dedupe_removes_hydrogens_before_deduping(atoms):
    duplicate_heavy = replace(atoms[1], atom_serial=99, occupancy=0.25, temp_factor=99.0)
    duplicate_hydrogen = replace(atoms[2], atom_serial=98, occupancy=0.25, temp_factor=99.0)
    heavy_deduped = dedupe_exact_atoms(heavy_atoms([atoms[0], atoms[1], duplicate_heavy, atoms[2], duplicate_hydrogen]))

    assert [atom.atom_serial for atom in heavy_deduped] == [1, 2]


def test_write_pdb_atoms_round_trip(tmp_path: Path, atoms):
    filtered = heavy_atoms(atoms)
    out_path = tmp_path / "written.pdb"
    write_pdb_atoms(filtered, out_path)

    reloaded = load_pdb_atoms(out_path)
    assert atom_count(reloaded) == 4
    assert residue_count(reloaded) == 3
    assert chain_ids(reloaded) == ["A", "B"]
    assert [atom.atom_serial for atom in reloaded] == [1, 2, 3, 4]
    assert centroid(reloaded) == pytest.approx((0.75, 0.5, 2.75))
