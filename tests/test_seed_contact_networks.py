import csv
import importlib.util
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


network = _load_script_module("analyze_seed_contact_networks", "scripts/analyze_seed_contact_networks.py")


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


def six_chain_seed(unit_count: int = 1) -> list[PDBAtom]:
    atoms: list[PDBAtom] = []
    serial = 1
    for chain_index, chain_id in enumerate("ABCDEF"):
        base_x = chain_index * 10.0
        for unit_index in range(unit_count):
            residue_number = unit_index * 2 + 1
            z = unit_index * 3.0
            atoms.append(atom(serial, "CA", "CYP" if chain_index % 2 == 0 else "MEP", residue_number, chain_id, base_x, 0.0, z, "C"))
            serial += 1
            atoms.append(atom(serial, "N", "GLU", residue_number + 1, chain_id, base_x, 1.0, z, "N"))
            serial += 1
    return atoms


def test_heavy_atom_contact_detection_on_small_synthetic_structure():
    atoms = [
        atom(1, "CA", "CYP", 1, "A", 0.0, 0.0, 0.0, "C"),
        atom(2, "H", "CYP", 1, "A", 0.1, 0.0, 0.0, "H"),
        atom(3, "CA", "MEP", 1, "B", 3.0, 0.0, 0.0, "C"),
        atom(4, "N", "GLU", 2, "B", 8.0, 0.0, 0.0, "N"),
    ]
    unit_lookup = {
        ("A", "CYP", 1, ""): 1,
        ("B", "MEP", 1, ""): 1,
        ("B", "GLU", 2, ""): 1,
    }

    contacts = network.find_contacts(atoms, unit_lookup, 4.5)

    assert len(contacts) == 1
    assert contacts[0].chain_pair == ("A", "B")
    assert contacts[0].category == "CYP_MEP_vs_CYP_MEP"


def test_chain_graph_connectivity_calculation():
    nodes = set("ABCDEF")
    connected_edges = {("A", "B"), ("B", "C"), ("C", "D"), ("D", "E"), ("E", "F")}
    disconnected_edges = {("A", "B"), ("C", "D"), ("E", "F")}

    assert len(network.connected_components(nodes, connected_edges)) == 1
    assert len(network.connected_components(nodes, disconnected_edges)) == 3


def test_unit_graph_largest_component_calculation(tmp_path):
    atoms = six_chain_seed(1)
    path = tmp_path / "mini_hexaplex_central1_units.pdb"
    write_pdb_atoms(atoms, path)
    variant = network.load_variant(path, 1, "central1_units", "central")
    contacts = [
        network.Contact(0, 2, atoms[0], atoms[2], ("A", "B"), (("A", 1), ("B", 1)), "CYP_MEP_vs_CYP_MEP", True),
        network.Contact(2, 4, atoms[2], atoms[4], ("B", "C"), (("B", 1), ("C", 1)), "CYP_MEP_vs_CYP_MEP", True),
    ]

    summary, _ = network.summarize_network(variant, contacts, 4.5, 2.0, 1.0, None)

    assert summary["unit_graph_node_count"] == "6"
    assert float(summary["unit_graph_largest_component_fraction"]) == pytest.approx(0.5)


def test_contact_categories_are_assigned():
    cyp = atom(1, "CA", "CYP", 1, "A", 0, 0, 0, "C")
    mep = atom(2, "CA", "MEP", 1, "B", 0, 0, 0, "C")
    glu = atom(3, "OE1", "GLU", 2, "C", 0, 0, 0, "O")

    assert network.contact_category(cyp, mep) == "CYP_MEP_vs_CYP_MEP"
    assert network.contact_category(cyp, glu) == "CYP_MEP_vs_GLU"
    assert network.contact_category(glu, glu) == "GLU_vs_GLU"
    assert network.is_backbone_like_contact(cyp, mep)
    assert not network.is_backbone_like_contact(glu, glu)


def test_output_csvs_are_written_with_expected_columns(tmp_path):
    structures_dir = tmp_path / "structures"
    structures_dir.mkdir()
    seed_atoms = six_chain_seed(1)
    seed_atoms[2] = atom(3, "CA", "MEP", 1, "B", 3.0, 0.0, 0.0, "C")
    write_pdb_atoms(seed_atoms, structures_dir / "mini_hexaplex_central1_units.pdb")
    args = type(
        "Args",
        (),
        {
            "structures_dir": structures_dir,
            "ensemble_dir": tmp_path / "ensembles",
            "unit_counts": "1",
            "cutoffs": "4.5",
            "include_lower_end": False,
            "skip_perturbation_robustness": True,
            "summary_csv": tmp_path / "summary.csv",
            "edges_csv": tmp_path / "edges.csv",
            "plot_dir": tmp_path / "plots",
            "out_report": tmp_path / "report.md",
        },
    )()

    summary_rows, edge_rows = network.run(args)

    assert len(summary_rows) == 1
    assert edge_rows
    with args.summary_csv.open("r", newline="", encoding="utf-8") as handle:
        assert csv.DictReader(handle).fieldnames == network.SUMMARY_COLUMNS
    with args.edges_csv.open("r", newline="", encoding="utf-8") as handle:
        assert csv.DictReader(handle).fieldnames == network.EDGE_COLUMNS
    assert args.out_report.exists()


def test_invalid_pdb_path_fails_clearly(tmp_path):
    with pytest.raises(FileNotFoundError, match="Mini-hexaplex PDB not found"):
        network.load_variant(tmp_path / "missing.pdb", 4, "missing", "central")

    invalid = tmp_path / "invalid.pdb"
    invalid.write_text("END\n", encoding="utf-8")
    with pytest.raises(ValueError, match="No ATOM/HETATM records"):
        network.load_variant(invalid, 4, "invalid", "central")
