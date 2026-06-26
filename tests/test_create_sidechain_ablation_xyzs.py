from pathlib import Path

import pytest

from scripts.create_sidechain_ablation_xyzs import (
    REMOVED_GLU_ATOMS,
    parse_pdb_for_ablation,
    process_manifest,
)


def pdb_line(serial: int, atom: str, residue: str, x: float, element: str) -> str:
    return (
        f"ATOM  {serial:5d} {atom:^4s} {residue:>3s} A{1:4d}    "
        f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00          {element:>2s}"
    )


def test_ablation_preserves_core_and_removes_only_distal_glu_atoms(tmp_path: Path):
    source = tmp_path / "model.pdb"
    names = ["N", "CA", "CB", "CG", "CD", "OE1", "OE2", "C", "O"]
    source.write_text(
        "\n".join(
            pdb_line(index, name, "GLU", float(index), "O" if name.startswith("O") else "C")
            for index, name in enumerate(names, start=1)
        )
        + "\n"
        + pdb_line(20, "C1", "CYP", 20.0, "C")
        + "\n",
        encoding="utf-8",
    )

    atoms, counts = parse_pdb_for_ablation(source)
    retained_x = {atom[1] for atom in atoms}

    assert {4.0, 5.0, 6.0, 7.0}.isdisjoint(retained_x)
    assert {1.0, 2.0, 3.0, 8.0, 9.0, 20.0}.issubset(retained_x)
    assert counts["removed_atom_count"] == 4
    assert counts["removed_heavy_atom_count"] == 4
    assert REMOVED_GLU_ATOMS == {"CG", "HG2", "HG3", "CD", "OE1", "OE2"}


def test_ablation_refuses_structure_without_removable_atoms(tmp_path: Path):
    source = tmp_path / "model.pdb"
    source.write_text(pdb_line(1, "CA", "GLU", 1.0, "C") + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="No removable GLU side-chain atoms"):
        parse_pdb_for_ablation(source)


def test_manifest_contains_counts_rule_and_output(tmp_path: Path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    source = Path("raw/model.pdb")
    source.parent.mkdir()
    source.write_text(
        pdb_line(1, "CA", "GLU", 1.0, "C")
        + "\n"
        + pdb_line(2, "CG", "GLU", 2.0, "C")
        + "\n",
        encoding="utf-8",
    )
    manifest = Path("ranked.csv")
    manifest.write_text(
        "model_id,twist_deg,rise_A,source_pdb\n"
        "angle30_cand0,30,3.4,raw/model.pdb\n",
        encoding="utf-8",
    )

    rows = process_manifest(manifest, Path("xyz"), Path("output_manifest.csv"))

    assert rows[0]["status"] == "complete"
    assert rows[0]["source_atom_count"] == 2
    assert rows[0]["removed_atom_count"] == 1
    assert rows[0]["retained_atom_count"] == 1
    assert "CG" in str(rows[0]["stripping_rule"])
    assert Path("xyz/angle30_cand0.xyz").is_file()
    assert "removed_heavy_atom_count" in Path("output_manifest.csv").read_text(encoding="utf-8")
