from scripts.inventory_asem_sidechain_candidates import inventory


def test_inventory_orders_angles_and_candidates_and_counts_structures(tmp_path):
    raw = tmp_path / "raw"
    cand10 = raw / "30" / "cand10"
    cand2 = raw / "30" / "cand2"
    cand1 = raw / "29" / "cand1"
    for path in (cand10, cand2, cand1):
        path.mkdir(parents=True)

    (cand1 / "initial.xyz").write_text(
        "2\nfixture\nC 0 0 0\nN 1 0 0\n",
        encoding="utf-8",
    )
    (cand2 / "initial.pdb").write_text(
        "ATOM      1  C   FIX A   1       0.000   0.000   0.000\n"
        "HETATM    2  O   FIX A   1       1.000   0.000   0.000\n",
        encoding="utf-8",
    )
    (cand10 / "model.mol2").write_text("@<TRIPOS>MOLECULE\n", encoding="utf-8")

    rows = inventory(raw)

    assert [(row["angle_deg"], row["candidate_id"]) for row in rows] == [
        ("29", "cand1"),
        ("30", "cand2"),
        ("30", "cand10"),
    ]
    assert rows[0]["structure_file_type"] == "xyz"
    assert rows[0]["atom_count"] == "2"
    assert rows[1]["structure_file_type"] == "pdb"
    assert rows[1]["atom_count"] == "2"
    assert rows[1]["has_pdb"] == "yes"
    assert rows[2]["has_mol2"] == "yes"
    assert rows[2]["atom_count"] == ""


def test_inventory_marks_candidate_without_recognized_structure(tmp_path):
    raw = tmp_path / "raw"
    candidate = raw / "27" / "cand0"
    candidate.mkdir(parents=True)
    (candidate / "notes.txt").write_text("fixture", encoding="utf-8")

    rows = inventory(raw)

    assert rows[0]["structure_file"] == ""
    assert "no recognized structure file" in rows[0]["notes"]
    assert "notes.txt" in rows[0]["notes"]
