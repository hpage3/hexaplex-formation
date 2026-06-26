from pathlib import Path

from scripts.extract_asem_candidate_overlay_pdf_metadata import parse_overlay_text


def test_parse_overlay_text_extracts_candidate_metrics():
    text = """
    29/cand0/initial.pdb
    r=0.341  Rwp=0.931  E=-134.8 kcal/mol
    """

    rows = parse_overlay_text(text, pdf_page=7, source_pdf=Path("overlay.pdf"))

    assert rows == [
        {
            "angle_deg": "29",
            "candidate_id": "cand0",
            "candidate_path": "29/cand0/initial.pdb",
            "pdf_page": "7",
            "r": "0.341",
            "rwp": "0.931",
            "energy_kcal_mol": "-134.8",
            "source_pdf": "overlay.pdf",
            "notes": "Asem PDF visual overlay; raw numeric profile data not present in this PDF",
        }
    ]


def test_parse_overlay_text_allows_nan_energy_and_multiple_entries():
    text = """
    18/cand0/initial.pdb
    r=0.083  Rwp=1.174  E=nan kcal/mol
    30/cand10/initial.pdb
    r=-0.205 Rwp=1.312 E=-99 kcal/mol
    """

    rows = parse_overlay_text(text, pdf_page=2, source_pdf=Path("pNAB.pdf"))

    assert [row["candidate_path"] for row in rows] == [
        "18/cand0/initial.pdb",
        "30/cand10/initial.pdb",
    ]
    assert rows[0]["energy_kcal_mol"] == "nan"
    assert rows[1]["r"] == "-0.205"
