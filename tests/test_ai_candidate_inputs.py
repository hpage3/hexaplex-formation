import importlib.util
import sys
from pathlib import Path

import pytest


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


build_ai = _load_script_module("build_ai_candidate_inputs", "scripts/build_ai_candidate_inputs.py")


def test_chain_assignment_helpers():
    assert build_ai.infer_chain_index(0, 180, 6) == 0
    assert build_ai.infer_chain_index(29, 180, 6) == 0
    assert build_ai.infer_chain_index(30, 180, 6) == 1
    assert build_ai.assign_chain_id(0, 180, 6) == "A"
    assert build_ai.assign_chain_id(150, 180, 6) == "F"


def test_proxy_residue_name_uses_x_for_nonstandard_residues():
    assert build_ai.proxy_residue_name("GLU") == "E"
    assert build_ai.proxy_residue_name("CYP") == "X"
    assert build_ai.proxy_residue_name("MEP") == "X"


def test_build_chain_rows_preserves_residue_patterns():
    residue_items = [
        (("", "CYP", 1, ""), []),
        (("", "GLU", 2, ""), []),
        (("", "MEP", 31, ""), []),
        (("", "GLU", 32, ""), []),
    ]

    residue_rows = build_ai.build_residue_rows(residue_items, 2)
    chain_rows = build_ai.build_chain_rows(residue_rows)

    assert chain_rows[0]["residue_pattern"] == "CYP-GLU"
    assert chain_rows[0]["proxy_sequence"] == "XE"
    assert chain_rows[1]["residue_pattern"] == "MEP-GLU"
    assert chain_rows[1]["proxy_sequence"] == "XE"


def test_full_pdb_builds_six_inferred_chains(tmp_path):
    out_dir = tmp_path / "ai"
    report = tmp_path / "report.md"

    argv = [
        "build_ai_candidate_inputs.py",
        "--pdb",
        str(REPO_ROOT / "inputs/structures/full_hexaplex_anti_parallel_30deg_ideal.pdb"),
        "--out-dir",
        str(out_dir),
        "--report",
        str(report),
    ]
    old_argv = sys.argv
    sys.argv = argv
    try:
        assert build_ai.main() == 0
    finally:
        sys.argv = old_argv

    chain_csv = out_dir / "full_hexaplex_anti_parallel_30deg_ideal_chain_pattern_summary.csv"
    fasta = out_dir / "full_hexaplex_anti_parallel_30deg_ideal_alphafold_esm_proxy.fasta"
    deduped_pdb = out_dir / "full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain.pdb"

    assert chain_csv.exists()
    assert fasta.exists()
    assert deduped_pdb.exists()

    text = chain_csv.read_text(encoding="utf-8")
    assert "chain_id,chain_index,residue_count,nonstandard_residue_count" in text
    assert "A,1,30,15" in text
    assert "F,6,30,15" in text

    fasta_text = fasta.read_text(encoding="utf-8")
    assert "# WARNING: CYP and MEP are nonstandard residues" in fasta_text
    assert "proxy" in fasta_text
    assert "XE" in fasta_text

    report_text = report.read_text(encoding="utf-8")
    assert "AlphaFold/ESM Candidate Inputs" in report_text
    assert "CYP and MEP are nonstandard residues" in report_text
    assert "proxy FASTA uses `X`" in report_text


def test_report_mentions_losssy_proxy_usage(tmp_path):
    residue_rows = [
        {
            "residue_index_in_pdb_order": "1",
            "chain_id": "A",
            "residue_name": "CYP",
            "residue_number": "1",
            "insertion_code": "",
            "residue_label": "A:CYP1",
            "proxy_residue": "X",
            "map_warning": "nonstandard_residue",
        }
    ]
    chain_rows = [
        {
            "chain_id": "A",
            "chain_index": "1",
            "residue_count": "1",
            "nonstandard_residue_count": "1",
            "residue_pattern": "CYP",
            "proxy_sequence": "X",
            "first_residue_label": "A:CYP1",
            "last_residue_label": "A:CYP1",
            "contains_nonstandard_residues": "yes",
        }
    ]
    out = tmp_path / "report.md"
    build_ai.write_report(Path("source.pdb"), Path("deduped.pdb"), residue_rows, chain_rows, out)
    text = out.read_text(encoding="utf-8")
    assert "candidate-only exploratory use" in text
    assert "does not make the structure a biologically standard AlphaFold/ESM sequence" in text
