import csv
import importlib.util
from pathlib import Path

from hexaplex_formation.pdb_utils import PDBAtom
from hexaplex_formation.strand_map import (
    block_for_atom,
    infer_component,
    load_strand_map,
    residue_label_from_atom,
)


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_script_module(name: str, relative_path: str):
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


decompose = _load_script_module("block_contact_decomposition", "scripts/block_contact_decomposition.py")
summarize = _load_script_module("summarize_block_contacts", "scripts/summarize_block_contacts.py")
report = _load_script_module("report_block_contact_analysis", "scripts/report_block_contact_analysis.py")


def atom(
    serial: int,
    atom_name: str,
    residue_name: str,
    residue_number: int,
    element: str,
    chain_id: str = "A",
) -> PDBAtom:
    return PDBAtom(
        record_type="ATOM",
        atom_serial=serial,
        atom_name=atom_name,
        alt_loc="",
        residue_name=residue_name,
        chain_id=chain_id,
        residue_number=residue_number,
        insertion_code="",
        x=float(serial),
        y=0.0,
        z=0.0,
        occupancy=1.0,
        temp_factor=0.0,
        element=element,
    )


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def test_strand_map_loading_and_atom_block_assignment(tmp_path):
    path = tmp_path / "strand_map_candidate.csv"
    _write_csv(
        path,
        ["block_id", "residue_index_in_pdb_order", "chain_id", "residue_name", "residue_number", "insertion_code", "residue_label"],
        [
            {
                "block_id": "2",
                "residue_index_in_pdb_order": "1",
                "chain_id": "A",
                "residue_name": "GLU",
                "residue_number": "10",
                "insertion_code": "",
                "residue_label": "A:GLU10",
            }
        ],
    )
    mapping = load_strand_map(path)
    scaffold_atom = atom(1, "CA", "GLU", 10, "C")
    hexad_atom = atom(3, "N1", "GLU", 10, "N")
    other_atom = atom(2, "CA", "ALA", 11, "C")

    assert block_for_atom(scaffold_atom, mapping) == "2"
    assert block_for_atom(other_atom, mapping) is None
    assert infer_component(scaffold_atom, mapping) == "scaffold"
    assert infer_component(hexad_atom, mapping) == "hexad_or_other"
    assert infer_component(other_atom, mapping) == "hexad_or_other"
    assert residue_label_from_atom(scaffold_atom) == "A:GLU10"


def test_contact_category_classification():
    assert decompose.classify_contact("scaffold", "scaffold", "1", "1") == "scaffold_within_block"
    assert decompose.classify_contact("scaffold", "scaffold", "1", "2") == "scaffold_between_blocks"
    assert decompose.classify_contact("hexad_or_other", "hexad_or_other", None, None) == "hexad_or_other_internal"
    assert decompose.classify_contact("scaffold", "hexad_or_other", "1", None) == "scaffold_hexad_or_other"


def test_glu_motif_boolean_flags():
    flags = decompose.motif_flags(
        atom(1, "O", "ALA", 1, "O"),
        atom(2, "OE1", "GLU", 2, "O"),
    )
    assert flags["is_GLU_involved"] == "yes"
    assert flags["is_GLU_GLU"] == "no"
    assert flags["is_backbone_O_to_GLU_sidechain_O"] == "yes"
    assert flags["is_GLU_sidechain_O_to_GLU_sidechain_O"] == "no"

    glu_glu_flags = decompose.motif_flags(
        atom(3, "OE1", "GLU", 3, "O"),
        atom(4, "OE2", "GLU", 4, "O"),
    )
    assert glu_glu_flags["is_GLU_GLU"] == "yes"
    assert glu_glu_flags["is_GLU_sidechain_O_to_GLU_sidechain_O"] == "yes"


def test_block_pair_summary_aggregation():
    rows = [
        {
            "component_i": "scaffold",
            "component_j": "scaffold",
            "block_i": "2",
            "block_j": "1",
            "is_GLU_involved": "yes",
            "is_GLU_GLU": "yes",
            "min_distance_A": "3.500000",
        },
        {
            "component_i": "scaffold",
            "component_j": "hexad_or_other",
            "block_i": "2",
            "block_j": "",
            "is_GLU_involved": "no",
            "is_GLU_GLU": "no",
            "min_distance_A": "4.000000",
        },
        {
            "component_i": "hexad_or_other",
            "component_j": "hexad_or_other",
            "block_i": "",
            "block_j": "",
            "is_GLU_involved": "no",
            "is_GLU_GLU": "no",
            "min_distance_A": "2.500000",
        },
    ]

    summary = {row["block_pair"]: row for row in summarize.summarize_by_block_pair(rows)}

    assert summary["1--2"]["contact_count"] == "1"
    assert summary["1--2"]["GLU_GLU_count"] == "1"
    assert summary["1--2"]["min_distance_A"] == "3.500000"
    assert summary["scaffold--hexad_or_other"]["contact_count"] == "1"
    assert summary["hexad_or_other--hexad_or_other"]["contact_count"] == "1"


def test_block_contact_report_includes_caution_language(tmp_path):
    out = tmp_path / "block_contact_analysis.md"
    summaries = {
        "hexaplex_scaffold_only_complement_heavy_deduped": {
            "scaffold_within_block": {
                "contact_count": "10",
                "GLU_involved_count": "8",
                "GLU_GLU_count": "4",
            },
            "scaffold_between_blocks": {
                "contact_count": "5",
                "GLU_involved_count": "3",
                "GLU_GLU_count": "1",
            },
        },
        "full_hexaplex_anti_parallel_30deg_ideal_heavy_deduped": {
            "scaffold_hexad_or_other": {
                "contact_count": "12",
                "GLU_involved_count": "0",
                "GLU_GLU_count": "0",
            }
        },
        "scaffold_blocks_1_heavy_deduped": {
            "scaffold_between_blocks": {
                "contact_count": "0",
                "GLU_involved_count": "0",
                "GLU_GLU_count": "0",
            }
        },
        "scaffold_blocks_1_2_heavy_deduped": {
            "scaffold_between_blocks": {
                "contact_count": "3",
                "GLU_involved_count": "2",
                "GLU_GLU_count": "1",
            }
        },
    }

    report.write_markdown_report(summaries, out)

    text = out.read_text(encoding="utf-8")
    assert "reciprocal-space scaffold signature" in text
    assert "candidate contiguous-residue mapping" in text
    assert "one candidate block already spans nearly full angular coverage" in text
    assert "does not prove temporal order" in text
