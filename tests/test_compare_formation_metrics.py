import csv
import importlib.util
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_script_module(name: str, relative_path: str):
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


compare = _load_script_module("compare_formation_metrics", "scripts/compare_formation_metrics.py")


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def test_structure_base_from_metric_names():
    assert (
        compare.structure_base_from_name("hexaplex_scaffold_only_complement_heavy_deduped_contacts_4p5A.csv")
        == "hexaplex_scaffold_only_complement"
    )
    assert (
        compare.structure_base_from_name("full_hexaplex_anti_parallel_30deg_ideal_allatom_deduped_hbond_candidates.csv")
        == "full_hexaplex_anti_parallel_30deg_ideal"
    )
    assert (
        compare.structure_base_from_name("outputs/intermediates/normalized_structures/alanine_beta_sheet_full_heavy.pdb")
        == "alanine_beta_sheet_full"
    )


def test_read_motif_csv_uses_zero_for_absent_motifs(tmp_path):
    path = tmp_path / "example_heavy_deduped_glu_motifs.csv"
    _write_csv(
        path,
        ["motif_type", "count"],
        [
            {"motif_type": "GLU-GLU", "count": "3"},
            {"motif_type": "GLU-any", "count": "5"},
        ],
    )

    counts = compare.read_motif_counts(path)

    assert counts["GLU-GLU"] == 3
    assert counts["GLU-any"] == 5
    assert counts["other_GLU_contact"] == 0
    assert counts["backbone_N_to_GLU_sidechain_O"] == 0


def test_aggregation_from_temporary_fixtures(tmp_path):
    formation_dir = tmp_path / "formation"
    normalization = tmp_path / "structure_normalization_summary.csv"
    base = "hexaplex_scaffold_only_complement"

    _write_csv(
        normalization,
        [
            "structure_file",
            "original_atom_count",
            "allatom_deduped_atom_count",
            "heavy_atom_count",
            "heavy_deduped_atom_count",
            "hydrogens_detected",
            "output_heavy_deduped_pdb",
        ],
        [
            {
                "structure_file": f"inputs/structures/{base}.pdb",
                "original_atom_count": "20",
                "allatom_deduped_atom_count": "18",
                "heavy_atom_count": "12",
                "heavy_deduped_atom_count": "10",
                "hydrogens_detected": "8",
                "output_heavy_deduped_pdb": f"outputs/intermediates/normalized_structures/{base}_heavy_deduped.pdb",
            }
        ],
    )
    _write_csv(
        formation_dir / f"{base}_heavy_deduped_contacts_4p5A.csv",
        ["residue_name_i", "residue_name_j"],
        [
            {"residue_name_i": "GLU", "residue_name_j": "ALA"},
            {"residue_name_i": "GLU", "residue_name_j": "GLU"},
            {"residue_name_i": "ALA", "residue_name_j": "ALA"},
        ],
    )
    _write_csv(
        formation_dir / f"{base}_heavy_deduped_glu_motifs.csv",
        ["motif_type", "count"],
        [
            {"motif_type": "GLU-GLU", "count": "1"},
            {"motif_type": "GLU-any", "count": "2"},
        ],
    )
    _write_csv(
        formation_dir / f"{base}_heavy_deduped_helical_order_summary.csv",
        ["residue_count", "mean_radius_xy", "min_radius_xy", "max_radius_xy", "z_span", "angular_coverage_rad"],
        [
            {
                "residue_count": "4",
                "mean_radius_xy": "2.5",
                "min_radius_xy": "2.0",
                "max_radius_xy": "3.0",
                "z_span": "8.0",
                "angular_coverage_rad": "5.5",
            }
        ],
    )
    _write_csv(
        formation_dir / f"{base}_allatom_deduped_hbond_candidates.csv",
        ["donor_residue", "acceptor_residue"],
        [
            {"donor_residue": "A:GLU1", "acceptor_residue": "B:ALA2"},
            {"donor_residue": "A:ALA3", "acceptor_residue": "B:ALA4"},
        ],
    )

    rows = compare.aggregate_metrics(formation_dir, normalization)

    assert len(rows) == 1
    row = rows[0]
    assert row["structure_base"] == base
    assert row["heavy_deduped_atom_count"] == "10"
    assert row["contact_count_4p5A"] == "3"
    assert row["GLU_contact_count"] == "2"
    assert row["GLU_GLU_contact_count"] == "1"
    assert row["contact_count_per_residue"] == "0.750000"
    assert row["motif_GLU_GLU"] == "1"
    assert row["motif_other_GLU_contact"] == "0"
    assert row["hbond_candidate_count"] == "2"
    assert row["GLU_involved_hbond_candidate_count"] == "1"
    assert row["hbond_candidates_per_residue"] == "0.500000"


def test_markdown_report_includes_caution_language(tmp_path):
    out = tmp_path / "formation_comparison.md"
    rows = [
        {
            "structure_base": "alanine_beta_sheet_full",
            "residue_count": "10",
            "heavy_deduped_atom_count": "100",
            "contact_count_4p5A": "4",
            "motif_GLU_GLU": "0",
            "motif_GLU_any": "0",
            "hbond_candidate_count": "0",
            "GLU_involved_hbond_candidate_count": "0",
            "mean_radius_xy": "1.0",
            "z_span": "2.0",
            "angular_coverage_rad": "1.0",
        },
        {
            "structure_base": "hexaplex_scaffold_only_complement",
            "residue_count": "20",
            "heavy_deduped_atom_count": "200",
            "contact_count_4p5A": "8",
            "motif_GLU_GLU": "3",
            "motif_GLU_any": "6",
            "hbond_candidate_count": "2",
            "GLU_involved_hbond_candidate_count": "2",
            "mean_radius_xy": "5.0",
            "z_span": "9.0",
            "angular_coverage_rad": "6.0",
        },
    ]

    compare.write_markdown_report(rows, out)

    text = out.read_text(encoding="utf-8")
    assert "reciprocal-space feature" in text
    assert "does not prove" in text
    assert "requires follow-up" in text
    assert "supports the working hypothesis" in text
