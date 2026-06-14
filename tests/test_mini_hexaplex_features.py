import csv
import importlib.util
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_script_module(name: str, relative_path: str):
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


features = _load_script_module("analyze_mini_hexaplex_features", "scripts/analyze_mini_hexaplex_features.py")


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def test_feature_window_extraction_works_on_synthetic_radial_data(tmp_path):
    profile_dir = tmp_path / "profiles"
    manifest = tmp_path / "manifest.csv"
    geometry = tmp_path / "geometry.csv"
    _write_csv(
        manifest,
        ["variant_id", "truncation_rule", "units_per_chain", "residues_per_chain", "total_residue_count", "warnings"],
        [
            {
                "variant_id": "central8_units",
                "truncation_rule": "central test",
                "units_per_chain": "8",
                "residues_per_chain": "A:16",
                "total_residue_count": "96",
                "warnings": "",
            },
            {
                "variant_id": "central12_units",
                "truncation_rule": "central test 12",
                "units_per_chain": "12",
                "residues_per_chain": "A:24",
                "total_residue_count": "144",
                "warnings": "",
            }
        ],
    )
    _write_csv(
        geometry,
        ["variant_id", "warnings", "structural_coherence_flag"],
        [
            {"variant_id": "central8_units", "warnings": "", "structural_coherence_flag": "coherent"},
            {"variant_id": "central12_units", "warnings": "", "structural_coherence_flag": "coherent"},
        ],
    )
    _write_csv(
        profile_dir / "full_length_baseline_radial.csv",
        ["q_Ainv", "d_A", "intensity"],
        [
            {"q_Ainv": "1.20", "d_A": "5.2", "intensity": "1.0"},
            {"q_Ainv": "1.30", "d_A": "4.85", "intensity": "10.0"},
            {"q_Ainv": "1.36", "d_A": "4.62", "intensity": "10.0"},
            {"q_Ainv": "1.85", "d_A": "3.4", "intensity": "4.0"},
            {"q_Ainv": "2.09", "d_A": "3.0", "intensity": "3.0"},
        ],
    )
    _write_csv(
        profile_dir / "central8_units_radial.csv",
        ["q_Ainv", "d_A", "intensity"],
        [
            {"q_Ainv": "1.18", "d_A": "5.3", "intensity": "1.0"},
            {"q_Ainv": "1.26", "d_A": "5.0", "intensity": "2.0"},
            {"q_Ainv": "1.32", "d_A": "4.75", "intensity": "8.0"},
            {"q_Ainv": "1.39", "d_A": "4.52", "intensity": "2.0"},
            {"q_Ainv": "1.85", "d_A": "3.4", "intensity": "5.0"},
            {"q_Ainv": "2.09", "d_A": "3.0", "intensity": "6.0"},
        ],
    )
    _write_csv(
        profile_dir / "central12_units_radial.csv",
        ["q_Ainv", "d_A", "intensity"],
        [
            {"q_Ainv": "1.18", "d_A": "5.3", "intensity": "1.0"},
            {"q_Ainv": "1.26", "d_A": "5.0", "intensity": "3.0"},
            {"q_Ainv": "1.32", "d_A": "4.75", "intensity": "12.0"},
            {"q_Ainv": "1.39", "d_A": "4.52", "intensity": "3.0"},
            {"q_Ainv": "1.85", "d_A": "3.4", "intensity": "7.0"},
            {"q_Ainv": "2.09", "d_A": "3.0", "intensity": "8.0"},
        ],
    )

    rows = features.build_feature_rows(profile_dir, manifest, geometry, {"central8_units", "central12_units"})

    assert len(rows) == 3
    row = {row["variant_id"]: row for row in rows}["central8_units"]
    assert row["variant_id"] == "central8_units"
    assert row["structural_coherence_flag"] == "coherent"
    assert row["residues_per_chain"] == "A:16"
    assert row["total_residue_count"] == "96"
    assert row["d_A_at_max_in_4p5_5A_window"] == "4.750000"
    assert row["integrated_intensity_4p5_5A"] == "12.000000"
    assert row["ratio_to_full_length_4p5_5A"] == "0.600000"
    assert row["ratio_to_matching_8unit_model_4p5_5A"] == ""
    assert row["integrated_intensity_3p4A"] == "5.000000"
    assert row["integrated_intensity_3p0A"] == "6.000000"
    assert row["has_local_maximum_4p5_5A"] == "yes"
    assert "ratio vs full length = 0.600" in row["comparison_to_full_length"]
    assert row["comparison_to_8unit_models"] == "8-unit reference"

    row12 = {row["variant_id"]: row for row in rows}["central12_units"]
    assert row12["units_per_chain"] == "12"
    assert row12["structural_coherence_flag"] == "coherent"
    assert row12["total_residue_count"] == "144"
    assert row12["integrated_intensity_4p5_5A"] == "18.000000"
    assert row12["ratio_to_full_length_4p5_5A"] == "0.900000"
    assert row12["ratio_to_matching_8unit_model_4p5_5A"] == "1.500000"
    assert "ratio vs 8-unit = 1.500" in row12["comparison_to_8unit_models"]

    baseline = {row["variant_id"]: row for row in rows}["full_length_baseline"]
    assert baseline["units_per_chain"] == "15"
    assert baseline["structural_coherence_flag"] == "coherent"
    assert baseline["ratio_to_full_length_4p5_5A"] == "1.000000"

    plot_dir = tmp_path / "plots"
    plot_paths = features.write_plots(profile_dir, rows, plot_dir)
    expected_plots = {
        plot_dir / "mini_hexaplex_units_vs_integrated_4p5_5p0.png",
        plot_dir / "mini_hexaplex_units_vs_ratio_to_full_4p5_5p0.png",
        plot_dir / "mini_hexaplex_units_vs_3p4A_intensity.png",
        plot_dir / "mini_hexaplex_units_vs_3p0A_intensity.png",
        plot_dir / "mini_hexaplex_central_d_profile_overlay.png",
        plot_dir / "mini_hexaplex_lower_end_d_profile_overlay.png",
    }
    assert expected_plots.issubset(set(plot_paths))
    for path in expected_plots:
        assert path.exists()


def test_report_contains_required_cautions(tmp_path):
    report = tmp_path / "report.md"
    rows = [
        {
            "variant_id": "central8_units",
            "units_per_chain": "8",
            "structural_coherence_flag": "coherent",
            "d_A_at_max_in_4p5_5A_window": "4.75",
            "integrated_intensity_4p5_5A": "12.0",
            "ratio_to_full_length_4p5_5A": "0.6",
            "ratio_to_matching_8unit_model_4p5_5A": "",
            "has_local_maximum_4p5_5A": "yes",
            "integrated_intensity_3p4A": "5.0",
            "integrated_intensity_3p0A": "6.0",
            "comparison_to_8unit_models": "8-unit reference",
        }
    ]

    features.write_report(
        rows,
        report,
        [],
        tmp_path / "manifest.csv",
        tmp_path / "geometry.csv",
        "baseline.pdb",
        "debye",
    )

    text = report.read_text(encoding="utf-8")
    assert "length response" in text.lower()
    assert "anti-parallel" in text
    assert "q = 2*pi/d" in text
    assert "coordinate truncations, not independently relaxed" in text
