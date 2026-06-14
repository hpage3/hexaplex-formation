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


transition = _load_script_module("analyze_seed_6_7_transition", "scripts/analyze_seed_6_7_transition.py")


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def seed_row(unit: int, ensemble: str, index: int, offset: float) -> dict[str, str]:
    return {
        "unit_count": str(unit),
        "sample_id": f"central{unit}_{ensemble}_{index:04d}",
        "ensemble_type": ensemble,
        "radius_of_gyration_A": str(10 + offset),
        "RMSD_to_formed_seed_A": str(5 - offset),
        "contact_fraction_vs_target": str(0.2 + offset),
        "CYP_MEP_contact_fraction_vs_target": str(0.3 + offset),
        "axial_register_score": str(0.4 + offset),
        "compactness_score": str(0.5 + offset),
        "seed_formation_score": str(0.6 + offset),
    }


def test_summary_and_comparison_are_written(tmp_path):
    seed_rows = []
    for unit, unit_offset in [(6, 0.0), (7, 0.1)]:
        for ensemble, ensemble_offset in [
            ("loose_initial", 0.0),
            ("angular_randomized_loose_initial", 0.02),
            ("formed_perturbed", 0.4),
        ]:
            for index in range(2):
                seed_rows.append(seed_row(unit, ensemble, index, unit_offset + ensemble_offset + index * 0.01))
    seed_csv = tmp_path / "seed.csv"
    write_csv(seed_csv, seed_rows)

    contact_csv = tmp_path / "contacts.csv"
    write_csv(
        contact_csv,
        [
            {
                "variant_id": "central6_units",
                "units_per_chain": "6",
                "CYP_MEP_contact_count": "100",
                "CYP_MEP_contact_fraction": "0.8",
                "GLU_contact_count": "20",
                "backbone_like_contact_count": "10",
                "total_interchain_contacts": "200",
                "contacts_per_unit": "33.3",
                "perturbation_contact_fraction_vs_reference_mean": "0.7",
                "perturbation_chain_graph_connected_probability": "1.0",
            },
            {
                "variant_id": "central7_units",
                "units_per_chain": "7",
                "CYP_MEP_contact_count": "150",
                "CYP_MEP_contact_fraction": "0.9",
                "GLU_contact_count": "30",
                "backbone_like_contact_count": "15",
                "total_interchain_contacts": "280",
                "contacts_per_unit": "40.0",
                "perturbation_contact_fraction_vs_reference_mean": "0.8",
                "perturbation_chain_graph_connected_probability": "1.0",
            },
        ],
    )

    activation_csv = tmp_path / "activation.csv"
    activation_rows = []
    for unit, time in [(6, 0.6), (7, 0.5)]:
        for feature in ["CYP_MEP_contact_fraction_vs_target", "contact_fraction_vs_target", "compactness_score"]:
            activation_rows.append(
                {
                    "unit_count": str(unit),
                    "feature": feature,
                    "threshold_fraction": "0.75",
                    "median_activation_time": str(time),
                }
            )
    write_csv(activation_csv, activation_rows)

    geometry_csv = tmp_path / "geometry.csv"
    write_csv(
        geometry_csv,
        [
            {"variant_id": "central6_units", "units_per_chain": "6", "axial_extent_A": "20", "radial_extent_mean_A": "5", "suspicious_overlap_count": "0"},
            {"variant_id": "central7_units", "units_per_chain": "7", "axial_extent_A": "24", "radial_extent_mean_A": "5.2", "suspicious_overlap_count": "0"},
        ],
    )
    helicity_csv = tmp_path / "helicity.csv"
    write_csv(
        helicity_csv,
        [
            {"variant_id": "central6_units", "units_per_chain": "6", "mean_helical_r2": "0.7", "mean_twist_per_unit_deg": "30", "helical_coherence_score": "0.8"},
            {"variant_id": "central7_units", "units_per_chain": "7", "mean_helical_r2": "0.75", "mean_twist_per_unit_deg": "30", "helical_coherence_score": "0.85"},
        ],
    )

    args = type(
        "Args",
        (),
        {
            "seed_order_csv": seed_csv,
            "contact_summary_csv": contact_csv,
            "bridge_activation_csv": activation_csv,
            "bridge_ordering_csv": tmp_path / "unused_ordering.csv",
            "angular_bridge_activation_csv": activation_csv,
            "angular_bridge_ordering_csv": tmp_path / "unused_angular_ordering.csv",
            "mini_geometry_csv": geometry_csv,
            "mini_helicity_csv": helicity_csv,
            "summary_csv": tmp_path / "summary.csv",
            "comparison_csv": tmp_path / "comparison.csv",
            "plot_dir": tmp_path / "plots",
            "report": tmp_path / "report.md",
        },
    )()

    result = transition.run(args)

    assert result["summary_rows"] == 36
    assert result["comparison_rows"] > 20
    assert args.summary_csv.exists()
    assert args.comparison_csv.exists()
    assert args.report.exists()

    comparison_text = args.comparison_csv.read_text(encoding="utf-8")
    assert "CYP_MEP_contact_count" in comparison_text
    assert "formed-minus-loose separation" in comparison_text
    assert "7-unit stronger" in comparison_text
