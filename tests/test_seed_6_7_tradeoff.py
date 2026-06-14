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


tradeoff = _load_script_module("analyze_seed_6_7_tradeoff", "scripts/analyze_seed_6_7_tradeoff.py")


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def test_classification_respects_orientation_and_negligible_threshold():
    assert tradeoff.classify(10.0, 12.0, "higher", 0.005) == "favors_7"
    assert tradeoff.classify(10.0, 12.0, "lower", 0.005) == "favors_6"
    assert tradeoff.classify(10.0, 10.01, "higher", 0.005) == "mixed_or_negligible"


def test_tradeoff_outputs_are_written(tmp_path):
    contact_summary = tmp_path / "contact_summary.csv"
    write_csv(
        contact_summary,
        [
            {
                "variant_id": "central6_units",
                "units_per_chain": "6",
                "total_interchain_contacts": "100",
                "CYP_MEP_contact_count": "80",
                "GLU_contact_count": "12",
                "backbone_like_contact_count": "8",
                "contacts_per_unit": "16.7",
                "contacts_per_chain_pair_mean": "20",
                "contacts_per_chain_pair_max": "30",
                "CYP_MEP_contact_fraction": "0.90",
                "perturbation_contact_fraction_vs_reference_mean": "0.82",
                "perturbation_contact_fraction_vs_reference_std": "0.02",
                "contact_redundancy_score": "0.8",
                "nucleation_network_score": "0.8",
            },
            {
                "variant_id": "central7_units",
                "units_per_chain": "7",
                "total_interchain_contacts": "150",
                "CYP_MEP_contact_count": "120",
                "GLU_contact_count": "20",
                "backbone_like_contact_count": "10",
                "contacts_per_unit": "21.4",
                "contacts_per_chain_pair_mean": "25",
                "contacts_per_chain_pair_max": "35",
                "CYP_MEP_contact_fraction": "0.88",
                "perturbation_contact_fraction_vs_reference_mean": "0.78",
                "perturbation_contact_fraction_vs_reference_std": "0.04",
                "contact_redundancy_score": "0.9",
                "nucleation_network_score": "0.85",
            },
        ],
    )

    contact_edges = tmp_path / "contact_edges.csv"
    write_csv(
        contact_edges,
        [
            {
                "variant_id": "central6_units",
                "units_per_chain": "6",
                "edge_scope": "chain",
                "unit_a": "",
                "unit_b": "",
                "contact_count": "100",
                "CYP_MEP_vs_CYP_MEP_count": "70",
                "CYP_MEP_vs_GLU_count": "10",
                "GLU_vs_GLU_count": "20",
            },
            {
                "variant_id": "central7_units",
                "units_per_chain": "7",
                "edge_scope": "chain",
                "unit_a": "",
                "unit_b": "",
                "contact_count": "150",
                "CYP_MEP_vs_CYP_MEP_count": "100",
                "CYP_MEP_vs_GLU_count": "20",
                "GLU_vs_GLU_count": "30",
            },
            {
                "variant_id": "central6_units",
                "units_per_chain": "6",
                "edge_scope": "unit",
                "unit_a": "1",
                "unit_b": "2",
                "contact_count": "10",
                "CYP_MEP_vs_CYP_MEP_count": "10",
                "CYP_MEP_vs_GLU_count": "0",
                "GLU_vs_GLU_count": "0",
            },
            {
                "variant_id": "central7_units",
                "units_per_chain": "7",
                "edge_scope": "unit",
                "unit_a": "7",
                "unit_b": "6",
                "contact_count": "30",
                "CYP_MEP_vs_CYP_MEP_count": "30",
                "CYP_MEP_vs_GLU_count": "0",
                "GLU_vs_GLU_count": "0",
            },
        ],
    )

    transition = tmp_path / "transition.csv"
    rows = []
    for unit, formed_contact, compactness in [(6, 0.80, 0.99), (7, 0.81, 0.98)]:
        for ensemble in ["formed_perturbed", "loose_initial", "angular_randomized_loose_initial"]:
            for feature, value in [
                ("CYP_MEP_contact_fraction_vs_target", formed_contact),
                ("contact_fraction_vs_target", formed_contact),
                ("compactness_score", compactness),
                ("axial_register_score", compactness),
                ("RMSD_to_formed_seed_A", 0.6 if unit == 6 else 0.7),
                ("seed_formation_score", compactness),
            ]:
                rows.append({"unit_count": str(unit), "ensemble_type": ensemble, "feature": feature, "mean": str(value)})
    write_csv(transition, rows)

    activation = tmp_path / "activation.csv"
    write_csv(
        activation,
        [
            {"unit_count": "6", "feature": "CYP_MEP_contact_fraction_vs_target", "threshold_fraction": "0.75", "median_activation_time": "0.8"},
            {"unit_count": "7", "feature": "CYP_MEP_contact_fraction_vs_target", "threshold_fraction": "0.75", "median_activation_time": "0.7"},
        ],
    )

    args = type(
        "Args",
        (),
        {
            "transition_summary": transition,
            "transition_comparison": tmp_path / "unused_transition_comparison.csv",
            "seed_order_csv": tmp_path / "unused_seed_order.csv",
            "contact_summary_csv": contact_summary,
            "contact_edges_csv": contact_edges,
            "endpoint_means_csv": tmp_path / "unused_endpoint.csv",
            "angular_activation_csv": activation,
            "angular_ordering_csv": tmp_path / "unused_ordering.csv",
            "summary_csv": tmp_path / "summary.csv",
            "feature_table_csv": tmp_path / "features.csv",
            "report": tmp_path / "report.md",
            "plot_dir": tmp_path / "plots",
            "negligible_relative": 0.005,
        },
    )()

    result = tradeoff.run(args)

    assert result["feature_rows"] > 20
    assert result["summary_rows"] >= 4
    assert args.summary_csv.exists()
    assert args.feature_table_csv.exists()
    assert args.report.exists()
    assert result["plots"]

    table = args.feature_table_csv.read_text(encoding="utf-8")
    assert "absolute_network_growth" in table
    assert "perturbation_retention" in table
    assert "favors_7" in table
    assert "favors_6" in table
