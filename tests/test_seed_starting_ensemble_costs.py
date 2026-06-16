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


costs = _load_script_module("analyze_seed_starting_ensemble_costs", "scripts/analyze_seed_starting_ensemble_costs.py")


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def synthetic_row(unit: int, ensemble: str, sample: int, offset: float) -> dict[str, str]:
    return {
        "unit_count": str(unit),
        "sample_id": f"{ensemble}_{unit}_{sample}",
        "ensemble_type": ensemble,
        "RMSD_to_formed_seed_A": str(1.0 + offset),
        "compactness_score": str(0.9 - offset),
        "contact_fraction_vs_target": str(0.8 - offset),
        "CYP_MEP_contact_fraction_vs_target": str(0.7 - offset),
        "axial_register_score": str(0.85 - offset),
        "angular_phase_order_score": str(0.75 - offset),
        "refined_angular_phase_score": str(0.65 - offset),
    }


def run_args(tmp_path: Path, input_csv: Path):
    return type(
        "Args",
        (),
        {
            "input_csv": input_csv,
            "components_csv": tmp_path / "components.csv",
            "summary_csv": tmp_path / "summary.csv",
            "report": tmp_path / "report.md",
            "plot_dir": tmp_path / "plots",
            "reference_class": "formed_perturbed",
            "unit_counts": "",
            "std_floor": 1e-3,
        },
    )()


def test_cost_audit_runs_on_small_synthetic_input(tmp_path):
    rows = []
    for unit in [6, 7]:
        for index in range(3):
            rows.append(synthetic_row(unit, "formed_perturbed", index, index * 0.01))
            rows.append(synthetic_row(unit, "loose_initial", index, 0.25 + index * 0.01))
            rows.append(synthetic_row(unit, "angular_randomized_loose_initial", index, 0.35 + index * 0.01))
    input_csv = tmp_path / "seed.csv"
    write_csv(input_csv, rows)

    result = costs.run(run_args(tmp_path, input_csv))

    assert result["component_rows"] > 0
    assert result["summary_rows"] == 4
    assert (tmp_path / "components.csv").exists()
    assert (tmp_path / "summary.csv").exists()
    assert (tmp_path / "report.md").exists()


def test_missing_optional_feature_columns_are_skipped(tmp_path):
    rows = [
        {
            "unit_count": "6",
            "sample_id": "formed_0",
            "ensemble_type": "formed_perturbed",
            "RMSD_to_formed_seed_A": "1.0",
            "contact_fraction_vs_target": "0.9",
        },
        {
            "unit_count": "6",
            "sample_id": "loose_0",
            "ensemble_type": "loose_initial",
            "RMSD_to_formed_seed_A": "2.0",
            "contact_fraction_vs_target": "0.4",
        },
    ]
    input_csv = tmp_path / "seed.csv"
    write_csv(input_csv, rows)

    costs.run(run_args(tmp_path, input_csv))

    report = (tmp_path / "report.md").read_text(encoding="utf-8")
    assert "missing column" in report
    with (tmp_path / "components.csv").open(newline="", encoding="utf-8") as handle:
        component_rows = list(csv.DictReader(handle))
    assert {row["feature_name"] for row in component_rows} == {"RMSD_to_formed_seed_A", "contact_fraction_vs_target"}


def test_tiny_formed_reference_std_uses_floor(tmp_path):
    rows = []
    for index in range(2):
        rows.append(synthetic_row(6, "formed_perturbed", index, 0.0))
        rows.append(synthetic_row(6, "loose_initial", index, 0.1))
    input_csv = tmp_path / "seed.csv"
    write_csv(input_csv, rows)

    costs.run(run_args(tmp_path, input_csv))

    with (tmp_path / "components.csv").open(newline="", encoding="utf-8") as handle:
        component_rows = list(csv.DictReader(handle))
    assert any(row["std_floor_used"] == "true" for row in component_rows)
    assert all(row["abs_z_deviation"] for row in component_rows)


def test_dominant_cost_group_is_assigned_from_available_groups():
    component_rows = [
        {
            "unit_count": "6",
            "start_class": "loose_initial",
            "feature_group": "geometric",
            "abs_z_deviation": "1.0",
        },
        {
            "unit_count": "6",
            "start_class": "loose_initial",
            "feature_group": "contact_recovery",
            "abs_z_deviation": "4.0",
        },
    ]

    summary = costs.build_summary_rows(component_rows)

    assert summary[0]["dominant_cost_group"] == "contact_recovery"


def test_output_csvs_contain_expected_required_columns(tmp_path):
    rows = []
    for index in range(2):
        rows.append(synthetic_row(6, "formed_perturbed", index, index * 0.01))
        rows.append(synthetic_row(6, "loose_initial", index, 0.2 + index * 0.01))
    input_csv = tmp_path / "seed.csv"
    write_csv(input_csv, rows)

    costs.run(run_args(tmp_path, input_csv))

    with (tmp_path / "components.csv").open(newline="", encoding="utf-8") as handle:
        assert csv.DictReader(handle).fieldnames == costs.COMPONENT_COLUMNS
    with (tmp_path / "summary.csv").open(newline="", encoding="utf-8") as handle:
        assert csv.DictReader(handle).fieldnames == costs.SUMMARY_COLUMNS


def test_report_contains_conservative_disclaimer_language(tmp_path):
    rows = []
    for index in range(2):
        rows.append(synthetic_row(6, "formed_perturbed", index, index * 0.01))
        rows.append(synthetic_row(6, "loose_initial", index, 0.2 + index * 0.01))
    input_csv = tmp_path / "seed.csv"
    write_csv(input_csv, rows)

    costs.run(run_args(tmp_path, input_csv))

    report = (tmp_path / "report.md").read_text(encoding="utf-8")
    assert "hypothesis-generating" in report
    assert "not an atomistic Schrodinger bridge" in report
    assert "not molecular dynamics" in report
    assert "not evidence of a physical nucleation pathway" in report
