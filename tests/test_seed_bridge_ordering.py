import csv
import importlib.util
import sys
from pathlib import Path

import pytest

from hexaplex_formation.pdb_utils import PDBAtom


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_script_module(name: str, relative_path: str):
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


bridge = _load_script_module("analyze_seed_bridge_ordering", "scripts/analyze_seed_bridge_ordering.py")


def atom(serial: int, chain_id: str, x: float, y: float, z: float = 0.0) -> PDBAtom:
    return PDBAtom(
        record_type="ATOM",
        atom_serial=serial,
        atom_name="CA",
        alt_loc="",
        residue_name="CYP",
        chain_id=chain_id,
        residue_number=1,
        insertion_code="",
        x=x,
        y=y,
        z=z,
        occupancy=1.0,
        temp_factor=0.0,
        element="C",
    )


def sample(sample_id: str, values: dict[str, float]) -> bridge.Sample:
    return bridge.Sample(sample_id=sample_id, raw={"sample_id": sample_id}, features=values)


def feature_values(base: float) -> dict[str, float]:
    return {
        "compactness_score": base,
        "contact_fraction_vs_target": base + 0.01,
        "CYP_MEP_contact_fraction_vs_target": base + 0.02,
        "axial_register_score": base + 0.03,
        "angular_phase_order_score": base + 0.04,
        "refined_angular_phase_score": base + 0.045,
        "rmsd_formedness_score": base + 0.05,
        "seed_formation_score": base + 0.06,
    }


def test_metric_orientation_larger_is_more_formed_transformations():
    row = {
        "compactness_score": "0.4",
        "RMSD_to_formed_seed_A": "5.0",
    }

    assert bridge.oriented_feature_value(row, "compactness_score", "compactness_score") == pytest.approx(0.4)
    assert bridge.oriented_feature_value(row, "rmsd_formedness_score", "RMSD_to_formed_seed_A") == pytest.approx(0.5)
    assert bridge.rmsd_formedness_score(0.0) > bridge.rmsd_formedness_score(10.0)


def test_threshold_calculation_works():
    assert bridge.threshold_value(0.2, 0.8, 0.75) == pytest.approx(0.65)
    assert bridge.threshold_value(1.0, 0.5, 0.5) == pytest.approx(0.75)


def test_activation_time_detection_for_monotonic_trajectory():
    times = [0.0, 0.25, 0.5, 0.75, 1.0]
    values = [0.0, 0.2, 0.6, 0.8, 1.0]

    assert bridge.activation_time(times, values, 0.6) == pytest.approx(0.5)


def test_activation_time_returns_nan_for_non_crossing_trajectory():
    value = bridge.activation_time([0.0, 0.5, 1.0], [0.0, 0.1, 0.2], 0.5)

    assert value != value


def test_refined_angular_phase_allows_global_rotation_but_detects_label_swap():
    reference = [
        atom(1, "A", 1.0, 0.0, -1.0),
        atom(2, "A", 1.0, 0.0, 1.0),
        atom(3, "B", 0.0, 1.0, -1.0),
        atom(4, "B", 0.0, 1.0, 1.0),
        atom(5, "C", -1.0, 0.0, -1.0),
        atom(6, "C", -1.0, 0.0, 1.0),
    ]
    rotated = [
        atom(1, "A", 0.0, 1.0, -1.0),
        atom(2, "A", 0.0, 1.0, 1.0),
        atom(3, "B", -1.0, 0.0, -1.0),
        atom(4, "B", -1.0, 0.0, 1.0),
        atom(5, "C", 0.0, -1.0, -1.0),
        atom(6, "C", 0.0, -1.0, 1.0),
    ]
    swapped = [
        atom(1, "A", -1.0, 0.0, -1.0),
        atom(2, "A", -1.0, 0.0, 1.0),
        atom(3, "B", 0.0, 1.0, -1.0),
        atom(4, "B", 0.0, 1.0, 1.0),
        atom(5, "C", 1.0, 0.0, -1.0),
        atom(6, "C", 1.0, 0.0, 1.0),
    ]

    assert bridge.refined_angular_phase_score(rotated, reference) > 0.97
    assert bridge.refined_angular_phase_score(swapped, reference) < 0.9


def test_feature_ordering_ranking_handles_ties():
    ranks = bridge.average_ranks({"a": 0.2, "b": 0.2, "c": 0.8})

    assert ranks["a"] == pytest.approx(1.5)
    assert ranks["b"] == pytest.approx(1.5)
    assert ranks["c"] == pytest.approx(3.0)


def test_endpoint_matching_returns_expected_pair_count_on_synthetic_data():
    features = [feature for feature, _, _ in bridge.selected_feature_specs("legacy")]
    loose = [sample("l0", feature_values(0.0)), sample("l1", feature_values(0.2))]
    formed = [sample("f0", feature_values(0.8)), sample("f1", feature_values(1.0))]

    pairs, method, note = bridge.endpoint_matching(loose, formed, features, epsilon=0.5, iterations=100)

    assert len(pairs) == 2
    assert {pair.loose.sample_id for pair in pairs} == {"l0", "l1"}
    assert {pair.formed.sample_id for pair in pairs} == {"f0", "f1"}
    assert method in {"sinkhorn_greedy_unique", "greedy_minimum_cost_fallback"}
    assert note


def test_output_csvs_are_written_with_expected_columns(tmp_path):
    input_csv = tmp_path / "seed.csv"
    fieldnames = [
        "unit_count",
        "sample_id",
        "ensemble_type",
        "compactness_score",
        "contact_fraction_vs_target",
        "CYP_MEP_contact_fraction_vs_target",
        "axial_register_score",
        "angular_phase_order_score",
        "RMSD_to_formed_seed_A",
        "seed_formation_score",
    ]
    with input_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for index in range(2):
            writer.writerow(
                {
                    "unit_count": "4",
                    "sample_id": f"loose_{index}",
                    "ensemble_type": "loose_initial",
                    "compactness_score": str(0.4 + index * 0.01),
                    "contact_fraction_vs_target": str(0.1 + index * 0.01),
                    "CYP_MEP_contact_fraction_vs_target": str(0.1 + index * 0.01),
                    "axial_register_score": str(0.5 + index * 0.01),
                    "angular_phase_order_score": str(0.5 + index * 0.01),
                    "RMSD_to_formed_seed_A": "15.0",
                    "seed_formation_score": str(0.4 + index * 0.01),
                }
            )
            writer.writerow(
                {
                    "unit_count": "4",
                    "sample_id": f"formed_{index}",
                    "ensemble_type": "formed_perturbed",
                    "compactness_score": str(0.9 + index * 0.01),
                    "contact_fraction_vs_target": str(0.8 + index * 0.01),
                    "CYP_MEP_contact_fraction_vs_target": str(0.8 + index * 0.01),
                    "axial_register_score": str(0.9 + index * 0.01),
                    "angular_phase_order_score": str(0.9 + index * 0.01),
                    "RMSD_to_formed_seed_A": "0.5",
                    "seed_formation_score": str(0.9 + index * 0.01),
                }
            )

    args = type(
        "Args",
        (),
        {
            "input_csv": input_csv,
            "contact_network_summary": tmp_path / "unused.csv",
            "structures_dir": tmp_path / "structures",
            "ensemble_dir": tmp_path / "ensembles",
            "unit_counts": "4",
            "loose_ensemble_type": "loose_initial",
            "angular_phase_mode": "legacy",
            "coordinate_backed_only": False,
            "n_time_points": 3,
            "time_points": "",
            "threshold_mode": "formed_fraction",
            "formed_fraction": 0.75,
            "formed_fractions": "",
            "sinkhorn_epsilon": 0.5,
            "sinkhorn_iterations": 100,
            "max_pairs": 0,
            "pairs_csv": tmp_path / "pairs.csv",
            "paths_csv": tmp_path / "paths.csv",
            "activation_csv": tmp_path / "activation.csv",
            "ordering_csv": tmp_path / "ordering.csv",
            "plot_dir": tmp_path / "plots",
            "out_report": tmp_path / "report.md",
        },
    )()

    bridge.run(args)

    with args.pairs_csv.open("r", newline="", encoding="utf-8") as handle:
        assert csv.DictReader(handle).fieldnames == bridge.PAIR_COLUMNS
    with args.paths_csv.open("r", newline="", encoding="utf-8") as handle:
        assert csv.DictReader(handle).fieldnames == bridge.PATH_COLUMNS
    with args.activation_csv.open("r", newline="", encoding="utf-8") as handle:
        assert csv.DictReader(handle).fieldnames == bridge.ACTIVATION_COLUMNS
    with args.ordering_csv.open("r", newline="", encoding="utf-8") as handle:
        assert csv.DictReader(handle).fieldnames == bridge.ORDERING_COLUMNS
    assert args.out_report.exists()
