import csv
import importlib.util
import math
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_script_module(name: str, relative_path: str):
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


workflow = _load_script_module(
    "analyze_length_twist_diffraction_sensitivity",
    "scripts/analyze_length_twist_diffraction_sensitivity.py",
)


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def test_target_window_q_conversion_uses_two_pi_over_d():
    rows = {row["feature"]: row for row in workflow.q_window_table()}

    assert float(rows["4.5-5.0 A"]["q_min_Ainv"]) == pytest.approx(2.0 * math.pi / 5.0)
    assert float(rows["4.5-5.0 A"]["q_max_Ainv"]) == pytest.approx(2.0 * math.pi / 4.5)
    assert float(rows["3.4 A"]["q_min_Ainv"]) == pytest.approx(2.0 * math.pi / 3.55)


def test_window_metrics_reports_fwhm_and_width_proxy_on_synthetic_peak():
    points = [
        {"q_Ainv": 1.20, "d_A": 5.236, "intensity": 0.0},
        {"q_Ainv": 1.26, "d_A": 4.987, "intensity": 5.0},
        {"q_Ainv": 1.30, "d_A": 4.833, "intensity": 10.0},
        {"q_Ainv": 1.34, "d_A": 4.684, "intensity": 5.0},
        {"q_Ainv": 1.40, "d_A": 4.488, "intensity": 0.0},
    ]

    metrics = workflow.window_metrics(points, 4.5, 5.0)

    assert metrics["point_count"] == 3
    assert metrics["d_at_max"] == pytest.approx(4.833)
    assert metrics["fwhm_q"] == pytest.approx(0.08)
    assert metrics["equivalent_width_q"] == pytest.approx(0.06)
    assert metrics["second_moment_width_q"] is not None


def test_feature_extraction_handles_synthetic_radial_profile(tmp_path):
    profile_dir = tmp_path / "profiles"
    diffraction_dir = tmp_path / "diffraction"
    model_path = tmp_path / "model.pdb"
    model_path.write_text(
        "\n".join(
            [
                "ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00  0.00           N",
                "ATOM      2  C   GLY B   1       1.000   0.000   0.000  1.00  0.00           C",
                "ATOM      3  C   GLY C   1       0.000   1.000   0.000  1.00  0.00           C",
                "ATOM      4  C   GLY D   1       0.000   0.000   1.000  1.00  0.00           C",
                "ATOM      5  C   GLY E   1       1.000   1.000   0.000  1.00  0.00           C",
                "ATOM      6  C   GLY F   1       1.000   0.000   1.000  1.00  0.00           C",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    _write_csv(
        profile_dir / "central4_units_radial.csv",
        ["q_Ainv", "d_A", "intensity"],
        [
            {"q_Ainv": "1.25", "d_A": "5.00", "intensity": "2.0"},
            {"q_Ainv": "1.32", "d_A": "4.75", "intensity": "8.0"},
            {"q_Ainv": "1.39", "d_A": "4.52", "intensity": "2.0"},
            {"q_Ainv": "1.85", "d_A": "3.40", "intensity": "4.0"},
            {"q_Ainv": "2.09", "d_A": "3.00", "intensity": "3.0"},
        ],
    )
    model = workflow.Model(
        model_id="length_4_central_twist_30",
        radial_id="central4_units",
        length_units=4,
        twist_deg=30.0,
        model_path=model_path,
        radial_profile_path=profile_dir / "central4_units_radial.csv",
        native_detector_npy=diffraction_dir / "central4_units.npy",
        status="available",
        source="test",
    )

    features, widths = workflow.feature_and_width_rows([model])

    assert features[0]["chain_count"] == "6"
    assert features[0]["integrated_intensity_4p5_5p0_A"] == "12.000000"
    assert features[0]["integrated_intensity_3p4_A"] == "4.000000"
    assert features[0]["ratio_integrated_4p5_5p0_A_to_3p4_A"] == "3.000000"
    assert {row["window_name"] for row in widths} >= {"4.5-5.0 A", "3.4 A", "3.0 A", "4.1-8.4 A"}


def test_manifest_marks_planned_twist_variants_without_builder(tmp_path):
    class Args:
        profile_dir = tmp_path / "profiles"
        diffraction_dir = tmp_path / "diffraction"
        include_planned_twist_rows = True

    models = workflow.build_models(Args(), builder_available=False)
    by_id = {model.model_id: model for model in models}

    assert by_id["full_length_twist_30"].twist_deg == 30.0
    assert by_id["full_length_twist_30"].status == "missing_radial_profile"
    assert by_id["full_length_twist_28_planned"].status == "pending_official_builder"
    assert "official builder not found locally" in by_id["full_length_twist_28_planned"].warnings


def test_builder_stub_documents_required_inputs(tmp_path):
    stub = tmp_path / "run_pnab_twist_builder_template.sh"

    workflow.write_builder_stub(stub)

    text = stub.read_text(encoding="utf-8")
    assert "official Proto-Nucleic Acids Building / PNAB program was not found" in text
    assert "--helical-twist-deg" in text
    assert "strand-count 6" in text
