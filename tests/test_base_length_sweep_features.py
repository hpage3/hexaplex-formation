import csv
import importlib.util
import math
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_script_module(name: str, relative_path: str):
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


sweep = _load_script_module("analyze_base_length_sweep_features", "scripts/analyze_base_length_sweep_features.py")


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def test_q_d_conversion_uses_two_pi_convention():
    q_value = sweep.q_from_d(4.5)

    assert q_value == pytest.approx(2.0 * math.pi / 4.5)


def test_feature_window_extraction_and_integrated_intensity():
    points = [
        {"q_Ainv": 1.20, "d_A": 5.2, "intensity": 2.0},
        {"q_Ainv": 1.30, "d_A": 4.85, "intensity": 5.0},
        {"q_Ainv": 1.34, "d_A": 4.70, "intensity": 7.0},
        {"q_Ainv": 1.50, "d_A": 4.0, "intensity": 11.0},
    ]

    stats = sweep.window_stats(points, 4.5, 5.0)

    assert stats["point_count"] == 2
    assert stats["max_intensity"] == 7.0
    assert stats["mean_intensity"] == 6.0
    assert stats["integrated_intensity"] == 12.0
    assert stats["d_at_max"] == 4.70


def test_local_maximum_detector_distinguishes_interior_from_boundary():
    interior_peak = [
        {"q_Ainv": 1.10, "d_A": 5.4, "intensity": 10.0},
        {"q_Ainv": 1.25, "d_A": 5.0, "intensity": 12.0},
        {"q_Ainv": 1.32, "d_A": 4.75, "intensity": 20.0},
        {"q_Ainv": 1.39, "d_A": 4.55, "intensity": 12.0},
        {"q_Ainv": 1.55, "d_A": 4.0, "intensity": 9.0},
    ]
    boundary_peak = [
        {"q_Ainv": 1.10, "d_A": 5.4, "intensity": 10.0},
        {"q_Ainv": 1.25, "d_A": 5.0, "intensity": 20.0},
        {"q_Ainv": 1.32, "d_A": 4.75, "intensity": 12.0},
        {"q_Ainv": 1.39, "d_A": 4.55, "intensity": 11.0},
        {"q_Ainv": 1.55, "d_A": 4.0, "intensity": 9.0},
    ]

    assert sweep.local_maximum_note(interior_peak, 4.5, 5.0)[0] is True
    has_peak, note = sweep.local_maximum_note(boundary_peak, 4.5, 5.0)
    assert has_peak is False
    assert "boundary" in note


def test_feature_summary_handles_missing_or_empty_window(tmp_path):
    manifest = tmp_path / "manifest.csv"
    geometry = tmp_path / "geometry.csv"
    profile_dir = tmp_path / "profiles"
    profile_dir.mkdir()
    variant_id = "hexaplex_base_length_scale_1p00"
    _write_csv(
        manifest,
        ["variant_id", "scale_factor"],
        [{"variant_id": variant_id, "scale_factor": "1.00"}],
    )
    _write_csv(
        geometry,
        ["variant_id", "warnings"],
        [{"variant_id": variant_id, "warnings": ""}],
    )
    _write_csv(
        profile_dir / f"{variant_id}_radial.csv",
        ["q_Ainv", "d_A", "mean_intensity", "pixel_count"],
        [
            {"q_Ainv": "0.5", "d_A": "12.5", "mean_intensity": "3.0", "pixel_count": "4"},
            {"q_Ainv": "0.6", "d_A": "10.5", "mean_intensity": "4.0", "pixel_count": "4"},
        ],
    )

    rows = sweep.build_feature_rows(profile_dir, manifest, geometry, [1.0])

    assert len(rows) == 1
    assert rows[0]["variant_id"] == variant_id
    assert rows[0]["integrated_intensity_4p5_5A"] == ""
    assert "empty" in rows[0]["interpretation_note"]


def test_report_includes_sensitivity_caution_and_limitations(tmp_path):
    report = tmp_path / "report.md"
    rows = [
        {
            "variant_id": "hexaplex_base_length_scale_1p00",
            "scale_factor": "1.00",
            "geometry_warning": "",
            "d_A_at_max_in_4p5_5A_window": "4.75",
            "integrated_intensity_4p5_5A": "10.0",
            "has_local_maximum_4p5_5A": "no",
            "integrated_intensity_3p4A": "5.0",
            "integrated_intensity_3p0A": "4.0",
        }
    ]
    manifest = tmp_path / "manifest.csv"
    geometry = tmp_path / "geometry.csv"
    _write_csv(
        manifest,
        ["variant_id", "transformed_atom_count", "fixed_atom_count"],
        [{"variant_id": "hexaplex_base_length_scale_1p00", "transformed_atom_count": "2", "fixed_atom_count": "3"}],
    )
    _write_csv(
        geometry,
        ["variant_id", "warnings"],
        [{"variant_id": "hexaplex_base_length_scale_1p00", "warnings": ""}],
    )

    sweep.write_report(rows, report, [], manifest, geometry, "baseline.pdb")

    text = report.read_text(encoding="utf-8")
    assert "computational sensitivity study only" in text
    assert "q = 2*pi/d" in text
    assert "Scale 1.20 is excluded" in text
    assert "does not determine the structure" in text
