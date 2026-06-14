import csv
import importlib.util
import json
import math
from pathlib import Path

import pytest

from hexaplex_formation.pdb_utils import PDBAtom
from hexaplex_formation.scattering import (
    debye_intensity,
    debye_intensity_from_distance_histogram,
    d_from_q,
    pair_distance_histogram_for_debye,
    q_from_d,
    stratified_sample_atoms,
)


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_script_module(name: str, relative_path: str):
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


score_windows = _load_script_module("score_radial_windows", "scripts/score_radial_windows.py")
compare_diffraction = _load_script_module("compare_ladder_diffraction", "scripts/compare_ladder_diffraction.py")
debye_profile = _load_script_module("debye_radial_profile", "scripts/debye_radial_profile.py")


def atom(serial: int, x: float, y: float, z: float) -> PDBAtom:
    return PDBAtom(
        record_type="ATOM",
        atom_serial=serial,
        atom_name="CA",
        alt_loc="",
        residue_name="GLU",
        chain_id="A",
        residue_number=serial,
        insertion_code="",
        x=x,
        y=y,
        z=z,
        occupancy=1.0,
        temp_factor=0.0,
        element="C",
    )


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def test_q_d_round_trip():
    d_value = 4.5
    q_value = q_from_d(d_value)

    assert d_from_q(q_value) == pytest.approx(d_value)


def test_debye_intensity_two_atoms_is_finite():
    atoms = [atom(1, 0.0, 0.0, 0.0), atom(2, 1.0, 0.0, 0.0)]

    intensities = debye_intensity(atoms, [0.0, math.pi])

    assert intensities[0] == pytest.approx(4.0)
    assert math.isfinite(intensities[1])


def test_histogram_debye_matches_direct_approximately_for_small_atom_set():
    atoms = [
        atom(1, 0.0, 0.0, 0.0),
        atom(2, 1.0, 0.0, 0.0),
        atom(3, 0.0, 2.0, 0.0),
    ]
    q_values = [0.2, 1.0, 2.0]

    direct = debye_intensity(atoms, q_values)
    histogram = pair_distance_histogram_for_debye(atoms, bin_width=0.001)
    from_histogram = debye_intensity_from_distance_histogram(histogram, q_values)

    assert from_histogram == pytest.approx(direct, abs=1e-3)


def test_histogram_self_term_and_pair_double_counting():
    atoms = [atom(1, 0.0, 0.0, 0.0), atom(2, 1.0, 0.0, 0.0)]
    histogram = pair_distance_histogram_for_debye(atoms, bin_width=0.05, element_weights={"C": 2.0})

    assert histogram["self_term"] == pytest.approx(8.0)
    assert sum(histogram["bins"].values()) == pytest.approx(8.0)
    assert debye_intensity_from_distance_histogram(histogram, [0.0])[0] == pytest.approx(16.0)


def test_stratified_sampling_does_not_take_first_n_atoms():
    atoms = [atom(serial, float(serial), 0.0, 0.0) for serial in range(1, 21)]

    sampled = stratified_sample_atoms(atoms, 6)

    assert len(sampled) == 6
    assert sampled != atoms[:6]
    assert sampled[0].atom_serial != 1
    assert sampled[-1].atom_serial > 6


def test_window_scoring_on_tiny_profile():
    profile_rows = [
        {"q_Ainv": "1.40", "d_A": "4.49", "intensity": "10.0"},
        {"q_Ainv": "1.41", "d_A": "4.45", "intensity": "20.0"},
        {"q_Ainv": "2.00", "d_A": "3.14", "intensity": "5.0"},
    ]

    scores = {row["window_name"]: row for row in score_windows.score_windows(profile_rows)}

    assert scores["d_4p5"]["point_count"] == "2"
    assert scores["d_4p5"]["mean_intensity"] == "15.000000"
    assert float(scores["d_4p5"]["integrated_intensity_fraction_of_total"]) == pytest.approx(30.0 / 35.0, abs=1e-6)


def test_ladder_diffraction_aggregation_from_fixture_csvs(tmp_path):
    summary = tmp_path / "intermediate_ladder_summary.csv"
    window_dir = tmp_path / "window_scores"
    _write_csv(
        summary,
        ["model_name", "included_blocks", "includes_hexads", "atom_mode", "atom_count", "residue_count"],
        [
            {
                "model_name": "scaffold_blocks_1_heavy_deduped",
                "included_blocks": "1",
                "includes_hexads": "no",
                "atom_mode": "heavy_deduped",
                "atom_count": "10",
                "residue_count": "2",
            }
        ],
    )
    _write_csv(
        window_dir / "scaffold_blocks_1_heavy_deduped_window_scores.csv",
        [
            "window_name",
            "mean_intensity",
            "integrated_intensity_fraction_of_total",
        ],
        [
            {"window_name": "d_4p5", "mean_intensity": "12.0", "integrated_intensity_fraction_of_total": "0.2"},
            {"window_name": "d_3p4", "mean_intensity": "4.0", "integrated_intensity_fraction_of_total": "0.1"},
        ],
    )

    rows = compare_diffraction.aggregate_scores(window_dir, summary)

    assert len(rows) == 1
    assert rows[0]["model_name"] == "scaffold_blocks_1_heavy_deduped"
    assert rows[0]["d_4p5_mean"] == "12.0"
    assert rows[0]["d_4p5_fraction"] == "0.2"
    assert rows[0]["d_3p4_fraction"] == "0.1"


def test_ladder_diffraction_report_includes_caution_language(tmp_path):
    out = tmp_path / "report.md"
    rows = [
        {
            "model_name": "scaffold_blocks_1_heavy_deduped",
            "included_blocks": "1",
            "includes_hexads": "no",
            "atom_mode": "heavy_deduped",
            "d_4p5_fraction": "0.1",
            "d_3p4_fraction": "0.2",
        }
    ]

    compare_diffraction.write_markdown_report(rows, out, tmp_path / "missing.csv", tmp_path)

    text = out.read_text(encoding="utf-8")
    assert "simplified Debye-style isotropic scattering approximation" in text
    assert "not a replacement for full fiber-diffraction simulations" in text
    assert "not a literal 4.5 A atom-contact assignment" in text


def test_debye_profile_metadata_json_writing(tmp_path):
    out = tmp_path / "profile.csv"
    metadata = {
        "pdb": "fixture.pdb",
        "atom_count_input": 3,
        "atom_count_used": 2,
        "method": "histogram",
        "distance_bin_width": 0.05,
        "sample_atoms": 2,
        "max_atoms": None,
        "q_min": 0.2,
        "q_max": 0.4,
        "q_step": 0.2,
        "caution": "fixture caution",
    }

    debye_profile.write_profile(out, [0.2, 0.4], [1.0, 2.0], metadata)

    assert out.read_text(encoding="utf-8").splitlines()[0] == "q_Ainv,d_A,intensity"
    written = json.loads((tmp_path / "profile.metadata.json").read_text(encoding="utf-8"))
    assert written["method"] == "histogram"
    assert written["atom_count_used"] == 2


def test_metadata_notes_distinguish_sampling_and_first_n(tmp_path):
    profile_dir = tmp_path / "profiles"
    profile_dir.mkdir()
    (profile_dir / "sample_debye_profile.metadata.json").write_text(
        json.dumps(
            {
                "method": "histogram",
                "distance_bin_width": 0.05,
                "sample_atoms": 800,
                "sampling_used": True,
                "max_atoms": None,
                "first_n_truncated": False,
            }
        ),
        encoding="utf-8",
    )
    (profile_dir / "legacy_debye_profile.metadata.json").write_text(
        json.dumps(
            {
                "method": "direct",
                "distance_bin_width": 0.05,
                "sample_atoms": None,
                "sampling_used": False,
                "max_atoms": 500,
                "first_n_truncated": True,
            }
        ),
        encoding="utf-8",
    )

    notes = "\n".join(compare_diffraction._metadata_notes(profile_dir))

    assert "deterministic stratified sampling" in notes
    assert "STRONG WARNING" in notes
    assert "legacy first-N atom truncation" in notes
