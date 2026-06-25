import argparse
import csv
from pathlib import Path

from scripts.run_twist_rise_profile_extraction_smoke import (
    apply_model_id_map,
    load_model_id_map,
    run_smoke,
    write_demo_profiles,
)


def write_csv(path: Path, rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerows(rows)


def test_write_demo_profiles_creates_profile_csv(tmp_path):
    profile_path = tmp_path / "profiles.csv"

    write_demo_profiles(profile_path)

    text = profile_path.read_text(encoding="utf-8")
    assert "model_id,d_A,intensity" in text
    assert "twist30p0_rise3p400" in text


def test_run_smoke_writes_pipeline_outputs(tmp_path):
    profile_path = tmp_path / "profiles.csv"
    target_path = tmp_path / "targets.csv"
    manifest_path = tmp_path / "grid_manifest.csv"
    out_dir = tmp_path / "smoke"

    write_csv(
        target_path,
        [
            ["target_label", "target_d_A", "target_group", "notes"],
            ["base", "3.38", "base_stacking", "base target"],
            ["A", "3.8", "backbone_associated", "A target"],
        ],
    )
    write_csv(
        manifest_path,
        [
            [
                "model_id",
                "twist_deg",
                "rise_A",
                "model_status",
                "coordinate_file",
                "diffraction_status",
                "diffraction_file",
                "score_status",
                "base_rmsd",
                "helical_rmsd",
                "combined_rmsd",
                "notes",
            ],
            [
                "twist30p0_rise3p400",
                "30.000",
                "3.400",
                "pending",
                "",
                "pending",
                "",
                "pending",
                "",
                "",
                "",
                "",
            ],
        ],
    )

    args = argparse.Namespace(
        profiles=profile_path,
        write_demo_profiles=True,
        manifest=manifest_path,
        targets=target_path,
        output_dir=out_dir,
        window_half_width=0.25,
        include_missing=True,
        target_order="base,A",
        model_id_map=None,
    )

    payload = run_smoke(args)

    assert payload["status"] == "success"
    assert payload["counts"]["models_in_profile"] == 3
    assert payload["counts"]["targets"] == 2
    assert (out_dir / "observed_peaks_long.csv").exists()
    assert (out_dir / "observed_peaks_wide.csv").exists()
    assert (out_dir / "twist_rise_scored_manifest.csv").exists()
    assert (out_dir / "smoke_manifest.json").exists()

    scored_text = (out_dir / "twist_rise_scored_manifest.csv").read_text(encoding="utf-8")
    assert "score_completeness" in scored_text
    assert "twist30p0_rise3p400" in scored_text


def test_load_model_id_map_reads_source_to_target_mapping(tmp_path):
    mapping_path = tmp_path / "model_id_map.csv"
    write_csv(
        mapping_path,
        [
            ["source_model_id", "target_model_id"],
            ["compact_hexaplex_twist_30", "twist30p0_rise3p400"],
        ],
    )

    mapping = load_model_id_map(mapping_path)

    assert mapping == {
        "compact_hexaplex_twist_30": "twist30p0_rise3p400",
    }


def test_apply_model_id_map_preserves_source_model_id():
    rows = [
        {
            "model_id": "compact_hexaplex_twist_30",
            "target_label": "base",
            "observed_d_A": "3.407",
        },
        {
            "model_id": "compact_hexaplex_twist_28",
            "target_label": "base",
            "observed_d_A": "3.396",
        },
    ]
    mapping = {
        "compact_hexaplex_twist_30": "twist30p0_rise3p400",
    }

    mapped = apply_model_id_map(rows, mapping)

    assert mapped[0]["model_id"] == "twist30p0_rise3p400"
    assert mapped[0]["source_model_id"] == "compact_hexaplex_twist_30"
    assert mapped[1]["model_id"] == "compact_hexaplex_twist_28"
    assert "source_model_id" not in mapped[1]


def test_run_smoke_with_model_id_map_scores_mapped_row(tmp_path):
    profile_path = tmp_path / "profiles.csv"
    target_path = tmp_path / "targets.csv"
    manifest_path = tmp_path / "grid_manifest.csv"
    map_path = tmp_path / "model_id_map.csv"
    out_dir = tmp_path / "smoke"

    write_csv(
        profile_path,
        [
            ["model_id", "d_A", "intensity"],
            ["compact_hexaplex_twist_30", "3.407", "100"],
            ["compact_hexaplex_twist_30", "3.853", "90"],
        ],
    )
    write_csv(
        target_path,
        [
            ["target_label", "target_d_A", "target_group", "notes"],
            ["base", "3.38", "base_stacking", "base target"],
            ["A", "3.8", "backbone_associated", "A target"],
        ],
    )
    write_csv(
        manifest_path,
        [
            [
                "model_id",
                "twist_deg",
                "rise_A",
                "model_status",
                "coordinate_file",
                "diffraction_status",
                "diffraction_file",
                "score_status",
                "base_rmsd",
                "helical_rmsd",
                "combined_rmsd",
                "notes",
            ],
            [
                "twist30p0_rise3p400",
                "30.000",
                "3.400",
                "pending",
                "",
                "pending",
                "",
                "pending",
                "",
                "",
                "",
                "",
            ],
        ],
    )
    write_csv(
        map_path,
        [
            ["source_model_id", "target_model_id"],
            ["compact_hexaplex_twist_30", "twist30p0_rise3p400"],
        ],
    )

    args = argparse.Namespace(
        profiles=profile_path,
        write_demo_profiles=False,
        manifest=manifest_path,
        targets=target_path,
        output_dir=out_dir,
        window_half_width=0.25,
        include_missing=True,
        target_order="base,A",
        model_id_map=map_path,
    )

    payload = run_smoke(args)

    assert payload["status"] == "success"
    assert payload["mapped_model_count"] == 1
    assert payload["counts"]["scored_rows"] == 1

    scored_text = (out_dir / "twist_rise_scored_manifest.csv").read_text(encoding="utf-8")
    assert "twist30p0_rise3p400" in scored_text
    assert "scored" in scored_text
