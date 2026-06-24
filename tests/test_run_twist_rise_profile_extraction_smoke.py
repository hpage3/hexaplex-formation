import argparse
import csv
from pathlib import Path

from scripts.run_twist_rise_profile_extraction_smoke import run_smoke, write_demo_profiles


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
