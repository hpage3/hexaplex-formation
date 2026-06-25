import argparse
import csv
from pathlib import Path

import pytest

from scripts.generate_twist_rise_pilot_coordinates import (
    build_fields,
    generate_coordinates,
    patch_options_yaml,
    require_columns,
)


def write_csv(path: Path, rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerows(rows)


def test_patch_options_yaml_updates_twist_and_rise(tmp_path):
    options = tmp_path / "options.yaml"
    options.write_text(
        "HelicalParameters:\n"
        "  h_twist: 30.000000\n"
        "  h_rise: 3.400000\n",
        encoding="utf-8",
    )

    patch_options_yaml(options, twist_deg=29.0, rise_A=3.35)

    text = options.read_text(encoding="utf-8")
    assert "h_twist: 29.000000" in text
    assert "h_rise: 3.350000" in text


def test_patch_options_yaml_requires_keys(tmp_path):
    options = tmp_path / "options.yaml"
    options.write_text("Other: true\n", encoding="utf-8")

    with pytest.raises(ValueError, match="h_twist"):
        patch_options_yaml(options, twist_deg=29.0, rise_A=3.35)


def test_require_columns_rejects_missing_values(tmp_path):
    path = tmp_path / "manifest.csv"
    rows = [{"model_id": "m1", "twist_deg": "30"}]

    with pytest.raises(ValueError, match="missing required columns"):
        require_columns(rows, path)


def test_generate_coordinates_dry_run_writes_coordinate_paths(tmp_path):
    manifest = tmp_path / "pilot.csv"
    template = tmp_path / "template"
    work_dir = tmp_path / "work"
    coord_dir = tmp_path / "coords"
    output_manifest = tmp_path / "out.csv"

    write_csv(
        manifest,
        [
            ["model_id", "twist_deg", "rise_A", "model_status", "coordinate_file"],
            ["twist30p0_rise3p400", "30.000", "3.400", "pending", ""],
        ],
    )
    template.mkdir()
    (template / "options.yaml").write_text(
        "HelicalParameters:\n"
        "  h_twist: 30.000000\n"
        "  h_rise: 3.400000\n",
        encoding="utf-8",
    )
    (template / "run.py").write_text("print('stub')\n", encoding="utf-8")

    args = argparse.Namespace(
        manifest=manifest,
        template_dir=template,
        work_dir=work_dir,
        coordinate_dir=coord_dir,
        output_manifest=output_manifest,
        conda_exe=Path("missing-conda.exe"),
        conda_env="pnab",
        timeout_seconds=10,
        max_models=0,
        dry_run=True,
    )

    rows = generate_coordinates(args)

    assert len(rows) == 1
    assert rows[0]["model_status"] == "dry_run"
    assert rows[0]["coordinate_file"].endswith("twist30p0_rise3p400.pdb")


def test_build_fields_includes_pnab_metadata():
    fields = build_fields(
        [
            {
                "model_id": "m1",
                "pnab_work_dir": "work",
                "pnab_returncode": "0",
            }
        ]
    )

    assert "model_id" in fields
    assert "pnab_work_dir" in fields
    assert "pnab_returncode" in fields
