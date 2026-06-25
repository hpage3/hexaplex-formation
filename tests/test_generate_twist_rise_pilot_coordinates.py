import argparse
import csv
import subprocess
import time
from pathlib import Path

import pytest

from scripts.generate_twist_rise_pilot_coordinates import (
    build_parser,
    build_fields,
    finish_completed_run,
    generate_coordinates,
    patch_options_yaml,
    process_selected_rows,
    require_columns,
    summarize_statuses,
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
        "  h_twist:\n"
        "  - 30\n"
        "  - 30\n"
        "  - 1\n"
        "  h_rise:\n"
        "  - 3.40\n"
        "  - 3.40\n"
        "  - 1\n",
        encoding="utf-8",
    )

    patch_options_yaml(options, twist_deg=29.0, rise_A=3.35)

    text = options.read_text(encoding="utf-8")
    assert "  h_twist:\n  - 29.000000\n  - 29.000000\n  - 1" in text
    assert "  h_rise:\n  - 3.350000\n  - 3.350000\n  - 1" in text


def test_patch_options_yaml_requires_keys(tmp_path):
    options = tmp_path / "options.yaml"
    options.write_text("Other: true\n", encoding="utf-8")

    with pytest.raises(ValueError, match="h_rise|h_twist"):
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
        workers=1,
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


def completed(returncode: int) -> subprocess.CompletedProcess[str]:
    return subprocess.CompletedProcess(
        args=["conda", "run", "python", "run.py"],
        returncode=returncode,
        stdout="stdout detail",
        stderr="stderr detail",
    )


def finish_run(tmp_path: Path, returncode: int, fixed_exists: bool) -> dict[str, str]:
    fixed_pdb = tmp_path / "work" / "fixed.pdb"
    coordinate_path = tmp_path / "coordinates" / "model.pdb"
    fixed_pdb.parent.mkdir(parents=True)
    coordinate_path.parent.mkdir(parents=True)
    if fixed_exists:
        fixed_pdb.write_text("MODEL\nEND\n", encoding="utf-8")

    out = {"started_at": "2026-06-25T00:00:00Z"}
    return finish_completed_run(
        out=out,
        result=completed(returncode),
        fixed_pdb=fixed_pdb,
        coordinate_path=coordinate_path,
        start_time=time.monotonic(),
    )


def test_nonzero_return_with_fixed_pdb_is_salvaged_cautiously(tmp_path):
    row = finish_run(tmp_path, returncode=3221225477, fixed_exists=True)

    assert row["model_status"] == "generated_with_pnab_nonzero"
    assert row["pnab_returncode"] == "3221225477"
    assert row["pnab_stdout_tail"] == "stdout detail"
    assert row["pnab_stderr_tail"] == "stderr detail"
    assert row["started_at"] == "2026-06-25T00:00:00Z"
    assert row["finished_at"].endswith("Z")
    assert float(row["elapsed_seconds"]) >= 0
    assert "fixed.pdb produced despite nonzero pNAB exit code" in row["notes"]
    assert (tmp_path / "coordinates" / "model.pdb").read_text(encoding="utf-8") == "MODEL\nEND\n"


def test_zero_return_with_fixed_pdb_remains_generated(tmp_path):
    row = finish_run(tmp_path, returncode=0, fixed_exists=True)

    assert row["model_status"] == "generated"
    assert row["notes"] == "pNAB coordinate generated"
    assert (tmp_path / "coordinates" / "model.pdb").is_file()


def test_nonzero_return_without_fixed_pdb_remains_failed(tmp_path):
    row = finish_run(tmp_path, returncode=1, fixed_exists=False)

    assert row["model_status"] == "failed"
    assert row["notes"] == "pNAB returned non-zero exit code and fixed.pdb was not found"
    assert not (tmp_path / "coordinates" / "model.pdb").exists()


def test_zero_return_without_fixed_pdb_remains_missing(tmp_path):
    row = finish_run(tmp_path, returncode=0, fixed_exists=False)

    assert row["model_status"] == "missing_fixed_pdb"
    assert row["notes"] == "pNAB completed but fixed.pdb was not found"
    assert not (tmp_path / "coordinates" / "model.pdb").exists()


def test_summary_counts_coordinate_present_without_calling_salvage_clean():
    rows = [
        {"model_status": "generated"},
        *[{"model_status": "generated_with_pnab_nonzero"} for _ in range(19)],
        {"model_status": "failed"},
        {"model_status": "dry_run"},
    ]

    assert summarize_statuses(rows) == {
        "coordinate_present": 20,
        "clean_generated": 1,
        "generated_with_pnab_nonzero": 19,
        "dry_run": 1,
        "failed_or_other": 1,
    }


def test_workers_default_is_one():
    args = build_parser().parse_args([])
    assert args.workers == 1


def parallel_args(tmp_path: Path, workers: int) -> argparse.Namespace:
    return argparse.Namespace(
        template_dir=tmp_path / "template",
        work_dir=tmp_path / "work",
        coordinate_dir=tmp_path / "coordinates",
        conda_exe=Path("missing-conda.exe"),
        conda_env="pnab",
        timeout_seconds=10,
        dry_run=True,
        workers=workers,
    )


def test_parallel_path_preserves_manifest_row_order(tmp_path):
    rows = [
        {"model_id": "first", "twist_deg": "28", "rise_A": "3.35"},
        {"model_id": "second", "twist_deg": "29", "rise_A": "3.40"},
        {"model_id": "third", "twist_deg": "30", "rise_A": "3.45"},
    ]
    delays = {"first": 0.03, "second": 0.02, "third": 0.01}

    def out_of_order_process(*, row, **unused):
        time.sleep(delays[row["model_id"]])
        return {
            **row,
            "model_status": "dry_run",
            "elapsed_seconds": f"{delays[row['model_id']]:.3f}",
        }

    output = process_selected_rows(
        rows,
        parallel_args(tmp_path, workers=3),
        process_fn=out_of_order_process,
    )

    assert [row["model_id"] for row in output] == ["first", "second", "third"]


def test_parallel_dry_run_produces_expected_rows(tmp_path):
    rows = [
        {"model_id": "m1", "twist_deg": "28", "rise_A": "3.35"},
        {"model_id": "m2", "twist_deg": "30", "rise_A": "3.40"},
    ]

    output = process_selected_rows(rows, parallel_args(tmp_path, workers=2))

    assert [row["model_id"] for row in output] == ["m1", "m2"]
    assert all(row["model_status"] == "dry_run" for row in output)
    assert output[0]["coordinate_file"].endswith("m1.pdb")
    assert output[1]["coordinate_file"].endswith("m2.pdb")


def test_parallel_worker_exception_does_not_abort_other_rows(tmp_path):
    rows = [
        {"model_id": "bad", "twist_deg": "28", "rise_A": "3.35"},
        {"model_id": "good", "twist_deg": "30", "rise_A": "3.40"},
    ]

    def sometimes_fails(*, row, **unused):
        if row["model_id"] == "bad":
            raise RuntimeError("synthetic worker error")
        return {**row, "model_status": "dry_run", "elapsed_seconds": "0.001"}

    output = process_selected_rows(
        rows,
        parallel_args(tmp_path, workers=2),
        process_fn=sometimes_fails,
    )

    assert [row["model_id"] for row in output] == ["bad", "good"]
    assert output[0]["model_status"] == "failed"
    assert "synthetic worker error" in output[0]["notes"]
    assert output[1]["model_status"] == "dry_run"
