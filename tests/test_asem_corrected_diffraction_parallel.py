import importlib.util
import json
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


runner = load_module(
    "run_asem_corrected_diffraction",
    "scripts/run_asem_corrected_diffraction.py",
)


def tiny_xyz(path):
    path.write_text(
        "4\nsynthetic\n"
        "C 0.0 0.0 0.0\n"
        "N 1.1 -0.4 0.7\n"
        "O -0.8 0.5 -1.2\n"
        "C 0.3 0.9 1.4\n",
        encoding="utf-8",
    )


def arguments(tmp_path, extra=None):
    argv = [
        "--coordinate-file", str(tmp_path / "tiny.xyz"),
        "--output-dir", str(tmp_path / "out"),
        "--grid-size", "9",
        "--theta-count", "2",
        "--phi-count", "2",
        "--psi-count", "1",
    ]
    return runner.parse_args(argv + (extra or []))


def test_worker_parsing_and_serial_override(tmp_path):
    parallel = arguments(tmp_path, ["--workers", "16", "--skip-existing"])
    assert runner.analysis_mode(parallel) == ("parallel", 16)
    assert parallel.skip_existing is True
    serial = arguments(tmp_path, ["--workers", "16", "--serial"])
    assert runner.analysis_mode(serial) == ("serial", 1)


def test_missing_input_is_skipped(tmp_path):
    args = arguments(tmp_path)
    row = runner.run_model(args.coordinate_file[0], args, "serial", 1)
    assert row["status"] == "skipped_missing_input"


def test_skip_existing_avoids_compute(tmp_path):
    path = tmp_path / "tiny.xyz"
    tiny_xyz(path)
    args = arguments(tmp_path, ["--skip-existing"])
    args.output_dir.mkdir()
    for output in runner.output_paths(args.output_dir, path.stem).values():
        output.write_bytes(b"x")

    def should_not_run(*unused):
        raise AssertionError("compute should have been skipped")

    row = runner.run_model(path, args, "serial", 1, compute_fn=should_not_run)
    assert row["status"] == "skipped_existing"


def test_failed_model_is_recorded_without_raising(tmp_path):
    path = tmp_path / "tiny.xyz"
    tiny_xyz(path)
    args = arguments(tmp_path)

    def fail(*unused):
        raise RuntimeError("synthetic worker failure")

    row = runner.run_model(path, args, "parallel", 4, compute_fn=fail)
    assert row["status"] == "failed"
    assert "synthetic worker failure" in row["error"]
    metadata = json.loads(
        runner.output_paths(args.output_dir, path.stem)["metadata"].read_text(encoding="utf-8")
    )
    assert metadata["mode"] == "parallel"
    assert metadata["worker_count"] == 4


def test_serial_parallel_equivalence_on_tiny_fixture(tmp_path):
    path = tmp_path / "tiny.xyz"
    tiny_xyz(path)
    args = arguments(tmp_path)
    atomic_numbers, coords = runner.load_xyz(path)
    serial = runner.compute_image(atomic_numbers, coords * 1e-7, args, "serial", 1)
    parallel = runner.compute_image(atomic_numbers, coords * 1e-7, args, "parallel", 2)
    assert serial.shape == parallel.shape == (9, 9)
    assert np.allclose(serial, parallel, rtol=1e-12, atol=1e-9)


def test_success_metadata_records_mode_workers_and_grid(tmp_path):
    path = tmp_path / "tiny.xyz"
    tiny_xyz(path)
    args = arguments(tmp_path)

    def fake_compute(*unused):
        return np.ones((9, 9), dtype=float)

    row = runner.run_model(path, args, "parallel", 16, compute_fn=fake_compute)
    assert row["status"] == "success"
    metadata = json.loads(
        runner.output_paths(args.output_dir, path.stem)["metadata"].read_text(encoding="utf-8")
    )
    assert metadata["mode"] == "parallel"
    assert metadata["worker_count"] == 16
    assert metadata["effective_worker_count"] == 4
    assert metadata["grid_size"] == 9


def test_script_has_windows_multiprocessing_guard():
    source = (ROOT / "scripts" / "run_asem_corrected_diffraction.py").read_text(encoding="utf-8")
    assert 'if __name__ == "__main__":' in source
    assert "average_powder_diffraction_parallel" in source
