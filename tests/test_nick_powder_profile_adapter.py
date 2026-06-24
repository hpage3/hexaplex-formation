import importlib.util
import json
import sys
import types
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


adapter = load_module(
    "run_nick_powderfast10_profile",
    "scripts/run_nick_powderfast10_profile.py",
)


def test_reference_script_imports_without_writing_outputs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    scipy = types.ModuleType("scipy")
    scipy_ndimage = types.ModuleType("scipy.ndimage")
    scipy_ndimage.rotate = lambda array, angle, reshape=False: array
    scipy.ndimage = scipy_ndimage
    monkeypatch.setitem(sys.modules, "scipy", scipy)
    monkeypatch.setitem(sys.modules, "scipy.ndimage", scipy_ndimage)
    reference = load_module("nick_powderfast10_test", "scripts/powderfast10.py")
    assert callable(reference.plot_powder_diffraction)
    assert not (tmp_path / "PowderPattern.png").exists()
    assert not (tmp_path / "middlerow.txt").exists()


def test_argument_parsing_and_output_paths(tmp_path):
    args = adapter.parse_args(
        [
            "--diffraction-npy",
            "pattern.npy",
            "--output-dir",
            str(tmp_path),
            "--z-min",
            "-100",
            "--z-max",
            "100",
            "--x-min",
            "-100",
            "--x-max",
            "100",
            "--max-intensity-scaling",
            "0.05",
            "--run-label",
            "small-test",
        ]
    )
    assert args.diffraction_npy == Path("pattern.npy")
    assert args.run_label == "small-test"
    paths = adapter.output_paths(args.output_dir)
    assert paths["powder_pattern"] == tmp_path / "PowderPattern.png"
    assert paths["middle_row_profile"] == tmp_path / "middlerow.txt"


def test_scipy_availability_detection(monkeypatch):
    monkeypatch.setattr(adapter.importlib.util, "find_spec", lambda name: object())
    assert adapter.scipy_available() is True
    monkeypatch.setattr(adapter.importlib.util, "find_spec", lambda name: None)
    assert adapter.scipy_available() is False


def test_adapter_confines_outputs_and_writes_manifest(tmp_path, monkeypatch):
    caller = tmp_path / "caller"
    output_dir = tmp_path / "controlled-output"
    caller.mkdir()
    monkeypatch.chdir(caller)
    input_path = tmp_path / "tiny.npy"
    np.save(input_path, np.arange(9, dtype=float).reshape(3, 3))
    monkeypatch.setattr(adapter, "scipy_available", lambda: True)

    def fake_plotter(data, z_limits, x_limits, scaling, z_grid_size):
        assert data.shape == (3, 3)
        assert z_limits == (-5.0, 5.0)
        assert x_limits == (-4.0, 4.0)
        assert scaling == 0.2
        assert z_grid_size == 3
        Path("plot2.png").write_bytes(b"fiber")
        Path("PowderPattern.png").write_bytes(b"powder")
        Path("middlerow.txt").write_text("1\n2\n3\n", encoding="utf-8")

    metadata = adapter.run_profile(
        input_path,
        output_dir,
        z_limits=(-5.0, 5.0),
        x_limits=(-4.0, 4.0),
        max_intensity_scaling=0.2,
        run_label="tiny-test",
        plotter=fake_plotter,
        generated_at_utc="2026-06-24T12:00:00+00:00",
    )

    assert not (caller / "PowderPattern.png").exists()
    assert not (caller / "middlerow.txt").exists()
    assert (output_dir / "PowderPattern.png").is_file()
    assert (output_dir / "middlerow.txt").is_file()
    manifest = json.loads(
        (output_dir / "nick_powder_profile_manifest.json").read_text(encoding="utf-8")
    )
    assert manifest == metadata
    assert manifest["status"] == "success"
    assert manifest["error"] == ""
    assert manifest["scipy_available"] is True
    assert manifest["run_label"] == "tiny-test"
    assert manifest["grid_settings"]["z_grid_size"] == 3
    assert manifest["powder_settings"]["rotation_count"] == 360
    assert manifest["reference"]["function"] == "plot_powder_diffraction"


def test_failed_plotter_records_failure_manifest(tmp_path, monkeypatch):
    input_path = tmp_path / "tiny.npy"
    output_dir = tmp_path / "controlled-output"
    np.save(input_path, np.ones((3, 3)))
    monkeypatch.setattr(adapter, "scipy_available", lambda: False)

    def failed_plotter(*args):
        raise RuntimeError("synthetic plotting failure")

    try:
        adapter.run_profile(input_path, output_dir, plotter=failed_plotter)
    except RuntimeError as exc:
        assert "synthetic plotting failure" in str(exc)
    else:
        raise AssertionError("Expected plotting failure")

    manifest = json.loads(
        (output_dir / "nick_powder_profile_manifest.json").read_text(encoding="utf-8")
    )
    assert manifest["status"] == "failed"
    assert manifest["scipy_available"] is False
    assert "synthetic plotting failure" in manifest["error"]


def test_adapter_rejects_non_2d_input(tmp_path):
    input_path = tmp_path / "vector.npy"
    np.save(input_path, np.arange(3))
    try:
        adapter.run_profile(input_path, tmp_path / "out", plotter=lambda *args: None)
    except ValueError as exc:
        assert "two-dimensional" in str(exc)
    else:
        raise AssertionError("Expected a one-dimensional array to be rejected")
