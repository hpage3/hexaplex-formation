import csv
import importlib.util
import sys
from argparse import Namespace
from pathlib import Path

import pytest
import yaml


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_script_module(name: str, relative_path: str):
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


pnab_twist = _load_script_module("generate_pnab_twist_variants", "scripts/generate_pnab_twist_variants.py")


def baseline_options() -> dict:
    return {
        "Backbone": {"file_path": "backbone.pdb", "interconnects": [1, 2], "linker": [3, 4]},
        "HelicalParameters": {
            "is_helical": True,
            "x_displacement": [0, 0, 1],
            "y_displacement": [0, 0, 1],
            "h_rise": [3.4, 3.4, 1],
            "inclination": [0, 0, 1],
            "tip": [0, 0, 1],
            "h_twist": [30, 30, 1],
            "twist": [99, 99, 1],
        },
        "RuntimeParameters": {
            "is_hexad": True,
            "strand": "OOOO",
            "build_strand": [True, True, True, True, True, True],
            "strand_orientation": [True, False, True, False, True, False],
        },
    }


def test_set_h_twist_preserves_other_yaml_fields():
    original = baseline_options()

    updated = pnab_twist.set_h_twist(original, 28.0)

    assert updated["HelicalParameters"]["h_twist"] == [28.0, 28.0, 1]
    assert updated["HelicalParameters"]["h_rise"] == [3.4, 3.4, 1]
    assert updated["HelicalParameters"]["twist"] == [99, 99, 1]
    assert updated["Backbone"] == original["Backbone"]
    assert original["HelicalParameters"]["h_twist"] == [30, 30, 1]


def test_parse_twists_and_model_ids():
    assert pnab_twist.parse_twists("28,30,32") == [28.0, 30.0, 32.0]
    assert pnab_twist.parse_twists("", narrow=True) == [27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0]
    assert pnab_twist.model_id_for_twist(30.0) == "pnab_twist_30"
    assert pnab_twist.model_id_for_twist(28.5) == "pnab_twist_28p50"


def test_missing_h_twist_is_rejected():
    options = baseline_options()
    del options["HelicalParameters"]["h_twist"]

    with pytest.raises(ValueError, match="missing h_twist"):
        pnab_twist.set_h_twist(options, 30.0)


def test_dry_run_manifest_has_expected_twist_values(tmp_path):
    baseline_yaml = tmp_path / "baseline.yaml"
    baseline_yaml.write_text(yaml.safe_dump(baseline_options(), sort_keys=False), encoding="utf-8")
    manifest = tmp_path / "manifest.csv"
    args = Namespace(
        baseline_yaml=baseline_yaml,
        out_dir=tmp_path / "pnab_twist_sweep",
        twists="28,30,32",
        narrow_sweep=False,
        dry_run=True,
        python_bin=sys.executable,
        number_of_cpus=1,
        manifest=manifest,
    )

    rows = pnab_twist.generate_variants(args)
    pnab_twist.write_manifest(manifest, rows)

    assert [row["twist_deg"] for row in rows] == ["28.00", "30.00", "32.00"]
    assert all(row["model_status"] == "dry_run_yaml_only" for row in rows)
    assert (tmp_path / "pnab_twist_sweep/yaml/pnab_twist_28.yaml").exists()

    with manifest.open("r", newline="", encoding="utf-8") as handle:
        loaded = list(csv.DictReader(handle))
    assert loaded[1]["model_id"] == "pnab_twist_30"
