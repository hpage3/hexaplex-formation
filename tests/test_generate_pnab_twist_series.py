import importlib.util
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


generate = load_module("generate_pnab_twist_series", "scripts/generate_pnab_twist_series.py")
analyze = load_module("analyze_pnab_twist_series", "scripts/analyze_pnab_twist_series.py")


def options():
    return {
        "Backbone": {"file_path": "backbone.pdb"},
        "HelicalParameters": {"h_rise": [3.4, 3.4, 1], "h_twist": [30, 30, 1], "tip": [0, 0, 1]},
        "RuntimeParameters": {"strand": "Y" * 15, "strand_orientation": [True, False] * 3},
    }


def test_twist_label_renders_safely():
    assert generate.number_token(29.5) == "29p5"
    assert generate.model_id(29.5, 3.38) == "pnab_twist29p5_rise3p38"


def test_yaml_update_changes_only_rise_and_twist():
    original = options()
    updated = generate.update_helical_parameters(original, 3.38, 29.5)
    assert updated["HelicalParameters"]["h_rise"] == [3.38, 3.38, 1]
    assert updated["HelicalParameters"]["h_twist"] == [29.5, 29.5, 1]
    assert updated["HelicalParameters"]["tip"] == original["HelicalParameters"]["tip"]
    assert updated["RuntimeParameters"] == original["RuntimeParameters"]
    assert original["HelicalParameters"]["h_rise"] == [3.4, 3.4, 1]


def test_results_parser_extracts_conformer_and_distance(tmp_path):
    path = tmp_path / "results.csv"
    path.write_text("# comment\n1,558214,0.099703,1,2,3,4,-130\n", encoding="utf-8")
    parsed = generate.parse_results(path)
    assert parsed["accepted_candidate_count"] == "1"
    assert parsed["conformer_index"] == "558214"
    assert parsed["best_distance_A"] == "0.099703"


def test_teardown_code_is_warning_when_outputs_exist(tmp_path):
    for name in ("fixed.pdb", "results.csv", "prefix.yaml"):
        (tmp_path / name).write_text("x", encoding="utf-8")
    status, error, warning = generate.classify_run(
        generate.TEARDOWN_CODE,
        tmp_path,
        {"conformer_index": "1"},
    )
    assert status == "success_with_warning"
    assert error == ""
    assert "teardown" in warning


def test_failed_run_is_recorded_when_outputs_missing(tmp_path):
    status, error, warning = generate.classify_run(1, tmp_path, {})
    assert status == "failed"
    assert "no accepted" in error
    assert warning == ""


def test_timeout_is_recorded_explicitly(tmp_path):
    status, error, warning = generate.classify_run(124, tmp_path, {})
    assert status == "timed_out"
    assert "timeout" in error
    assert warning == ""


def test_manifest_writer_supports_preserved_resume_rows(tmp_path):
    path = tmp_path / "manifest.csv"
    rows = [
        {field: "" for field in generate.MANIFEST_FIELDS},
        {field: "" for field in generate.MANIFEST_FIELDS},
    ]
    rows[0].update({"model_id": "pnab_twist28_rise3p38", "twist_deg": "28.0"})
    rows[1].update({"model_id": "pnab_twist29_rise3p38", "twist_deg": "29.0"})
    generate.write_manifest(path, rows)
    assert "pnab_twist28_rise3p38" in path.read_text(encoding="utf-8")


def test_local_peak_requires_interior_maximum():
    rows = [
        {"d_A": 3.6, "intensity": 1.0},
        {"d_A": 3.7, "intensity": 3.0},
        {"d_A": 3.8, "intensity": 2.0},
        {"d_A": 3.9, "intensity": 1.0},
    ]
    peak, status = analyze.local_peak(rows, 3.6, 3.9)
    assert status == "local_maximum"
    assert peak["d_A"] == 3.7


def test_parallel_argument_parsing_and_serial_override():
    args = analyze.parse_args(["--workers", "16", "--skip-existing"])
    assert analyze.analysis_mode(args) == ("parallel", 16)
    assert args.skip_existing is True

    serial = analyze.parse_args(["--workers", "16", "--serial"])
    assert analyze.analysis_mode(serial) == ("serial", 1)


def generation_row(model_id="pnab_twist30_rise3p38", structure_path=""):
    return {
        "model_id": model_id,
        "twist_deg": "30.0",
        "rise_A": "3.38",
        "status": "success_with_warning" if structure_path else "timed_out",
        "structure_path": structure_path,
    }


def test_missing_generation_is_skipped_cleanly(tmp_path):
    jobs, rows = analyze.plan_jobs(
        [generation_row()],
        tmp_path,
        "parallel",
        16,
        False,
        None,
    )
    assert jobs == []
    assert rows[0]["analysis_status"] == "skipped_missing_generation"
    assert rows[0]["worker_count"] == "16"
    assert rows[0]["analysis_mode"] == "parallel"
    assert rows[0]["radial_profile_path"] == ""
    assert rows[0]["peak_list_path"] == ""


def test_skip_existing_avoids_recomputation(tmp_path):
    pdb = ROOT / "tests" / "fixtures" / "temporary_parallel_test.pdb"
    pdb.parent.mkdir(parents=True, exist_ok=True)
    pdb.write_text("END\n", encoding="utf-8")
    try:
        item = generation_row(structure_path=str(pdb.relative_to(ROOT)))
        profile, peaks = analyze.output_paths(tmp_path, item["model_id"])
        profile.parent.mkdir(parents=True)
        peaks.parent.mkdir(parents=True)
        profile.write_text("q_Ainv,d_A,intensity,intensity_norm\n", encoding="utf-8")
        peaks.write_text(",".join(analyze.PEAK_FIELDS) + "\n", encoding="utf-8")
        jobs, rows = analyze.plan_jobs([item], tmp_path, "parallel", 16, True, None)
        assert jobs == []
        assert rows[0]["analysis_status"] == "skipped_existing"
    finally:
        pdb.unlink(missing_ok=True)


def failing_worker(_payload):
    raise RuntimeError("synthetic worker failure")


def test_failed_serial_worker_is_recorded_without_stopping_batch(tmp_path):
    payloads = [
        {"item": generation_row("model-a", "tests/test_score_peak_position_fit.py"), "output_root": str(tmp_path)},
        {"item": generation_row("model-b", "tests/test_score_peak_position_fit.py"), "output_root": str(tmp_path)},
    ]
    rows = analyze.run_jobs(payloads, "serial", 1, None, worker_fn=failing_worker)
    assert len(rows) == 2
    assert all(row["analysis_status"] == "failed" for row in rows)
    assert all("synthetic worker failure" in row["error_summary"] for row in rows)


def test_analysis_script_has_windows_multiprocessing_guard():
    source = (ROOT / "scripts" / "analyze_pnab_twist_series.py").read_text(encoding="utf-8")
    assert 'if __name__ == "__main__":' in source
    assert "ProcessPoolExecutor" in source
