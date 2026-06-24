from pathlib import Path

import pytest

from scripts.score_twist_rise_scan import (
    load_observed_peaks,
    load_targets,
    merge_observed_peaks,
    score_rows,
)


def test_load_targets_reads_combined_target_file(tmp_path):
    target_file = tmp_path / "targets.csv"
    target_file.write_text(
        "target_label,target_d_A,target_group,notes\n"
        "base,3.38,base_stacking,base target\n"
        "A,3.8,backbone_associated,A target\n",
        encoding="utf-8",
    )

    targets = load_targets(target_file)

    assert targets["base"]["target_group"] == "base_stacking"
    assert targets["A"]["target_d_A"] == "3.8"


def test_score_rows_leaves_pending_without_observed_peaks():
    rows = [
        {
            "model_id": "twist30p0_rise3p380",
            "twist_deg": "30.000",
            "rise_A": "3.380",
            "score_status": "pending",
        }
    ]
    targets = {
        "base": {
            "target_label": "base",
            "target_d_A": "3.38",
            "target_group": "base_stacking",
        }
    }

    scored = score_rows(rows, targets)

    assert scored[0]["score_status"] == "pending"
    assert scored[0]["observed_peak_count"] == "0"
    assert scored[0]["expected_peak_count"] == "1"
    assert scored[0]["missing_peak_count"] == "1"
    assert scored[0]["score_completeness"] == "0.000000000"
    assert scored[0]["base_rmsd"] == ""
    assert scored[0]["combined_rmsd"] == ""


def test_load_observed_peaks_keyed_by_model_id(tmp_path):
    observed_file = tmp_path / "observed.csv"
    observed_file.write_text(
        "model_id,observed_base_d_A,observed_A_d_A\n"
        "twist30p0_rise3p380,3.40,3.90\n",
        encoding="utf-8",
    )

    observed = load_observed_peaks(observed_file)

    assert observed["twist30p0_rise3p380"]["observed_base_d_A"] == "3.40"


def test_load_observed_peaks_rejects_duplicates(tmp_path):
    observed_file = tmp_path / "observed.csv"
    observed_file.write_text(
        "model_id,observed_base_d_A\n"
        "twist30p0_rise3p380,3.40\n"
        "twist30p0_rise3p380,3.41\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="Duplicate model_id"):
        load_observed_peaks(observed_file)


def test_merge_observed_peaks_only_adds_matching_model_rows():
    rows = [
        {"model_id": "twist30p0_rise3p380", "score_status": "pending"},
        {"model_id": "twist29p5_rise3p380", "score_status": "pending"},
    ]
    observed = {
        "twist30p0_rise3p380": {
            "model_id": "twist30p0_rise3p380",
            "observed_base_d_A": "3.40",
            "notes": "ignored because it is not an observed_ column",
        }
    }

    merged = merge_observed_peaks(rows, observed)

    assert merged[0]["observed_base_d_A"] == "3.40"
    assert "notes" not in merged[0]
    assert "observed_base_d_A" not in merged[1]


def test_score_rows_computes_base_and_helical_rmsd_from_observed_table():
    rows = [
        {
            "model_id": "twist30p0_rise3p380",
            "twist_deg": "30.000",
            "rise_A": "3.380",
        }
    ]
    observed = {
        "twist30p0_rise3p380": {
            "model_id": "twist30p0_rise3p380",
            "observed_base_d_A": "3.40",
            "observed_A_d_A": "3.90",
            "observed_B_d_A": "4.50",
        }
    }
    targets = {
        "base": {
            "target_label": "base",
            "target_d_A": "3.38",
            "target_group": "base_stacking",
        },
        "A": {
            "target_label": "A",
            "target_d_A": "3.8",
            "target_group": "backbone_associated",
        },
        "B": {
            "target_label": "B",
            "target_d_A": "4.4",
            "target_group": "backbone_associated",
        },
    }

    scored = score_rows(rows, targets, observed)

    assert scored[0]["score_status"] == "scored"
    assert scored[0]["observed_peak_count"] == "3"
    assert scored[0]["expected_peak_count"] == "3"
    assert scored[0]["missing_peak_count"] == "0"
    assert scored[0]["score_completeness"] == "1.000000000"
    assert scored[0]["observed_base_d_A"] == "3.40"
    assert scored[0]["base_rmsd"] == "0.020000000"
    assert scored[0]["helical_rmsd"] == "0.100000000"
    assert scored[0]["combined_rmsd"] != ""

def test_score_rows_tracks_partial_peak_coverage():
    rows = [
        {
            "model_id": "twist30p0_rise3p400",
            "twist_deg": "30.000",
            "rise_A": "3.400",
        }
    ]
    observed = {
        "twist30p0_rise3p400": {
            "model_id": "twist30p0_rise3p400",
            "observed_base_d_A": "3.40",
            "observed_A_d_A": "3.90",
        }
    }
    targets = {
        "base": {
            "target_label": "base",
            "target_d_A": "3.38",
            "target_group": "base_stacking",
        },
        "A": {
            "target_label": "A",
            "target_d_A": "3.8",
            "target_group": "backbone_associated",
        },
        "D": {
            "target_label": "D",
            "target_d_A": "7.3",
            "target_group": "backbone_associated",
        },
    }

    scored = score_rows(rows, targets, observed)

    assert scored[0]["score_status"] == "scored"
    assert scored[0]["observed_peak_count"] == "2"
    assert scored[0]["expected_peak_count"] == "3"
    assert scored[0]["missing_peak_count"] == "1"
    assert scored[0]["score_completeness"] == "0.666666667"
    assert scored[0]["combined_rmsd"] != ""
