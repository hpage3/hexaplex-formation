import csv
from pathlib import Path

from scripts.score_twist_rise_scan import load_targets, score_rows


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
    assert scored[0]["base_rmsd"] == ""
    assert scored[0]["combined_rmsd"] == ""


def test_score_rows_computes_base_and_helical_rmsd():
    rows = [
        {
            "model_id": "twist30p0_rise3p380",
            "twist_deg": "30.000",
            "rise_A": "3.380",
            "observed_base_d_A": "3.40",
            "observed_A_d_A": "3.90",
            "observed_B_d_A": "4.50",
        }
    ]
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

    scored = score_rows(rows, targets)

    assert scored[0]["score_status"] == "scored"
    assert scored[0]["base_rmsd"] == "0.020000000"
    assert scored[0]["helical_rmsd"] == "0.100000000"
    assert scored[0]["combined_rmsd"] != ""
