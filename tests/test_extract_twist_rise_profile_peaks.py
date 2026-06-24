import pytest

from scripts.extract_twist_rise_profile_peaks import (
    choose_peak,
    extract_peaks,
    load_profiles,
    load_targets,
)


def test_choose_peak_selects_max_intensity_inside_window():
    profile_rows = [
        {"d_A": "3.30", "intensity": "10"},
        {"d_A": "3.38", "intensity": "40"},
        {"d_A": "3.41", "intensity": "35"},
        {"d_A": "3.80", "intensity": "100"},
    ]

    peak = choose_peak(profile_rows, target_d_A=3.38, window_half_width=0.10)

    assert peak is not None
    assert peak["d_A"] == "3.38"
    assert peak["intensity"] == "40"


def test_choose_peak_returns_none_when_window_empty():
    profile_rows = [
        {"d_A": "4.00", "intensity": "10"},
        {"d_A": "4.10", "intensity": "40"},
    ]

    peak = choose_peak(profile_rows, target_d_A=3.38, window_half_width=0.10)

    assert peak is None


def test_extract_peaks_outputs_long_form_rows():
    profiles_by_model = {
        "twist30p0_rise3p400": [
            {"d_A": "3.39", "intensity": "10"},
            {"d_A": "3.41", "intensity": "50"},
            {"d_A": "3.83", "intensity": "20"},
        ]
    }
    targets = [
        {"target_label": "base", "target_d_A": "3.38"},
        {"target_label": "A", "target_d_A": "3.80"},
    ]

    rows = extract_peaks(
        profiles_by_model=profiles_by_model,
        targets=targets,
        window_half_width=0.10,
    )

    assert len(rows) == 2
    assert rows[0]["model_id"] == "twist30p0_rise3p400"
    assert rows[0]["target_label"] == "base"
    assert rows[0]["observed_d_A"] == "3.410000000"
    assert rows[0]["observed_intensity"] == "50.000000000"
    assert rows[0]["peak_status"] == "found"

    assert rows[1]["target_label"] == "A"
    assert rows[1]["observed_d_A"] == "3.830000000"


def test_extract_peaks_can_include_missing_rows():
    profiles_by_model = {
        "twist30p0_rise3p400": [
            {"d_A": "5.00", "intensity": "10"},
        ]
    }
    targets = [
        {"target_label": "base", "target_d_A": "3.38"},
    ]

    rows = extract_peaks(
        profiles_by_model=profiles_by_model,
        targets=targets,
        window_half_width=0.10,
        include_missing=True,
    )

    assert len(rows) == 1
    assert rows[0]["target_label"] == "base"
    assert rows[0]["observed_d_A"] == ""
    assert rows[0]["peak_status"] == "missing"


def test_load_targets_rejects_duplicate_labels(tmp_path):
    targets = tmp_path / "targets.csv"
    targets.write_text(
        "target_label,target_d_A,target_group,notes\n"
        "base,3.38,base_stacking,one\n"
        "base,3.39,base_stacking,two\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="Duplicate target_label"):
        load_targets(targets)


def test_load_profiles_groups_by_model_id(tmp_path):
    profiles = tmp_path / "profiles.csv"
    profiles.write_text(
        "model_id,d_A,intensity\n"
        "m1,3.38,10\n"
        "m1,3.40,20\n"
        "m2,3.38,30\n",
        encoding="utf-8",
    )

    grouped = load_profiles(profiles)

    assert sorted(grouped) == ["m1", "m2"]
    assert len(grouped["m1"]) == 2
    assert grouped["m1"][0]["d_A"] == "3.380000000"


def test_load_profiles_rejects_missing_columns(tmp_path):
    profiles = tmp_path / "bad_profiles.csv"
    profiles.write_text(
        "model_id,d_A\n"
        "m1,3.38\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="missing required columns"):
        load_profiles(profiles)
