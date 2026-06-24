import pytest

from scripts.build_twist_rise_observed_peaks_table import (
    long_to_wide,
    parse_target_order,
    validate_long_rows,
)


def test_parse_target_order_default():
    assert parse_target_order(None) == ["base", "A", "B", "C", "D"]


def test_parse_target_order_rejects_duplicates():
    with pytest.raises(ValueError, match="duplicates"):
        parse_target_order("base,A,A")


def test_long_to_wide_converts_peak_assignments():
    rows = [
        {"model_id": "twist30p0_rise3p400", "target_label": "base", "observed_d_A": "3.407"},
        {"model_id": "twist30p0_rise3p400", "target_label": "A", "observed_d_A": "3.853"},
        {"model_id": "twist29p0_rise3p400", "target_label": "base", "observed_d_A": "3.396"},
        {"model_id": "twist29p0_rise3p400", "target_label": "D", "observed_d_A": "7.465"},
    ]

    wide = long_to_wide(rows, ["base", "A", "B", "C", "D"])

    assert wide[0]["model_id"] == "twist29p0_rise3p400"
    assert wide[0]["observed_base_d_A"] == "3.396"
    assert wide[0]["observed_D_d_A"] == "7.465"
    assert wide[0]["observed_A_d_A"] == ""

    assert wide[1]["model_id"] == "twist30p0_rise3p400"
    assert wide[1]["observed_base_d_A"] == "3.407"
    assert wide[1]["observed_A_d_A"] == "3.853"


def test_long_to_wide_rejects_duplicate_model_target_pair():
    rows = [
        {"model_id": "twist30p0_rise3p400", "target_label": "base", "observed_d_A": "3.407"},
        {"model_id": "twist30p0_rise3p400", "target_label": "base", "observed_d_A": "3.408"},
    ]

    with pytest.raises(ValueError, match="Duplicate observed peak"):
        long_to_wide(rows, ["base"])


def test_validate_long_rows_rejects_missing_required_columns(tmp_path):
    source = tmp_path / "bad.csv"
    rows = [{"model_id": "m1", "target_label": "base"}]

    with pytest.raises(ValueError, match="missing required columns"):
        validate_long_rows(rows, source)


def test_validate_long_rows_rejects_non_numeric_observed_d(tmp_path):
    source = tmp_path / "bad.csv"
    rows = [{"model_id": "m1", "target_label": "base", "observed_d_A": "not-a-number"}]

    with pytest.raises(ValueError):
        validate_long_rows(rows, source)

def test_long_to_wide_skips_missing_peak_rows():
    rows = [
        {
            "model_id": "twist30p0_rise3p400",
            "target_label": "base",
            "observed_d_A": "3.407",
            "peak_status": "found",
        },
        {
            "model_id": "twist30p0_rise3p400",
            "target_label": "D",
            "observed_d_A": "",
            "peak_status": "missing",
        },
    ]

    wide = long_to_wide(rows, ["base", "A", "B", "C", "D"])

    assert len(wide) == 1
    assert wide[0]["observed_base_d_A"] == "3.407"
    assert wide[0]["observed_D_d_A"] == ""


def test_validate_long_rows_allows_missing_peak_rows(tmp_path):
    source = tmp_path / "observed_long.csv"
    rows = [
        {
            "model_id": "twist30p0_rise3p400",
            "target_label": "D",
            "observed_d_A": "",
            "peak_status": "missing",
        }
    ]

    validate_long_rows(rows, source)
