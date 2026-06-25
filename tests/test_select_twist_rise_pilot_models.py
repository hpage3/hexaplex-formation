import csv
from pathlib import Path

import pytest

from scripts.select_twist_rise_pilot_models import (
    parse_float_set,
    require_columns,
    select_pilot_rows,
)


def test_parse_float_set_handles_commas_and_spacing():
    assert parse_float_set("28, 29,30") == {28.0, 29.0, 30.0}


def test_parse_float_set_rejects_empty_selection():
    with pytest.raises(ValueError, match="cannot be empty"):
        parse_float_set(" , ")


def test_select_pilot_rows_selects_twist_rise_cross_product():
    rows = [
        {"model_id": "twist28p0_rise3p350", "twist_deg": "28.000", "rise_A": "3.350"},
        {"model_id": "twist28p0_rise3p380", "twist_deg": "28.000", "rise_A": "3.380"},
        {"model_id": "twist29p0_rise3p350", "twist_deg": "29.000", "rise_A": "3.350"},
        {"model_id": "twist30p0_rise3p400", "twist_deg": "30.000", "rise_A": "3.400"},
        {"model_id": "twist31p0_rise3p450", "twist_deg": "31.000", "rise_A": "3.450"},
        {"model_id": "twist32p0_rise3p450", "twist_deg": "32.000", "rise_A": "3.450"},
        {"model_id": "twist33p0_rise3p450", "twist_deg": "33.000", "rise_A": "3.450"},
    ]

    selected = select_pilot_rows(
        rows,
        twists={28.0, 29.0, 30.0, 31.0, 32.0},
        rises={3.35, 3.38, 3.40, 3.45},
    )

    assert [row["model_id"] for row in selected] == [
        "twist28p0_rise3p350",
        "twist28p0_rise3p380",
        "twist29p0_rise3p350",
        "twist30p0_rise3p400",
        "twist31p0_rise3p450",
        "twist32p0_rise3p450",
    ]
    assert all(row["pilot_selection"] == "true" for row in selected)


def test_require_columns_rejects_missing_manifest_columns(tmp_path):
    path = tmp_path / "manifest.csv"
    rows = [{"model_id": "m1", "twist_deg": "30"}]

    with pytest.raises(ValueError, match="missing required columns"):
        require_columns(rows, path)
