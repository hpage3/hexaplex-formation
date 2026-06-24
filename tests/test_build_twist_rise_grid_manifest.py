from scripts.build_twist_rise_grid_manifest import build_rows, decimal_range, model_id


def test_decimal_range_includes_stop():
    values = decimal_range("20.0", "21.0", "0.5")
    assert [str(v) for v in values] == ["20.0", "20.5", "21.0"]


def test_model_id_is_stable():
    assert model_id(values_twist := decimal_range("29.5", "29.5", "0.5")[0], decimal_range("3.38", "3.38", "0.01")[0]) == "twist29p5_rise3p380"


def test_build_rows_small_grid():
    rows = build_rows(
        twist_min="29.0",
        twist_max="30.0",
        twist_step="0.5",
        rise_min="3.35",
        rise_max="3.40",
        rise_step="0.05",
    )

    assert len(rows) == 6
    assert rows[0]["model_id"] == "twist29p0_rise3p350"
    assert rows[-1]["model_id"] == "twist30p0_rise3p400"
    assert rows[0]["model_status"] == "pending"
    assert rows[0]["diffraction_status"] == "pending"
    assert rows[0]["score_status"] == "pending"
