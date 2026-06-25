import pytest

from scripts.adapt_asem_radial_profiles import (
    adapt_profile_rows,
    collect_input_files,
    derive_model_id,
)


def test_derive_model_id_removes_radial_suffix():
    assert derive_model_id(__import__("pathlib").Path("compact_hexaplex_twist_30_radial.csv")) == "compact_hexaplex_twist_30"


def test_adapt_profile_rows_maps_d_and_intensity_columns(tmp_path):
    source = tmp_path / "compact_hexaplex_twist_30_radial.csv"
    rows = [
        {"d_A": "3.407", "intensity_mean": "100.5"},
        {"d_A": "3.853", "intensity_mean": "90.25"},
    ]

    adapted = adapt_profile_rows(rows, "compact_hexaplex_twist_30", source)

    assert adapted == [
        {
            "model_id": "compact_hexaplex_twist_30",
            "d_A": "3.407000000",
            "intensity": "100.500000000",
        },
        {
            "model_id": "compact_hexaplex_twist_30",
            "d_A": "3.853000000",
            "intensity": "90.250000000",
        },
    ]


def test_adapt_profile_rows_rejects_missing_columns(tmp_path):
    source = tmp_path / "bad.csv"
    rows = [{"d_A": "3.407"}]

    with pytest.raises(ValueError, match="missing required columns"):
        adapt_profile_rows(rows, "m1", source)


def test_collect_input_files_accepts_single_file(tmp_path):
    source = tmp_path / "one_radial.csv"
    source.write_text("d_A,intensity_mean\n3.4,10\n", encoding="utf-8")

    files = collect_input_files(source, "*_radial.csv")

    assert files == [source]


def test_collect_input_files_uses_directory_pattern(tmp_path):
    a = tmp_path / "a_radial.csv"
    b = tmp_path / "b_radial.csv"
    ignored = tmp_path / "ignored.csv"
    a.write_text("d_A,intensity_mean\n3.4,10\n", encoding="utf-8")
    b.write_text("d_A,intensity_mean\n3.5,20\n", encoding="utf-8")
    ignored.write_text("d_A,intensity_mean\n3.6,30\n", encoding="utf-8")

    files = collect_input_files(tmp_path, "*_radial.csv")

    assert files == [a, b]
