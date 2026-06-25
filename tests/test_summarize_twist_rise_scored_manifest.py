import csv
from pathlib import Path

import pytest

from scripts.summarize_twist_rise_scored_manifest import (
    rank_rows,
    scored_rows,
    summarize_manifest,
)


def write_csv(path: Path, rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerows(rows)


def test_scored_rows_filters_only_scored_status():
    rows = [
        {"model_id": "m1", "score_status": "pending"},
        {"model_id": "m2", "score_status": "scored"},
    ]

    assert scored_rows(rows) == [{"model_id": "m2", "score_status": "scored"}]


def test_rank_rows_orders_by_completeness_then_combined_then_helical():
    rows = [
        {
            "model_id": "incomplete_low_rmsd",
            "score_completeness": "0.800000000",
            "combined_rmsd": "0.010",
            "helical_rmsd": "0.010",
        },
        {
            "model_id": "complete_higher_rmsd",
            "score_completeness": "1.000000000",
            "combined_rmsd": "0.090",
            "helical_rmsd": "0.080",
        },
        {
            "model_id": "complete_best",
            "score_completeness": "1.000000000",
            "combined_rmsd": "0.050",
            "helical_rmsd": "0.060",
        },
        {
            "model_id": "complete_tie_better_helical",
            "score_completeness": "1.000000000",
            "combined_rmsd": "0.050",
            "helical_rmsd": "0.040",
        },
    ]

    ranked = rank_rows(rows)

    assert [row["model_id"] for row in ranked] == [
        "complete_tie_better_helical",
        "complete_best",
        "complete_higher_rmsd",
        "incomplete_low_rmsd",
    ]
    assert [row["rank"] for row in ranked] == ["1", "2", "3", "4"]


def test_summarize_manifest_writes_ranked_csv_and_markdown(tmp_path):
    input_path = tmp_path / "scored.csv"
    ranked_csv = tmp_path / "ranked.csv"
    summary_md = tmp_path / "summary.md"

    write_csv(
        input_path,
        [
            [
                "model_id",
                "twist_deg",
                "rise_A",
                "score_status",
                "observed_peak_count",
                "expected_peak_count",
                "missing_peak_count",
                "score_completeness",
                "base_rmsd",
                "helical_rmsd",
                "combined_rmsd",
                "notes",
            ],
            ["m_pending", "30", "3.4", "pending", "0", "5", "5", "0", "", "", "", ""],
            ["m_good", "30", "3.4", "scored", "5", "5", "0", "1", "0.02", "0.04", "0.05", ""],
            ["m_less_good", "28", "3.4", "scored", "5", "5", "0", "1", "0.03", "0.08", "0.09", ""],
        ],
    )

    ranked = summarize_manifest(input_path, ranked_csv, summary_md, top_n=2)

    assert len(ranked) == 2
    assert ranked_csv.exists()
    assert summary_md.exists()

    ranked_text = ranked_csv.read_text(encoding="utf-8")
    assert "rank,model_id" in ranked_text
    assert "m_good" in ranked_text
    assert "m_pending" not in ranked_text

    summary_text = summary_md.read_text(encoding="utf-8")
    assert "Twist/Rise Scored Manifest Summary" in summary_text
    assert "Interpretation guardrail" in summary_text


def test_summarize_manifest_rejects_missing_required_columns(tmp_path):
    input_path = tmp_path / "bad.csv"
    ranked_csv = tmp_path / "ranked.csv"
    summary_md = tmp_path / "summary.md"

    write_csv(
        input_path,
        [
            ["model_id", "score_status"],
            ["m1", "scored"],
        ],
    )

    with pytest.raises(ValueError, match="missing required columns"):
        summarize_manifest(input_path, ranked_csv, summary_md, top_n=5)


def test_summarize_manifest_writes_header_when_no_scored_rows(tmp_path):
    input_path = tmp_path / "scored.csv"
    ranked_csv = tmp_path / "ranked.csv"
    summary_md = tmp_path / "summary.md"

    write_csv(
        input_path,
        [
            [
                "model_id",
                "twist_deg",
                "rise_A",
                "score_status",
                "score_completeness",
                "helical_rmsd",
                "combined_rmsd",
            ],
            ["m_pending", "30", "3.4", "pending", "0", "", ""],
        ],
    )

    ranked = summarize_manifest(input_path, ranked_csv, summary_md, top_n=5)

    assert ranked == []
    ranked_text = ranked_csv.read_text(encoding="utf-8")
    assert ranked_text.startswith("rank,model_id")
    assert "m_pending" not in ranked_text

    summary_text = summary_md.read_text(encoding="utf-8")
    assert "No scored rows were found." in summary_text
