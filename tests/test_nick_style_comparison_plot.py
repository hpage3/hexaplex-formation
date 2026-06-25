from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def write_profile(path: Path, intensity_column: str) -> None:
    rows = [
        (3.0, 0.10),
        (3.4, 1.00),
        (3.8, 0.45),
        (4.4, 0.65),
        (5.6, 0.40),
        (7.3, 0.55),
        (8.0, 0.05),
    ]
    path.write_text(
        "d_A," + intensity_column + "\n"
        + "\n".join(f"{d},{i}" for d, i in rows)
        + "\n",
        encoding="utf-8",
    )


def test_nick_style_comparison_accepts_intensity_mean_label_and_title(tmp_path: Path) -> None:
    experimental = tmp_path / "experimental.csv"
    simulated = tmp_path / "simulated.csv"
    output = tmp_path / "comparison.png"

    write_profile(experimental, "intensity_normalized")
    write_profile(simulated, "intensity_mean")

    result = subprocess.run(
        [
            sys.executable,
            "scripts/plot_nick_ideal_16mer_style_comparison.py",
            "--experimental-profile",
            str(experimental),
            "--simulated-profile",
            str(simulated),
            "--simulated-intensity-column",
            "intensity_mean",
            "--simulated-label",
            "30 deg / 3.40 A smoke",
            "--title",
            "Smoke comparison",
            "--output",
            str(output),
        ],
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    assert result.returncode == 0, result.stderr
    assert output.exists()
    assert output.stat().st_size > 0
