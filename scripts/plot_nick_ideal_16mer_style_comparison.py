#!/usr/bin/env python3
"""Create a Nick-style visual comparison from existing ideal 16-mer CSV outputs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_EXPERIMENTAL = ROOT / "inputs" / "experimental" / "nick_powder_profile_corrected_emory.csv"
DEFAULT_SIMULATED = (
    ROOT
    / "outputs"
    / "nick_ideal_16mer_corrected_emory_profile"
    / "radial_profiles"
    / "ideal_16mer_antiparallel_30deg_radial_profile.csv"
)
DEFAULT_OUTPUT = (
    ROOT
    / "outputs"
    / "nick_ideal_16mer_corrected_emory_profile"
    / "plots"
    / "nick_ideal_16mer_vs_corrected_experimental_profile_nick_style.png"
)

FEATURE_WINDOWS = [
    ("3.4 A", 3.3, 3.5),
    ("3.7 A", 3.6, 3.8),
    ("4.4 A", 4.3, 4.5),
    ("5.6 A", 5.5, 5.7),
    ("7.25 A", 7.0, 7.5),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--experimental-profile", type=Path, default=DEFAULT_EXPERIMENTAL)
    parser.add_argument("--simulated-profile", type=Path, default=DEFAULT_SIMULATED)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--d-min", type=float, default=3.0)
    parser.add_argument("--d-max", type=float, default=8.0)
    parser.add_argument("--simulated-offset", type=float, default=1.2)
    return parser.parse_args()


def read_profile(path: Path, intensity_column: str) -> tuple[np.ndarray, np.ndarray]:
    d_values = []
    intensities = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            d_values.append(float(row["d_A"]))
            intensities.append(float(row[intensity_column]))
    return np.asarray(d_values, dtype=np.float64), np.asarray(intensities, dtype=np.float64)


def windowed_normalize(
    d_values: np.ndarray,
    intensities: np.ndarray,
    d_min: float,
    d_max: float,
) -> tuple[np.ndarray, np.ndarray]:
    mask = np.isfinite(d_values) & np.isfinite(intensities) & (d_values >= d_min) & (d_values <= d_max)
    x = d_values[mask]
    y = intensities[mask]
    order = np.argsort(x)
    x = x[order]
    y = y[order]
    y = y - np.min(y)
    peak = np.max(y)
    if peak > 0:
        y = y / peak
    return x, y


def main() -> int:
    args = parse_args()
    exp_d, exp_i = read_profile(args.experimental_profile, "intensity_normalized")
    sim_d, sim_i = read_profile(args.simulated_profile, "intensity_norm")
    exp_x, exp_y = windowed_normalize(exp_d, exp_i, args.d_min, args.d_max)
    sim_x, sim_y = windowed_normalize(sim_d, sim_i, args.d_min, args.d_max)

    fig, ax = plt.subplots(figsize=(9.5, 5.25), dpi=180)
    for label, d_lo, d_hi in FEATURE_WINDOWS:
        ax.axvspan(d_lo, d_hi, color="#9aa0a6", alpha=0.10, linewidth=0)
        ax.text(
            (d_lo + d_hi) / 2.0,
            args.simulated_offset + 1.04,
            label,
            ha="center",
            va="bottom",
            fontsize=7.5,
            color="#555555",
        )

    ax.plot(exp_x, exp_y, color="black", linewidth=1.65, label="Corrected Experimental")
    ax.plot(
        sim_x,
        sim_y + args.simulated_offset,
        color="#1f77b4",
        linewidth=1.55,
        label="Ideal 16-mer AntiParallel 30deg",
    )

    ax.set_xlim(args.d_max, args.d_min)
    ax.set_ylim(-0.08, args.simulated_offset + 1.18)
    ax.set_xlabel("d-spacing (A)")
    ax.set_ylabel("independently normalized intensity; theoretical trace offset")
    ax.set_title("Nick-style ideal 16-mer comparison")
    ax.grid(True, alpha=0.25, linewidth=0.6)
    ax.legend(loc="upper right", frameon=False)
    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output)
    plt.close(fig)
    print(f"Wrote Nick-style plot: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
