#!/usr/bin/env python3
"""Broaden existing ideal 16-mer radial profiles for Nick-style visual comparison."""

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
DEFAULT_PLOT_DIR = ROOT / "outputs" / "nick_ideal_16mer_corrected_emory_profile" / "plots"
DEFAULT_AUDIT = ROOT / "outputs" / "metrics" / "nick_ideal_16mer_broadening_audit.csv"

FEATURE_WINDOWS = [
    ("7.25 A", 7.25, 7.0, 7.5),
    ("5.6 A", 5.6, 5.5, 5.7),
    ("4.4 A", 4.4, 4.3, 4.5),
    ("3.7 A", 3.7, 3.6, 3.8),
    ("3.4 A", 3.4, 3.3, 3.5),
]

D_SIGMAS = [0.03, 0.05, 0.10]
Q_SIGMAS = [0.010]
RECOMMENDED_D_SIGMA = 0.10


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--experimental-profile", type=Path, default=DEFAULT_EXPERIMENTAL)
    parser.add_argument("--simulated-profile", type=Path, default=DEFAULT_SIMULATED)
    parser.add_argument("--plot-dir", type=Path, default=DEFAULT_PLOT_DIR)
    parser.add_argument("--audit-output", type=Path, default=DEFAULT_AUDIT)
    parser.add_argument("--d-min", type=float, default=3.0)
    parser.add_argument("--d-max", type=float, default=8.0)
    parser.add_argument("--grid-points", type=int, default=1600)
    return parser.parse_args()


def read_profile(path: Path, intensity_column: str) -> list[dict[str, float]]:
    rows = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            d_value = float(row["d_A"])
            q_value = float(row.get("q_Ainv", 2.0 * np.pi / d_value))
            rows.append(
                {
                    "d_A": d_value,
                    "q_Ainv": q_value,
                    "intensity": float(row[intensity_column]),
                }
            )
    return rows


def grid_profile(
    rows: list[dict[str, float]],
    x_key: str,
    d_min: float,
    d_max: float,
    points: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    selected = [
        row
        for row in rows
        if np.isfinite(row["d_A"]) and np.isfinite(row[x_key]) and d_min <= row["d_A"] <= d_max
    ]
    selected.sort(key=lambda row: row[x_key])
    x = np.asarray([row[x_key] for row in selected], dtype=np.float64)
    d = np.asarray([row["d_A"] for row in selected], dtype=np.float64)
    y = np.asarray([row["intensity"] for row in selected], dtype=np.float64)
    unique_x, unique_indices = np.unique(x, return_index=True)
    x = unique_x
    d = d[unique_indices]
    y = y[unique_indices]
    x_grid = np.linspace(float(np.min(x)), float(np.max(x)), points)
    y_grid = np.interp(x_grid, x, y)
    d_grid = np.interp(x_grid, x, d)
    return x_grid, d_grid, y_grid


def gaussian_kernel(sigma: float, step: float) -> np.ndarray:
    radius = max(3, int(np.ceil(4.0 * sigma / step)))
    offsets = np.arange(-radius, radius + 1, dtype=np.float64) * step
    kernel = np.exp(-0.5 * (offsets / sigma) ** 2)
    kernel /= np.sum(kernel)
    return kernel


def smooth_uniform(y_values: np.ndarray, sigma: float, step: float) -> np.ndarray:
    kernel = gaussian_kernel(sigma, step)
    padded = np.pad(y_values, (len(kernel) // 2, len(kernel) // 2), mode="edge")
    return np.convolve(padded, kernel, mode="valid")


def normalize(y_values: np.ndarray) -> np.ndarray:
    y = np.asarray(y_values, dtype=np.float64)
    y = y - np.nanmin(y)
    peak = np.nanmax(y)
    if peak > 0:
        y = y / peak
    return y


def profile_peak(rows: list[dict[str, float]], d_min: float, d_max: float) -> tuple[float | None, float | None]:
    points = [row for row in rows if d_min <= row["d_A"] <= d_max]
    if not points:
        return None, None
    peak = max(points, key=lambda row: row["intensity"])
    return peak["d_A"], peak["intensity"]


def smoothed_rows(d_values: np.ndarray, y_values: np.ndarray) -> list[dict[str, float]]:
    return [
        {"d_A": float(d_value), "intensity": float(intensity)}
        for d_value, intensity in zip(d_values, y_values)
        if np.isfinite(d_value) and np.isfinite(intensity)
    ]


def slug_sigma(sigma: float, suffix: str) -> str:
    return f"{sigma:.2f}{suffix}".replace(".", "p")


def plot_stacked(
    path: Path,
    exp_d: np.ndarray,
    exp_y: np.ndarray,
    sim_d: np.ndarray,
    sim_y: np.ndarray,
    title: str,
    d_min: float,
    d_max: float,
) -> None:
    fig, ax = plt.subplots(figsize=(6.6, 5.0), dpi=180)
    fig.patch.set_facecolor("#8a8a8a")
    ax.set_facecolor("#8a8a8a")

    offsets = {"experimental": 3.0, "simulated": 2.0}
    for label, _center, d_lo, d_hi in FEATURE_WINDOWS:
        ax.axvspan(d_lo, d_hi, color="white", alpha=0.08, linewidth=0)
        ax.text(
            (d_lo + d_hi) / 2.0,
            1.84,
            label,
            ha="center",
            va="bottom",
            fontsize=6.8,
            color="black",
            alpha=0.75,
        )

    ax.plot(exp_d, normalize(exp_y) + offsets["experimental"], color="black", linewidth=1.35, label="Corrected Experimental")
    ax.plot(
        sim_d,
        normalize(sim_y) + offsets["simulated"],
        color="red",
        linewidth=1.25,
        label="Ideal 16-mer AntiParallel 30deg",
    )

    ax.set_xlim(d_max, d_min)
    ax.set_ylim(1.78, 4.15)
    ax.set_xlabel("A")
    ax.set_ylabel("offset normalized intensity")
    ax.set_title(title, fontsize=10)
    ax.legend(loc="upper left", fontsize=6.5, frameon=True, facecolor="white", edgecolor="none")
    ax.tick_params(axis="both", colors="black", labelsize=7, direction="out", length=3)
    for spine in ax.spines.values():
        spine.set_color("black")
    ax.grid(False)
    fig.tight_layout(pad=1.0)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)


def write_audit(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames = [
        "smoothing_domain",
        "sigma",
        "sigma_units",
        "feature_window",
        "simulated_peak_d_A_after_smoothing",
        "experimental_peak_d_A",
        "peak_offset_d_A",
        "notes",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    experimental = read_profile(args.experimental_profile, "intensity_normalized")
    simulated = read_profile(args.simulated_profile, "intensity_norm")
    exp_d_grid = np.linspace(args.d_min, args.d_max, args.grid_points)
    exp_source = sorted(
        [row for row in experimental if args.d_min <= row["d_A"] <= args.d_max],
        key=lambda row: row["d_A"],
    )
    exp_d = np.asarray([row["d_A"] for row in exp_source], dtype=np.float64)
    exp_y = np.asarray([row["intensity"] for row in exp_source], dtype=np.float64)
    exp_y_grid = np.interp(exp_d_grid, exp_d, exp_y)
    audit_rows: list[dict[str, object]] = []

    d_x, d_grid, d_y = grid_profile(simulated, "d_A", args.d_min, args.d_max, args.grid_points)
    d_step = float(np.mean(np.diff(d_x)))
    for sigma in D_SIGMAS:
        smooth_y = smooth_uniform(d_y, sigma, d_step)
        output = args.plot_dir / f"nick_ideal_16mer_broadened_sigma_{slug_sigma(sigma, 'A')}.png"
        plot_stacked(
            output,
            exp_d_grid,
            exp_y_grid,
            d_grid,
            smooth_y,
            f"Nick-style broadened ideal 16-mer (sigma {sigma:.2f} A)",
            args.d_min,
            args.d_max,
        )
        rows = smoothed_rows(d_grid, smooth_y)
        for feature, _center, d_lo, d_hi in FEATURE_WINDOWS:
            sim_peak_d, _sim_peak_i = profile_peak(rows, d_lo, d_hi)
            exp_peak_d, _exp_peak_i = profile_peak(experimental, d_lo, d_hi)
            audit_rows.append(
                {
                    "smoothing_domain": "d",
                    "sigma": f"{sigma:.6g}",
                    "sigma_units": "A",
                    "feature_window": feature,
                    "simulated_peak_d_A_after_smoothing": f"{sim_peak_d:.12g}" if sim_peak_d is not None else "",
                    "experimental_peak_d_A": f"{exp_peak_d:.12g}" if exp_peak_d is not None else "",
                    "peak_offset_d_A": f"{(sim_peak_d - exp_peak_d):.12g}"
                    if sim_peak_d is not None and exp_peak_d is not None
                    else "",
                    "notes": "Gaussian smoothing on uniform d-spacing grid for visualization only.",
                }
            )
        if sigma == RECOMMENDED_D_SIGMA:
            recommended = args.plot_dir / "nick_ideal_16mer_vs_corrected_experimental_profile_broadened_nick_style.png"
            plot_stacked(
                recommended,
                exp_d_grid,
                exp_y_grid,
                d_grid,
                smooth_y,
                f"Nick-style broadened ideal 16-mer (sigma {sigma:.2f} A)",
                args.d_min,
                args.d_max,
            )

    q_x, q_d_grid, q_y = grid_profile(simulated, "q_Ainv", args.d_min, args.d_max, args.grid_points)
    q_step = float(np.mean(np.diff(q_x)))
    for sigma in Q_SIGMAS:
        smooth_y = smooth_uniform(q_y, sigma, q_step)
        rows = smoothed_rows(q_d_grid, smooth_y)
        for feature, _center, d_lo, d_hi in FEATURE_WINDOWS:
            sim_peak_d, _sim_peak_i = profile_peak(rows, d_lo, d_hi)
            exp_peak_d, _exp_peak_i = profile_peak(experimental, d_lo, d_hi)
            audit_rows.append(
                {
                    "smoothing_domain": "q",
                    "sigma": f"{sigma:.6g}",
                    "sigma_units": "A^-1",
                    "feature_window": feature,
                    "simulated_peak_d_A_after_smoothing": f"{sim_peak_d:.12g}" if sim_peak_d is not None else "",
                    "experimental_peak_d_A": f"{exp_peak_d:.12g}" if exp_peak_d is not None else "",
                    "peak_offset_d_A": f"{(sim_peak_d - exp_peak_d):.12g}"
                    if sim_peak_d is not None and exp_peak_d is not None
                    else "",
                    "notes": "Gaussian smoothing on uniform q grid for audit comparison only.",
                }
            )

    write_audit(args.audit_output, audit_rows)
    print(f"Wrote broadening audit: {args.audit_output}")
    print(f"Wrote broadened plots under: {args.plot_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
