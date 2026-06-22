#!/usr/bin/env python3
"""Run corrected diffraction for synthetic extended Asem twist stacks."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BUILD_AUDIT = ROOT / "outputs" / "metrics" / "asem_twist_extended_stack_build_audit.csv"
EXPERIMENTAL = ROOT / "inputs" / "experimental" / "nick_powder_profile_corrected_emory.csv"
OUTPUT_ROOT = ROOT / "outputs" / "asem_twist_extended_stack_corrected_emory_profile"
RADIAL_DIR = OUTPUT_ROOT / "radial_profiles"
PLOT_DIR = OUTPUT_ROOT / "plots"
TABLE_DIR = OUTPUT_ROOT / "tables"
METRICS = ROOT / "outputs" / "metrics"
FEATURE_CSV = METRICS / "asem_twist_extended_stack_feature_summary.csv"
RANKING_CSV = METRICS / "asem_twist_extended_stack_ranking.csv"
REPORT = ROOT / "outputs" / "reports" / "asem_twist_extended_stack_diffraction_report.md"
NICK_RUNNER = ROOT / "scripts" / "run_nick_ideal_16mer_corrected.py"
ORIGINAL_FEATURES = METRICS / "asem_twist_series_feature_summary_corrected_emory_profile.csv"
ORIGINAL_RANKING = METRICS / "asem_twist_series_ranking_corrected_emory_profile.csv"

FEATURE_WINDOWS = [
    ("3.38/3.4 A", 3.4, 3.30, 3.50),
    ("3.77 A", 3.77, 3.67, 3.87),
    ("4.4 A", 4.4, 4.30, 4.50),
    ("5.6 A", 5.6, 5.50, 5.70),
    ("7.3 A", 7.3, 7.10, 7.50),
]
RANK_FIELDS = {
    "3.38/3.4 A": "base_stack_offset_A",
    "3.77 A": "feature_3p77_offset_A",
    "4.4 A": "feature_4p4_offset_A",
    "5.6 A": "feature_5p6_offset_A",
    "7.3 A": "feature_7p3_offset_A",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid-size", type=int, default=129)
    parser.add_argument("--grid-limit-mm", type=float, default=100.0)
    parser.add_argument("--radial-bins", type=int, default=420)
    parser.add_argument("--detector-distance-mm", type=float, default=338.4)
    parser.add_argument("--wavelength-angstrom", type=float, default=0.7749)
    parser.add_argument("--normalize-q-min", type=float, default=0.15)
    return parser.parse_args()


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load {name} from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def load_nick_runner():
    return load_module("nick_ideal_16mer_corrected_extended_helpers", NICK_RUNNER)


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT)).replace("\\", "/")
    except ValueError:
        return str(path)


def f(value: object) -> float | None:
    if value in (None, ""):
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def fmt(value: float | None) -> str:
    return "" if value is None or not math.isfinite(value) else f"{value:.12g}"


def profile_peak_and_area(nick, profile, experimental, d_min, d_max):
    sim_peak_d, sim_peak_i, sim_area = nick.profile_peak_and_area(profile, d_min, d_max)
    exp_peak_d, _exp_peak_i, exp_area = nick.profile_peak_and_area(experimental, d_min, d_max)
    offset = sim_peak_d - exp_peak_d if not (math.isnan(sim_peak_d) or math.isnan(exp_peak_d)) else math.nan
    return sim_peak_d, sim_peak_i, sim_area, exp_peak_d, exp_area, offset


def run_profile(args, nick, scripts_module, orientation_module, radial_module, xyz_path: Path):
    _atoms, coords_angstrom, atomic_numbers = nick.load_xyz(xyz_path, scripts_module.atomic_number)
    wavelength_mm = args.wavelength_angstrom * 1e-7
    coords_mm = coords_angstrom * 1e-7
    grid_limits = (-args.grid_limit_mm, args.grid_limit_mm)
    image = orientation_module.average_fiber_diffraction(
        atomic_numbers,
        coords_mm,
        wavelength_mm,
        args.detector_distance_mm,
        grid_limits,
        grid_limits,
        args.grid_size,
        args.grid_size,
        [0],
        list(range(0, 181, 5)),
    )
    return nick.radial_rows_from_image(image, args, radial_module)


def feature_rows_for_profile(nick, profile, experimental, twist, label, target) -> list[dict[str, object]]:
    rows = []
    for feature, _center, d_min, d_max in FEATURE_WINDOWS:
        sim_peak_d, sim_peak_i, sim_area, exp_peak_d, exp_area, offset = profile_peak_and_area(
            nick, profile, experimental, d_min, d_max
        )
        rows.append(
            {
                "twist_deg": twist,
                "stack_length_label": label,
                "target_repeat_count": target,
                "feature_window": feature,
                "experimental_peak_d_A": fmt(None if math.isnan(exp_peak_d) else exp_peak_d),
                "simulated_peak_d_A": fmt(None if math.isnan(sim_peak_d) else sim_peak_d),
                "peak_offset_A": fmt(None if math.isnan(offset) else offset),
                "abs_peak_offset_A": fmt(None if math.isnan(offset) else abs(offset)),
                "simulated_peak_intensity_norm": fmt(None if math.isnan(sim_peak_i) else sim_peak_i),
                "window_area_simulated": fmt(None if math.isnan(sim_area) else sim_area),
                "window_area_experimental": fmt(None if math.isnan(exp_area) else exp_area),
            }
        )
    return rows


def add_original_feature_rows() -> list[dict[str, object]]:
    rows = []
    for row in read_csv(ORIGINAL_FEATURES):
        if row["feature_window"] not in RANK_FIELDS:
            continue
        rows.append(
            {
                "twist_deg": row["twist_label_deg"],
                "stack_length_label": "original",
                "target_repeat_count": "",
                "feature_window": row["feature_window"],
                "experimental_peak_d_A": row["experimental_peak_d_A"],
                "simulated_peak_d_A": row["simulated_peak_d_A"],
                "peak_offset_A": row["peak_offset_d_A"],
                "abs_peak_offset_A": row["abs_peak_offset_d_A"],
                "simulated_peak_intensity_norm": row["simulated_peak_intensity_norm"],
                "window_area_simulated": row["window_area_simulated"],
                "window_area_experimental": row["window_area_experimental"],
            }
        )
    return rows


def build_ranking(feature_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    by_model: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in feature_rows:
        by_model[(str(row["twist_deg"]), str(row["stack_length_label"]), str(row["target_repeat_count"]))].append(row)
    ranking = []
    original_31_mean = None
    for row in read_csv(ORIGINAL_RANKING):
        if row.get("twist_label_deg") == "31":
            original_31_mean = f(row.get("mean_abs_primary_peak_offset_A"))
    for (twist, label, target), rows in sorted(by_model.items(), key=lambda item: (str(item[0][2]).zfill(4), item[0][0])):
        offsets = {RANK_FIELDS[row["feature_window"]]: f(row["peak_offset_A"]) for row in rows}
        abs_values = [abs(value) for value in offsets.values() if value is not None]
        mean = float(np.mean(abs_values)) if abs_values else math.nan
        max_v = float(np.max(abs_values)) if abs_values else math.nan
        ranking.append(
            {
                "twist_deg": twist,
                "stack_length_label": label,
                "target_repeat_count": target,
                "mean_abs_primary_peak_offset_A": fmt(None if math.isnan(mean) else mean),
                "max_abs_primary_peak_offset_A": fmt(None if math.isnan(max_v) else max_v),
                "base_stack_offset_A": fmt(offsets.get("base_stack_offset_A")),
                "feature_3p77_offset_A": fmt(offsets.get("feature_3p77_offset_A")),
                "feature_4p4_offset_A": fmt(offsets.get("feature_4p4_offset_A")),
                "feature_5p6_offset_A": fmt(offsets.get("feature_5p6_offset_A")),
                "feature_7p3_offset_A": fmt(offsets.get("feature_7p3_offset_A")),
                "rank_within_length": "",
                "status_within_length": "",
                "does_29_30_tie_persist": "",
                "is_31_more_disfavored_than_original": "",
                "qualitative_note": "",
            }
        )
    by_length: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in ranking:
        by_length[str(row["stack_length_label"])].append(row)
    for label, rows in by_length.items():
        rows.sort(key=lambda row: f(row["mean_abs_primary_peak_offset_A"]) or math.inf)
        best = f(rows[0]["mean_abs_primary_peak_offset_A"]) if rows else None
        cutoff = None if best is None else best + 0.005
        mean_by_twist = {row["twist_deg"]: f(row["mean_abs_primary_peak_offset_A"]) for row in rows}
        tie_29_30 = (
            mean_by_twist.get("29") is not None
            and mean_by_twist.get("30") is not None
            and abs(mean_by_twist["29"] - mean_by_twist["30"]) <= 0.005
        )
        for idx, row in enumerate(rows, start=1):
            mean = f(row["mean_abs_primary_peak_offset_A"])
            row["rank_within_length"] = idx
            row["status_within_length"] = (
                "survives_near_best_filter"
                if mean is not None and cutoff is not None and mean <= cutoff
                else "disfavored_by_peak_offsets"
            )
            row["does_29_30_tie_persist"] = "yes" if tie_29_30 else "no"
            if row["twist_deg"] == "31" and mean is not None and original_31_mean is not None:
                row["is_31_more_disfavored_than_original"] = "yes" if mean > original_31_mean + 0.005 else "no"
            row["qualitative_note"] = f"Near-best cutoff for {label}: {fmt(cutoff)} A."
    ranking.sort(key=lambda row: (str(row["stack_length_label"]), int(row["rank_within_length"])))
    return ranking


def plot_mean_offsets(ranking: list[dict[str, object]]) -> None:
    labels = ["original", "extended_32layer_equiv", "extended_64layer_equiv", "extended_100layer_equiv"]
    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(8, 4.8), dpi=180)
    for twist, color in [("29", "#4c78a8"), ("30", "#59a14f"), ("31", "#e15759")]:
        vals = []
        for label in labels:
            row = next((r for r in ranking if r["twist_deg"] == twist and r["stack_length_label"] == label), None)
            vals.append(f(row["mean_abs_primary_peak_offset_A"]) if row else np.nan)
        ax.plot(x, vals, marker="o", label=f"{twist} deg", color=color)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=25, ha="right")
    ax.set_ylabel("Mean abs primary peak offset (A)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "mean_primary_offset_by_twist_and_length.png")
    plt.close(fig)


def plot_feature_offsets(feature_rows: list[dict[str, object]]) -> None:
    rows = [r for r in feature_rows if str(r["stack_length_label"]) != "original"]
    labels = sorted({str(r["stack_length_label"]) for r in rows})
    fig, axes = plt.subplots(len(FEATURE_WINDOWS), 1, figsize=(8.5, 10), dpi=180, sharex=True)
    for ax, (feature, *_rest) in zip(axes, FEATURE_WINDOWS):
        for twist, color in [("29", "#4c78a8"), ("30", "#59a14f"), ("31", "#e15759")]:
            vals = []
            for label in labels:
                row = next((r for r in rows if r["twist_deg"] == twist and r["stack_length_label"] == label and r["feature_window"] == feature), None)
                vals.append(f(row["peak_offset_A"]) if row else np.nan)
            ax.plot(labels, vals, marker="o", label=f"{twist} deg", color=color)
        ax.axhline(0, color="black", lw=0.6)
        ax.set_ylabel(feature)
    axes[0].legend(fontsize=7)
    axes[-1].tick_params(axis="x", rotation=25)
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "feature_offsets_by_twist_and_length.png")
    plt.close(fig)


def plot_profiles(profile_map: dict[tuple[str, str], list[dict[str, float]]], experimental: list[dict[str, float]]) -> None:
    exp = sorted([r for r in experimental if 2.8 <= r["d_A"] <= 8.4], key=lambda r: r["d_A"])
    for twist in ("29", "30", "31"):
        fig, ax = plt.subplots(figsize=(8.5, 4.8), dpi=180)
        ax.plot([r["d_A"] for r in exp], [r["intensity_norm"] for r in exp], color="black", label="experimental")
        for idx, label in enumerate(["extended_32layer_equiv", "extended_64layer_equiv", "extended_100layer_equiv"]):
            rows = sorted([r for r in profile_map[(twist, label)] if 2.8 <= r["d_A"] <= 8.4], key=lambda r: r["d_A"])
            ax.plot([r["d_A"] for r in rows], [r["intensity_norm"] + 0.9 * (idx + 1) for r in rows], label=label)
        ax.set_xlim(8.4, 2.8)
        ax.set_title(f"{twist} degree profiles by length")
        ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(PLOT_DIR / f"extended_stack_profiles_by_length_{twist}deg.png")
        plt.close(fig)
    for target, label in [(32, "extended_32layer_equiv"), (64, "extended_64layer_equiv"), (100, "extended_100layer_equiv")]:
        fig, ax = plt.subplots(figsize=(8.5, 4.8), dpi=180)
        ax.plot([r["d_A"] for r in exp], [r["intensity_norm"] for r in exp], color="black", label="experimental")
        for idx, twist in enumerate(["29", "30", "31"]):
            rows = sorted([r for r in profile_map[(twist, label)] if 2.8 <= r["d_A"] <= 8.4], key=lambda r: r["d_A"])
            ax.plot([r["d_A"] for r in rows], [r["intensity_norm"] + 0.9 * (idx + 1) for r in rows], label=f"{twist} deg")
        ax.set_xlim(8.4, 2.8)
        ax.set_title(f"Same-length overlay: {target} layer-equivalent")
        ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(PLOT_DIR / f"same_length_overlay_{target}.png")
        plt.close(fig)


def plot_differences(ranking: list[dict[str, object]], twist_a: str, twist_b: str, filename: str) -> None:
    labels = ["original", "extended_32layer_equiv", "extended_64layer_equiv", "extended_100layer_equiv"]
    diffs = []
    for label in labels:
        a = next((r for r in ranking if r["twist_deg"] == twist_a and r["stack_length_label"] == label), None)
        b = next((r for r in ranking if r["twist_deg"] == twist_b and r["stack_length_label"] == label), None)
        diffs.append((f(a["mean_abs_primary_peak_offset_A"]) or np.nan) - (f(b["mean_abs_primary_peak_offset_A"]) or np.nan))
    fig, ax = plt.subplots(figsize=(7.2, 4.2), dpi=180)
    ax.plot(labels, diffs, marker="o")
    ax.axhline(0, color="black", lw=0.7)
    ax.set_ylabel(f"Mean offset difference: {twist_a} - {twist_b} (A)")
    ax.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    fig.savefig(PLOT_DIR / filename)
    plt.close(fig)


def write_report(build_rows, feature_rows, ranking) -> None:
    by_label = defaultdict(list)
    for row in ranking:
        by_label[str(row["stack_length_label"])].append(row)
    lines = [
        "# Asem Twist Extended Stack Diffraction Report",
        "",
        "## Purpose",
        "",
        "Test whether longer synthetic stack lengths improve corrected powder-diffraction discrimination between Asem 29, 30, and 31 degree twist models.",
        "",
        "## Inputs",
        "",
        "Existing Asem 29/30/31 models and the Emory-corrected experimental profile were used. The synthetic extensions built by `scripts/build_extended_asem_twist_stacks.py` are layer-equivalent repeats, not new chemically searched candidates.",
        "",
        "## Build Method",
        "",
        "Complete six-residue layer-equivalent units were repeated with the source model's mean screw transform inferred from layer-centroid PCA axis, rise, and twist. New residue IDs were assigned and heavy-atom deduped XYZ files were generated.",
        "",
        "## Diffraction Method",
        "",
        "Corrected Asem non-accumulating/vectorized diffraction was run with `tilts = [0]`, `rotations = range(0, 181, 5)`, hydrogens excluded, exact heavy-atom deduplication, and the same grid/radial settings as the corrected baseline.",
        "",
        "## Ranking By Length",
        "",
        "| length | twist | mean abs offset A | max abs offset A | 3.77 offset A | rank | status |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for label in sorted(by_label):
        for row in sorted(by_label[label], key=lambda r: int(r["rank_within_length"])):
            lines.append(
                f"| {label} | {row['twist_deg']} | {row['mean_abs_primary_peak_offset_A']} | {row['max_abs_primary_peak_offset_A']} | {row['feature_3p77_offset_A']} | {row['rank_within_length']} | {row['status_within_length']} |"
            )
    lines.extend(["", "## Key Questions", ""])
    labels = ["extended_32layer_equiv", "extended_64layer_equiv", "extended_100layer_equiv"]
    for label in labels:
        rows = {r["twist_deg"]: r for r in by_label[label]}
        tie = rows.get("29", {}).get("does_29_30_tie_persist", "")
        status31 = rows.get("31", {}).get("status_within_length", "")
        lines.append(f"- {label}: 29/30 tie persists? {tie}; 31 status: {status31}.")
    row_lookup = {(r["stack_length_label"], r["twist_deg"]): r for r in ranking}
    def mean_for(label: str, twist: str) -> float | None:
        row = row_lookup.get((label, twist))
        return f(row.get("mean_abs_primary_peak_offset_A")) if row else None
    def offset_for(label: str, twist: str, field: str) -> float | None:
        row = row_lookup.get((label, twist))
        return f(row.get(field)) if row else None
    original_31 = mean_for("original", "31")
    lines.extend(
        [
            "- Does extending stack length separate 29 from 30? Partly. They remain within the near-best tolerance at 32 and 100 layer-equivalent lengths, but 64 layer-equivalent stacks favor 30 over 29 by the current mean-offset metric.",
            "- Does 31 become more disfavored as length increases? Mostly yes relative to the original at 64 and 100 layers, while 32 layers approximately preserves the original 31-degree mean offset.",
            "- Does the 3.77 A window remain the main discriminator? At original and 32/64 layers, yes: 31 carries the large 3.77 A offset. At 100 layers, 31's 3.77 A peak returns near the experimental position and the largest 31-degree error moves to the 7.3 A window, which is a warning that synthetic coherent repetition can change the apparent discriminator.",
            "- Are any new windows sensitive to length? Yes. The 7.3 A and 4.4 A windows become more length-sensitive in the 100-layer synthetic stacks, especially for 31 degrees.",
            "- Are peak positions stable but intensities changing? Not entirely. Some peak positions remain stable, but several extended-stack windows jump between nearby radial bins, so this is not only an intensity effect.",
            "- Do longer synthetic stacks sharpen or shift radial features? They can do both. The 64/100 layer-equivalent results show enough peak-position changes to treat the longest stacks as coherence sensitivity tests, not literal refined models.",
            "- Does the result strengthen near-30 family or exact-30 framing? Overall it still supports a near-30 family. The 64-layer result hints at possible 30-degree preference, but 100 layers returns 29 to the near-best set, so exact 30 degrees is not uniquely established.",
            "",
            "## Interpretation",
            "",
            "The length-extension test is mixed. It does not produce a clean monotonic separation of 29 from 30. It does keep 31 disfavored at every tested length, but the feature responsible for that disfavoring changes in the 100-layer synthetic stack. That pattern is consistent with finite-length and coherence sensitivity, and it cautions against overclaiming either exact 30-degree selection or a literal long-stack structural prediction from these synthetic repeats.",
            "",
            "## Synthetic-Length Artifacts",
            "",
            "The inferred centroid screw transform is a layer-equivalent proxy, not a chemically generated periodic model. The 100-layer stacks in particular show feature-window jumps that are likely affected by ideal coherent repetition and radial-bin sampling. These artifacts are useful for sensitivity testing but should not be interpreted as minimized molecular behavior.",
            "",
            "## Limitations",
            "",
            "- Synthetic repeated stacks are not energy-minimized or independently generated.",
            "- No AmberTools/tleap, pNAB, minimization, notebooks, or new candidate search was run.",
            "- Longer ideal repeats may exaggerate coherent diffraction relative to real disordered fibers.",
            "- This model-length test does not replace experimental structural refinement.",
            "",
            "## Recommendation",
            "",
            "Use these results as a sensitivity screen. The decisive next chemistry-facing test remains multiple completed full 30 degree candidates from Asem, ideally with carboxylates built, screened with the same workflow.",
        ]
    )
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    for directory in (RADIAL_DIR, PLOT_DIR, TABLE_DIR, METRICS, REPORT.parent):
        directory.mkdir(parents=True, exist_ok=True)
    nick = load_nick_runner()
    scripts_module, orientation_module, radial_module = nick.load_engine()
    experimental = nick.read_experimental_profile(EXPERIMENTAL)
    build_rows = [row for row in read_csv(BUILD_AUDIT) if row.get("build_status") == "built"]
    feature_rows = add_original_feature_rows()
    profile_map = {}
    metadata = {
        "diffraction_rerun": "extended synthetic stacks only",
        "tilts": [0],
        "rotations": "range(0, 181, 5)",
        "hydrogens_excluded": True,
        "exact_heavy_atom_deduplication": True,
    }
    for row in build_rows:
        twist = row["source_twist_deg"]
        target = int(row["target_repeat_count"])
        label = f"extended_{target}layer_equiv"
        xyz_path = ROOT / row["output_xyz"]
        profile = run_profile(args, nick, scripts_module, orientation_module, radial_module, xyz_path)
        profile_map[(twist, label)] = profile
        radial_path = RADIAL_DIR / f"asem_{twist}deg_extended_{target}layer_equiv_radial_profile.csv"
        nick.write_profile(radial_path, profile)
        feature_rows.extend(feature_rows_for_profile(nick, profile, experimental, twist, label, str(target)))
    feature_fields = [
        "twist_deg",
        "stack_length_label",
        "target_repeat_count",
        "feature_window",
        "experimental_peak_d_A",
        "simulated_peak_d_A",
        "peak_offset_A",
        "abs_peak_offset_A",
        "simulated_peak_intensity_norm",
        "window_area_simulated",
        "window_area_experimental",
    ]
    write_csv(FEATURE_CSV, feature_rows, feature_fields)
    write_csv(TABLE_DIR / FEATURE_CSV.name, feature_rows, feature_fields)
    ranking = build_ranking(feature_rows)
    ranking_fields = [
        "twist_deg",
        "stack_length_label",
        "target_repeat_count",
        "mean_abs_primary_peak_offset_A",
        "max_abs_primary_peak_offset_A",
        "base_stack_offset_A",
        "feature_3p77_offset_A",
        "feature_4p4_offset_A",
        "feature_5p6_offset_A",
        "feature_7p3_offset_A",
        "rank_within_length",
        "status_within_length",
        "does_29_30_tie_persist",
        "is_31_more_disfavored_than_original",
        "qualitative_note",
    ]
    write_csv(RANKING_CSV, ranking, ranking_fields)
    write_csv(TABLE_DIR / RANKING_CSV.name, ranking, ranking_fields)
    (TABLE_DIR / "run_metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    plot_mean_offsets(ranking)
    plot_feature_offsets(feature_rows)
    plot_profiles(profile_map, experimental)
    plot_differences(ranking, "29", "30", "difference_29_vs_30_by_length.png")
    plot_differences(ranking, "31", "30", "difference_31_vs_30_by_length.png")
    write_report(build_rows, feature_rows, ranking)
    print(f"Wrote {FEATURE_CSV}")
    print(f"Wrote {RANKING_CSV}")
    print(f"Wrote {REPORT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

