#!/usr/bin/env python3
"""Validate and benchmark Asem's 29/30/31 twist-series models."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from hexaplex_formation.pdb_utils import (  # noqa: E402
    dedupe_exact_atoms,
    heavy_atoms,
    load_pdb_atoms,
    residue_count,
    residue_names,
)

INPUT_ROOT = ROOT / "inputs" / "asem_twist_series_29_30_31"
RAW_ROOT = INPUT_ROOT / "raw"
XYZ_ROOT = INPUT_ROOT / "xyz"
BASELINE_PDB = ROOT / "inputs" / "nick_ideal_models" / "Hexaplex_AntiParallel_30deg_Ideal.pdb"
EXPERIMENTAL = ROOT / "inputs" / "experimental" / "nick_powder_profile_corrected_emory.csv"
OUTPUT_ROOT = ROOT / "outputs" / "asem_twist_series_corrected_emory_profile"
RADIAL_DIR = OUTPUT_ROOT / "radial_profiles"
PLOT_DIR = OUTPUT_ROOT / "plots"
TABLE_DIR = OUTPUT_ROOT / "tables"
METRICS_DIR = ROOT / "outputs" / "metrics"
REPORT = ROOT / "outputs" / "reports" / "asem_twist_series_corrected_emory_profile_report.md"
STRUCTURAL_CSV = METRICS_DIR / "asem_twist_series_structural_summary.csv"
FEATURE_CSV = METRICS_DIR / "asem_twist_series_feature_summary_corrected_emory_profile.csv"
RANKING_CSV = METRICS_DIR / "asem_twist_series_ranking_corrected_emory_profile.csv"
NICK_RUNNER = ROOT / "scripts" / "run_nick_ideal_16mer_corrected.py"

FEATURE_WINDOWS = [
    ("3.38/3.4 A", 3.4, 3.30, 3.50),
    ("3.77 A", 3.77, 3.67, 3.87),
    ("4.4 A", 4.4, 4.30, 4.50),
    ("5.6 A", 5.6, 5.50, 5.70),
    ("7.3 A", 7.3, 7.10, 7.50),
]
EXPECTED_RESIDUES = {"CYP", "MEP", "GLU"}
PDB_SUFFIXES = (".pdb", ".pdb.txt", ".ent")


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
    return load_module("nick_ideal_16mer_corrected_helpers", NICK_RUNNER)


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT)).replace("\\", "/")
    except ValueError:
        return str(path)


def safe_label(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_").lower()


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def is_pdb_path(path: Path) -> bool:
    return path.name.lower().endswith(PDB_SUFFIXES)


def residue_key(atom) -> tuple[str, int | None, str, str]:
    return (atom.chain_id, atom.residue_number, atom.insertion_code, atom.residue_name)


def coord_array(atoms) -> np.ndarray:
    return np.asarray([[atom.x, atom.y, atom.z] for atom in atoms], dtype=float)


def candidate_paths() -> list[Path]:
    if not RAW_ROOT.exists():
        return []
    return sorted(path for path in RAW_ROOT.rglob("*") if path.is_file() and is_pdb_path(path))


def twist_from_path(path: Path) -> str:
    return path.relative_to(RAW_ROOT).parts[0]


def is_diffraction_ready(atoms: list) -> tuple[bool, str]:
    baseline_atoms = len(load_pdb_atoms(BASELINE_PDB))
    residues = set(residue_names(atoms))
    reasons = []
    if len(atoms) < 0.45 * baseline_atoms:
        reasons.append("atom count below full-model scale")
    if residue_count(atoms) < 0.8 * 180:
        reasons.append("residue count below full-model scale")
    if not residues.issuperset(EXPECTED_RESIDUES):
        reasons.append("missing expected CYP/MEP/GLU scaffold residues")
    if reasons:
        return False, "; ".join(reasons)
    return True, "full-scale CYP/MEP/GLU scaffold; comparable to ideal 16-mer heavy-atom XYZ convention"


def layer_stats(atoms: list) -> dict[str, object]:
    residues: dict[tuple[str, int | None, str, str], list] = defaultdict(list)
    for atom in atoms:
        residues[residue_key(atom)].append(atom)
    centroids = []
    for key, group in residues.items():
        centroids.append((key, coord_array(group).mean(axis=0)))
    if len(centroids) < 2:
        return {"rise_median_A": "", "rise_mean_A": "", "rough_twist_median_deg": "", "rough_twist_mean_deg": ""}
    centroids.sort(key=lambda item: (item[0][0], item[0][1] if item[0][1] is not None else -9999, item[0][3]))
    z = np.asarray([c[1][2] for c in centroids])
    dz = np.diff(np.sort(z))
    dz = dz[np.abs(dz) > 1e-6]
    xy = np.asarray([[c[1][0], c[1][1]] for c in centroids])
    center = xy.mean(axis=0)
    angles = np.degrees(np.unwrap(np.arctan2(xy[:, 1] - center[1], xy[:, 0] - center[0])))
    da = np.abs(np.diff(angles))
    da = da[(da > 1e-6) & (da < 180)]
    return {
        "rise_median_A": "" if len(dz) == 0 else f"{float(np.median(np.abs(dz))):.6g}",
        "rise_mean_A": "" if len(dz) == 0 else f"{float(np.mean(np.abs(dz))):.6g}",
        "rough_twist_median_deg": "" if len(da) == 0 else f"{float(np.median(da)):.6g}",
        "rough_twist_mean_deg": "" if len(da) == 0 else f"{float(np.mean(da)):.6g}",
    }


def hbond_proxy(atoms: list) -> dict[str, object]:
    n_atoms = [atom for atom in atoms if atom.element.upper() == "N" or atom.atom_name.upper().startswith("N")]
    o_atoms = [atom for atom in atoms if atom.element.upper() == "O" or atom.atom_name.upper().startswith("O")]
    h_atoms = [atom for atom in atoms if atom.element.upper() == "H" or atom.atom_name.upper().startswith("H")]
    best_no = math.inf
    best_h_o = math.inf
    best_angle = ""
    n_coord = coord_array(n_atoms) if n_atoms else np.empty((0, 3))
    o_coord = coord_array(o_atoms) if o_atoms else np.empty((0, 3))
    if len(n_coord) and len(o_coord):
        for i, n_atom in enumerate(n_atoms):
            d = np.linalg.norm(o_coord - n_coord[i], axis=1)
            mask = np.asarray([residue_key(o) != residue_key(n_atom) for o in o_atoms])
            if mask.any():
                local = d[mask]
                if len(local):
                    best_no = min(best_no, float(local.min()))
    if h_atoms and o_atoms:
        h_coord = coord_array(h_atoms)
        for h_i, h_atom in enumerate(h_atoms):
            distances = np.linalg.norm(o_coord - h_coord[h_i], axis=1)
            for o_i in np.argsort(distances)[:3]:
                if residue_key(o_atoms[o_i]) == residue_key(h_atom):
                    continue
                ho = float(distances[o_i])
                if ho < best_h_o:
                    best_h_o = ho
                    nearest_n = min(n_atoms, key=lambda n: np.linalg.norm(np.asarray([n.x, n.y, n.z]) - h_coord[h_i])) if n_atoms else None
                    if nearest_n is not None:
                        n_vec = np.asarray([nearest_n.x, nearest_n.y, nearest_n.z]) - h_coord[h_i]
                        o_vec = o_coord[o_i] - h_coord[h_i]
                        denom = np.linalg.norm(n_vec) * np.linalg.norm(o_vec)
                        if denom:
                            angle = math.degrees(math.acos(float(np.clip(np.dot(n_vec, o_vec) / denom, -1.0, 1.0))))
                            best_angle = f"{angle:.6g}"
    return {
        "hbond_proxy_best_n_o_distance_A": "" if math.isinf(best_no) else f"{best_no:.6g}",
        "hbond_proxy_best_h_o_distance_A": "" if math.isinf(best_h_o) else f"{best_h_o:.6g}",
        "hbond_proxy_best_n_h_o_angle_deg": best_angle,
    }


def clash_proxy(heavy: list) -> dict[str, object]:
    coords = coord_array(heavy)
    min_dist = math.inf
    counts = {2.0: 0, 2.2: 0, 2.5: 0}
    residue_keys = [residue_key(atom) for atom in heavy]
    for i in range(len(coords)):
        diff = coords[i + 1 :] - coords[i]
        if len(diff) == 0:
            continue
        d = np.linalg.norm(diff, axis=1)
        mask = np.asarray([key != residue_keys[i] for key in residue_keys[i + 1 :]])
        d = d[mask]
        if len(d) == 0:
            continue
        min_dist = min(min_dist, float(d.min()))
        for threshold in counts:
            counts[threshold] += int(np.count_nonzero(d < threshold))
    return {
        "min_heavy_heavy_nonbonded_distance_A": "" if math.isinf(min_dist) else f"{min_dist:.6g}",
        "clash_count_lt_2p0": counts[2.0],
        "clash_count_lt_2p2": counts[2.2],
        "clash_count_lt_2p5": counts[2.5],
    }


def write_xyz(path: Path, heavy_deduped: list, source: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f"{len(heavy_deduped)}\n")
        handle.write(f"Heavy atoms only, exact deduplication applied from {display_path(source)}.\n")
        for atom in heavy_deduped:
            handle.write(f"{atom.element.upper():<6}{atom.x:12.6f}{atom.y:12.6f}{atom.z:12.6f}\n")


def validate_candidates() -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    structural_rows = []
    candidates = []
    per_twist_counter: Counter[str] = Counter()
    for path in candidate_paths():
        twist = twist_from_path(path)
        per_twist_counter[twist] += 1
        label = f"twist_{twist}_candidate_{per_twist_counter[twist]:02d}_{safe_label(path.name)}"
        atoms = load_pdb_atoms(path)
        heavy = heavy_atoms(atoms)
        deduped = dedupe_exact_atoms(heavy)
        ready, note = is_diffraction_ready(atoms)
        element_counts = Counter(atom.element.upper() for atom in atoms)
        row = {
            "source_folder": twist,
            "candidate_file": display_path(path),
            "twist_label_deg": twist,
            "candidate_label": label,
            "diffraction_ready": "yes" if ready else "no",
            "readiness_note": note,
            "atom_count": len(atoms),
            "residue_count": residue_count(atoms),
            "residue_names": ";".join(residue_names(atoms)),
            "element_counts": ";".join(f"{k}:{v}" for k, v in sorted(element_counts.items())),
            "heavy_atom_count": len(heavy),
            "exact_heavy_duplicate_count": len(heavy) - len(deduped),
            "limitations": "Layer/twist statistics are rough centroid proxies; H-bond screen is geometric only; no minimization was run.",
        }
        row.update(layer_stats(atoms))
        row.update(hbond_proxy(atoms))
        row.update(clash_proxy(deduped))
        xyz_path = XYZ_ROOT / f"{label}.xyz"
        if ready:
            write_xyz(xyz_path, deduped, path)
            candidates.append({"path": path, "xyz": xyz_path, "label": label, "twist": twist, "structural": row})
        row["xyz_output"] = display_path(xyz_path) if ready else ""
        structural_rows.append(row)
    return structural_rows, candidates


def run_diffraction(args: argparse.Namespace, candidates: list[dict[str, object]], nick) -> dict[str, list[dict[str, float]]]:
    scripts_module, orientation_module, radial_module = nick.load_engine()
    rotations = list(range(0, 181, 5))
    profiles = {}
    helper_args = SimpleNamespace(
        grid_size=args.grid_size,
        grid_limit_mm=args.grid_limit_mm,
        radial_bins=args.radial_bins,
        detector_distance_mm=args.detector_distance_mm,
        wavelength_angstrom=args.wavelength_angstrom,
        normalize_q_min=args.normalize_q_min,
    )
    for cand in candidates:
        _atoms, coords_angstrom, atomic_numbers = nick.load_xyz(cand["xyz"], scripts_module.atomic_number)
        wavelength_mm = args.wavelength_angstrom * 1e-7
        coords_mm = coords_angstrom * 1e-7
        z_grid_limits = (-args.grid_limit_mm, args.grid_limit_mm)
        x_grid_limits = (-args.grid_limit_mm, args.grid_limit_mm)
        image = orientation_module.average_fiber_diffraction(
            atomic_numbers,
            coords_mm,
            wavelength_mm,
            args.detector_distance_mm,
            z_grid_limits,
            x_grid_limits,
            args.grid_size,
            args.grid_size,
            [0],
            rotations,
        )
        profile = nick.radial_rows_from_image(image, helper_args, radial_module)
        radial_path = RADIAL_DIR / f"{cand['label']}_radial_profile.csv"
        nick.write_profile(radial_path, profile)
        cand["radial_profile"] = radial_path
        profiles[cand["label"]] = profile
    return profiles


def feature_rows(candidates: list[dict[str, object]], profiles: dict[str, list[dict[str, float]]], experimental: list[dict[str, float]], nick) -> list[dict[str, object]]:
    rows = []
    for cand in candidates:
        profile = profiles[cand["label"]]
        for feature, _center, d_min, d_max in FEATURE_WINDOWS:
            sim_peak_d, sim_peak_i, sim_area = nick.profile_peak_and_area(profile, d_min, d_max)
            exp_peak_d, _exp_peak_i, exp_area = nick.profile_peak_and_area(experimental, d_min, d_max)
            offset = sim_peak_d - exp_peak_d if not (math.isnan(sim_peak_d) or math.isnan(exp_peak_d)) else math.nan
            rows.append({
                "source_folder": cand["twist"],
                "candidate_file": display_path(cand["path"]),
                "twist_label_deg": cand["twist"],
                "candidate_label": cand["label"],
                "feature_window": feature,
                "experimental_peak_d_A": "" if math.isnan(exp_peak_d) else f"{exp_peak_d:.12g}",
                "simulated_peak_d_A": "" if math.isnan(sim_peak_d) else f"{sim_peak_d:.12g}",
                "peak_offset_d_A": "" if math.isnan(offset) else f"{offset:.12g}",
                "abs_peak_offset_d_A": "" if math.isnan(offset) else f"{abs(offset):.12g}",
                "simulated_peak_intensity_norm": "" if math.isnan(sim_peak_i) else f"{sim_peak_i:.12g}",
                "window_area_simulated": "" if math.isnan(sim_area) else f"{sim_area:.12g}",
                "window_area_experimental": "" if math.isnan(exp_area) else f"{exp_area:.12g}",
            })
    return rows


def ranking_rows(candidates: list[dict[str, object]], features: list[dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    by_label: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in features:
        by_label[str(row["candidate_label"])].append(row)
    for cand in candidates:
        offsets = []
        feature_offsets = {}
        for row in by_label[cand["label"]]:
            value = row["peak_offset_d_A"]
            if value != "":
                offset = float(value)
                offsets.append(abs(offset))
                feature_offsets[row["feature_window"]] = offset
        structural = cand["structural"]
        mean_offset = float(np.mean(offsets)) if offsets else math.nan
        max_offset = float(np.max(offsets)) if offsets else math.nan
        clash_22 = int(structural["clash_count_lt_2p2"])
        hbond = structural["hbond_proxy_best_n_o_distance_A"]
        note_parts = []
        if not math.isnan(mean_offset) and mean_offset < 0.08:
            note_parts.append("strong primary-window match")
        elif not math.isnan(mean_offset) and mean_offset < 0.15:
            note_parts.append("moderate primary-window match")
        else:
            note_parts.append("weaker primary-window match")
        if hbond != "" and float(hbond) <= 3.2:
            note_parts.append("plausible N-O H-bond proxy")
        if clash_22:
            note_parts.append(f"{clash_22} heavy-heavy contacts below 2.2 A")
        rows.append({
            "source_folder": cand["twist"],
            "candidate_file": display_path(cand["path"]),
            "twist_label_deg": cand["twist"],
            "candidate_label": cand["label"],
            "mean_abs_primary_peak_offset_A": "" if math.isnan(mean_offset) else f"{mean_offset:.12g}",
            "max_abs_primary_peak_offset_A": "" if math.isnan(max_offset) else f"{max_offset:.12g}",
            "base_stack_peak_offset_A": f"{feature_offsets.get('3.38/3.4 A', math.nan):.12g}",
            "feature_3p77_offset_A": f"{feature_offsets.get('3.77 A', math.nan):.12g}",
            "feature_4p4_offset_A": f"{feature_offsets.get('4.4 A', math.nan):.12g}",
            "feature_5p6_offset_A": f"{feature_offsets.get('5.6 A', math.nan):.12g}",
            "feature_7p3_offset_A": f"{feature_offsets.get('7.3 A', math.nan):.12g}",
            "hbond_proxy_best_distance_A": hbond,
            "clash_count_lt_2p2": clash_22,
            "qualitative_note": "; ".join(note_parts),
        })
    rows.sort(key=lambda row: float(row["mean_abs_primary_peak_offset_A"]) if row["mean_abs_primary_peak_offset_A"] else math.inf)
    return rows


def plot_profiles(path: Path, title: str, profiles: dict[str, list[dict[str, float]]], experimental: list[dict[str, float]], labels: list[str]) -> None:
    fig, ax = plt.subplots(figsize=(9.5, 5.5), dpi=180)
    exp_rows = sorted([row for row in experimental if 2.8 <= row["d_A"] <= 8.4], key=lambda row: row["d_A"])
    ax.plot([row["d_A"] for row in exp_rows], [row["intensity_norm"] for row in exp_rows], color="black", lw=1.4, label="Corrected experimental")
    colors = ["#2166ac", "#b2182b", "#1b7837", "#762a83", "#e08214", "#35978f"]
    for idx, label in enumerate(labels):
        rows = sorted([row for row in profiles[label] if 2.8 <= row["d_A"] <= 8.4], key=lambda row: row["d_A"])
        ax.plot([row["d_A"] for row in rows], [row["intensity_norm"] + 1.05 + idx * 0.82 for row in rows], lw=1.1, color=colors[idx % len(colors)], label=label)
    for _feature, _center, d_min, d_max in FEATURE_WINDOWS:
        ax.axvspan(d_min, d_max, color="#d8d8d8", alpha=0.12)
    ax.set_xlim(8.4, 2.8)
    ax.set_xlabel("d spacing (A)")
    ax.set_ylabel("Normalized intensity, offset for display")
    ax.set_title(title)
    ax.legend(fontsize=7)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)


def gaussian_smooth_d(rows: list[dict[str, float]], sigma: float = 0.10) -> tuple[np.ndarray, np.ndarray]:
    subset = sorted([row for row in rows if 2.8 <= row["d_A"] <= 8.4], key=lambda row: row["d_A"])
    d = np.asarray([row["d_A"] for row in subset], dtype=float)
    y = np.asarray([row["intensity_norm"] for row in subset], dtype=float)
    out = np.zeros_like(y)
    for i, x in enumerate(d):
        w = np.exp(-0.5 * ((d - x) / sigma) ** 2)
        out[i] = float(np.sum(w * y) / np.sum(w)) if np.sum(w) else y[i]
    if len(out) and out.max() > 0:
        out = out / out.max()
    return d, out


def make_plots(candidates: list[dict[str, object]], profiles: dict[str, list[dict[str, float]]], experimental: list[dict[str, float]], ranking: list[dict[str, object]]) -> None:
    labels = [cand["label"] for cand in candidates]
    plot_profiles(PLOT_DIR / "corrected_experimental_vs_29_30_31_overlay.png", "Corrected experimental vs Asem twist candidates", profiles, experimental, labels)
    labels_30 = [cand["label"] for cand in candidates if cand["twist"] == "30"]
    if len(labels_30) > 1:
        plot_profiles(PLOT_DIR / "corrected_experimental_vs_all_30deg_candidates_overlay.png", "Corrected experimental vs all 30 degree candidates", profiles, experimental, labels_30)
    elif len(labels_30) == 1:
        plot_profiles(PLOT_DIR / "corrected_experimental_vs_30deg_candidate_overlay.png", "Corrected experimental vs 30 degree candidate", profiles, experimental, labels_30)

    fig, ax = plt.subplots(figsize=(7.4, 4.6), dpi=180)
    ax.bar([row["candidate_label"] for row in ranking], [float(row["mean_abs_primary_peak_offset_A"]) for row in ranking], color="#4d9221")
    ax.set_ylabel("Mean abs primary peak offset (A)")
    ax.tick_params(axis="x", rotation=35, labelsize=7)
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "mean_peak_offset_by_twist_candidate.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.5, 4.9), dpi=180)
    x = np.arange(len(ranking))
    width = 0.15
    fields = ["base_stack_peak_offset_A", "feature_3p77_offset_A", "feature_4p4_offset_A", "feature_5p6_offset_A", "feature_7p3_offset_A"]
    for i, field in enumerate(fields):
        ax.bar(x + (i - 2) * width, [float(row[field]) for row in ranking], width=width, label=field.replace("_offset_A", ""))
    ax.axhline(0, color="black", lw=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels([row["candidate_label"] for row in ranking], rotation=35, ha="right", fontsize=7)
    ax.set_ylabel("Peak offset (A)")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "feature_specific_offsets_by_twist_candidate.png")
    plt.close(fig)

    best_by_twist = []
    seen = set()
    for row in ranking:
        if row["twist_label_deg"] not in seen:
            best_by_twist.append(row["candidate_label"])
            seen.add(row["twist_label_deg"])
    fig, ax = plt.subplots(figsize=(9.5, 5.5), dpi=180)
    exp_d, exp_y = gaussian_smooth_d(experimental, sigma=0.10)
    ax.plot(exp_d, exp_y, color="black", lw=1.5, label="Corrected experimental, smoothed")
    colors = ["#2166ac", "#b2182b", "#1b7837"]
    for idx, label in enumerate(best_by_twist):
        d, y = gaussian_smooth_d(profiles[label], sigma=0.10)
        ax.plot(d, y + 1.05 + idx * 0.85, color=colors[idx % len(colors)], lw=1.2, label=f"{label}, smoothed")
    ax.set_xlim(8.4, 2.8)
    ax.set_xlabel("d spacing (A)")
    ax.set_ylabel("Normalized intensity, offset for display")
    ax.set_title("Nick-style d-space Gaussian smoothing, sigma 0.10 A, visualization only")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "nick_style_smoothed_best_by_twist_visualization_only.png")
    plt.close(fig)


def markdown_table(rows: list[dict[str, object]], fields: list[str]) -> str:
    if not rows:
        return "_No rows._"
    lines = ["| " + " | ".join(fields) + " |", "| " + " | ".join(["---"] * len(fields)) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(str(row.get(field, "")) for field in fields) + " |")
    return "\n".join(lines)


def write_report(structural: list[dict[str, object]], candidates: list[dict[str, object]], ranking: list[dict[str, object]], args: argparse.Namespace) -> None:
    best = ranking[0] if ranking else None
    best_twist = str(best["twist_label_deg"]) if best else "none"
    tested = sorted({cand["twist"] for cand in candidates})
    candidates_30 = [cand for cand in candidates if cand["twist"] == "30"]
    stable_30 = "Only one 30 degree candidate was present, so cross-candidate stability could not be tested."
    if len(candidates_30) > 1:
        offsets = [float(row["mean_abs_primary_peak_offset_A"]) for row in ranking if row["twist_label_deg"] == "30"]
        stable_30 = f"Multiple 30 degree candidates were present; mean offsets span {min(offsets):.4g} to {max(offsets):.4g} A."
    hbond_lines = [row for row in structural if row["diffraction_ready"] == "yes"]
    best_hbond = min((float(row["hbond_proxy_best_n_o_distance_A"]) for row in hbond_lines if row["hbond_proxy_best_n_o_distance_A"]), default=math.nan)
    excluded = [row for row in structural if row["diffraction_ready"] != "yes"]
    lines = [
        "# Asem Twist Series Corrected Emory Profile Benchmark",
        "",
        "## Purpose",
        "",
        "Benchmark Asem's provided 29, 30, and 31 degree twist models against the Emory-corrected powder profile using the corrected non-accumulating/vectorized diffraction workflow.",
        "",
        "## Asem's Note",
        "",
        "Asem attached models for 29, 30, and 31 degree twist, plus generation inputs for handedness and parallel/antiparallel variants when present. Future alternative structures can be generated by editing `options.yaml`, running `run.py`, and then running `tleap -f initial.in` with AmberTools to build missing carboxylate groups. This first pass does not execute those generation steps.",
        "",
        "## Imported Source Folders",
        "",
        "- `inputs/asem_twist_series_29_30_31/raw/29`",
        "- `inputs/asem_twist_series_29_30_31/raw/30`",
        "- `inputs/asem_twist_series_29_30_31/raw/31`",
        "",
        "## Inventory Summary",
        "",
        markdown_table(structural, ["source_folder", "candidate_file", "atom_count", "residue_count", "residue_names", "diffraction_ready", "readiness_note"]),
        "",
        "## Diffraction-Ready Selection",
        "",
        f"Tested twist folders: {', '.join(tested) if tested else 'none'}.",
        "",
        "All diffraction-ready files were full-scale PDB-like files with CYP/MEP/GLU residues and were converted using heavy-atom exclusion of hydrogens plus exact heavy-atom deduplication. Excluded files are listed below if present.",
        "",
        markdown_table(excluded, ["source_folder", "candidate_file", "readiness_note"]) if excluded else "_No PDB candidates were excluded as local/intermediate._",
        "",
        "## Structural Validation Summary",
        "",
        markdown_table(structural, ["candidate_label", "heavy_atom_count", "exact_heavy_duplicate_count", "hbond_proxy_best_n_o_distance_A", "hbond_proxy_best_h_o_distance_A", "min_heavy_heavy_nonbonded_distance_A", "clash_count_lt_2p2", "rise_median_A", "rough_twist_median_deg"]),
        "",
        "## Corrected Diffraction Settings",
        "",
        "- Corrected Asem non-accumulating/vectorized diffraction path from `reference/asem_corrected_diffraction_engine`.",
        "- Tilts: `[0]`.",
        "- Rotations: `range(0, 181, 5)`.",
        "- Hydrogens excluded; exact heavy-atom deduplication applied.",
        f"- Grid size: {args.grid_size} x {args.grid_size}; grid limit: {args.grid_limit_mm:g} mm; radial bins: {args.radial_bins}.",
        f"- Detector distance: {args.detector_distance_mm:g} mm; wavelength: {args.wavelength_angstrom:g} A.",
        "",
        "## Ranking",
        "",
        markdown_table(ranking, ["candidate_label", "twist_label_deg", "mean_abs_primary_peak_offset_A", "max_abs_primary_peak_offset_A", "base_stack_peak_offset_A", "feature_3p77_offset_A", "feature_4p4_offset_A", "feature_5p6_offset_A", "feature_7p3_offset_A", "hbond_proxy_best_distance_A", "clash_count_lt_2p2", "qualitative_note"]),
        "",
        "## Conclusions",
        "",
        f"- Best candidate by mean primary-window peak offset: `{best['candidate_label']}` ({best['mean_abs_primary_peak_offset_A']} A) from twist {best_twist}." if best else "- No diffraction-ready candidates were available.",
        f"- Does 30 degree remain best relative to 29 and 31 degree? {'yes' if best_twist == '30' else 'no'}.",
        f"- 30 degree candidate stability: {stable_30}",
        f"- Best backbone-backbone H-bond proxy N-O distance among validated candidates: {best_hbond:.4g} A." if not math.isnan(best_hbond) else "- H-bond proxy could not be evaluated.",
        "- Candidates with plausible N-O distance proxies preserve the diffraction comparison only to the extent shown in the ranking table; this is not energetic validation.",
        "- Recommendation: use the top-ranked provided full model as the next structure to discuss with Nick/Asem, while treating clash and H-bond proxy results as screening flags.",
        "",
        "## Plots",
        "",
        "- `outputs/asem_twist_series_corrected_emory_profile/plots/corrected_experimental_vs_29_30_31_overlay.png`",
        "- `outputs/asem_twist_series_corrected_emory_profile/plots/mean_peak_offset_by_twist_candidate.png`",
        "- `outputs/asem_twist_series_corrected_emory_profile/plots/feature_specific_offsets_by_twist_candidate.png`",
        "- `outputs/asem_twist_series_corrected_emory_profile/plots/nick_style_smoothed_best_by_twist_visualization_only.png`",
        "",
        "## Limitations",
        "",
        "- This does not generate new structures.",
        "- This does not relax energy filters.",
        "- This does not run AmberTools/tleap.",
        "- This does not run pNAB, execute notebooks, or run minimization.",
        "- Powder/radial comparison is not full structural refinement.",
        "",
        "## Future Work",
        "",
        "After benchmarking the provided 29/30/31 models, a later task can relax energy thresholds in `options.yaml`, generate multiple 30 degree candidates, run `tleap -f initial.in` if AmberTools is available, and rerun this same structural/diffraction screen.",
    ]
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    for directory in (XYZ_ROOT, RADIAL_DIR, PLOT_DIR, TABLE_DIR, METRICS_DIR, REPORT.parent):
        directory.mkdir(parents=True, exist_ok=True)
    nick = load_nick_runner()
    structural, candidates = validate_candidates()
    structural_fields = [
        "source_folder", "candidate_file", "twist_label_deg", "candidate_label", "diffraction_ready", "readiness_note", "atom_count", "residue_count", "residue_names", "element_counts", "heavy_atom_count", "exact_heavy_duplicate_count", "hbond_proxy_best_n_o_distance_A", "hbond_proxy_best_h_o_distance_A", "hbond_proxy_best_n_h_o_angle_deg", "min_heavy_heavy_nonbonded_distance_A", "clash_count_lt_2p0", "clash_count_lt_2p2", "clash_count_lt_2p5", "rise_median_A", "rise_mean_A", "rough_twist_median_deg", "rough_twist_mean_deg", "xyz_output", "limitations",
    ]
    write_csv(STRUCTURAL_CSV, structural, structural_fields)
    profiles = run_diffraction(args, candidates, nick) if candidates else {}
    experimental = nick.read_experimental_profile(EXPERIMENTAL)
    features = feature_rows(candidates, profiles, experimental, nick)
    feature_fields = ["source_folder", "candidate_file", "twist_label_deg", "candidate_label", "feature_window", "experimental_peak_d_A", "simulated_peak_d_A", "peak_offset_d_A", "abs_peak_offset_d_A", "simulated_peak_intensity_norm", "window_area_simulated", "window_area_experimental"]
    write_csv(FEATURE_CSV, features, feature_fields)
    write_csv(TABLE_DIR / FEATURE_CSV.name, features, feature_fields)
    ranking = ranking_rows(candidates, features)
    ranking_fields = ["source_folder", "candidate_file", "twist_label_deg", "candidate_label", "mean_abs_primary_peak_offset_A", "max_abs_primary_peak_offset_A", "base_stack_peak_offset_A", "feature_3p77_offset_A", "feature_4p4_offset_A", "feature_5p6_offset_A", "feature_7p3_offset_A", "hbond_proxy_best_distance_A", "clash_count_lt_2p2", "qualitative_note"]
    write_csv(RANKING_CSV, ranking, ranking_fields)
    write_csv(TABLE_DIR / RANKING_CSV.name, ranking, ranking_fields)
    make_plots(candidates, profiles, experimental, ranking)
    metadata = {
        "candidate_count": len(candidates),
        "experimental_profile": display_path(EXPERIMENTAL),
        "rotation_sampling": "range(0, 181, 5)",
        "tilts": [0],
        "hydrogens_excluded": True,
        "exact_heavy_atom_deduplication": True,
        "notebooks_executed": False,
        "pnab_run": False,
        "ambertools_tleap_run": False,
        "minimization_run": False,
    }
    (TABLE_DIR / "run_metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_report(structural, candidates, ranking, args)
    print(f"Wrote structural summary: {STRUCTURAL_CSV}")
    print(f"Wrote feature summary: {FEATURE_CSV}")
    print(f"Wrote ranking: {RANKING_CSV}")
    print(f"Wrote report: {REPORT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

