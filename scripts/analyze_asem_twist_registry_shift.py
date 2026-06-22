#!/usr/bin/env python3
"""Audit whether Asem 31 degree shifts inter-hexad registry versus 29/30."""

from __future__ import annotations

import csv
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from hexaplex_formation.pdb_utils import heavy_atoms, is_hydrogen, load_pdb_atoms  # noqa: E402

INPUTS = {
    "29": ROOT / "inputs" / "asem_twist_series_29_30_31" / "raw" / "29" / "29.pdb.txt",
    "30": ROOT / "inputs" / "asem_twist_series_29_30_31" / "raw" / "30" / "30.pdb.txt",
    "31": ROOT / "inputs" / "asem_twist_series_29_30_31" / "raw" / "31" / "31.pdb.txt",
}
METRICS = ROOT / "outputs" / "metrics"
PLOTS = ROOT / "outputs" / "asem_twist_registry_shift" / "plots"
REPORT = ROOT / "outputs" / "reports" / "asem_twist_registry_shift_hypothesis_report.md"
FEATURE_CSV = METRICS / "asem_twist_series_feature_summary_corrected_emory_profile.csv"
FEATURE_DELTA_CSV = METRICS / "asem_twist_registry_feature_delta.csv"
PAIR_SUMMARY_CSV = METRICS / "asem_twist_registry_pair_distance_summary.csv"
SHIFT_CSV = METRICS / "asem_twist_registry_31_specific_shifts.csv"
SHORT_CONTACT_CSV = METRICS / "asem_twist_registry_short_contact_localization.csv"

BANDS = [
    ("base_stack_3p2_3p6", 3.2, 3.6),
    ("feature_3p77_3p6_3p95", 3.6, 3.95),
    ("feature_4p4_4p15_4p65", 4.15, 4.65),
    ("feature_5p6_5p35_5p85", 5.35, 5.85),
    ("feature_7p3_7p0_7p6", 7.0, 7.6),
]
FEATURE_ORDER = ["3.38/3.4 A", "3.77 A", "4.4 A", "5.6 A", "7.3 A"]
CLASS_ORDER = ["base_core", "backbone", "side_chain_or_carboxylate", "unknown"]
BACKBONE_NAMES = {"N", "CA", "C", "O", "OXT"}
BASE_CORE_NAMES = {
    "N1",
    "C2",
    "OC2",
    "N3",
    "C4",
    "N5",
    "C6",
    "OC6",
    "OC4",
    "N7",
    "C8",
    "N9",
    "C5",
}
SIDE_CHAIN_NAMES = {"CB", "CG", "CD", "OE1", "OE2", "OD1", "OD2", "CC"}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


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


def atom_class(atom) -> str:
    name = atom.atom_name.strip().upper()
    residue = atom.residue_name.strip().upper()
    if is_hydrogen(atom):
        return "hydrogen"
    if name in BACKBONE_NAMES:
        return "backbone"
    if residue == "GLU" or name in SIDE_CHAIN_NAMES or name.startswith(("OE", "OD")):
        return "side_chain_or_carboxylate"
    if residue in {"CYP", "MEP"} and (name in BASE_CORE_NAMES or any(char.isdigit() for char in name)):
        return "base_core"
    return "unknown"


def residue_key(atom) -> tuple[str, int | None, str, str]:
    return (atom.chain_id, atom.residue_number, atom.insertion_code, atom.residue_name)


def atom_id(atom) -> str:
    chain = atom.chain_id or "."
    resnum = "" if atom.residue_number is None else str(atom.residue_number)
    return f"{atom.record_type}:{atom.atom_serial}:{chain}:{atom.residue_name}{resnum}:{atom.atom_name}"


def coords(atoms: list) -> np.ndarray:
    return np.asarray([[atom.x, atom.y, atom.z] for atom in atoms], dtype=float)


def layer_groups(atoms: list) -> dict[int, list]:
    residues: dict[tuple[str, int | None, str, str], list] = defaultdict(list)
    for atom in atoms:
        residues[residue_key(atom)].append(atom)
    ordered = sorted(residues, key=lambda key: (key[1] if key[1] is not None else -1, key[3], key[0], key[2]))
    groups: dict[int, list] = defaultdict(list)
    for idx, key in enumerate(ordered):
        layer = idx // 6
        groups[layer].extend(residues[key])
    return dict(groups)


def class_pair(a: str, b: str) -> tuple[str, str]:
    return tuple(sorted((a, b)))


def hist_peak(values: list[float], lo: float, hi: float) -> float | None:
    if not values:
        return None
    bins = np.arange(lo, hi + 0.0201, 0.02)
    counts, edges = np.histogram(values, bins=bins)
    if counts.size == 0 or counts.max() == 0:
        return None
    idx = int(np.argmax(counts))
    return float((edges[idx] + edges[idx + 1]) / 2.0)


def pair_distance_summary_for_model(twist: str, atoms: list) -> tuple[list[dict[str, object]], dict[tuple[str, str, str], list[float]]]:
    layers = layer_groups(heavy_atoms(atoms))
    values: dict[tuple[str, str, str], list[float]] = defaultdict(list)
    for layer in range(max(layers) if layers else -1):
        a_atoms = layers.get(layer, [])
        b_atoms = layers.get(layer + 1, [])
        if not a_atoms or not b_atoms:
            continue
        a_coords = coords(a_atoms)
        b_coords = coords(b_atoms)
        a_classes = [atom_class(atom) for atom in a_atoms]
        b_classes = [atom_class(atom) for atom in b_atoms]
        for i, atom_a in enumerate(a_atoms):
            d = np.linalg.norm(b_coords - a_coords[i], axis=1)
            for band, lo, hi in BANDS:
                mask = (d >= lo) & (d <= hi)
                if not np.any(mask):
                    continue
                for j in np.where(mask)[0]:
                    ca, cb = class_pair(a_classes[i], b_classes[j])
                    values[("adjacent_layers", band, f"{ca}|{cb}")].append(float(d[j]))
    rows = []
    for (scope, band, pair), distances in sorted(values.items()):
        arr = np.asarray(distances, dtype=float)
        lo, hi = next((lo, hi) for name, lo, hi in BANDS if name == band)
        a_class, b_class = pair.split("|")
        rows.append(
            {
                "twist_deg": twist,
                "pair_scope": scope,
                "distance_band": band,
                "atom_class_a": a_class,
                "atom_class_b": b_class,
                "pair_count": len(distances),
                "mean_distance_A": fmt(float(np.mean(arr))),
                "median_distance_A": fmt(float(np.median(arr))),
                "p10_distance_A": fmt(float(np.percentile(arr, 10))),
                "p90_distance_A": fmt(float(np.percentile(arr, 90))),
                "hist_peak_distance_A": fmt(hist_peak(distances, lo, hi)),
                "delta_hist_peak_vs_30_A": "",
                "delta_pair_count_vs_30": "",
                "note": "Adjacent-layer heavy-atom pairs; atom classes are heuristic.",
            }
        )
    return rows, values


def add_pair_deltas(rows: list[dict[str, object]]) -> None:
    by_key: dict[tuple[str, str, str, str], dict[str, dict[str, object]]] = defaultdict(dict)
    for row in rows:
        key = (str(row["pair_scope"]), str(row["distance_band"]), str(row["atom_class_a"]), str(row["atom_class_b"]))
        by_key[key][str(row["twist_deg"])] = row
    for group in by_key.values():
        ref = group.get("30")
        if not ref:
            continue
        ref_peak = f(ref.get("hist_peak_distance_A"))
        ref_count = int(ref.get("pair_count", 0))
        for twist, row in group.items():
            peak = f(row.get("hist_peak_distance_A"))
            row["delta_hist_peak_vs_30_A"] = fmt(None if peak is None or ref_peak is None else peak - ref_peak)
            row["delta_pair_count_vs_30"] = int(row.get("pair_count", 0)) - ref_count


def feature_delta_rows() -> list[dict[str, object]]:
    by_feature: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    for row in read_csv(FEATURE_CSV):
        by_feature[row["feature_window"]][row["twist_label_deg"]] = row
    out = []
    for feature in FEATURE_ORDER:
        rows = by_feature.get(feature, {})
        offsets = {twist: f(rows.get(twist, {}).get("peak_offset_d_A")) for twist in ("29", "30", "31")}
        abs_offsets = {twist: abs(value) if value is not None else None for twist, value in offsets.items()}
        best_abs = min(value for value in abs_offsets.values() if value is not None)
        winner = ",".join(twist for twist, value in abs_offsets.items() if value is not None and abs(value - best_abs) < 1e-9)
        delta_31 = None if abs_offsets["31"] is None else abs_offsets["31"] - best_abs
        drives = delta_31 is not None and delta_31 > 0.005
        note = "31 differs from the best/tied basin in this window." if drives else "31 does not drive the disfavoring in this window."
        out.append(
            {
                "feature_window": feature,
                "offset_29_A": fmt(offsets["29"]),
                "offset_30_A": fmt(offsets["30"]),
                "offset_31_A": fmt(offsets["31"]),
                "abs_offset_29_A": fmt(abs_offsets["29"]),
                "abs_offset_30_A": fmt(abs_offsets["30"]),
                "abs_offset_31_A": fmt(abs_offsets["31"]),
                "delta_31_minus_best_abs_offset_A": fmt(delta_31),
                "winner": winner,
                "does_31_drive_failure": "yes" if drives else "no",
                "note": note,
            }
        )
    return out


def shifted_family_rows(summary_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    by_key: dict[tuple[str, str, str], dict[str, dict[str, object]]] = defaultdict(dict)
    for row in summary_rows:
        key = (str(row["distance_band"]), f"{row['atom_class_a']}|{row['atom_class_b']}", str(row["pair_scope"]))
        by_key[key][str(row["twist_deg"])] = row
    shifts = []
    for (band, pair, _scope), group in sorted(by_key.items()):
        if not all(twist in group for twist in ("29", "30", "31")):
            continue
        for metric in ("hist_peak_distance_A", "pair_count", "median_distance_A"):
            v29 = f(group["29"].get(metric))
            v30 = f(group["30"].get(metric))
            v31 = f(group["31"].get(metric))
            if v29 is None or v30 is None or v31 is None:
                continue
            d31_30 = v31 - v30
            d31_29 = v31 - v29
            similar_29_30 = abs(v29 - v30) <= (0.03 if metric != "pair_count" else max(5.0, 0.10 * max(v29, v30)))
            shifted_31 = abs(d31_30) >= (0.06 if metric != "pair_count" else max(10.0, 0.20 * max(v30, 1.0)))
            supports = similar_29_30 and shifted_31
            if supports or (metric == "pair_count" and abs(d31_30) >= max(20.0, 0.30 * max(v30, 1.0))):
                interp = (
                    "29 and 30 are similar while 31 changes, consistent with a 31-specific registry shift."
                    if supports
                    else "31 changes strongly, but 29/30 similarity is weaker; treat as contextual."
                )
                shifts.append(
                    {
                        "distance_band": band,
                        "atom_class_pair": pair,
                        "metric": metric,
                        "value_29": fmt(v29),
                        "value_30": fmt(v30),
                        "value_31": fmt(v31),
                        "delta_31_vs_30": fmt(d31_30),
                        "delta_31_vs_29": fmt(d31_29),
                        "supports_overrotation_hypothesis": "yes" if supports else "contextual",
                        "interpretation": interp,
                    }
                )
    shifts.sort(key=lambda row: (row["supports_overrotation_hypothesis"] != "yes", row["distance_band"], row["atom_class_pair"], row["metric"]))
    return shifts


def short_contacts_for_model(twist: str, atoms: list) -> list[dict[str, object]]:
    heavy = heavy_atoms(atoms)
    layer_by_residue = {}
    for layer, layer_atoms in layer_groups(heavy).items():
        for atom in layer_atoms:
            layer_by_residue[residue_key(atom)] = layer
    c = coords(heavy)
    classes = [atom_class(atom) for atom in heavy]
    rows = []
    for i, atom_a in enumerate(heavy):
        d = np.linalg.norm(c[i + 1 :] - c[i], axis=1)
        for local_j in np.where(d < 2.5)[0]:
            j = i + 1 + int(local_j)
            atom_b = heavy[j]
            same_res = residue_key(atom_a) == residue_key(atom_b)
            layer_a = layer_by_residue.get(residue_key(atom_a))
            layer_b = layer_by_residue.get(residue_key(atom_b))
            if same_res:
                continue
            if layer_a == layer_b:
                relation = "same_layer"
            elif layer_a is not None and layer_b is not None and abs(layer_a - layer_b) == 1:
                relation = "adjacent_layer"
            else:
                relation = "other"
            ca, cb = class_pair(classes[i], classes[j])
            dist = float(d[local_j])
            rows.append(
                {
                    "twist_deg": twist,
                    "threshold_A": "2.0" if dist < 2.0 else "2.2" if dist < 2.2 else "2.5",
                    "atom_a": atom_id(atom_a),
                    "atom_b": atom_id(atom_b),
                    "atom_class_a": ca,
                    "atom_class_b": cb,
                    "atom_class_pair": f"{ca}|{cb}",
                    "distance_A": fmt(dist),
                    "same_residue": "yes" if same_res else "no",
                    "layer_a": "" if layer_a is None else layer_a,
                    "layer_b": "" if layer_b is None else layer_b,
                    "layer_relation": relation,
                    "enrichment_note": "",
                }
            )
    return rows


def annotate_short_contact_enrichment(rows: list[dict[str, object]]) -> None:
    counts = Counter((row["twist_deg"], row["threshold_A"], row["atom_class_pair"], row["layer_relation"]) for row in rows)
    for row in rows:
        key_base = (row["threshold_A"], row["atom_class_pair"], row["layer_relation"])
        c29 = counts.get(("29", *key_base), 0)
        c30 = counts.get(("30", *key_base), 0)
        c31 = counts.get(("31", *key_base), 0)
        if row["twist_deg"] == "31" and c31 > max(c29, c30) * 1.2 + 5:
            row["enrichment_note"] = f"31 enriched for this contact class/relation: 29={c29}, 30={c30}, 31={c31}"


def plot_feature_delta(rows: list[dict[str, object]]) -> None:
    x = np.arange(len(rows))
    fig, ax = plt.subplots(figsize=(8.2, 4.7), dpi=180)
    for twist, color in [("29", "#4c78a8"), ("30", "#59a14f"), ("31", "#e15759")]:
        ax.plot(x, [float(row[f"abs_offset_{twist}_A"]) for row in rows], marker="o", label=f"{twist} deg", color=color)
    ax.set_xticks(x)
    ax.set_xticklabels([row["feature_window"] for row in rows], rotation=25, ha="right")
    ax.set_ylabel("Absolute peak offset (A)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(PLOTS / "feature_delta_29_30_31.png")
    plt.close(fig)


def plot_pair_histograms(raw_values: dict[str, dict[tuple[str, str, str], list[float]]], shifts: list[dict[str, object]]) -> None:
    diagnostic = next((row for row in shifts if row["supports_overrotation_hypothesis"] == "yes"), None)
    if diagnostic:
        pair = str(diagnostic["atom_class_pair"])
        key_suffix = ("adjacent_layers", str(diagnostic["distance_band"]), pair)
    else:
        key_suffix = ("adjacent_layers", "feature_3p77_3p6_3p95", "base_core|backbone")
    fig, ax = plt.subplots(figsize=(8, 4.7), dpi=180)
    lo, hi = next((lo, hi) for name, lo, hi in BANDS if name == key_suffix[1])
    for twist, color in [("29", "#4c78a8"), ("30", "#59a14f"), ("31", "#e15759")]:
        vals = raw_values.get(twist, {}).get(key_suffix, [])
        if vals:
            ax.hist(vals, bins=np.arange(lo, hi + 0.025, 0.025), histtype="step", lw=1.4, density=True, label=f"{twist} deg", color=color)
    ax.set_xlabel("Distance (A)")
    ax.set_ylabel("Density")
    ax.set_title(f"{key_suffix[1]} {key_suffix[2]}")
    ax.legend()
    fig.tight_layout()
    fig.savefig(PLOTS / "pair_distance_histograms_by_band.png")
    plt.close(fig)


def plot_shift_heatmap(shifts: list[dict[str, object]]) -> None:
    rows = [row for row in shifts if row["metric"] in {"hist_peak_distance_A", "pair_count"}][:30]
    if not rows:
        rows = [{"distance_band": band, "atom_class_pair": "none", "metric": "none", "delta_31_vs_30": "0"} for band, *_ in BANDS]
    labels = [f"{row['distance_band']}\n{row['atom_class_pair']}\n{row['metric']}" for row in rows]
    values = np.asarray([[float(row["delta_31_vs_30"])] for row in rows], dtype=float)
    fig, ax = plt.subplots(figsize=(7, max(4, len(rows) * 0.22)), dpi=180)
    im = ax.imshow(values, cmap="coolwarm", aspect="auto")
    ax.set_yticks(np.arange(len(rows)))
    ax.set_yticklabels(labels, fontsize=6)
    ax.set_xticks([0])
    ax.set_xticklabels(["31 - 30"])
    fig.colorbar(im, ax=ax)
    fig.tight_layout()
    fig.savefig(PLOTS / "atom_class_pair_shift_heatmap.png")
    plt.close(fig)


def plot_short_contacts(rows: list[dict[str, object]]) -> None:
    counts = Counter((row["twist_deg"], row["atom_class_pair"]) for row in rows if row["threshold_A"] in {"2.0", "2.2"})
    pairs = sorted({pair for _twist, pair in counts})
    x = np.arange(len(pairs))
    fig, ax = plt.subplots(figsize=(8.5, 4.8), dpi=180)
    width = 0.25
    for idx, twist in enumerate(["29", "30", "31"]):
        ax.bar(x + (idx - 1) * width, [counts.get((twist, pair), 0) for pair in pairs], width=width, label=f"{twist} deg")
    ax.set_xticks(x)
    ax.set_xticklabels(pairs, rotation=35, ha="right", fontsize=7)
    ax.set_ylabel("Contacts < 2.2 A")
    ax.legend()
    fig.tight_layout()
    fig.savefig(PLOTS / "short_contact_class_counts.png")
    plt.close(fig)


def plot_registry_summary(feature_rows: list[dict[str, object]], shifts: list[dict[str, object]]) -> None:
    fig, ax = plt.subplots(figsize=(7, 4.5), dpi=180)
    feature_penalty = sum(float(row["delta_31_minus_best_abs_offset_A"]) for row in feature_rows)
    supported = sum(1 for row in shifts if row["supports_overrotation_hypothesis"] == "yes")
    contextual = sum(1 for row in shifts if row["supports_overrotation_hypothesis"] == "contextual")
    ax.bar(["31 feature penalty", "supported shifts", "contextual shifts"], [feature_penalty, supported, contextual], color=["#e15759", "#4c78a8", "#f28e2b"])
    ax.set_title("Registry-shift summary")
    fig.tight_layout()
    fig.savefig(PLOTS / "registry_shift_summary.png")
    plt.close(fig)


def write_report(feature_rows, pair_rows, shift_rows, short_rows) -> None:
    drivers = [row for row in feature_rows if row["does_31_drive_failure"] == "yes"]
    supported = [row for row in shift_rows if row["supports_overrotation_hypothesis"] == "yes"]
    contextual = [row for row in shift_rows if row["supports_overrotation_hypothesis"] == "contextual"]
    verdict = "partially_supported" if supported else "inconclusive" if contextual else "not_supported"
    short_counts = Counter((row["twist_deg"], row["threshold_A"]) for row in short_rows)
    top_supported = supported[:8]
    lines = [
        "# Asem Twist Registry-Shift Hypothesis Report",
        "",
        "## Purpose",
        "",
        "Test whether the 31 degree Asem model is disfavored because over-twist shifts medium-range inter-hexad registry, while 29 and 30 degrees remain in the same radial-distance basin.",
        "",
        "## Inputs",
        "",
        "- `inputs/asem_twist_series_29_30_31/raw/29/29.pdb.txt`",
        "- `inputs/asem_twist_series_29_30_31/raw/30/30.pdb.txt`",
        "- `inputs/asem_twist_series_29_30_31/raw/31/31.pdb.txt`",
        "- Existing corrected diffraction metrics: `outputs/metrics/asem_twist_series_feature_summary_corrected_emory_profile.csv`.",
        "",
        "## Methods",
        "",
        "No diffraction was rerun. The analysis uses existing feature metrics plus direct structural pair-distance audits of the imported PDBs. Layers are inferred as consecutive groups of six residues, giving 30 layers for 29/30 degree models and 32 layers for the 31 degree model. Atom classes are heuristic: backbone (`N`, `CA`, `C`, `O`), base core for CYP/MEP ring-like atoms, side-chain/carboxylate for GLU and carboxylate-like atoms, hydrogen, and unknown.",
        "",
        "Distance bands: 3.2-3.6 A, 3.6-3.95 A, 4.15-4.65 A, 5.35-5.85 A, and 7.0-7.6 A. Pair-distance summaries focus on adjacent layer heavy-atom pairs.",
        "",
        "## Feature Failure",
        "",
        "31 degree is worse mainly where its absolute peak offset exceeds the 29/30 tied/best basin:",
        "",
        "| feature | offset 29 | offset 30 | offset 31 | 31 penalty | drives? |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for row in feature_rows:
        lines.append(
            f"| {row['feature_window']} | {row['offset_29_A']} | {row['offset_30_A']} | {row['offset_31_A']} | {row['delta_31_minus_best_abs_offset_A']} | {row['does_31_drive_failure']} |"
        )
    lines.extend(
        [
            "",
            "## Pair-Distance Findings",
            "",
            f"Flagged 31-specific shifted families: {len(supported)} supported and {len(contextual)} contextual. Top supported examples:",
            "",
            "| band | atom-class pair | metric | 29 | 30 | 31 | delta 31 vs 30 |",
            "| --- | --- | --- | --- | --- | --- | --- |",
        ]
    )
    for row in top_supported:
        lines.append(
            f"| {row['distance_band']} | {row['atom_class_pair']} | {row['metric']} | {row['value_29']} | {row['value_30']} | {row['value_31']} | {row['delta_31_vs_30']} |"
        )
    if not top_supported:
        lines.append("| none | none | none |  |  |  |  |")
    lines.extend(
        [
            "",
            "## Short-Contact Findings",
            "",
            f"Non-same-residue contacts below 2.2 A: 29 degree = {short_counts.get(('29', '2.0'), 0) + short_counts.get(('29', '2.2'), 0)}, 30 degree = {short_counts.get(('30', '2.0'), 0) + short_counts.get(('30', '2.2'), 0)}, 31 degree = {short_counts.get(('31', '2.0'), 0) + short_counts.get(('31', '2.2'), 0)}. These follow the existing structural-warning convention by excluding same-residue covalent-neighbor pairs.",
            "",
            "## Hypothesis Verdict",
            "",
            f"Verdict: `{verdict}`.",
            "",
        ]
    )
    if verdict == "partially_supported":
        lines.append(
            "The feature deltas and pair-distance audit support a cautious interpretation that 31 degree moves one or more medium-range adjacent-layer pair-distance families out of the 29/30 basin. This is consistent with over-rotation disrupting inter-hexad registry, but it does not prove a literal leading/trailing end effect."
        )
    elif verdict == "not_supported":
        lines.append("The audit did not find a clear 31-specific shifted pair-distance family, so the over-rotation/registry-shift hypothesis is not supported by this proxy analysis.")
    else:
        lines.append("The audit finds contextual 31-specific differences but not enough clean 29/30 similarity plus 31 displacement to call the hypothesis supported.")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The results should be interpreted as a proxy-level registry analysis. If supported or partially supported, they suggest asymmetric tolerance: under-twist to 29 degree is tolerated, but over-twist to 31 degree begins to degrade registry. This does not replace Asem/Nick visual chemistry review.",
            "",
            "## Limitations",
            "",
            "- Atom-class assignment is heuristic.",
            "- Pair-distance counts are proxies for diffraction contributions.",
            "- Powder/radial features aggregate many atom pairs.",
            "- This is not an energy calculation.",
            "- No AmberTools/tleap, minimization, notebooks, pNAB, diffraction rerun, or model generation was performed.",
        ]
    )
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    METRICS.mkdir(parents=True, exist_ok=True)
    PLOTS.mkdir(parents=True, exist_ok=True)
    atoms_by_twist = {twist: load_pdb_atoms(path) for twist, path in INPUTS.items()}
    feature_rows = feature_delta_rows()
    write_csv(
        FEATURE_DELTA_CSV,
        feature_rows,
        [
            "feature_window",
            "offset_29_A",
            "offset_30_A",
            "offset_31_A",
            "abs_offset_29_A",
            "abs_offset_30_A",
            "abs_offset_31_A",
            "delta_31_minus_best_abs_offset_A",
            "winner",
            "does_31_drive_failure",
            "note",
        ],
    )
    pair_rows = []
    raw_values = {}
    for twist, atoms in atoms_by_twist.items():
        rows, values = pair_distance_summary_for_model(twist, atoms)
        pair_rows.extend(rows)
        raw_values[twist] = values
    add_pair_deltas(pair_rows)
    write_csv(
        PAIR_SUMMARY_CSV,
        pair_rows,
        [
            "twist_deg",
            "pair_scope",
            "distance_band",
            "atom_class_a",
            "atom_class_b",
            "pair_count",
            "mean_distance_A",
            "median_distance_A",
            "p10_distance_A",
            "p90_distance_A",
            "hist_peak_distance_A",
            "delta_hist_peak_vs_30_A",
            "delta_pair_count_vs_30",
            "note",
        ],
    )
    shifts = shifted_family_rows(pair_rows)
    write_csv(
        SHIFT_CSV,
        shifts,
        [
            "distance_band",
            "atom_class_pair",
            "metric",
            "value_29",
            "value_30",
            "value_31",
            "delta_31_vs_30",
            "delta_31_vs_29",
            "supports_overrotation_hypothesis",
            "interpretation",
        ],
    )
    short_rows = []
    for twist, atoms in atoms_by_twist.items():
        short_rows.extend(short_contacts_for_model(twist, atoms))
    annotate_short_contact_enrichment(short_rows)
    write_csv(
        SHORT_CONTACT_CSV,
        short_rows,
        [
            "twist_deg",
            "threshold_A",
            "atom_a",
            "atom_b",
            "atom_class_a",
            "atom_class_b",
            "atom_class_pair",
            "distance_A",
            "same_residue",
            "layer_a",
            "layer_b",
            "layer_relation",
            "enrichment_note",
        ],
    )
    plot_feature_delta(feature_rows)
    plot_pair_histograms(raw_values, shifts)
    plot_shift_heatmap(shifts)
    plot_short_contacts(short_rows)
    plot_registry_summary(feature_rows, shifts)
    write_report(feature_rows, pair_rows, shifts, short_rows)
    print(f"Wrote {FEATURE_DELTA_CSV}")
    print(f"Wrote {PAIR_SUMMARY_CSV}")
    print(f"Wrote {SHIFT_CSV}")
    print(f"Wrote {SHORT_CONTACT_CSV}")
    print(f"Wrote {REPORT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


