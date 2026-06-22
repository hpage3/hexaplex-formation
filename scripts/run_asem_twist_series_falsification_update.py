#!/usr/bin/env python3
"""Apply a falsification-style screen to Asem's twist-series benchmark."""

from __future__ import annotations

import csv
import math
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "outputs" / "metrics"
PLOTS = ROOT / "outputs" / "asem_twist_series_falsification" / "plots"
REPORT = ROOT / "outputs" / "reports" / "asem_twist_series_falsification_update_report.md"
SUMMARY_CSV = METRICS / "asem_twist_series_falsification_summary.csv"
RANKING_CSV = METRICS / "asem_twist_series_falsification_ranking.csv"

PRIMARY = {
    "3.38/3.4 A": "base_stack_offset_A",
    "3.4 A": "base_stack_offset_A",
    "3.77 A": "feature_3p77_offset_A",
    "3.7 A": "feature_3p77_offset_A",
    "4.4 A": "feature_4p4_offset_A",
    "5.6 A": "feature_5p6_offset_A",
    "7.3 A": "feature_7p3_offset_A",
    "7.25 A": "feature_7p3_offset_A",
}
ORDERED_FIELDS = [
    "base_stack_offset_A",
    "feature_3p77_offset_A",
    "feature_4p4_offset_A",
    "feature_5p6_offset_A",
    "feature_7p3_offset_A",
]


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


def offsets_from_feature_rows(rows: list[dict[str, str]], offset_col: str = "peak_offset_d_A") -> dict[str, float | None]:
    offsets = {field: None for field in ORDERED_FIELDS}
    for row in rows:
        field = PRIMARY.get(row.get("feature_window", ""))
        if field:
            offsets[field] = f(row.get(offset_col))
    return offsets


def summarize_offsets(offsets: dict[str, float | None]) -> tuple[float | None, float | None]:
    values = [abs(v) for v in offsets.values() if v is not None]
    if not values:
        return None, None
    return float(np.mean(values)), float(np.max(values))


def consistency_filter(offsets: dict[str, float | None]) -> bool:
    values = [abs(v) for v in offsets.values() if v is not None]
    if len(values) < 5:
        return False
    mean = float(np.mean(values))
    max_v = float(np.max(values))
    # Reject single-feature chasing: one bad window cannot dominate an otherwise small mean.
    return max_v <= 0.12 and max_v <= max(0.08, 3.0 * mean)


def structural_warning(clash_count: float | None, min_dist: float | None) -> bool:
    return (clash_count is not None and clash_count > 0) or (min_dist is not None and min_dist < 2.0)


def status_for(row: dict[str, object], offsets: dict[str, float | None], role: str, near_best_pass: bool | None = None) -> str:
    if row["is_diffraction_ready"] != "yes":
        return "context_only" if role == "structural_only_not_diffraction_ready" else "not_diffraction_ready"
    mean = f(row["mean_abs_primary_peak_offset_A"])
    max_v = f(row["max_abs_primary_peak_offset_A"])
    warning = row["structural_warning_flag"] == "yes"
    primary_pass = mean is not None and max_v is not None and mean <= 0.05 and max_v <= 0.09
    multi_pass = consistency_filter(offsets)
    if warning and role == "new_asem_candidate":
        if primary_pass and multi_pass:
            return "survives_but_not_unique"
        return "disfavored_by_structural_warning"
    if primary_pass and multi_pass:
        return "survives_current_filters" if role == "current_positive_baseline" else "survives_but_not_unique"
    return "disfavored_by_peak_offsets"


def make_row(
    model_group: str,
    model_label: str,
    source_file: str,
    twist: str,
    ready: str,
    role: str,
    offsets: dict[str, float | None],
    hbond_no: float | None = None,
    hbond_ho: float | None = None,
    clash_22: float | None = None,
    min_dist: float | None = None,
    note: str = "",
) -> dict[str, object]:
    mean, max_v = summarize_offsets(offsets)
    row: dict[str, object] = {
        "model_group": model_group,
        "model_label": model_label,
        "source_file": source_file,
        "twist_label_deg": twist,
        "is_diffraction_ready": ready,
        "mean_abs_primary_peak_offset_A": fmt(mean),
        "max_abs_primary_peak_offset_A": fmt(max_v),
        "passes_primary_peak_position_filter": "yes" if mean is not None and max_v is not None and mean <= 0.05 and max_v <= 0.09 else "no",
        "passes_multi_window_consistency_filter": "yes" if consistency_filter(offsets) else "no",
        "hbond_proxy_best_N_O_A": fmt(hbond_no),
        "hbond_proxy_best_H_O_A": fmt(hbond_ho),
        "clash_count_lt_2p2": "" if clash_22 is None else str(int(clash_22)),
        "structural_warning_flag": "yes" if structural_warning(clash_22, min_dist) else "no",
        "control_role": role,
        "qualitative_note": note,
    }
    for field in ORDERED_FIELDS:
        row[field] = fmt(offsets.get(field))
    row["falsification_status"] = status_for(row, offsets, role)
    return row


def build_summary() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []

    asem_struct = {r["candidate_label"]: r for r in read_csv(METRICS / "asem_twist_series_structural_summary.csv")}
    for rank in read_csv(METRICS / "asem_twist_series_ranking_corrected_emory_profile.csv"):
        label = rank["candidate_label"]
        s = asem_struct.get(label, {})
        offsets = {
            "base_stack_offset_A": f(rank.get("base_stack_peak_offset_A")),
            "feature_3p77_offset_A": f(rank.get("feature_3p77_offset_A")),
            "feature_4p4_offset_A": f(rank.get("feature_4p4_offset_A")),
            "feature_5p6_offset_A": f(rank.get("feature_5p6_offset_A")),
            "feature_7p3_offset_A": f(rank.get("feature_7p3_offset_A")),
        }
        rows.append(make_row(
            "Asem twist series",
            label,
            rank.get("candidate_file", ""),
            rank.get("twist_label_deg", ""),
            "yes",
            "new_asem_candidate",
            offsets,
            f(s.get("hbond_proxy_best_n_o_distance_A")),
            f(s.get("hbond_proxy_best_h_o_distance_A")),
            f(s.get("clash_count_lt_2p2")),
            f(s.get("min_heavy_heavy_nonbonded_distance_A")),
            "New full Asem candidate; diffraction screen is strong for 29/30, but short-contact proxy flags remain.",
        ))

    baseline_rows = read_csv(METRICS / "nick_ideal_16mer_feature_summary_corrected_emory_profile.csv")
    rows.append(make_row(
        "Ideal baseline",
        "ideal_16mer_antiparallel_30deg",
        "inputs/nick_ideal_models/Hexaplex_AntiParallel_30deg_Ideal.pdb",
        "30",
        "yes",
        "current_positive_baseline",
        offsets_from_feature_rows(baseline_rows),
        note="Current positive corrected-Emory baseline; same peak positions as Asem 30 in the primary windows.",
    ))

    parallel_rows = read_csv(METRICS / "parallel_control_feature_comparison_corrected_emory_profile.csv")
    rows.append(make_row(
        "Parallel control",
        "parallel_sheet_control",
        "inputs/nick_ideal_models/parallel_control/TwoBetaSheetBackbones180_180.pdb",
        "",
        "yes",
        "parallel_control",
        offsets_from_feature_rows(parallel_rows, "parallel_offset_d_A"),
        note="Negative/falsification control; generally worse than antiparallel in prior comparison.",
    ))

    for rank in read_csv(METRICS / "rise_sensitivity_ranking_corrected_emory_profile.csv"):
        rise = rank.get("rise_A", "")
        if rise not in {"3.36", "3.37", "3.38", "3.39", "3.40"}:
            continue
        offsets = {
            "base_stack_offset_A": f(rank.get("base_stack_peak_offset_A")),
            "feature_3p77_offset_A": f(rank.get("feature_3p77_offset_A")),
            "feature_4p4_offset_A": f(rank.get("feature_4p4_offset_A")),
            "feature_5p6_offset_A": f(rank.get("feature_5p6_offset_A")),
            "feature_7p3_offset_A": f(rank.get("feature_7p3_offset_A")),
        }
        rows.append(make_row(
            "Rise sensitivity",
            f"rise_{rise}_control",
            "outputs/metrics/rise_sensitivity_ranking_corrected_emory_profile.csv",
            "",
            "yes",
            "rise_variant_control",
            offsets,
            note="Existing rise-sensitivity control reused; no diffraction rerun.",
        ))

    atom_rows = read_csv(METRICS / "nick_atom_contribution_feature_summary_corrected_emory_profile.csv")
    by_model: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in atom_rows:
        by_model[row.get("model", "")].append(row)
    for model, model_rows in sorted(by_model.items()):
        label = model_rows[0].get("model_label", model) if model_rows else model
        rows.append(make_row(
            "Atom contribution",
            f"atom_contribution_{model}",
            "outputs/metrics/nick_atom_contribution_feature_summary_corrected_emory_profile.csv",
            "",
            "yes",
            "atom_contribution_control",
            offsets_from_feature_rows(model_rows),
            note=f"Atom-contribution context: {label}; useful for peak-assignment interpretation, not a new full candidate.",
        ))

    rows.append(make_row(
        "Omega structural context",
        "omega172_structural_context",
        "outputs/reports/full_omega172_transfer_stop_report.md",
        "",
        "no",
        "structural_only_not_diffraction_ready",
        {field: None for field in ORDERED_FIELDS},
        note="Omega172 remains local structural context; no full diffraction-ready transferred model exists.",
    ))
    rows.append(make_row(
        "Asem omega series",
        "asem_struct2_omega_series_context",
        "outputs/metrics/asem_omega_series_structural_audit.csv",
        "",
        "no",
        "structural_only_not_diffraction_ready",
        {field: None for field in ORDERED_FIELDS},
        note="Struct2 omega-series audit is structural context only for this comparison.",
    ))
    return rows


def apply_asem_near_best_filter(summary: list[dict[str, object]], tolerance_A: float = 0.005) -> None:
    """Apply the relative Asem-family survival rule requested for this update."""
    asem_rows = [row for row in summary if row.get("model_group") == "Asem twist series"]
    means = [f(row.get("mean_abs_primary_peak_offset_A")) for row in asem_rows]
    finite_means = [value for value in means if value is not None]
    if not finite_means:
        return
    best_mean = min(finite_means)
    cutoff = best_mean + tolerance_A
    for row in asem_rows:
        mean = f(row.get("mean_abs_primary_peak_offset_A"))
        max_v = f(row.get("max_abs_primary_peak_offset_A"))
        multi_pass = row.get("passes_multi_window_consistency_filter") == "yes"
        near_best = mean is not None and mean <= cutoff
        absolute_peak_pass = mean is not None and max_v is not None and max_v <= 0.09
        if near_best and absolute_peak_pass and multi_pass:
            row["passes_primary_peak_position_filter"] = "yes"
            row["falsification_status"] = "survives_current_filters"
            row["qualitative_note"] = (
                str(row.get("qualitative_note", ""))
                + f" Relative near-best Asem filter passed: mean offset <= {cutoff:.6g} A."
            ).strip()
        else:
            row["passes_primary_peak_position_filter"] = "no"
            row["falsification_status"] = "disfavored_by_peak_offsets"
            row["qualitative_note"] = (
                str(row.get("qualitative_note", ""))
                + f" Relative near-best Asem filter failed: best mean {best_mean:.6g} A, cutoff {cutoff:.6g} A."
            ).strip()
def rank_rows(summary: list[dict[str, object]]) -> list[dict[str, object]]:
    ranked = []
    role_bonus = {
        "new_asem_candidate": 0.0,
        "current_positive_baseline": 0.005,
        "rise_variant_control": 0.02,
        "parallel_control": 0.04,
        "atom_contribution_control": 0.05,
    }
    for row in summary:
        if row["is_diffraction_ready"] != "yes":
            continue
        mean = f(row["mean_abs_primary_peak_offset_A"]) or 9.0
        max_v = f(row["max_abs_primary_peak_offset_A"]) or 9.0
        structural_penalty = 0.03 if row["structural_warning_flag"] == "yes" else 0.0
        hbond_no = f(row["hbond_proxy_best_N_O_A"])
        hbond_penalty = 0.0 if hbond_no is None or 2.0 <= hbond_no <= 3.2 else 0.02
        score = mean + 0.25 * max_v + structural_penalty + hbond_penalty + role_bonus.get(str(row["control_role"]), 0.08)
        out = dict(row)
        out["combined_screen_score"] = fmt(score)
        out["rank_formula"] = "mean_abs_offset + 0.25*max_abs_offset + structural_warning(0.03) + Hbond_penalty(0.02) + role/context penalty"
        ranked.append(out)
    ranked.sort(key=lambda r: f(r["combined_screen_score"]) or 9.0)
    for idx, row in enumerate(ranked, start=1):
        row["rank"] = idx
    return ranked


def plot_bar(summary: list[dict[str, object]]) -> None:
    rows = [r for r in summary if r["is_diffraction_ready"] == "yes"]
    labels = [str(r["model_label"]) for r in rows]
    values = [f(r["mean_abs_primary_peak_offset_A"]) or 0 for r in rows]
    fig, ax = plt.subplots(figsize=(10, 5), dpi=180)
    ax.bar(labels, values, color="#4c78a8")
    ax.set_ylabel("Mean abs primary peak offset (A)")
    ax.tick_params(axis="x", rotation=45, labelsize=7)
    fig.tight_layout()
    fig.savefig(PLOTS / "primary_window_mean_offset_by_model.png")
    plt.close(fig)


def plot_heatmap(summary: list[dict[str, object]]) -> None:
    rows = [r for r in summary if r["is_diffraction_ready"] == "yes"]
    data = np.asarray([[f(r[field]) or np.nan for field in ORDERED_FIELDS] for r in rows], dtype=float)
    fig, ax = plt.subplots(figsize=(8.5, 6), dpi=180)
    im = ax.imshow(data, cmap="coolwarm", aspect="auto", vmin=-0.15, vmax=0.15)
    ax.set_yticks(np.arange(len(rows)))
    ax.set_yticklabels([str(r["model_label"]) for r in rows], fontsize=7)
    ax.set_xticks(np.arange(len(ORDERED_FIELDS)))
    ax.set_xticklabels(ORDERED_FIELDS, rotation=35, ha="right", fontsize=7)
    fig.colorbar(im, ax=ax, label="Peak offset (A)")
    fig.tight_layout()
    fig.savefig(PLOTS / "feature_by_feature_offset_heatmap.png")
    plt.close(fig)


def plot_asem_vs_baseline(summary: list[dict[str, object]]) -> None:
    keep = [r for r in summary if r["model_group"] == "Asem twist series" or r["control_role"] == "current_positive_baseline"]
    x = np.arange(len(ORDERED_FIELDS))
    fig, ax = plt.subplots(figsize=(8.5, 5), dpi=180)
    for row in keep:
        ax.plot(x, [f(row[field]) or np.nan for field in ORDERED_FIELDS], marker="o", label=str(row["model_label"]))
    ax.axhline(0, color="black", lw=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(ORDERED_FIELDS, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("Peak offset (A)")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(PLOTS / "asem_29_30_31_vs_current_ideal_baseline.png")
    plt.close(fig)


def plot_scatter(summary: list[dict[str, object]], x_field: str, output: str, xlabel: str) -> None:
    rows = [r for r in summary if r["model_group"] == "Asem twist series"]
    fig, ax = plt.subplots(figsize=(6.5, 4.5), dpi=180)
    for row in rows:
        x = f(row[x_field])
        y = f(row["mean_abs_primary_peak_offset_A"])
        if x is None or y is None:
            continue
        ax.scatter(x, y, s=55)
        ax.annotate(str(row["twist_label_deg"]), (x, y), textcoords="offset points", xytext=(5, 3), fontsize=8)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Mean abs primary peak offset (A)")
    fig.tight_layout()
    fig.savefig(PLOTS / output)
    plt.close(fig)


def plot_status(summary: list[dict[str, object]]) -> None:
    counts = Counter(str(r["falsification_status"]) for r in summary)
    fig, ax = plt.subplots(figsize=(7, 4), dpi=180)
    ax.bar(list(counts), list(counts.values()), color="#59a14f")
    ax.tick_params(axis="x", rotation=30, labelsize=8)
    ax.set_ylabel("Model/control count")
    fig.tight_layout()
    fig.savefig(PLOTS / "falsification_status_summary.png")
    plt.close(fig)


def make_plots(summary: list[dict[str, object]]) -> None:
    PLOTS.mkdir(parents=True, exist_ok=True)
    plot_bar(summary)
    plot_heatmap(summary)
    plot_asem_vs_baseline(summary)
    plot_scatter(summary, "hbond_proxy_best_N_O_A", "hbond_proxy_vs_mean_diffraction_offset.png", "Best N-O H-bond proxy (A)")
    plot_scatter(summary, "clash_count_lt_2p2", "short_contact_count_vs_mean_diffraction_offset.png", "Heavy-heavy contacts < 2.2 A")
    plot_status(summary)


def md_table(rows: list[dict[str, object]], fields: list[str]) -> str:
    lines = ["| " + " | ".join(fields) + " |", "| " + " | ".join(["---"] * len(fields)) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(str(row.get(field, "")) for field in fields) + " |")
    return "\n".join(lines)


def write_report(summary: list[dict[str, object]], ranking: list[dict[str, object]]) -> None:
    asem = [r for r in summary if r["model_group"] == "Asem twist series"]
    rows_29_30 = [r for r in asem if str(r["twist_label_deg"]) in {"29", "30"}]
    row_31 = next((r for r in asem if str(r["twist_label_deg"]) == "31"), None)
    lines = [
        "# Asem Twist Series Falsification Update",
        "",
        "## Purpose",
        "",
        "Apply the prior falsification-style framework to Asem's new 29, 30, and 31 degree twist-series models.",
        "",
        "## Inputs",
        "",
        "- Imported Asem source folders: `inputs/asem_twist_series_29_30_31/raw/29`, `raw/30`, `raw/31`.",
        "- Corrected experimental profile: `inputs/experimental/nick_powder_profile_corrected_emory.csv`.",
        "- Reused benchmark metrics from `outputs/metrics/asem_twist_series_*`.",
        "- Reused current ideal baseline, parallel control, rise-sensitivity, atom-contribution, omega audit, and omega172 stop-context outputs.",
        "",
        "## Methods",
        "",
        "No diffraction was rerun. The script reuses existing corrected-diffraction CSVs produced with the Asem-corrected non-accumulating/vectorized path, `tilts = [0]`, `rotations = range(0, 181, 5)`, hydrogen exclusion, exact heavy-atom deduplication, and the current baseline grid/radial settings. Primary windows are 3.38/3.4 A, 3.77 A, 4.4 A, 5.6 A, and 7.3 A.",
        "",
        "The combined ranking score is documented in the CSV as: mean absolute primary-window offset + 0.25 times max offset + structural warning penalty + H-bond plausibility penalty + role/context penalty. This is a transparent screen, not a fitted physical score.",
        "",
        "## Main Findings",
        "",
        "- 29 degree and 30 degree tie by the current mean primary-window peak offset.",
        "- 31 degree is worse by the same metric and is comparatively disfavored.",
        "- 30 degree remains supported, but it is not uniquely selected by this metric.",
        "- The result supports a near-30 degree twist family rather than exact unique 30 degree from this run alone.",
        "- All Asem candidates show plausible H-bond proxies, but all also carry short-contact flags.",
        "",
        "## Falsification Summary",
        "",
        md_table(summary, ["model_label", "control_role", "mean_abs_primary_peak_offset_A", "max_abs_primary_peak_offset_A", "passes_primary_peak_position_filter", "passes_multi_window_consistency_filter", "structural_warning_flag", "falsification_status"]),
        "",
        "## Ranking",
        "",
        md_table(ranking, ["rank", "model_label", "control_role", "combined_screen_score", "mean_abs_primary_peak_offset_A", "max_abs_primary_peak_offset_A", "falsification_status"]),
        "",
        "## Interpretation",
        "",
        "These results strengthen the case that the corrected experimental profile is selecting a narrow helical/twist family rather than arbitrary models. The current screen should not be read as proof of exact unique 30 degree twist because 29 and 30 tie on peak positions and only one full 30 degree candidate is available.",
        "",
        "Prior control framing remains intact: the parallel control is comparatively disfavored by peak offsets, rise 3.38 did not beat rise 3.40 overall in the previous controlled test, and omega172 remains a local structural recommendation because no full diffraction-ready omega172 model exists yet.",
        "",
        "## Recommended Next Step",
        "",
        "Ask Asem to send several completed full 30 degree PDB candidates with carboxylates built. Then rerun this same falsification screen as an ensemble-stability test.",
        "",
        "## Limitations",
        "",
        "- No AmberTools/tleap generation was run; AmberTools/tleap is not available here.",
        "- No minimization was run.",
        "- No notebooks were executed.",
        "- Powder/radial scores are screening metrics, not full structural refinement.",
        "- H-bond and steric/contact metrics are proxies.",
        "- The 29/30 tie may depend on the current radial-window metric and should not be overclaimed.",
        "",
        "## Appendix: Short Nick/Asem Update Draft",
        "",
        "Asem's provided 29, 30, and 31 degree full PDBs all benchmark cleanly through the corrected Emory-profile screen. The 29 and 30 degree models tie by the current mean primary-window peak-offset metric, while 31 degree is worse, so the data support a near-30 degree twist family rather than uniquely selecting exact 30 degree from this single-candidate run. All three show plausible backbone H-bond proxy distances, but they also retain short-contact flags, so I would treat them as screening candidates. The most decisive next test would be several independent completed 30 degree full PDB candidates with carboxylates built, then rerunning this same falsification screen for ensemble stability.",
    ]
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    summary = build_summary()
    apply_asem_near_best_filter(summary)
    fields = [
        "model_group", "model_label", "source_file", "twist_label_deg", "is_diffraction_ready",
        "mean_abs_primary_peak_offset_A", "max_abs_primary_peak_offset_A", "base_stack_offset_A",
        "feature_3p77_offset_A", "feature_4p4_offset_A", "feature_5p6_offset_A", "feature_7p3_offset_A",
        "passes_primary_peak_position_filter", "passes_multi_window_consistency_filter",
        "hbond_proxy_best_N_O_A", "hbond_proxy_best_H_O_A", "clash_count_lt_2p2",
        "structural_warning_flag", "control_role", "falsification_status", "qualitative_note",
    ]
    write_csv(SUMMARY_CSV, summary, fields)
    ranking = rank_rows(summary)
    rank_fields = ["rank", *fields, "combined_screen_score", "rank_formula"]
    write_csv(RANKING_CSV, ranking, rank_fields)
    make_plots(summary)
    write_report(summary, ranking)
    print(f"Wrote {SUMMARY_CSV}")
    print(f"Wrote {RANKING_CSV}")
    print(f"Wrote {REPORT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


