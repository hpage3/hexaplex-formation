#!/usr/bin/env python3
"""Analyze generated pNAB twists with serial or model-level parallel workers."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import math
import statistics
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from hexaplex_formation.pdb_utils import heavy_atoms, load_pdb_atoms  # noqa: E402
from hexaplex_formation.scattering import (  # noqa: E402
    debye_intensity_from_distance_histogram,
    d_from_q,
    make_q_grid,
    pair_distance_histogram_for_debye,
)


TARGETS = [
    ("base", "base_stacking", 3.38, 3.25, 3.52),
    ("A", "backbone_associated", 3.80, 3.60, 4.00),
    ("B", "backbone_associated", 4.40, 4.15, 4.65),
    ("C", "backbone_associated", 5.65, 5.30, 6.00),
    ("D", "backbone_associated", 7.30, 6.80, 7.80),
]
PROFILE_FIELDS = ["q_Ainv", "d_A", "intensity", "intensity_norm"]
PEAK_FIELDS = [
    "model_id", "twist_deg", "rise_A", "target_label", "target_group",
    "target_d_A", "window_min_A", "window_max_A", "peak_d_A",
    "peak_intensity_optional", "assignment_status", "notes",
]
ANALYSIS_FIELDS = [
    "model_id", "twist_deg", "rise_A", "generation_status", "analysis_status",
    "analysis_mode", "worker_count", "structure_path", "radial_profile_path",
    "peak_list_path", "elapsed_seconds", "error_summary",
]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=ROOT / "outputs/pnab_twist_series_rise3p38/pnab_twist_series_manifest.csv",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=ROOT / "outputs/pnab_twist_series_rise3p38",
    )
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--serial", action="store_true")
    parser.add_argument("--max-models", type=int)
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument("--timeout-seconds", type=int, default=None)
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    args = build_parser().parse_args(argv)
    if args.workers < 1:
        raise SystemExit("--workers must be at least 1")
    if args.max_models is not None and args.max_models < 1:
        raise SystemExit("--max-models must be at least 1")
    if args.timeout_seconds is not None and args.timeout_seconds < 1:
        raise SystemExit("--timeout-seconds must be at least 1")
    return args


def analysis_mode(args: argparse.Namespace) -> tuple[str, int]:
    if args.serial or args.workers == 1:
        return "serial", 1
    return "parallel", args.workers


def load_score_module():
    path = ROOT / "scripts" / "score_peak_position_fit.py"
    spec = importlib.util.spec_from_file_location("pnab_score_peak_position_fit", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def radial_profile(pdb_path: Path) -> list[dict[str, float]]:
    atoms = heavy_atoms(load_pdb_atoms(pdb_path))
    histogram = pair_distance_histogram_for_debye(atoms, bin_width=0.05, element_weights=None)
    q_values = make_q_grid(0.65, 2.15, 0.005)
    intensities = debye_intensity_from_distance_histogram(histogram, q_values)
    maximum = max(intensities) if intensities else 1.0
    return [
        {"q_Ainv": q, "d_A": d_from_q(q), "intensity": value, "intensity_norm": value / maximum}
        for q, value in zip(q_values, intensities)
    ]


def local_peak(rows: list[dict[str, float]], low: float, high: float) -> tuple[dict[str, float] | None, str]:
    window = sorted((row for row in rows if low <= row["d_A"] <= high), key=lambda row: row["d_A"])
    if len(window) < 3:
        return None, "missing_window"
    candidates = []
    for index in range(1, len(window) - 1):
        if (
            window[index]["intensity"] >= window[index - 1]["intensity"]
            and window[index]["intensity"] >= window[index + 1]["intensity"]
        ):
            candidates.append(window[index])
    if not candidates:
        return max(window, key=lambda row: row["intensity"]), "ambiguous_boundary_or_monotonic"
    return max(candidates, key=lambda row: row["intensity"]), "local_maximum"


def output_paths(output_root: Path, model_id: str) -> tuple[Path, Path]:
    return (
        output_root / "radial_profiles" / f"{model_id}_radial.csv",
        output_root / "peak_lists" / f"{model_id}_peaks.csv",
    )


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def make_peak_rows(item: dict[str, str], profile: list[dict[str, float]]) -> list[dict[str, object]]:
    rows = []
    for label, group, center, low, high in TARGETS:
        peak, status = local_peak(profile, low, high)
        peak_d = peak["d_A"] if peak else math.nan
        intensity = peak["intensity_norm"] if peak else math.nan
        rows.append(
            {
                "model_id": item["model_id"],
                "twist_deg": item["twist_deg"],
                "rise_A": item["rise_A"],
                "target_label": label,
                "target_group": group,
                "target_d_A": f"{center:.6g}",
                "window_min_A": f"{low:.6g}",
                "window_max_A": f"{high:.6g}",
                "peak_d_A": "" if math.isnan(peak_d) else f"{peak_d:.8g}",
                "peak_intensity_optional": "" if math.isnan(intensity) else f"{intensity:.8g}",
                "assignment_status": status,
                "notes": (
                    "Simplified isotropic Debye radial profile; assignment requires "
                    "an interior local maximum."
                ),
            }
        )
    return rows


def analyze_model_job(payload: dict[str, str]) -> dict[str, str]:
    started = time.perf_counter()
    item = payload["item"]
    output_root = Path(payload["output_root"])
    profile_path, peak_path = output_paths(output_root, item["model_id"])
    profile = radial_profile(ROOT / item["structure_path"])
    peaks = make_peak_rows(item, profile)
    write_csv(profile_path, PROFILE_FIELDS, profile)
    write_csv(peak_path, PEAK_FIELDS, peaks)
    return {
        "model_id": item["model_id"],
        "radial_profile_path": display_path(profile_path),
        "peak_list_path": display_path(peak_path),
        "elapsed_seconds": f"{time.perf_counter() - started:.6f}",
    }


def base_manifest_row(
    item: dict[str, str],
    mode: str,
    workers: int,
    status: str,
    error: str = "",
) -> dict[str, str]:
    return {
        "model_id": item["model_id"],
        "twist_deg": item["twist_deg"],
        "rise_A": item["rise_A"],
        "generation_status": item["status"],
        "analysis_status": status,
        "analysis_mode": mode,
        "worker_count": str(workers),
        "structure_path": item.get("structure_path", ""),
        "radial_profile_path": "",
        "peak_list_path": "",
        "elapsed_seconds": "",
        "error_summary": error,
    }


def plan_jobs(
    generation_rows: list[dict[str, str]],
    output_root: Path,
    mode: str,
    workers: int,
    skip_existing: bool,
    max_models: int | None,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    jobs = []
    manifest_rows = []
    considered = 0
    for source in generation_rows:
        item = dict(source)
        item["_output_root"] = str(output_root.resolve())
        if not item.get("structure_path") or not (ROOT / item.get("structure_path", "")).exists():
            manifest_rows.append(base_manifest_row(item, mode, workers, "skipped_missing_generation"))
            continue
        if max_models is not None and considered >= max_models:
            manifest_rows.append(base_manifest_row(item, mode, workers, "skipped_max_models"))
            continue
        considered += 1
        profile_path, peak_path = output_paths(output_root, item["model_id"])
        if skip_existing and profile_path.exists() and peak_path.exists():
            row = base_manifest_row(item, mode, workers, "skipped_existing")
            row["radial_profile_path"] = display_path(profile_path)
            row["peak_list_path"] = display_path(peak_path)
            manifest_rows.append(row)
            continue
        jobs.append({"item": item, "output_root": str(output_root.resolve())})
    return jobs, manifest_rows


def run_jobs(
    jobs: list[dict[str, str]],
    mode: str,
    workers: int,
    timeout_seconds: int | None,
    worker_fn=analyze_model_job,
) -> list[dict[str, str]]:
    rows = []
    if mode == "serial":
        for payload in jobs:
            item = payload["item"]
            started = time.perf_counter()
            try:
                result = worker_fn(payload)
                row = base_manifest_row(item, mode, workers, "success")
                row.update(result)
            except Exception as exc:  # noqa: BLE001 - batch must continue
                row = base_manifest_row(item, mode, workers, "failed", f"{type(exc).__name__}: {exc}")
                row["elapsed_seconds"] = f"{time.perf_counter() - started:.6f}"
            rows.append(row)
        return rows

    with ProcessPoolExecutor(max_workers=workers) as executor:
        future_to_payload = {executor.submit(worker_fn, payload): payload for payload in jobs}
        try:
            iterator = as_completed(future_to_payload, timeout=timeout_seconds)
            for future in iterator:
                payload = future_to_payload[future]
                item = payload["item"]
                try:
                    result = future.result()
                    row = base_manifest_row(item, mode, workers, "success")
                    row.update(result)
                except Exception as exc:  # noqa: BLE001 - isolate worker failure
                    row = base_manifest_row(item, mode, workers, "failed", f"{type(exc).__name__}: {exc}")
                rows.append(row)
        except TimeoutError:
            finished = {row["model_id"] for row in rows}
            for future, payload in future_to_payload.items():
                item = payload["item"]
                if item["model_id"] not in finished:
                    future.cancel()
                    rows.append(base_manifest_row(item, mode, workers, "timed_out", "analysis timeout exceeded"))
    return rows


def score_group(rows: list[dict[str, object]], group: str) -> float:
    values = [
        float(row["abs_relative_error"])
        for row in rows
        if row["target_group"] == group and row["is_missing"] == "false"
    ]
    return math.sqrt(statistics.mean(value * value for value in values)) if values else math.nan


def aggregate_outputs(analysis_rows: list[dict[str, str]], output_root: Path) -> list[dict[str, str]]:
    score = load_score_module()
    targets = [score.TargetPeak(label, center, 1.0, group) for label, group, center, _low, _high in TARGETS]
    all_peak_rows = []
    assignment_rows = []
    window_rows = []
    for item in analysis_rows:
        if item["analysis_status"] not in {"success", "skipped_existing"}:
            continue
        peak_path = Path(item["peak_list_path"])
        if not peak_path.is_absolute():
            peak_path = ROOT / peak_path
        if not peak_path.exists():
            continue
        model_peaks = read_csv(peak_path)
        all_peak_rows.extend(model_peaks)
        for row in model_peaks:
            peak_d = float(row["peak_d_A"]) if row["peak_d_A"] else math.nan
            center = float(row["target_d_A"])
            relative = abs((center - peak_d) / center) if not math.isnan(peak_d) else math.nan
            window_rows.append(
                {
                    **row,
                    "abs_relative_error": "" if math.isnan(relative) else f"{relative:.12g}",
                    "accepted_for_peak_scoring": "yes" if row["assignment_status"] == "local_maximum" else "no",
                }
            )
            if row["assignment_status"] == "local_maximum":
                assignment_rows.append(
                    {
                        "model_id": row["model_id"],
                        "twist_deg": row["twist_deg"],
                        "rise_A": row["rise_A"],
                        "target_label": row["target_label"],
                        "theoretical_d_A": row["peak_d_A"],
                        "assignment_method": (
                            f"interior_local_maximum_in_{row['window_min_A']}_{row['window_max_A']}A"
                        ),
                        "notes": "Simplified isotropic Debye radial profile.",
                    }
                )

    assignment_path = ROOT / "inputs/theoretical_peak_assignments/pnab_rise3p38_twist_series_peak_assignments.csv"
    write_csv(assignment_path, score.ASSIGNMENT_FIELDS, assignment_rows)
    if all_peak_rows:
        write_csv(output_root / "peak_lists/pnab_rise3p38_twist_series_all_peaks.csv", PEAK_FIELDS, all_peak_rows)
        write_csv(
            output_root / "window_scores/pnab_rise3p38_twist_series_window_details.csv",
            list(window_rows[0]),
            window_rows,
        )

    assignments = score.read_assignments(assignment_path)
    summary, per_peak = score.score_models(targets, assignments, missing_peak_penalty=0.25)
    score.write_csv(
        ROOT / "outputs/metrics/pnab_rise3p38_twist_series_peak_position_fit_summary.csv",
        score.SUMMARY_FIELDS,
        summary,
    )
    score.write_csv(
        ROOT / "outputs/metrics/pnab_rise3p38_twist_series_peak_position_fit_per_peak_errors.csv",
        score.PER_PEAK_FIELDS,
        per_peak,
    )

    by_model = {}
    for row in per_peak:
        by_model.setdefault(row["model_id"], []).append(row)
    fit_rows = []
    for model, rows in by_model.items():
        base = score_group(rows, "base_stacking")
        backbone = score_group(rows, "backbone_associated")
        combined_values = [float(row["abs_relative_error"]) for row in rows]
        combined = math.sqrt(statistics.mean(value * value for value in combined_values))
        fit_rows.append(
            {
                "model_id": model,
                "twist_deg": rows[0]["twist_deg"],
                "rise_A": rows[0]["rise_A"],
                "base_stacking_score": "" if math.isnan(base) else f"{base:.12g}",
                "backbone_associated_score": "" if math.isnan(backbone) else f"{backbone:.12g}",
                "combined_score": f"{combined:.12g}",
                "missing_peak_count": str(sum(row["is_missing"] == "true" for row in rows)),
                "metric_note": (
                    "RMS relative local-peak-position error; missing assignments use the scorer penalty."
                ),
            }
        )
    fit_rows.sort(key=lambda row: float(row["combined_score"]))
    if fit_rows:
        write_csv(
            ROOT / "outputs/metrics/pnab_rise3p38_twist_series_window_fit_summary.csv",
            list(fit_rows[0]),
            fit_rows,
        )
    return fit_rows


def write_report(
    path: Path,
    analysis_rows: list[dict[str, str]],
    fit_rows: list[dict[str, str]],
    mode: str,
    workers: int,
) -> None:
    counts = {status: sum(row["analysis_status"] == status for row in analysis_rows) for status in {
        "success", "skipped_existing", "skipped_missing_generation", "failed", "timed_out", "skipped_max_models"
    }}
    available_twists = [
        row["twist_deg"]
        for row in analysis_rows
        if row["analysis_status"] in {"success", "skipped_existing"}
    ]
    missing_twists = [
        row["twist_deg"]
        for row in analysis_rows
        if row["analysis_status"] == "skipped_missing_generation"
    ]
    best = fit_rows[0] if fit_rows else None
    ranked = {row["twist_deg"]: row for row in fit_rows}
    lines = [
        "# pNAB Fixed-Rise Twist-Series Analysis",
        "",
        "- Source workflow: `inputs/pnab_asem_30deg`",
        "- Original source: `C:\\Users\\hpage3\\OneDrive - Georgia Institute of Technology\\Documents\\GitHub\\research\\30_A`",
        "- Fixed rise: 3.38 A",
        "- Attempted twist grid: 28.0 to 32.0 degrees in 0.5-degree increments",
        "- Generated structure directory: `outputs/pnab_twist_series_rise3p38/structures`",
        f"- Analysis mode: {mode}",
        f"- Worker count: {workers}",
        "",
        "## Analysis Status",
        "",
        f"- Successful computations: {counts['success']}",
        f"- Existing outputs reused: {counts['skipped_existing']}",
        f"- Available analyzed model outputs: {len(available_twists)} ({', '.join(available_twists)})",
        f"- Missing generation outputs skipped: {counts['skipped_missing_generation']}",
        f"- Generation timeouts without promoted structures: {', '.join(missing_twists) or 'none'}",
        f"- Analysis failures: {counts['failed']}",
        f"- Analysis timeouts: {counts['timed_out']}",
        "",
        "## Comparative Readout",
        "",
    ]
    if best:
        lines.append(
            f"- Best available penalized combined local-window score: {best['twist_deg']} degrees "
            f"({best['combined_score']})."
        )
    else:
        lines.append("- No successful model profiles were available for comparative scoring.")
    if "31.0" in ranked:
        rank = next(index for index, row in enumerate(fit_rows, 1) if row["twist_deg"] == "31.0")
        lines.append(f"- The available 31.0-degree model ranks {rank} of {len(fit_rows)} by combined score.")
        other_missing = [
            row for row in fit_rows
            if row["twist_deg"] != "31.0" and int(row["missing_peak_count"]) > 0
        ]
        if rank == 1 and other_missing:
            lines.append(
                "- This does not establish that 31.0 degrees is preferred or refute the earlier disfavored "
                "assessment: it is the only available model with an accepted interior B-window maximum, "
                "while the other available models receive a missing-peak penalty."
            )
    if "29.5" in ranked and "30.0" in ranked:
        relation = "lower" if float(ranked["29.5"]["combined_score"]) < float(ranked["30.0"]["combined_score"]) else "higher"
        lines.append(
            f"- The 29.5-degree combined score is {relation} than the 30.0-degree score, "
            "but both lack an accepted B-window assignment and the numerical separation is small."
        )
    if "30.5" in ranked and "30.0" in ranked:
        relation = "lower" if float(ranked["30.5"]["combined_score"]) < float(ranked["30.0"]["combined_score"]) else "higher"
        lines.append(f"- The 30.5-degree combined score is {relation} than the 30.0-degree score.")
    lines.extend(
        [
            "",
            "These rankings are comparative and are not proof of a unique structure. Missing pNAB models are "
            "not excluded scientifically; they were unavailable because generation timed out.",
            "",
            "## Caveats",
            "",
            "- pNAB structures are not MD- or QM-relaxed in this workflow.",
            "- The radial profiles use a simplified isotropic Debye approximation.",
            "- Powder/radial profiles discard directional two-dimensional information.",
            "- Peak assignments require an interior local maximum in a predefined target window; ambiguous windows are left unassigned.",
            "- Parallel execution distributes independent model jobs. Each model calculation remains serial, avoiding nested worker pools.",
            "- Generated PDB structures and per-twist work directories are retained locally but are not intended for this commit; the generator and manifests reproduce their locations.",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    mode, workers = analysis_mode(args)
    generation_rows = read_csv(args.manifest)
    jobs, analysis_rows = plan_jobs(
        generation_rows,
        args.output_root,
        mode,
        workers,
        args.skip_existing,
        args.max_models,
    )
    analysis_rows.extend(run_jobs(jobs, mode, workers, args.timeout_seconds))
    analysis_rows.sort(key=lambda row: float(row["twist_deg"]))
    analysis_manifest = args.output_root / "analysis_manifest.csv"
    write_csv(analysis_manifest, ANALYSIS_FIELDS, analysis_rows)
    fit_rows = aggregate_outputs(analysis_rows, args.output_root)
    write_report(
        ROOT / "outputs/reports/pnab_rise3p38_twist_series_report.md",
        analysis_rows,
        fit_rows,
        mode,
        workers,
    )
    print(f"Wrote {analysis_manifest}")
    print(f"Analyzed {sum(row['analysis_status'] == 'success' for row in analysis_rows)} model(s)")
    return 0 if not any(row["analysis_status"] == "failed" for row in analysis_rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
