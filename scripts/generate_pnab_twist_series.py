#!/usr/bin/env python3
"""Generate a fixed-rise pNAB twist series from Asem's copied workflow."""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_WORKFLOW = ROOT / "inputs" / "pnab_asem_30deg"
DEFAULT_OUTPUT = ROOT / "outputs" / "pnab_twist_series_rise3p38"
DEFAULT_CONDA = Path(r"C:\Users\Public\hexaplex-tools\miniforge3\Scripts\conda.exe")
DEFAULT_TWISTS = [28.0, 28.5, 29.0, 29.5, 30.0, 30.5, 31.0, 31.5, 32.0]
TEARDOWN_CODE = -1073741819

MANIFEST_FIELDS = [
    "model_id",
    "twist_deg",
    "rise_A",
    "status",
    "work_dir",
    "structure_path",
    "results_csv_path",
    "prefix_yaml_path",
    "log_path",
    "accepted_candidate_count",
    "best_distance_A",
    "conformer_index",
    "atom_count",
    "residue_count",
    "chain_count",
    "return_code",
    "error_summary",
    "warnings",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workflow-dir", type=Path, default=DEFAULT_WORKFLOW)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--rise", type=float, default=3.38)
    parser.add_argument("--twists", nargs="+", type=float, default=DEFAULT_TWISTS)
    parser.add_argument("--conda-exe", type=Path, default=DEFAULT_CONDA)
    parser.add_argument("--conda-env", default="pnab")
    parser.add_argument("--timeout-seconds", type=int, default=420)
    parser.add_argument("--resume", action="store_true")
    return parser.parse_args()


def number_token(value: float) -> str:
    text = f"{value:.6f}".rstrip("0").rstrip(".")
    return text.replace("-", "m").replace(".", "p")


def model_id(twist: float, rise: float) -> str:
    return f"pnab_twist{number_token(twist)}_rise{number_token(rise)}"


def update_helical_parameters(options: dict, rise: float, twist: float) -> dict:
    updated = yaml.safe_load(yaml.safe_dump(options, sort_keys=False))
    helical = updated.get("HelicalParameters")
    if not isinstance(helical, dict):
        raise ValueError("options.yaml is missing HelicalParameters")
    if "h_rise" not in helical or "h_twist" not in helical:
        raise ValueError("options.yaml must define h_rise and h_twist")
    helical["h_rise"] = [float(rise), float(rise), 1]
    helical["h_twist"] = [float(twist), float(twist), 1]
    return updated


def parse_results(path: Path) -> dict[str, str]:
    rows = []
    if path.exists():
        for line in path.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                rows.append([item.strip() for item in line.split(",")])
    if not rows:
        return {}
    rows.sort(key=lambda row: float(row[2]))
    best = rows[0]
    return {
        "accepted_candidate_count": str(len(rows)),
        "prefix": str(int(float(best[0]))),
        "conformer_index": str(int(float(best[1]))),
        "best_distance_A": f"{float(best[2]):.6f}",
    }


def pdb_counts(path: Path) -> tuple[int, int, int]:
    atoms = 0
    residues = set()
    chains = set()
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        atoms += 1
        chain = line[21:22].strip()
        residues.add((chain, line[22:26].strip(), line[17:20].strip()))
        chains.add(chain)
    return atoms, len(residues), len(chains)


def completed_outputs(work_dir: Path) -> bool:
    return all((work_dir / name).exists() for name in ("fixed.pdb", "results.csv", "prefix.yaml"))


def classify_run(return_code: int, work_dir: Path, result: dict[str, str]) -> tuple[str, str, str]:
    outputs_complete = completed_outputs(work_dir) and bool(result)
    if return_code == 0 and outputs_complete:
        return "success", "", ""
    if outputs_complete and return_code in {TEARDOWN_CODE, 3221225477}:
        return (
            "success_with_warning",
            "",
            "Windows teardown access violation occurred after complete outputs were written",
        )
    if outputs_complete:
        return (
            "success_with_warning",
            "",
            f"non-zero return code {return_code} after complete outputs were written",
        )
    if return_code == 124:
        return "timed_out", "pNAB exceeded the configured timeout", ""
    if not result:
        return "failed", "no accepted pNAB conformer found", ""
    return "failed", f"pNAB return code {return_code}; required outputs incomplete", ""


def copy_workflow(source: Path, destination: Path) -> None:
    if destination.exists():
        shutil.rmtree(destination)
    shutil.copytree(
        source,
        destination,
        ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "fixed.pdb", "results.csv", "prefix.yaml", "[0-9]*_[0-9]*.pdb"),
    )


def run_one(
    workflow_dir: Path,
    output_root: Path,
    rise: float,
    twist: float,
    conda_exe: Path,
    conda_env: str,
    timeout_seconds: int,
) -> dict[str, str]:
    token = number_token(twist)
    current_id = model_id(twist, rise)
    work_dir = output_root / "work" / f"twist_{token}"
    log_path = output_root / "logs" / f"twist_{token}.log"
    structure_path = output_root / "structures" / f"{current_id}.pdb"
    metadata_dir = output_root / "metadata"
    results_out = metadata_dir / f"{current_id}_results.csv"
    prefix_out = metadata_dir / f"{current_id}_prefix.yaml"

    copy_workflow(workflow_dir, work_dir)
    options_path = work_dir / "options.yaml"
    options = yaml.safe_load(options_path.read_text(encoding="utf-8"))
    updated = update_helical_parameters(options, rise, twist)
    options_path.write_text(yaml.safe_dump(updated, sort_keys=False), encoding="utf-8")

    log_path.parent.mkdir(parents=True, exist_ok=True)
    command = [str(conda_exe), "run", "-n", conda_env, "python", "run.py"]
    return_code = 124
    output = ""
    creationflags = subprocess.CREATE_NEW_PROCESS_GROUP if os.name == "nt" else 0
    process = subprocess.Popen(
        command,
        cwd=work_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        creationflags=creationflags,
    )
    try:
        output, _ = process.communicate(timeout=timeout_seconds)
        return_code = process.returncode
    except subprocess.TimeoutExpired as exc:
        if os.name == "nt":
            subprocess.run(
                ["taskkill", "/PID", str(process.pid), "/T", "/F"],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                check=False,
            )
        else:
            process.kill()
        trailing, _ = process.communicate()
        output = (exc.stdout or "") + (trailing or "") + "\nTIMEOUT: process tree terminated\n"
    log_path.write_text(output, encoding="utf-8")

    result = parse_results(work_dir / "results.csv")
    status, error, warnings = classify_run(return_code, work_dir, result)
    atoms = residues = chains = 0
    if status.startswith("success"):
        structure_path.parent.mkdir(parents=True, exist_ok=True)
        metadata_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(work_dir / "fixed.pdb", structure_path)
        shutil.copy2(work_dir / "results.csv", results_out)
        shutil.copy2(work_dir / "prefix.yaml", prefix_out)
        atoms, residues, chains = pdb_counts(structure_path)

    return {
        "model_id": current_id,
        "twist_deg": f"{twist:.1f}",
        "rise_A": f"{rise:.2f}",
        "status": status,
        "work_dir": str(work_dir.relative_to(ROOT)),
        "structure_path": str(structure_path.relative_to(ROOT)) if structure_path.exists() else "",
        "results_csv_path": str(results_out.relative_to(ROOT)) if results_out.exists() else "",
        "prefix_yaml_path": str(prefix_out.relative_to(ROOT)) if prefix_out.exists() else "",
        "log_path": str(log_path.relative_to(ROOT)),
        "accepted_candidate_count": result.get("accepted_candidate_count", "0"),
        "best_distance_A": result.get("best_distance_A", ""),
        "conformer_index": result.get("conformer_index", ""),
        "atom_count": str(atoms) if atoms else "",
        "residue_count": str(residues) if residues else "",
        "chain_count": str(chains) if chains else "",
        "return_code": str(return_code),
        "error_summary": error,
        "warnings": warnings,
    }


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    manifest = args.output_root / "pnab_twist_series_manifest.csv"
    previous = {}
    if args.resume and manifest.exists():
        with manifest.open(newline="", encoding="utf-8") as handle:
            previous = {row["model_id"]: row for row in csv.DictReader(handle)}

    requested_ids = {model_id(twist, args.rise) for twist in args.twists}
    rows = [row for key, row in previous.items() if key not in requested_ids]
    for twist in args.twists:
        current_id = model_id(twist, args.rise)
        old = previous.get(current_id)
        if old and old["status"].startswith("success") and (ROOT / old["structure_path"]).exists():
            rows.append(old)
            print(f"{current_id}: resumed")
            continue
        row = run_one(
            args.workflow_dir,
            args.output_root,
            args.rise,
            twist,
            args.conda_exe,
            args.conda_env,
            args.timeout_seconds,
        )
        rows.append(row)
        rows.sort(key=lambda row: float(row["twist_deg"]))
        write_manifest(manifest, rows)
        print(f"{current_id}: {row['status']}")
    rows.sort(key=lambda row: float(row["twist_deg"]))
    write_manifest(manifest, rows)
    print(f"Wrote {manifest}")
    return 0 if all(row["status"].startswith("success") for row in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
