#!/usr/bin/env python3
"""Report pNAB availability and baseline-input readiness for twist sweeps."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import subprocess
import sys
from pathlib import Path

import yaml


REPO_ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pnab-repo", type=Path, default=Path("external/pnab"))
    parser.add_argument("--baseline-yaml", type=Path, default=None)
    parser.add_argument("--manifest", type=Path, default=Path("outputs/metrics/pnab_twist_model_manifest.csv"))
    parser.add_argument("--report", type=Path, default=Path("outputs/reports/pnab_twist_diffraction_sensitivity_report.md"))
    return parser.parse_args()


def command_output(command: list[str], cwd: Path | None = None) -> tuple[int, str]:
    try:
        completed = subprocess.run(command, cwd=cwd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
    except FileNotFoundError:
        return 127, f"{command[0]} not found"
    return completed.returncode, completed.stdout.strip()


def pnab_import_status() -> tuple[bool, str]:
    if importlib.util.find_spec("pnab") is None:
        return False, "python cannot import pnab"
    code = "import pnab; print(pnab.__file__)"
    status, output = command_output([sys.executable, "-c", code])
    return status == 0, output


def repo_commit(path: Path) -> str:
    if not (path / ".git").exists():
        return "not cloned"
    status, output = command_output(["git", "rev-parse", "HEAD"], cwd=path)
    return output if status == 0 else "unknown"


def read_manifest(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def helical_summary(path: Path | None) -> tuple[str, list[str]]:
    if path is None:
        return "no baseline YAML supplied", []
    if not path.exists():
        return f"baseline YAML not found: {path}", []
    options = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(options, dict):
        return "baseline YAML is not a top-level mapping", []
    helical = options.get("HelicalParameters", {})
    runtime = options.get("RuntimeParameters", {})
    backbone = options.get("Backbone", {})
    lines = [
        f"- YAML: {path}",
        f"- Backbone file: {backbone.get('file_path', '')}",
        f"- Strand sequence: {runtime.get('strand', '')}",
        f"- is_hexad: {runtime.get('is_hexad', '')}",
        f"- h_twist: {helical.get('h_twist', '')}",
        f"- h_rise: {helical.get('h_rise', '')}",
        f"- x_displacement: {helical.get('x_displacement', '')}",
        f"- y_displacement: {helical.get('y_displacement', '')}",
        f"- inclination: {helical.get('inclination', '')}",
        f"- tip: {helical.get('tip', '')}",
    ]
    return "baseline YAML supplied", lines


def markdown_table(rows: list[dict[str, str]], columns: list[str]) -> str:
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join("---" for _ in columns) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(row.get(column, "") for column in columns) + " |")
    return "\n".join(lines)


def write_report(args: argparse.Namespace) -> None:
    importable, import_detail = pnab_import_status()
    which_status, which_output = command_output(["which", "pnab"])
    conda_status, conda_output = command_output(["conda", "env", "list"])
    commit = repo_commit(args.pnab_repo)
    baseline_status, baseline_lines = helical_summary(args.baseline_yaml)
    manifest_rows = read_manifest(args.manifest)
    official_hexad = args.pnab_repo / "pnab/data/Hexad_Antiparallel.yaml"
    official_hexad_status, official_hexad_lines = helical_summary(official_hexad if official_hexad.exists() else None)

    lines = [
        "# pNAB Twist Diffraction Sensitivity Report",
        "",
        "## Purpose",
        "",
        "This report stages the official pNAB path for helical-twist variants around the current 30 degree hexaplex baseline. It is a forward-model sensitivity workflow, not structure determination.",
        "",
        "## pNAB Source And Availability",
        "",
        f"- Upstream checkout: {args.pnab_repo}",
        f"- Upstream commit: {commit}",
        f"- Python import status: {'available' if importable else 'not available'} ({import_detail})",
        f"- `which pnab`: {which_output if which_status == 0 else 'not on PATH'}",
        f"- `conda env list`: {'available' if conda_status == 0 else 'conda not installed/on PATH'}",
        "",
        "Install command documented by pNAB README:",
        "",
        "```bash",
        "conda create -n pnab -c conda-forge pnab",
        "conda activate pnab",
        "```",
        "",
        "No environment changes were made by this report.",
        "",
        "## pNAB YAML Format",
        "",
        "- pNAB input files are YAML mappings with `Backbone`, `HelicalParameters`, and `RuntimeParameters` sections.",
        "- `HelicalParameters.h_twist` is in degrees.",
        "- A single fixed value can be represented as `[30, 30, 1]`.",
        "- A range can be represented as `[28, 32, 5]`, which pNAB expands into uniformly spaced values including endpoints.",
        "- pNAB writes `results.csv`, `prefix.yaml`, and accepted conformer PDBs named `<Prefix>_<Conformer Index>.pdb`.",
        "",
        "## Current Baseline Input Status",
        "",
        f"- Status: {baseline_status}",
        "- Current baseline PDB: outputs/intermediates/ai_candidate_inputs/full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain.pdb",
        "- Current baseline model has CYP/MEP/GLU residues, six chains, and 15 base/GLU units per chain.",
    ]
    lines.extend(baseline_lines or ["- No current-project baseline pNAB YAML was found or supplied."])
    lines.extend(
        [
            "",
            "## Official pNAB Hexad Example",
            "",
            f"- Status: {official_hexad_status}",
        ]
    )
    lines.extend(official_hexad_lines)
    lines.extend(
        [
            "- This official example is schema-relevant, but it is not automatically assumed to be the current 30 degree baseline because its sequence/components differ from the current CYP/MEP/GLU model.",
            "",
            "## Missing Inputs Before Generating Current Twist Variants",
            "",
            "- Current 30 degree pNAB YAML or equivalent builder parameter file.",
            "- Backbone component file used for the current CYP/MEP/GLU model.",
            "- Base/nucleobase component definitions for CYP/MEP or their pNAB one-letter codes.",
            "- Sequence corresponding to 15 units per strand.",
            "- Runtime parameters, energy filters, search algorithm, random seed, strand orientation, and build_strand flags used for the current model.",
            "- Confirmation that only `h_twist` should vary while h_rise, x/y displacement, inclination, tip, sequence, and component files stay fixed.",
            "",
            "## Generated Twist Manifest",
            "",
        ]
    )
    if manifest_rows:
        lines.append(markdown_table(manifest_rows, ["model_id", "twist_deg", "input_yaml", "output_structure", "model_status", "warnings"]))
    else:
        lines.append("- No pNAB twist manifest was found yet.")
    lines.extend(
        [
            "",
            "## Diffraction Workflow",
            "",
            "- Once selected pNAB PDBs exist, run the existing normalization/heavy-atom conversion, detector simulation, and `radial_average.py` path.",
            "- Native `.npy` detector/image-plate outputs remain unchanged.",
            "- Comparisons to John Bacsa's experimental features must use radial `q_Ainv` or `d_A`, with q = 2*pi/d.",
            "",
            "## Conservative Interpretation",
            "",
            "- No twist-dependent structural or diffraction conclusion is made until pNAB-generated structures are available and validated.",
            "- Width metrics should remain simulated width proxies, not Scherrer/domain-size estimates.",
            "- The Emory data are fiber-like/oriented while current simulations are powder-averaged, limiting direct interpretation.",
            "",
            "## Next Steps",
            "",
            "1. Provide or reconstruct the exact baseline pNAB YAML and component files for the current 30 degree model.",
            "2. Run `python3 scripts/generate_pnab_twist_variants.py --baseline-yaml <baseline.yaml>` in a pNAB-enabled environment.",
            "3. Validate the regenerated 30 degree pNAB model against the current baseline PDB before interpreting twist variants.",
            "4. Run diffraction/radial averaging and feature extraction on accepted pNAB structures.",
            "5. Only after the single-parameter twist sweep works, run a small length-by-twist grid.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    write_report(args)
    print(f"Wrote {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
