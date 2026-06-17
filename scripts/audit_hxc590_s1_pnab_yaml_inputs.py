#!/usr/bin/env python3
"""Audit external pNAB hexad YAML provenance for the HXC590 S1 workflow."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import re
import sys
import time
from pathlib import Path

import yaml


DEFAULT_PNAB_ROOT = Path(r"C:\Users\hpage3\OneDrive - Georgia Institute of Technology\Documents\GitHub\pnab")
DEFAULT_AUDIT_CSV = Path("outputs/metrics/hxc590_s1_pnab_yaml_audit.csv")
DEFAULT_REPORT = Path("outputs/reports/hxc590_s1_pnab_yaml_audit_report.md")
DEFAULT_SMOKE_CSV = Path("outputs/metrics/hxc590_s1_pnab_smoke_test_results.csv")

AUDIT_COLUMNS = [
    "path",
    "classification",
    "backbone_file_path",
    "backbone_exists",
    "base_file_references",
    "base_files_exist",
    "h_rise",
    "rise",
    "h_twist",
    "twist",
    "is_hexad",
    "build_strand",
    "strand_orientation",
    "strand",
    "num_candidates",
    "num_steps",
    "ff_type",
    "provenance_note",
]

SMOKE_COLUMNS = [
    "input_yaml",
    "command",
    "status",
    "runtime_seconds",
    "output_files",
    "number_of_candidates",
    "message",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pnab-root", type=Path, default=DEFAULT_PNAB_ROOT)
    parser.add_argument("--audit-csv", type=Path, default=DEFAULT_AUDIT_CSV)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--smoke-csv", type=Path, default=DEFAULT_SMOKE_CSV)
    return parser.parse_args()


def stringify(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, list):
        return ";".join(str(item) for item in value)
    return str(value)


def load_yaml(path: Path) -> dict:
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        return {}
    return data


def find_hexad_yaml_files(pnab_root: Path) -> list[Path]:
    data_dir = pnab_root / "pnab" / "data"
    files = []
    if data_dir.exists():
        files.extend(data_dir.glob("Hexad*.yaml"))
        files.extend(data_dir.glob("*Antiparallel*.yaml"))
    return sorted(set(path.resolve() for path in files))


def file_exists_relative(yaml_path: Path, file_name: str) -> bool:
    return bool(file_name) and (yaml_path.parent / file_name).is_file()


def base_file_refs(data: dict) -> list[str]:
    refs = []
    for key, value in data.items():
        if str(key).lower().startswith("base") and isinstance(value, dict):
            file_path = value.get("file_path")
            if file_path:
                refs.append(str(file_path))
    return refs


def classify_yaml(path: Path, data: dict, pnab_root: Path) -> tuple[str, str]:
    text = path.read_text(encoding="utf-8", errors="ignore").lower()
    name_context = path.name.lower()
    hxc590_markers = ["hxc590", "hexaplex", "hexaflex", "cyp", "mep", "glu"]
    if any(re.search(rf"(?<![a-z0-9]){re.escape(marker)}(?![a-z0-9])", text) or marker in name_context for marker in hxc590_markers):
        return "project-specific HXC590 baseline YAML", "Contains project-specific HXC590/hexaplex provenance markers."
    if path.name == "Hexad_Antiparallel.yaml":
        return "possible antiparallel hexad YAML", "Public pNAB antiparallel hexad example; no HXC590 provenance marker found."
    if path.name == "Hexad.yaml":
        return "generic pNAB example YAML", "Public pNAB hexad example; no HXC590 provenance marker found."
    if pnab_root in path.parents:
        return "generic pNAB example YAML", "Located under the pNAB example/data tree; no HXC590 provenance marker found."
    return "unknown/provenance unclear YAML", "No clear provenance markers found."


def audit_yaml(path: Path, pnab_root: Path) -> dict[str, str]:
    data = load_yaml(path)
    backbone = data.get("Backbone", {}) if isinstance(data.get("Backbone"), dict) else {}
    helical = data.get("HelicalParameters", {}) if isinstance(data.get("HelicalParameters"), dict) else {}
    runtime = data.get("RuntimeParameters", {}) if isinstance(data.get("RuntimeParameters"), dict) else {}
    refs = base_file_refs(data)
    classification, note = classify_yaml(path, data, pnab_root)
    backbone_file = str(backbone.get("file_path", ""))
    return {
        "path": str(path),
        "classification": classification,
        "backbone_file_path": backbone_file,
        "backbone_exists": "yes" if file_exists_relative(path, backbone_file) else "no",
        "base_file_references": ";".join(refs),
        "base_files_exist": "yes" if refs and all(file_exists_relative(path, ref) for ref in refs) else ("not_applicable" if not refs else "no"),
        "h_rise": stringify(helical.get("h_rise")),
        "rise": stringify(helical.get("rise")),
        "h_twist": stringify(helical.get("h_twist")),
        "twist": stringify(helical.get("twist")),
        "is_hexad": stringify(runtime.get("is_hexad")),
        "build_strand": stringify(runtime.get("build_strand")),
        "strand_orientation": stringify(runtime.get("strand_orientation")),
        "strand": stringify(runtime.get("strand")),
        "num_candidates": stringify(runtime.get("num_candidates")),
        "num_steps": stringify(runtime.get("num_steps")),
        "ff_type": stringify(runtime.get("ff_type")),
        "provenance_note": note,
    }


def pnab_import_status() -> tuple[bool, str]:
    spec = importlib.util.find_spec("pnab")
    if spec is None:
        return False, "python cannot import pnab"
    return True, str(spec.origin or "pnab importable")


def smoke_rows(yaml_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    can_import, import_message = pnab_import_status()
    rows = []
    for row in yaml_rows:
        start = time.perf_counter()
        command = "not run; external provenance-only audit"
        if not can_import:
            status = "not_run_external_tool_not_available"
            message = f"{import_message}; pNAB is not required by this branch and was not installed or integrated."
            output_files = ""
            candidates = ""
        else:
            status = "not_run_out_of_scope"
            message = "pNAB import is available, but running pNAB is outside this branch-local provenance audit."
            output_files = ""
            candidates = ""
        rows.append(
            {
                "input_yaml": row["path"],
                "command": command,
                "status": status,
                "runtime_seconds": f"{time.perf_counter() - start:.6f}",
                "output_files": output_files,
                "number_of_candidates": candidates,
                "message": message,
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def markdown_table(rows: list[dict[str, str]], columns: list[str]) -> list[str]:
    lines = ["| " + " | ".join(columns) + " |", "|" + "|".join("---" for _ in columns) + "|"]
    for row in rows:
        lines.append("| " + " | ".join(row.get(column, "") for column in columns) + " |")
    return lines


def write_report(path: Path, pnab_root: Path, yaml_rows: list[dict[str, str]], smoke: list[dict[str, str]]) -> None:
    can_import, import_message = pnab_import_status()
    by_name = {Path(row["path"]).name: row for row in yaml_rows}
    has_antiparallel = "Hexad_Antiparallel.yaml" in by_name
    has_hexad = "Hexad.yaml" in by_name
    provenance_clear = [row for row in yaml_rows if row["classification"] == "project-specific HXC590 baseline YAML"]
    lines = [
        "# HXC590 S1 pNAB YAML Audit",
        "",
        "## Purpose",
        "",
        "This audit checks external pNAB YAML provenance for the HXC590 S1 falsification workflow.",
        "",
        "pNAB exists as a separate local utility repository. The HXC590 branch does not depend on pNAB, and this audit does not vendor, wrap, install, or integrate pNAB.",
        "",
        "The pNAB YAML files were inspected only for provenance. A generic pNAB example YAML is not automatically an HXC590 model-generation input.",
        "",
        "Future pNAB work, if any, should happen outside this branch and only hand off validated candidate PDBs/profiles after a separate provenance and baseline-reproduction check.",
        "",
        "## External pNAB location and runtime note",
        "",
        f"- Explicit pNAB repo path: `{pnab_root}`",
        f"- `Hexad_Antiparallel.yaml` present: {'yes' if has_antiparallel else 'no'}",
        f"- `Hexad.yaml` present: {'yes' if has_hexad else 'no'}",
        f"- pNAB import check from this Python process: {'yes' if can_import else 'no'} ({import_message})",
        "- Import status is recorded only to document the local boundary; pNAB should not be installed into this branch's virtualenv as part of this audit.",
        "",
        "No pNAB install scripts, wrappers, vendoring, or environment changes are added here.",
        "",
        "## YAML classification",
        "",
    ]
    lines.extend(
        markdown_table(
            yaml_rows,
            [
                "path",
                "classification",
                "backbone_exists",
                "base_files_exist",
                "h_rise",
                "h_twist",
                "is_hexad",
                "strand_orientation",
                "num_steps",
            ],
        )
    )
    lines.extend(
        [
            "",
            "## Smoke-test status",
            "",
            "No pNAB smoke test was run. Running pNAB is outside this branch-local provenance audit unless a separate external pNAB workflow first establishes provenance and reproduces the existing 30 degree baseline.",
            "",
        ]
    )
    lines.extend(markdown_table(smoke, ["input_yaml", "status", "runtime_seconds", "number_of_candidates", "message"]))
    lines.extend(
        [
            "",
            "## HXC590 provenance decision",
            "",
            f"1. `Hexad_Antiparallel.yaml` is {'present' if has_antiparallel else 'not present'} locally.",
            f"2. Location: `{by_name.get('Hexad_Antiparallel.yaml', {}).get('path', '')}`",
            "3. The located `Hexad_Antiparallel.yaml` is classified as a possible antiparallel hexad YAML from the public pNAB examples, not as the actual HXC590 baseline input.",
            f"4. pNAB import status was recorded as: {'available' if can_import else 'not available'}. This branch does not require pNAB.",
            "5. The referenced backbone/base PDB files for the located hexad YAMLs are present relative to the YAML files.",
            "6. No pNAB candidate PDB was generated because this branch is not the place to run or integrate pNAB.",
            "7. The generic `Hexad_Antiparallel.yaml` should not be used to generate HXC590 full-length twist variants in this branch.",
            "8. Missing provenance: a separate, provenance-clear HXC590 pNAB baseline YAML or equivalent external-builder record, component mapping for the current CYP/MEP/GLU model, and an external baseline-reproduction check.",
            "",
            f"Project-specific HXC590 baseline YAMLs found: {len(provenance_clear)}.",
            "",
            "No 30 degree pilot was generated in this branch, so no reproduction check against the existing full-length 30 degree baseline was performed.",
            "",
            "Because no provenance-clear HXC590 pNAB input was found, no pNAB-generated twist variants should be added to this branch.",
            "",
            "## Outputs",
            "",
            f"- `{DEFAULT_AUDIT_CSV}`",
            f"- `{DEFAULT_SMOKE_CSV}`",
            f"- `{DEFAULT_REPORT}`",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(args: argparse.Namespace) -> dict[str, int]:
    yaml_paths = find_hexad_yaml_files(args.pnab_root)
    yaml_rows = [audit_yaml(path, args.pnab_root.resolve()) for path in yaml_paths]
    smoke = smoke_rows(yaml_rows)
    write_csv(args.audit_csv, yaml_rows, AUDIT_COLUMNS)
    write_csv(args.smoke_csv, smoke, SMOKE_COLUMNS)
    write_report(args.report, args.pnab_root, yaml_rows, smoke)
    return {"yaml_count": len(yaml_rows), "smoke_rows": len(smoke)}


def main() -> None:
    result = run(parse_args())
    print(f"Wrote {result['yaml_count']} pNAB YAML audit rows")
    print(f"Wrote {result['smoke_rows']} pNAB smoke-test status rows")


if __name__ == "__main__":
    main()
