#!/usr/bin/env python3
"""Generate pNAB helical-twist YAML variants and optionally run pNAB."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import os
import shutil
import subprocess
import sys
from pathlib import Path

import yaml


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from hexaplex_formation.pdb_utils import chain_ids, load_pdb_atoms  # noqa: E402


BROAD_TWISTS = [24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0]
NARROW_TWISTS = [27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0]

MANIFEST_FIELDNAMES = [
    "model_id",
    "twist_deg",
    "input_yaml",
    "output_structure",
    "pnab_results_dir",
    "pnab_results_csv",
    "pnab_prefix",
    "pnab_conformer_index",
    "total_energy_kcal_mol_per_nt",
    "atom_count",
    "chain_count",
    "model_status",
    "notes",
    "warnings",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-yaml", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/pnab_twist_sweep"))
    parser.add_argument("--twists", default=",".join(str(value).rstrip("0").rstrip(".") for value in BROAD_TWISTS))
    parser.add_argument("--narrow-sweep", action="store_true", help="Use 27,28,29,30,31,32,33 instead of --twists.")
    parser.add_argument("--dry-run", action="store_true", help="Write YAML/manifest but do not invoke pNAB.")
    parser.add_argument("--python-bin", default=sys.executable)
    parser.add_argument("--number-of-cpus", type=int, default=1)
    parser.add_argument("--manifest", type=Path, default=Path("outputs/metrics/pnab_twist_model_manifest.csv"))
    return parser.parse_args()


def parse_twists(text: str, narrow: bool = False) -> list[float]:
    if narrow:
        return NARROW_TWISTS
    twists = []
    for item in text.split(","):
        stripped = item.strip()
        if stripped:
            twists.append(float(stripped))
    if not twists:
        raise ValueError("at least one twist value is required")
    return twists


def twist_token(twist: float) -> str:
    if float(twist).is_integer():
        return f"{int(twist):02d}"
    return f"{twist:.2f}".replace(".", "p").replace("-", "m")


def model_id_for_twist(twist: float) -> str:
    return f"pnab_twist_{twist_token(twist)}"


def set_h_twist(options: dict, twist: float) -> dict:
    updated = yaml.safe_load(yaml.safe_dump(options, sort_keys=False))
    helical = updated.get("HelicalParameters")
    if not isinstance(helical, dict):
        raise ValueError("baseline YAML is missing HelicalParameters")
    if "h_twist" not in helical:
        raise ValueError("baseline YAML HelicalParameters is missing h_twist")
    helical["h_twist"] = [float(twist), float(twist), 1]
    return updated


def write_yaml_variant(options: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(options, handle, sort_keys=False)


def pnab_importable() -> bool:
    return importlib.util.find_spec("pnab") is not None


def run_pnab(python_bin: str, yaml_path: Path, work_dir: Path, number_of_cpus: int) -> tuple[int, str]:
    work_dir.mkdir(parents=True, exist_ok=True)
    local_yaml = work_dir / yaml_path.name
    shutil.copy2(yaml_path, local_yaml)
    code = (
        "import pnab\n"
        f"run = pnab.pNAB({local_yaml.name!r})\n"
        f"run.run(number_of_cpus={number_of_cpus}, verbose=False)\n"
    )
    completed = subprocess.run(
        [python_bin, "-c", code],
        cwd=work_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    return completed.returncode, completed.stdout


def parse_results_csv(path: Path) -> dict[str, str] | None:
    if not path.exists():
        return None
    rows: list[list[str]] = []
    with path.open("r", newline="", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            rows.append([part.strip() for part in line.split(",")])
    if not rows:
        return None
    rows.sort(key=lambda row: float(row[7]))
    row = rows[0]
    return {
        "prefix": str(int(float(row[0]))),
        "conformer_index": str(int(float(row[1]))),
        "total_energy": f"{float(row[7]):.6f}",
    }


def copy_selected_structure(result: dict[str, str], results_dir: Path, structure_out: Path) -> tuple[str, str]:
    source = results_dir / f"{result['prefix']}_{result['conformer_index']}.pdb"
    if not source.exists():
        return "", f"selected pNAB PDB not found: {source}"
    structure_out.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, structure_out)
    return str(structure_out), ""


def structure_counts(path_text: str) -> tuple[str, str, str]:
    if not path_text:
        return "", "", ""
    path = Path(path_text)
    if not path.exists():
        return "", "", f"missing selected structure: {path}"
    atoms = load_pdb_atoms(path)
    return str(len(atoms)), str(len(chain_ids(atoms))), ""


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDNAMES, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def generate_variants(args: argparse.Namespace) -> list[dict[str, str]]:
    if not args.baseline_yaml.exists():
        raise FileNotFoundError(f"baseline pNAB YAML not found: {args.baseline_yaml}")
    baseline = yaml.safe_load(args.baseline_yaml.read_text(encoding="utf-8"))
    if not isinstance(baseline, dict):
        raise ValueError("baseline YAML must contain a mapping at top level")

    twists = parse_twists(args.twists, args.narrow_sweep)
    yaml_dir = args.out_dir / "yaml"
    results_root = args.out_dir / "pnab_results"
    structures_dir = args.out_dir / "structures"
    can_run = pnab_importable() and not args.dry_run
    rows: list[dict[str, str]] = []
    for twist in twists:
        model_id = model_id_for_twist(twist)
        yaml_path = yaml_dir / f"{model_id}.yaml"
        result_dir = results_root / model_id
        structure_path = structures_dir / f"{model_id}.pdb"
        write_yaml_variant(set_h_twist(baseline, twist), yaml_path)
        warnings: list[str] = []
        notes: list[str] = ["h_twist changed; other YAML fields preserved"]
        status = "yaml_written"
        result = None
        output_structure = ""
        if can_run:
            returncode, output = run_pnab(args.python_bin, yaml_path, result_dir, args.number_of_cpus)
            (result_dir / "pnab_stdout_stderr.txt").write_text(output, encoding="utf-8")
            if returncode != 0:
                status = "pnab_failed"
                warnings.append(f"pNAB exited with status {returncode}")
            else:
                result = parse_results_csv(result_dir / "results.csv")
                if result is None:
                    status = "no_accepted_candidates"
                    warnings.append("results.csv contains no accepted candidates")
                else:
                    output_structure, warning = copy_selected_structure(result, result_dir, structure_path)
                    if warning:
                        status = "missing_selected_structure"
                        warnings.append(warning)
                    else:
                        status = "selected_structure"
        elif args.dry_run:
            status = "dry_run_yaml_only"
            notes.append("pNAB was not invoked")
        else:
            status = "pnab_not_importable"
            warnings.append("python cannot import pnab; install with conda create -n pnab -c conda-forge pnab")

        atom_count, chain_count, count_warning = structure_counts(output_structure)
        if count_warning:
            warnings.append(count_warning)
        rows.append(
            {
                "model_id": model_id,
                "twist_deg": f"{twist:.2f}",
                "input_yaml": str(yaml_path),
                "output_structure": output_structure,
                "pnab_results_dir": str(result_dir),
                "pnab_results_csv": str(result_dir / "results.csv"),
                "pnab_prefix": result["prefix"] if result else "",
                "pnab_conformer_index": result["conformer_index"] if result else "",
                "total_energy_kcal_mol_per_nt": result["total_energy"] if result else "",
                "atom_count": atom_count,
                "chain_count": chain_count,
                "model_status": status,
                "notes": "; ".join(notes),
                "warnings": "; ".join(warnings),
            }
        )
    return rows


def main() -> int:
    args = parse_args()
    rows = generate_variants(args)
    write_manifest(args.manifest, rows)
    print(f"Wrote {args.manifest}")
    for row in rows:
        print(f"{row['model_id']}: {row['model_status']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
