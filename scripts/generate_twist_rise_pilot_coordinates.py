"""Generate pNAB coordinates for the twist/rise pilot manifest.

This runner reads a pilot manifest with model_id, twist_deg, and rise_A columns.
For each row it:

  1. copies the Asem pNAB template workflow,
  2. patches options.yaml with the requested h_twist and h_rise,
  3. runs pNAB through conda,
  4. copies fixed.pdb to a stable coordinate path,
  5. writes an updated manifest with coordinate_file and model_status.

Use --dry-run first to validate paths and manifest wiring without running pNAB.
"""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
from pathlib import Path


DEFAULT_CONDA_EXE = Path(r"C:\Users\Public\hexaplex-tools\miniforge3\Scripts\conda.exe")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def require_columns(rows: list[dict[str, str]], path: Path) -> None:
    if not rows:
        raise ValueError(f"No rows found in {path}")

    required = {"model_id", "twist_deg", "rise_A"}
    missing = required - set(rows[0])
    if missing:
        raise ValueError(f"{path} is missing required columns: {sorted(missing)}")


def patch_options_yaml(options_path: Path, twist_deg: float, rise_A: float) -> None:
    """Patch h_twist and h_rise in an Asem-style options.yaml file."""

    text = options_path.read_text(encoding="utf-8")

    lines = []
    found_twist = False
    found_rise = False

    for line in text.splitlines():
        stripped = line.strip()

        if stripped.startswith("h_twist:"):
            indent = line[: len(line) - len(line.lstrip())]
            lines.append(f"{indent}h_twist: {twist_deg:.6f}")
            found_twist = True
        elif stripped.startswith("h_rise:"):
            indent = line[: len(line) - len(line.lstrip())]
            lines.append(f"{indent}h_rise: {rise_A:.6f}")
            found_rise = True
        else:
            lines.append(line)

    if not found_twist:
        raise ValueError(f"Could not find h_twist in {options_path}")
    if not found_rise:
        raise ValueError(f"Could not find h_rise in {options_path}")

    options_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def copy_template(template_dir: Path, work_dir: Path) -> None:
    if work_dir.exists():
        shutil.rmtree(work_dir)
    shutil.copytree(template_dir, work_dir)


def run_pnab(
    work_dir: Path,
    conda_exe: Path,
    conda_env: str,
    timeout_seconds: int,
) -> subprocess.CompletedProcess[str]:
    command = [
        str(conda_exe),
        "run",
        "-n",
        conda_env,
        "python",
        "run.py",
    ]

    return subprocess.run(
        command,
        cwd=work_dir,
        text=True,
        capture_output=True,
        timeout=timeout_seconds,
        check=False,
    )


def build_fields(rows: list[dict[str, str]]) -> list[str]:
    preferred = [
        "model_id",
        "twist_deg",
        "rise_A",
        "model_status",
        "coordinate_file",
        "diffraction_status",
        "diffraction_file",
        "score_status",
        "notes",
        "pilot_selection",
        "pilot_notes",
        "pnab_work_dir",
        "pnab_returncode",
        "pnab_stdout_tail",
        "pnab_stderr_tail",
    ]

    fields: list[str] = []
    seen: set[str] = set()

    for field in preferred:
        if any(field in row for row in rows):
            fields.append(field)
            seen.add(field)

    for row in rows:
        for field in row:
            if field not in seen:
                fields.append(field)
                seen.add(field)

    return fields


def tail(text: str, max_chars: int = 1000) -> str:
    cleaned = text.replace("\r\n", "\n").replace("\r", "\n").strip()
    if len(cleaned) <= max_chars:
        return cleaned
    return cleaned[-max_chars:]


def process_row(
    row: dict[str, str],
    template_dir: Path,
    work_root: Path,
    coordinate_dir: Path,
    conda_exe: Path,
    conda_env: str,
    timeout_seconds: int,
    dry_run: bool,
) -> dict[str, str]:
    out = dict(row)

    model_id = out["model_id"]
    twist = float(out["twist_deg"])
    rise = float(out["rise_A"])

    work_dir = work_root / model_id
    coordinate_path = coordinate_dir / f"{model_id}.pdb"

    out["pnab_work_dir"] = str(work_dir)
    out["coordinate_file"] = str(coordinate_path)

    if dry_run:
        out["model_status"] = "dry_run"
        out["notes"] = "pNAB coordinate generation dry run"
        return out

    copy_template(template_dir, work_dir)
    patch_options_yaml(work_dir / "options.yaml", twist, rise)

    coordinate_dir.mkdir(parents=True, exist_ok=True)

    try:
        result = run_pnab(
            work_dir=work_dir,
            conda_exe=conda_exe,
            conda_env=conda_env,
            timeout_seconds=timeout_seconds,
        )
    except subprocess.TimeoutExpired as exc:
        out["model_status"] = "timeout"
        out["pnab_returncode"] = ""
        out["pnab_stdout_tail"] = tail(exc.stdout or "")
        out["pnab_stderr_tail"] = tail(exc.stderr or "")
        out["notes"] = f"pNAB timed out after {timeout_seconds} seconds"
        return out

    out["pnab_returncode"] = str(result.returncode)
    out["pnab_stdout_tail"] = tail(result.stdout)
    out["pnab_stderr_tail"] = tail(result.stderr)

    fixed_pdb = work_dir / "fixed.pdb"
    if result.returncode != 0:
        out["model_status"] = "failed"
        out["notes"] = "pNAB returned non-zero exit code"
        return out

    if not fixed_pdb.exists():
        out["model_status"] = "missing_fixed_pdb"
        out["notes"] = "pNAB completed but fixed.pdb was not found"
        return out

    shutil.copy2(fixed_pdb, coordinate_path)
    out["model_status"] = "generated"
    out["notes"] = "pNAB coordinate generated"
    return out


def generate_coordinates(args: argparse.Namespace) -> list[dict[str, str]]:
    rows = read_csv(args.manifest)
    require_columns(rows, args.manifest)

    if not args.template_dir.exists():
        raise ValueError(f"Template directory does not exist: {args.template_dir}")
    if not (args.template_dir / "options.yaml").exists():
        raise ValueError(f"Template directory is missing options.yaml: {args.template_dir}")
    if not (args.template_dir / "run.py").exists():
        raise ValueError(f"Template directory is missing run.py: {args.template_dir}")

    if not args.dry_run and not args.conda_exe.exists():
        raise ValueError(f"Conda executable does not exist: {args.conda_exe}")

    selected_rows = rows
    if args.max_models:
        selected_rows = rows[: args.max_models]

    output_rows = []
    for row in selected_rows:
        output_rows.append(
            process_row(
                row=row,
                template_dir=args.template_dir,
                work_root=args.work_dir,
                coordinate_dir=args.coordinate_dir,
                conda_exe=args.conda_exe,
                conda_env=args.conda_env,
                timeout_seconds=args.timeout_seconds,
                dry_run=args.dry_run,
            )
        )

    return output_rows


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("outputs/twist_rise_scan/twist_rise_pilot_manifest.csv"),
    )
    parser.add_argument(
        "--template-dir",
        type=Path,
        default=Path("inputs/pnab_asem_30deg"),
    )
    parser.add_argument(
        "--work-dir",
        type=Path,
        default=Path("outputs/twist_rise_scan/pilot_pnab_work"),
    )
    parser.add_argument(
        "--coordinate-dir",
        type=Path,
        default=Path("outputs/twist_rise_scan/pilot_coordinates"),
    )
    parser.add_argument(
        "--output-manifest",
        type=Path,
        default=Path("outputs/twist_rise_scan/twist_rise_pilot_coordinate_manifest.csv"),
    )
    parser.add_argument(
        "--conda-exe",
        type=Path,
        default=DEFAULT_CONDA_EXE,
    )
    parser.add_argument(
        "--conda-env",
        default="pnab",
    )
    parser.add_argument(
        "--timeout-seconds",
        type=int,
        default=900,
    )
    parser.add_argument(
        "--max-models",
        type=int,
        default=0,
        help="Limit number of manifest rows processed. Use 0 for all rows.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()

    if args.timeout_seconds <= 0:
        raise ValueError("--timeout-seconds must be positive")
    if args.max_models < 0:
        raise ValueError("--max-models cannot be negative")

    rows = generate_coordinates(args)
    fields = build_fields(rows)
    write_csv(args.output_manifest, rows, fields)

    generated = sum(1 for row in rows if row.get("model_status") == "generated")
    dry_run = sum(1 for row in rows if row.get("model_status") == "dry_run")
    failed = len(rows) - generated - dry_run

    print(f"Wrote {len(rows)} rows to {args.output_manifest}")
    print(f"generated={generated} dry_run={dry_run} failed_or_other={failed}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
