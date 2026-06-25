"""Inventory Asem side-chain candidate folders into an auditable CSV manifest."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


STRUCTURE_EXTENSIONS = (".pdb", ".xyz", ".mol2", ".prmtop", ".inpcrd")
FIELDS = [
    "angle_deg",
    "candidate_id",
    "candidate_dir",
    "structure_file",
    "structure_file_type",
    "structure_file_size",
    "atom_count",
    "has_pdb",
    "has_xyz",
    "has_mol2",
    "has_prmtop",
    "has_inpcrd",
    "notes",
]


def candidate_number(path: Path) -> tuple[int, str]:
    match = re.fullmatch(r"cand(\d+)", path.name, flags=re.IGNORECASE)
    if match:
        return int(match.group(1)), path.name
    return 10**9, path.name


def count_atoms(path: Path) -> int | None:
    suffix = path.suffix.lower()
    if suffix == ".pdb":
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            return sum(
                1
                for line in handle
                if line.startswith("ATOM  ") or line.startswith("HETATM")
            )
    if suffix == ".xyz":
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            first_line = handle.readline().strip()
        try:
            return int(first_line)
        except ValueError:
            return None
    return None


def relative_text(path: Path, root: Path) -> str:
    return str(path.relative_to(root)).replace("/", "\\")


def inventory_candidate(candidate_dir: Path, raw_root: Path, angle: int) -> dict[str, str]:
    files = sorted(path for path in candidate_dir.iterdir() if path.is_file())
    by_suffix = {
        suffix: [path for path in files if path.suffix.lower() == suffix]
        for suffix in STRUCTURE_EXTENSIONS
    }
    structure_files = [
        path
        for suffix in STRUCTURE_EXTENSIONS
        for path in by_suffix[suffix]
    ]
    primary = structure_files[0] if structure_files else None
    atom_count = count_atoms(primary) if primary else None
    notes = []
    if not primary:
        notes.append("no recognized structure file")
    elif len(structure_files) > 1:
        notes.append(f"multiple recognized structure files: {len(structure_files)}")
    unrecognized = [path.name for path in files if path not in structure_files]
    if unrecognized:
        notes.append("other files: " + ";".join(unrecognized))

    return {
        "angle_deg": str(angle),
        "candidate_id": candidate_dir.name,
        "candidate_dir": relative_text(candidate_dir, raw_root),
        "structure_file": relative_text(primary, raw_root) if primary else "",
        "structure_file_type": primary.suffix.lower().lstrip(".") if primary else "",
        "structure_file_size": str(primary.stat().st_size) if primary else "",
        "atom_count": str(atom_count) if atom_count is not None else "",
        "has_pdb": "yes" if by_suffix[".pdb"] else "no",
        "has_xyz": "yes" if by_suffix[".xyz"] else "no",
        "has_mol2": "yes" if by_suffix[".mol2"] else "no",
        "has_prmtop": "yes" if by_suffix[".prmtop"] else "no",
        "has_inpcrd": "yes" if by_suffix[".inpcrd"] else "no",
        "notes": "; ".join(notes),
    }


def inventory(raw_root: Path) -> list[dict[str, str]]:
    if not raw_root.is_dir():
        raise ValueError(f"Raw import directory does not exist: {raw_root}")
    rows = []
    angle_dirs = sorted(
        (path for path in raw_root.iterdir() if path.is_dir() and path.name.isdigit()),
        key=lambda path: int(path.name),
    )
    for angle_dir in angle_dirs:
        candidates = sorted(
            (path for path in angle_dir.iterdir() if path.is_dir()),
            key=candidate_number,
        )
        for candidate_dir in candidates:
            rows.append(inventory_candidate(candidate_dir, raw_root, int(angle_dir.name)))
    return rows


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--raw-root",
        type=Path,
        default=Path("inputs/asem_sidechains_20260625/raw"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "outputs/asem_sidechains_20260625/asem_sidechain_candidate_manifest.csv"
        ),
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    rows = inventory(args.raw_root)
    if not rows:
        raise ValueError(f"No candidate folders found under {args.raw_root}")
    write_manifest(args.output, rows)
    print(f"Wrote {len(rows)} candidate rows to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
