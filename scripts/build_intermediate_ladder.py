#!/usr/bin/env python3
"""Build candidate Hexaplex assembly ladder structures."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from hexaplex_formation.geometry import residue_key  # noqa: E402
from hexaplex_formation.pdb_utils import PDBAtom, load_pdb_atoms, residue_count, write_pdb_atoms  # noqa: E402


SUMMARY_FIELDNAMES = [
    "model_name",
    "output_pdb",
    "atom_count",
    "residue_count",
    "included_blocks",
    "includes_hexads",
    "atom_mode",
    "source_scaffold_pdb",
    "source_hexads_pdb",
    "strand_map",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--scaffold-pdb",
        type=Path,
        default=Path("outputs/intermediates/normalized_structures/hexaplex_scaffold_only_complement_heavy_deduped.pdb"),
    )
    parser.add_argument(
        "--hexads-pdb",
        type=Path,
        default=Path("outputs/intermediates/normalized_structures/hexaplex_hexads_only_heavy_deduped.pdb"),
    )
    parser.add_argument("--strand-map", type=Path, default=Path("inputs/metadata/strand_map_candidate.csv"))
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/intermediates/ladder_structures"))
    parser.add_argument(
        "--allatom-scaffold-pdb",
        type=Path,
        default=Path("outputs/intermediates/normalized_structures/hexaplex_scaffold_only_complement_allatom_deduped.pdb"),
    )
    parser.add_argument(
        "--allatom-hexads-pdb",
        type=Path,
        default=Path("outputs/intermediates/normalized_structures/hexaplex_hexads_only_allatom_deduped.pdb"),
    )
    parser.add_argument(
        "--full-pdb",
        type=Path,
        default=Path("outputs/intermediates/normalized_structures/full_hexaplex_anti_parallel_30deg_ideal_heavy_deduped.pdb"),
    )
    parser.add_argument("--summary", type=Path, default=Path("outputs/metrics/intermediate_ladder_summary.csv"))
    return parser.parse_args()


def residue_identity_from_row(row: dict[str, str]) -> tuple[str, str, int | None, str]:
    residue_number = row["residue_number"].strip()
    return (
        row["chain_id"],
        row["residue_name"],
        int(residue_number) if residue_number else None,
        row["insertion_code"],
    )


def read_strand_map(path: Path) -> dict[int, set[tuple[str, str, int | None, str]]]:
    blocks: dict[int, set[tuple[str, str, int | None, str]]] = {}
    with path.open("r", newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            block_id = int(row["block_id"])
            blocks.setdefault(block_id, set()).add(residue_identity_from_row(row))
    return blocks


def model_blocks_name(prefix: str, block_ids: list[int], atom_mode: str) -> str:
    return f"{prefix}_{'_'.join(str(block_id) for block_id in block_ids)}_{atom_mode}"


def select_atoms_by_blocks(
    atoms: list[PDBAtom],
    blocks_by_id: dict[int, set[tuple[str, str, int | None, str]]],
    included_blocks: list[int],
) -> list[PDBAtom]:
    selected_residues: set[tuple[str, str, int | None, str]] = set()
    for block_id in included_blocks:
        selected_residues.update(blocks_by_id[block_id])
    return [atom for atom in atoms if residue_key(atom) in selected_residues]


def _write_model(
    atoms: list[PDBAtom],
    out_path: Path,
    model_name: str,
    included_blocks: str,
    includes_hexads: bool,
    atom_mode: str,
    source_scaffold_pdb: Path,
    source_hexads_pdb: Path | str,
    strand_map: Path,
) -> dict[str, str]:
    write_pdb_atoms(atoms, out_path)
    return {
        "model_name": model_name,
        "output_pdb": str(out_path),
        "atom_count": str(len(atoms)),
        "residue_count": str(residue_count(atoms)),
        "included_blocks": included_blocks,
        "includes_hexads": "yes" if includes_hexads else "no",
        "atom_mode": atom_mode,
        "source_scaffold_pdb": str(source_scaffold_pdb),
        "source_hexads_pdb": str(source_hexads_pdb),
        "strand_map": str(strand_map),
    }


def build_ladder(args: argparse.Namespace) -> list[dict[str, str]]:
    blocks_by_id = read_strand_map(args.strand_map)
    block_ids = sorted(blocks_by_id)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str]] = []

    heavy_scaffold_atoms = load_pdb_atoms(args.scaffold_pdb)
    heavy_hexad_atoms = load_pdb_atoms(args.hexads_pdb)

    for block_count in range(1, len(block_ids) + 1):
        included = block_ids[:block_count]
        included_text = ",".join(str(block_id) for block_id in included)
        scaffold_atoms = select_atoms_by_blocks(heavy_scaffold_atoms, blocks_by_id, included)

        model_name = model_blocks_name("scaffold_blocks", included, "heavy_deduped")
        rows.append(
            _write_model(
                scaffold_atoms,
                args.out_dir / f"{model_name}.pdb",
                model_name,
                included_text,
                False,
                "heavy_deduped",
                args.scaffold_pdb,
                "",
                args.strand_map,
            )
        )

        model_name = model_blocks_name("hexads_plus_scaffold_blocks", included, "heavy_deduped")
        rows.append(
            _write_model(
                [*heavy_hexad_atoms, *scaffold_atoms],
                args.out_dir / f"{model_name}.pdb",
                model_name,
                included_text,
                True,
                "heavy_deduped",
                args.scaffold_pdb,
                args.hexads_pdb,
                args.strand_map,
            )
        )

    if args.allatom_scaffold_pdb.exists() and args.allatom_hexads_pdb.exists():
        allatom_scaffold_atoms = load_pdb_atoms(args.allatom_scaffold_pdb)
        allatom_hexad_atoms = load_pdb_atoms(args.allatom_hexads_pdb)
        for block_count in range(1, len(block_ids) + 1):
            included = block_ids[:block_count]
            included_text = ",".join(str(block_id) for block_id in included)
            scaffold_atoms = select_atoms_by_blocks(allatom_scaffold_atoms, blocks_by_id, included)

            model_name = model_blocks_name("scaffold_blocks", included, "allatom_deduped")
            rows.append(
                _write_model(
                    scaffold_atoms,
                    args.out_dir / f"{model_name}.pdb",
                    model_name,
                    included_text,
                    False,
                    "allatom_deduped",
                    args.allatom_scaffold_pdb,
                    "",
                    args.strand_map,
                )
            )

            model_name = model_blocks_name("hexads_plus_scaffold_blocks", included, "allatom_deduped")
            rows.append(
                _write_model(
                    [*allatom_hexad_atoms, *scaffold_atoms],
                    args.out_dir / f"{model_name}.pdb",
                    model_name,
                    included_text,
                    True,
                    "allatom_deduped",
                    args.allatom_scaffold_pdb,
                    args.allatom_hexads_pdb,
                    args.strand_map,
                )
            )

    rows.append(
        _write_model(
            heavy_hexad_atoms,
            args.out_dir / "reference_hexads_only_heavy_deduped.pdb",
            "reference_hexads_only_heavy_deduped",
            "",
            True,
            "heavy_deduped",
            "",
            args.hexads_pdb,
            args.strand_map,
        )
    )
    rows.append(
        _write_model(
            heavy_scaffold_atoms,
            args.out_dir / "reference_scaffold_only_complement_heavy_deduped.pdb",
            "reference_scaffold_only_complement_heavy_deduped",
            ",".join(str(block_id) for block_id in block_ids),
            False,
            "heavy_deduped",
            args.scaffold_pdb,
            "",
            args.strand_map,
        )
    )
    if args.full_pdb.exists():
        full_atoms = load_pdb_atoms(args.full_pdb)
        rows.append(
            _write_model(
                full_atoms,
                args.out_dir / "reference_full_hexaplex_heavy_deduped.pdb",
                "reference_full_hexaplex_heavy_deduped",
                "",
                True,
                "heavy_deduped",
                args.scaffold_pdb,
                args.hexads_pdb,
                args.strand_map,
            )
        )

    return rows


def write_summary(rows: list[dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDNAMES, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    rows = build_ladder(args)
    write_summary(rows, args.summary)
    print(f"Wrote {len(rows)} ladder structure row(s) to {args.summary}")
    print(f"Wrote PDBs under {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
