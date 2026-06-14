#!/usr/bin/env python3
"""Build chain-inferred, deduplicated candidate inputs for AlphaFold/ESM-only use."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import replace
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from hexaplex_formation.geometry import group_atoms_by_residue  # noqa: E402
from hexaplex_formation.pdb_utils import PDBAtom, dedupe_exact_atoms, load_pdb_atoms, write_pdb_atoms  # noqa: E402


CHAIN_IDS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
NONSTANDARD_RESIDUES = {"CYP", "MEP"}
PROXY_SEQUENCE_MAP = {"GLU": "E"}


RESIDUE_FIELDNAMES = [
    "residue_index_in_pdb_order",
    "chain_id",
    "residue_name",
    "residue_number",
    "insertion_code",
    "residue_label",
    "proxy_residue",
    "map_warning",
]

CHAIN_FIELDNAMES = [
    "chain_id",
    "chain_index",
    "residue_count",
    "nonstandard_residue_count",
    "residue_pattern",
    "proxy_sequence",
    "first_residue_label",
    "last_residue_label",
    "contains_nonstandard_residues",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pdb",
        type=Path,
        default=Path("inputs/structures/full_hexaplex_anti_parallel_30deg_ideal.pdb"),
    )
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/intermediates/ai_candidate_inputs"))
    parser.add_argument("--chains", type=int, default=6)
    parser.add_argument("--report", type=Path, default=Path("outputs/reports/ai_candidate_inputs_report.md"))
    return parser.parse_args()


def residue_label(chain_id: str, residue_name: str, residue_number: int | None, insertion_code: str) -> str:
    number = "" if residue_number is None else str(residue_number)
    residue = f"{residue_name}{number}{insertion_code}"
    return f"{chain_id}:{residue}" if chain_id else residue


def infer_chain_index(residue_index_zero_based: int, residue_count: int, chain_count: int) -> int:
    return min(chain_count - 1, (residue_index_zero_based * chain_count) // residue_count)


def assign_chain_id(residue_index_zero_based: int, residue_count: int, chain_count: int) -> str:
    return CHAIN_IDS[infer_chain_index(residue_index_zero_based, residue_count, chain_count)]


def proxy_residue_name(residue_name: str) -> str:
    return PROXY_SEQUENCE_MAP.get(residue_name.upper(), "X")


def load_deduped_residues(pdb_path: Path) -> list[tuple[tuple[str, str, int | None, str], list[PDBAtom]]]:
    atoms = dedupe_exact_atoms(load_pdb_atoms(pdb_path))
    return list(group_atoms_by_residue(atoms).items())


def build_chain_assigned_atoms(
    residue_items: list[tuple[tuple[str, str, int | None, str], list[PDBAtom]]],
    chain_count: int,
) -> list[PDBAtom]:
    residue_count = len(residue_items)
    assigned: list[PDBAtom] = []
    for residue_index, (key, residue_atoms) in enumerate(residue_items):
        chain_id = assign_chain_id(residue_index, residue_count, chain_count)
        for atom in residue_atoms:
            assigned.append(replace(atom, chain_id=chain_id))
    return assigned


def build_residue_rows(
    residue_items: list[tuple[tuple[str, str, int | None, str], list[PDBAtom]]],
    chain_count: int,
) -> list[dict[str, str]]:
    residue_count = len(residue_items)
    rows: list[dict[str, str]] = []
    for residue_index, (key, _atoms) in enumerate(residue_items):
        chain_id = assign_chain_id(residue_index, residue_count, chain_count)
        chain_index = infer_chain_index(residue_index, residue_count, chain_count) + 1
        chain_key, residue_name, residue_number, insertion_code = key
        residue_number_text = "" if residue_number is None else str(residue_number)
        residue_name = residue_name.upper()
        row_warning = ""
        if residue_name in NONSTANDARD_RESIDUES:
            row_warning = "nonstandard_residue"
        rows.append(
            {
                "residue_index_in_pdb_order": str(residue_index + 1),
                "chain_id": chain_id,
                "residue_name": residue_name,
                "residue_number": residue_number_text,
                "insertion_code": insertion_code,
                "residue_label": residue_label(chain_id, residue_name, residue_number, insertion_code),
                "proxy_residue": proxy_residue_name(residue_name),
                "map_warning": row_warning,
            }
        )
    return rows


def build_chain_rows(residue_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    by_chain: dict[str, list[dict[str, str]]] = {}
    for row in residue_rows:
        by_chain.setdefault(row["chain_id"], []).append(row)

    for chain_index, chain_id in enumerate(sorted(by_chain), start=1):
        chain_rows = by_chain[chain_id]
        residue_pattern = "-".join(row["residue_name"] for row in chain_rows)
        proxy_sequence = "".join(row["proxy_residue"] for row in chain_rows)
        nonstandard_count = sum(1 for row in chain_rows if row["residue_name"] in NONSTANDARD_RESIDUES)
        rows.append(
            {
                "chain_id": chain_id,
                "chain_index": str(chain_index),
                "residue_count": str(len(chain_rows)),
                "nonstandard_residue_count": str(nonstandard_count),
                "residue_pattern": residue_pattern,
                "proxy_sequence": proxy_sequence,
                "first_residue_label": chain_rows[0]["residue_label"],
                "last_residue_label": chain_rows[-1]["residue_label"],
                "contains_nonstandard_residues": "yes" if nonstandard_count else "no",
            }
        )
    return rows


def write_csv(rows: list[dict[str, str]], path: Path, fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_fasta(chain_rows: list[dict[str, str]], path: Path, source_pdb: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# AlphaFold/ESM proxy FASTA for chain-inferred Hexaplex scaffold.",
        "# WARNING: CYP and MEP are nonstandard residues and are represented as X in the proxy sequence.",
        "# This file is for candidate-only exploratory use, not for biological interpretation.",
    ]
    for row in chain_rows:
        header = (
            f">chain_{row['chain_id']}|chain_index={row['chain_index']}|"
            f"residues={row['residue_count']}|nonstandard={row['nonstandard_residue_count']}|"
            f"source={source_pdb.name}"
        )
        lines.append(header)
        lines.append(row["proxy_sequence"])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_report(
    source_pdb: Path,
    deduped_pdb: Path,
    residue_rows: list[dict[str, str]],
    chain_rows: list[dict[str, str]],
    out_path: Path,
) -> None:
    residue_count = len(residue_rows)
    nonstandard_total = sum(1 for row in residue_rows if row["residue_name"] in NONSTANDARD_RESIDUES)
    notes = [
        "# AlphaFold/ESM Candidate Inputs",
        "",
        "## Scientific cautions",
        "",
        "- CYP and MEP are nonstandard residues.",
        "- The proxy FASTA uses `X` for nonstandard residues and should be treated as an exploratory placeholder only.",
        "- These inputs preserve residue identity in CSV/PDB form, but they do not make the structure a biologically standard AlphaFold/ESM sequence.",
        "- This workflow validates chain inference and residue-pattern extraction; it does not prove the true assembly pathway or biological fold.",
        "",
        "## Input summary",
        "",
        f"- Source PDB: {source_pdb}",
        f"- Deduplicated inferred-chain PDB: {deduped_pdb}",
        f"- Residues after exact deduplication: {residue_count}",
        f"- Nonstandard residues retained: {nonstandard_total}",
        f"- Inferred chains: {len(chain_rows)}",
        "",
        "## Chain patterns",
        "",
        "| chain_id | residue_count | nonstandard_residue_count | first_residue_label | last_residue_label | residue_pattern | proxy_sequence |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in chain_rows:
        notes.append(
            "| {chain_id} | {residue_count} | {nonstandard_residue_count} | {first_residue_label} | {last_residue_label} | {residue_pattern} | {proxy_sequence} |".format(
                **row
            )
        )
    notes.extend(
        [
            "",
            "## Usage warning",
            "",
            "Use these outputs only as candidate-only exploratory use for AlphaFold/ESM-oriented exploration. The nonstandard residues remain explicit in the CSV and PDB outputs, and the proxy FASTA is intentionally lossy. This does not make the structure a biologically standard AlphaFold/ESM sequence.",
            "",
        ]
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(notes), encoding="utf-8")


def write_metadata(
    source_pdb: Path,
    deduped_pdb: Path,
    residue_rows: list[dict[str, str]],
    chain_rows: list[dict[str, str]],
    out_dir: Path,
) -> None:
    metadata = {
        "source_pdb": str(source_pdb),
        "deduped_pdb": str(deduped_pdb),
        "residue_count": len(residue_rows),
        "chain_count": len(chain_rows),
        "nonstandard_residues": sorted(NONSTANDARD_RESIDUES),
        "proxy_sequence_map": PROXY_SEQUENCE_MAP,
        "warning": (
            "CYP and MEP are nonstandard residues; the proxy FASTA uses X for them and should only be used for "
            "AlphaFold/ESM candidate exploration."
        ),
    }
    (out_dir / "ai_candidate_inputs.metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    args = parse_args()
    if args.chains <= 0:
        raise SystemExit("--chains must be greater than zero")
    if args.chains > len(CHAIN_IDS):
        raise SystemExit(f"--chains must be <= {len(CHAIN_IDS)}")

    residue_items = load_deduped_residues(args.pdb)
    if len(residue_items) % args.chains != 0:
        print(
            f"WARNING: {len(residue_items)} residues are not evenly divisible by {args.chains} chains; "
            "chain assignment will use proportional residue-order partitioning.",
            file=sys.stderr,
        )

    residue_rows = build_residue_rows(residue_items, args.chains)
    chain_rows = build_chain_rows(residue_rows)
    assigned_atoms = build_chain_assigned_atoms(residue_items, args.chains)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    base_name = args.pdb.stem
    deduped_pdb = args.out_dir / f"{base_name}_deduped_{args.chains}chain.pdb"
    residue_csv = args.out_dir / f"{base_name}_chain_residue_table.csv"
    chain_csv = args.out_dir / f"{base_name}_chain_pattern_summary.csv"
    proxy_fasta = args.out_dir / f"{base_name}_alphafold_esm_proxy.fasta"

    write_pdb_atoms(assigned_atoms, deduped_pdb)
    write_csv(residue_rows, residue_csv, RESIDUE_FIELDNAMES)
    write_csv(chain_rows, chain_csv, CHAIN_FIELDNAMES)
    write_fasta(chain_rows, proxy_fasta, args.pdb)
    write_report(args.pdb, deduped_pdb, residue_rows, chain_rows, args.report)
    write_metadata(args.pdb, deduped_pdb, residue_rows, chain_rows, args.out_dir)

    print(f"Wrote {deduped_pdb}")
    print(f"Wrote {residue_csv}")
    print(f"Wrote {chain_csv}")
    print(f"Wrote {proxy_fasta}")
    print(f"Wrote {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
