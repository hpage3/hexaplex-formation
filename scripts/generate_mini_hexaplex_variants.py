#!/usr/bin/env python3
"""Generate geometry-aware N-unit mini-hexaplex truncation variants."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import OrderedDict, defaultdict
from dataclasses import dataclass
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from hexaplex_formation.geometry import (  # noqa: E402
    covariance_matrix_3d,
    distance,
    group_atoms_by_residue,
    mean_point,
    power_iteration_principal_axis,
    project_point_to_axis,
)
from hexaplex_formation.pdb_utils import (  # noqa: E402
    PDBAtom,
    atom_identity_key,
    chain_ids,
    heavy_atoms,
    load_pdb_atoms,
    write_pdb_atoms,
)


BASE_RESIDUES = {"CYP", "MEP"}
BACKBONE_ATOM_NAMES = {"N", "CA", "C", "O", "OXT"}
DEFAULT_UNITS_PER_CHAIN = 8
DEFAULT_INPUT_PDB = Path(
    "outputs/intermediates/ai_candidate_inputs/full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain.pdb"
)

MANIFEST_FIELDNAMES = [
    "variant_id",
    "truncation_rule",
    "units_per_chain",
    "residues_per_chain",
    "total_residue_count",
    "total_atom_count",
    "chains_included",
    "residue_ranges_by_chain",
    "axial_selection_rule",
    "notes",
    "warnings",
]

GEOMETRY_FIELDNAMES = [
    "variant_id",
    "units_per_chain",
    "truncation_rule",
    "chain_count",
    "residues_per_chain",
    "total_residue_count",
    "total_atom_count",
    "axial_extent_A",
    "radial_extent_min_A",
    "radial_extent_max_A",
    "radial_extent_mean_A",
    "min_heavy_atom_distance_A",
    "suspicious_overlap_count",
    "structural_coherence_flag",
    "compactness_coherence_note",
    "chain_direction_notes",
    "warnings",
]


@dataclass(frozen=True)
class ResidueRecord:
    chain_id: str
    residue_name: str
    residue_number: int | None
    insertion_code: str
    atoms: tuple[PDBAtom, ...]
    order_index: int
    axial_t: float

    @property
    def key(self) -> tuple[str, str, int | None, str]:
        return self.chain_id, self.residue_name, self.residue_number, self.insertion_code


@dataclass(frozen=True)
class RepeatUnit:
    chain_id: str
    unit_index: int
    base_residue: ResidueRecord
    glu_residue: ResidueRecord
    axial_t: float

    @property
    def residue_keys(self) -> tuple[tuple[str, str, int | None, str], tuple[str, str, int | None, str]]:
        return self.base_residue.key, self.glu_residue.key

    @property
    def residue_numbers(self) -> tuple[int | None, int | None]:
        return self.base_residue.residue_number, self.glu_residue.residue_number


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb", type=Path, default=DEFAULT_INPUT_PDB)
    parser.add_argument("--units", "--units-per-chain", dest="units_per_chain", type=int, default=DEFAULT_UNITS_PER_CHAIN)
    parser.add_argument(
        "--unit-counts",
        default="",
        help="Optional comma-separated unit counts to generate, for example 8,12. Overrides --units.",
    )
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/mini_hexaplex/structures"))
    parser.add_argument("--manifest", type=Path, default=Path("outputs/mini_hexaplex/mini_hexaplex_variant_manifest.csv"))
    parser.add_argument("--geometry-out", type=Path, default=Path("outputs/metrics/mini_hexaplex_geometry_summary.csv"))
    parser.add_argument(
        "--diagnostic-report",
        type=Path,
        default=Path("outputs/reports/mini_hexaplex_truncation_diagnostics.md"),
    )
    parser.add_argument("--overlap-threshold", type=float, default=1.0)
    parser.add_argument(
        "--include-upper-end",
        action="store_true",
        help="Also generate the upper physical end N-unit truncation.",
    )
    return parser.parse_args()


def parse_unit_counts(text: str, fallback: int) -> list[int]:
    if not text.strip():
        return [fallback]
    counts: list[int] = []
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        count = int(item)
        if count not in counts:
            counts.append(count)
    return counts


def atom_xyz(atom: PDBAtom) -> tuple[float, float, float]:
    return atom.x, atom.y, atom.z


def format_float(value: float | None, digits: int = 6) -> str:
    if value is None or not math.isfinite(value):
        return ""
    return f"{value:.{digits}f}"


def residue_label(residue: ResidueRecord) -> str:
    suffix = residue.insertion_code or ""
    return f"{residue.residue_name}{residue.residue_number if residue.residue_number is not None else '?'}{suffix}"


def infer_axis(atoms: list[PDBAtom]) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    heavy = heavy_atoms(atoms)
    backbone = [atom for atom in heavy if atom.atom_name.strip().upper() in BACKBONE_ATOM_NAMES]
    axis_atoms = backbone if len(backbone) >= 3 else heavy
    if len(axis_atoms) < 3:
        raise ValueError("At least three atoms are required to infer a fitted axis")
    origin = mean_point(atom_xyz(atom) for atom in axis_atoms)
    axis = power_iteration_principal_axis(covariance_matrix_3d(atom_xyz(atom) for atom in axis_atoms))
    return origin, axis


def residue_centroid(atoms: list[PDBAtom]) -> tuple[float, float, float]:
    return mean_point(atom_xyz(atom) for atom in atoms)


def build_residues(
    atoms: list[PDBAtom],
    origin: tuple[float, float, float],
    axis: tuple[float, float, float],
) -> OrderedDict[str, list[ResidueRecord]]:
    by_chain: OrderedDict[str, list[ResidueRecord]] = OrderedDict()
    for order_index, ((chain_id, residue_name, residue_number, insertion_code), residue_atoms) in enumerate(
        group_atoms_by_residue(atoms).items()
    ):
        centroid = residue_centroid(residue_atoms)
        axial_t, _projected = project_point_to_axis(centroid, origin, axis)
        by_chain.setdefault(chain_id, []).append(
            ResidueRecord(
                chain_id=chain_id,
                residue_name=residue_name,
                residue_number=residue_number,
                insertion_code=insertion_code,
                atoms=tuple(residue_atoms),
                order_index=order_index,
                axial_t=axial_t,
            )
        )
    return by_chain


def confirm_six_chains(atoms: list[PDBAtom]) -> list[str]:
    ids = chain_ids(atoms)
    if len(ids) != 6:
        raise ValueError(f"Expected six chains, found {len(ids)}: {','.join(ids)}")
    return ids


def group_repeat_units(by_chain: OrderedDict[str, list[ResidueRecord]]) -> tuple[dict[str, list[RepeatUnit]], list[str]]:
    units_by_chain: dict[str, list[RepeatUnit]] = {}
    warnings: list[str] = []
    for chain_id, residues in by_chain.items():
        if len(residues) % 2:
            warnings.append(f"chain {chain_id} has odd residue count {len(residues)}")
        units: list[RepeatUnit] = []
        for unit_index, start in enumerate(range(0, len(residues) - 1, 2), start=1):
            base = residues[start]
            glu = residues[start + 1]
            if base.residue_name not in BASE_RESIDUES or glu.residue_name != "GLU":
                warnings.append(
                    f"chain {chain_id} unit {unit_index} is {base.residue_name}/{glu.residue_name}, expected CYP|MEP/GLU"
                )
            units.append(
                RepeatUnit(
                    chain_id=chain_id,
                    unit_index=unit_index,
                    base_residue=base,
                    glu_residue=glu,
                    axial_t=(base.axial_t + glu.axial_t) / 2.0,
                )
            )
        units_by_chain[chain_id] = units
    return units_by_chain, warnings


def chain_direction_notes(units_by_chain: dict[str, list[RepeatUnit]]) -> dict[str, str]:
    notes: dict[str, str] = {}
    for chain_id, units in units_by_chain.items():
        if len(units) < 2:
            notes[chain_id] = "insufficient units for direction call"
            continue
        delta = units[-1].axial_t - units[0].axial_t
        if abs(delta) < 1e-6:
            direction = "flat/ambiguous"
        elif delta > 0:
            direction = "residue order aligned with increasing axial coordinate"
        else:
            direction = "residue order anti-aligned with increasing axial coordinate"
        notes[chain_id] = f"{direction}; unit1_t={units[0].axial_t:.3f}, unit{units[-1].unit_index}_t={units[-1].axial_t:.3f}"
    return notes


def select_literal_first(units_by_chain: dict[str, list[RepeatUnit]], count: int) -> dict[str, list[RepeatUnit]]:
    return {chain_id: units[:count] for chain_id, units in units_by_chain.items()}


def select_physical_end(units_by_chain: dict[str, list[RepeatUnit]], count: int, end: str = "lower") -> dict[str, list[RepeatUnit]]:
    reverse = end == "upper"
    selected: dict[str, list[RepeatUnit]] = {}
    for chain_id, units in units_by_chain.items():
        chosen = sorted(units, key=lambda unit: unit.axial_t, reverse=reverse)[:count]
        selected[chain_id] = sorted(chosen, key=lambda unit: unit.unit_index)
    return selected


def select_central(units_by_chain: dict[str, list[RepeatUnit]], count: int) -> dict[str, list[RepeatUnit]]:
    all_unit_t = [unit.axial_t for units in units_by_chain.values() for unit in units]
    center = sum(all_unit_t) / len(all_unit_t)
    selected: dict[str, list[RepeatUnit]] = {}
    for chain_id, units in units_by_chain.items():
        chosen = sorted(units, key=lambda unit: (abs(unit.axial_t - center), unit.unit_index))[:count]
        selected[chain_id] = sorted(chosen, key=lambda unit: unit.unit_index)
    return selected


def selected_residue_keys(selection: dict[str, list[RepeatUnit]]) -> set[tuple[str, str, int | None, str]]:
    keys: set[tuple[str, str, int | None, str]] = set()
    for units in selection.values():
        for unit in units:
            keys.update(unit.residue_keys)
    return keys


def filter_atoms_for_selection(atoms: list[PDBAtom], selection: dict[str, list[RepeatUnit]]) -> list[PDBAtom]:
    keys = selected_residue_keys(selection)
    return [
        atom
        for atom in atoms
        if (atom.chain_id, atom.residue_name, atom.residue_number, atom.insertion_code) in keys
    ]


def residue_ranges_by_chain(selection: dict[str, list[RepeatUnit]]) -> str:
    parts: list[str] = []
    for chain_id in sorted(selection):
        numbers = [
            number
            for unit in selection[chain_id]
            for number in unit.residue_numbers
            if number is not None
        ]
        if not numbers:
            parts.append(f"{chain_id}:")
        else:
            parts.append(f"{chain_id}:{min(numbers)}-{max(numbers)}")
    return ";".join(parts)


def residues_per_chain(selection: dict[str, list[RepeatUnit]]) -> str:
    return ";".join(f"{chain_id}:{len(units) * 2}" for chain_id, units in sorted(selection.items()))


def validate_selection(selection: dict[str, list[RepeatUnit]], units_per_chain: int) -> list[str]:
    warnings: list[str] = []
    for chain_id, units in sorted(selection.items()):
        if len(units) != units_per_chain:
            warnings.append(f"chain {chain_id} selected {len(units)} units, expected {units_per_chain}")
        for unit in units:
            if unit.base_residue.residue_name not in BASE_RESIDUES or unit.glu_residue.residue_name != "GLU":
                warnings.append(
                    f"chain {chain_id} unit {unit.unit_index} is {unit.base_residue.residue_name}/{unit.glu_residue.residue_name}"
                )
    return warnings


def structural_coherence_assessment(
    variant_id: str,
    units_per_chain: int,
    axial_extent_A: float,
    suspicious_overlap_count: int,
) -> tuple[str, str]:
    if suspicious_overlap_count > 0:
        return "not_compact/control_only", "suspicious heavy-atom overlap detected"
    if variant_id.startswith("literal"):
        return "not_compact/control_only", "sequence-order control; geometry is less physically meaningful in anti-parallel strands"
    if units_per_chain >= 7:
        return "coherent", f"six-strand geometry retained with axial extent {axial_extent_A:.2f} A"
    if units_per_chain == 6:
        return "borderline", f"six-strand geometry retained but short with axial extent {axial_extent_A:.2f} A"
    return "not_compact/control_only", f"very short truncation with axial extent {axial_extent_A:.2f} A"


def radial_and_axial_extents(
    atoms: list[PDBAtom],
    origin: tuple[float, float, float],
    axis: tuple[float, float, float],
) -> tuple[float, float, float, float]:
    axial_values: list[float] = []
    radial_values: list[float] = []
    for atom in atoms:
        axial_t, projected = project_point_to_axis(atom_xyz(atom), origin, axis)
        axial_values.append(axial_t)
        radial_values.append(distance(atom_xyz(atom), projected))
    return (
        max(axial_values) - min(axial_values),
        min(radial_values),
        max(radial_values),
        sum(radial_values) / len(radial_values),
    )


def min_heavy_distance_and_overlaps(atoms: list[PDBAtom], overlap_threshold: float) -> tuple[float, int]:
    heavy = heavy_atoms(atoms)
    min_distance = math.inf
    overlap_count = 0
    for index, atom_i in enumerate(heavy):
        for atom_j in heavy[index + 1 :]:
            atom_distance = distance(atom_i, atom_j)
            if atom_distance < min_distance:
                min_distance = atom_distance
            if atom_distance < overlap_threshold:
                overlap_count += 1
    return min_distance, overlap_count


def duplicate_atom_count(atoms: list[PDBAtom]) -> int:
    seen = set()
    duplicates = 0
    for atom in atoms:
        key = atom_identity_key(atom)
        if key in seen:
            duplicates += 1
        else:
            seen.add(key)
    return duplicates


def write_csv(rows: list[dict[str, str]], path: Path, fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def variant_definitions(units_per_chain: int, include_upper_end: bool) -> list[tuple[str, str, str]]:
    variants = [
        (
            f"literal_first{units_per_chain}_units",
            f"first {units_per_chain} base/GLU repeat units by residue order in each chain",
            f"residue-order units 1-{units_per_chain} per chain",
        ),
        (
            f"lower_end_first{units_per_chain}_units",
            f"{units_per_chain} base/GLU repeat units closest to the lower fitted-axis coordinate end in each chain",
            f"sort units by axial coordinate ascending within each chain; take first {units_per_chain}",
        ),
        (
            f"central{units_per_chain}_units",
            f"{units_per_chain} base/GLU repeat units nearest the fitted-axis center of the full hexaplex",
            f"rank units in each chain by absolute distance from full-hexaplex axial center; take nearest {units_per_chain}",
        ),
    ]
    if include_upper_end:
        variants.insert(
            2,
            (
                f"upper_end_first{units_per_chain}_units",
                f"{units_per_chain} base/GLU repeat units closest to the upper fitted-axis coordinate end in each chain",
                f"sort units by axial coordinate descending within each chain; take first {units_per_chain}",
            ),
        )
    return variants


def selection_for_variant(
    variant_id: str,
    units_by_chain: dict[str, list[RepeatUnit]],
    units_per_chain: int,
) -> dict[str, list[RepeatUnit]]:
    if variant_id.startswith("literal_first") and variant_id.endswith("_units"):
        return select_literal_first(units_by_chain, units_per_chain)
    if variant_id.startswith("lower_end_first") and variant_id.endswith("_units"):
        return select_physical_end(units_by_chain, units_per_chain, "lower")
    if variant_id.startswith("upper_end_first") and variant_id.endswith("_units"):
        return select_physical_end(units_by_chain, units_per_chain, "upper")
    if variant_id.startswith("central") and variant_id.endswith("_units"):
        return select_central(units_by_chain, units_per_chain)
    raise ValueError(f"Unknown variant id: {variant_id}")


def build_variant_rows(
    input_pdb: Path,
    output_dir: Path,
    atoms: list[PDBAtom],
    units_by_chain: dict[str, list[RepeatUnit]],
    origin: tuple[float, float, float],
    axis: tuple[float, float, float],
    units_per_chain: int,
    overlap_threshold: float,
    include_upper_end: bool = False,
) -> tuple[list[dict[str, str]], list[dict[str, str]], dict[str, dict[str, list[RepeatUnit]]]]:
    output_dir.mkdir(parents=True, exist_ok=True)
    direction_notes = chain_direction_notes(units_by_chain)
    direction_text = "; ".join(f"{chain}:{note}" for chain, note in sorted(direction_notes.items()))
    manifest_rows: list[dict[str, str]] = []
    geometry_rows: list[dict[str, str]] = []
    selections: dict[str, dict[str, list[RepeatUnit]]] = {}

    for variant_id, truncation_rule, axial_rule in variant_definitions(units_per_chain, include_upper_end):
        selection = selection_for_variant(variant_id, units_by_chain, units_per_chain)
        selections[variant_id] = selection
        variant_atoms = filter_atoms_for_selection(atoms, selection)
        output_pdb = output_dir / f"mini_hexaplex_{variant_id}.pdb"
        write_pdb_atoms(variant_atoms, output_pdb)

        validation_warnings = validate_selection(selection, units_per_chain)
        duplicate_count = duplicate_atom_count(variant_atoms)
        if duplicate_count:
            validation_warnings.append(f"{duplicate_count} duplicate atom identity records in selected output")
        min_distance, overlap_count = min_heavy_distance_and_overlaps(variant_atoms, overlap_threshold)
        if overlap_count:
            validation_warnings.append(f"{overlap_count} heavy-atom pair(s) below {overlap_threshold:.2f} A")
        axial_extent, radial_min, radial_max, radial_mean = radial_and_axial_extents(variant_atoms, origin, axis)
        total_residues = sum(len(units) * 2 for units in selection.values())
        coherence_flag, coherence_note = structural_coherence_assessment(
            variant_id,
            units_per_chain,
            axial_extent,
            overlap_count,
        )

        notes = [
            "coordinates preserved from baseline",
            "original residue numbering preserved; output atom serials regenerated by write_pdb_atoms",
            f"input_pdb={input_pdb}",
        ]
        manifest_rows.append(
            {
                "variant_id": variant_id,
                "truncation_rule": truncation_rule,
                "units_per_chain": str(units_per_chain),
                "residues_per_chain": residues_per_chain(selection),
                "total_residue_count": str(total_residues),
                "total_atom_count": str(len(variant_atoms)),
                "chains_included": ",".join(sorted(selection)),
                "residue_ranges_by_chain": residue_ranges_by_chain(selection),
                "axial_selection_rule": axial_rule,
                "notes": "; ".join(notes),
                "warnings": "; ".join(sorted(set(validation_warnings))),
            }
        )
        geometry_rows.append(
            {
                "variant_id": variant_id,
                "units_per_chain": str(units_per_chain),
                "truncation_rule": truncation_rule,
                "chain_count": str(len({atom.chain_id for atom in variant_atoms})),
                "residues_per_chain": residues_per_chain(selection),
                "total_residue_count": str(total_residues),
                "total_atom_count": str(len(variant_atoms)),
                "axial_extent_A": format_float(axial_extent),
                "radial_extent_min_A": format_float(radial_min),
                "radial_extent_max_A": format_float(radial_max),
                "radial_extent_mean_A": format_float(radial_mean),
                "min_heavy_atom_distance_A": format_float(min_distance),
                "suspicious_overlap_count": str(overlap_count),
                "structural_coherence_flag": coherence_flag,
                "compactness_coherence_note": coherence_note,
                "chain_direction_notes": direction_text,
                "warnings": "; ".join(sorted(set(validation_warnings))),
            }
        )
    return manifest_rows, geometry_rows, selections


def write_diagnostic_report(
    path: Path,
    input_pdb: Path,
    atoms: list[PDBAtom],
    by_chain: OrderedDict[str, list[ResidueRecord]],
    units_by_chain: dict[str, list[RepeatUnit]],
    selections: dict[str, dict[str, list[RepeatUnit]]],
    axis: tuple[float, float, float],
    grouping_warnings: list[str],
) -> None:
    direction_notes = chain_direction_notes(units_by_chain)
    residue_names_by_chain = {
        chain_id: "/".join(residue.residue_name for residue in residues[:6]) + ("..." if len(residues) > 6 else "")
        for chain_id, residues in by_chain.items()
    }
    lines = [
        "# Mini-Hexaplex Truncation Diagnostic Report",
        "",
        "## Baseline",
        "",
        f"- Input PDB: {input_pdb}",
        f"- Atom count: {len(atoms)}",
        f"- Chains: {', '.join(chain_ids(atoms))}",
        f"- Fitted axis direction: ({axis[0]:.6f}, {axis[1]:.6f}, {axis[2]:.6f})",
        "",
        "## Residue and unit pattern",
        "",
    ]
    for chain_id, units in sorted(units_by_chain.items()):
        first_units = ", ".join(
            f"{unit.unit_index}:{residue_label(unit.base_residue)}/{residue_label(unit.glu_residue)}@t={unit.axial_t:.3f}"
            for unit in units[:4]
        )
        last_units = ", ".join(
            f"{unit.unit_index}:{residue_label(unit.base_residue)}/{residue_label(unit.glu_residue)}@t={unit.axial_t:.3f}"
            for unit in units[-4:]
        )
        lines.extend(
            [
                f"### Chain {chain_id}",
                "",
                f"- Residues: {len(by_chain[chain_id])}",
                f"- Repeat units: {len(units)}",
                f"- Pattern start: {residue_names_by_chain[chain_id]}",
                f"- Direction: {direction_notes[chain_id]}",
                f"- First units: {first_units}",
                f"- Last units: {last_units}",
                "",
            ]
        )
    lines.extend(["## Variant selections", ""])
    for variant_id, selection in selections.items():
        lines.append(f"### {variant_id}")
        lines.append("")
        for chain_id, units in sorted(selection.items()):
            labels = ", ".join(
                f"{unit.unit_index}:{residue_label(unit.base_residue)}/{residue_label(unit.glu_residue)}@t={unit.axial_t:.3f}"
                for unit in units
            )
            lines.append(f"- Chain {chain_id}: {labels}")
        lines.append("")
    lines.extend(
        [
            "## Warnings",
            "",
            "- " + "; ".join(grouping_warnings) if grouping_warnings else "- None.",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    unit_counts = parse_unit_counts(args.unit_counts, args.units_per_chain)
    if any(count <= 0 for count in unit_counts):
        raise SystemExit("--units and --unit-counts values must be positive")
    atoms = load_pdb_atoms(args.pdb)
    confirm_six_chains(atoms)
    origin, axis = infer_axis(atoms)
    by_chain = build_residues(atoms, origin, axis)
    units_by_chain, grouping_warnings = group_repeat_units(by_chain)
    for requested_count in unit_counts:
        for chain_id, units in units_by_chain.items():
            if len(units) < requested_count:
                raise SystemExit(f"Chain {chain_id} only has {len(units)} units; cannot select {requested_count}")

    manifest_rows: list[dict[str, str]] = []
    geometry_rows: list[dict[str, str]] = []
    selections: dict[str, dict[str, list[RepeatUnit]]] = {}
    for requested_count in unit_counts:
        count_manifest, count_geometry, count_selections = build_variant_rows(
            args.pdb,
            args.out_dir,
            atoms,
            units_by_chain,
            origin,
            axis,
            requested_count,
            args.overlap_threshold,
            args.include_upper_end,
        )
        manifest_rows.extend(count_manifest)
        geometry_rows.extend(count_geometry)
        selections.update(count_selections)
    write_csv(manifest_rows, args.manifest, MANIFEST_FIELDNAMES)
    write_csv(geometry_rows, args.geometry_out, GEOMETRY_FIELDNAMES)
    write_diagnostic_report(
        args.diagnostic_report,
        args.pdb,
        atoms,
        by_chain,
        units_by_chain,
        selections,
        axis,
        grouping_warnings,
    )
    print(f"Wrote {len(manifest_rows)} mini-hexaplex PDB(s) to {args.out_dir}")
    print(f"Wrote {args.manifest}")
    print(f"Wrote {args.geometry_out}")
    print(f"Wrote {args.diagnostic_report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
