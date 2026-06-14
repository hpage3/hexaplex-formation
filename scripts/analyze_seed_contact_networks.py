#!/usr/bin/env python3
"""Analyze formed seed interchain contact-network growth by mini-hexaplex length."""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from collections import OrderedDict, defaultdict
from dataclasses import dataclass
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from hexaplex_formation.geometry import (  # noqa: E402
    distance,
)
from hexaplex_formation.pdb_utils import (  # noqa: E402
    PDBAtom,
    chain_ids,
    dedupe_exact_atoms,
    heavy_atoms,
    is_hydrogen,
    load_pdb_atoms,
)


BASE_RESIDUES = {"CYP", "MEP"}
BACKBONE_LIKE_ATOMS = {"N", "CA", "C", "O", "OXT"}
DEFAULT_UNIT_COUNTS = [4, 5, 6, 7, 8]
DEFAULT_CUTOFFS = [4.5]

SUMMARY_COLUMNS = [
    "variant_id",
    "units_per_chain",
    "truncation_rule",
    "contact_cutoff_A",
    "chain_count",
    "total_atom_count",
    "total_interchain_contacts",
    "contacts_per_unit",
    "contacts_per_chain_pair_mean",
    "contacts_per_chain_pair_min",
    "contacts_per_chain_pair_max",
    "number_of_chain_pair_edges",
    "chain_graph_connected",
    "chain_graph_average_degree",
    "chain_graph_min_degree",
    "chain_graph_edge_density",
    "chain_graph_weighted_density",
    "number_of_connected_components",
    "largest_component_chain_count",
    "CYP_MEP_contact_count",
    "CYP_MEP_contact_fraction",
    "GLU_contact_count",
    "backbone_like_contact_count",
    "unit_graph_node_count",
    "unit_graph_edge_count",
    "unit_graph_connected_component_count",
    "unit_graph_largest_component_fraction",
    "axial_contact_span_A",
    "contact_redundancy_score",
    "nucleation_network_score",
    "perturbation_sample_count",
    "perturbation_contact_fraction_vs_reference_mean",
    "perturbation_contact_fraction_vs_reference_std",
    "perturbation_chain_graph_connected_probability",
    "perturbation_all_six_one_component_probability",
    "notes",
    "warnings",
]

EDGE_COLUMNS = [
    "variant_id",
    "units_per_chain",
    "truncation_rule",
    "contact_cutoff_A",
    "edge_scope",
    "node_a",
    "node_b",
    "chain_a",
    "chain_b",
    "unit_a",
    "unit_b",
    "contact_count",
    "CYP_MEP_vs_CYP_MEP_count",
    "CYP_MEP_vs_GLU_count",
    "GLU_vs_GLU_count",
    "backbone_like_count",
    "non_backbone_count",
]


@dataclass(frozen=True)
class Contact:
    atom_index_a: int
    atom_index_b: int
    atom_a: PDBAtom
    atom_b: PDBAtom
    chain_pair: tuple[str, str]
    unit_pair: tuple[tuple[str, int], tuple[str, int]]
    category: str
    backbone_like: bool


@dataclass(frozen=True)
class Variant:
    variant_id: str
    unit_count: int
    truncation_rule: str
    path: Path
    atoms: tuple[PDBAtom, ...]
    unit_by_residue: dict[tuple[str, str, int | None, str], int]
    warnings: tuple[str, ...]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--structures-dir", type=Path, default=Path("outputs/mini_hexaplex/structures"))
    parser.add_argument("--ensemble-dir", type=Path, default=Path("outputs/seed_formation/ensembles"))
    parser.add_argument("--unit-counts", default=",".join(str(value) for value in DEFAULT_UNIT_COUNTS))
    parser.add_argument("--cutoffs", default=",".join(str(value) for value in DEFAULT_CUTOFFS))
    parser.add_argument("--include-lower-end", action="store_true")
    parser.add_argument("--skip-perturbation-robustness", action="store_true")
    parser.add_argument("--summary-csv", type=Path, default=Path("outputs/metrics/seed_contact_network_summary.csv"))
    parser.add_argument("--edges-csv", type=Path, default=Path("outputs/metrics/seed_contact_network_edges.csv"))
    parser.add_argument("--plot-dir", type=Path, default=Path("outputs/seed_formation/plots"))
    parser.add_argument("--out-report", type=Path, default=Path("outputs/reports/seed_contact_network_report.md"))
    return parser.parse_args()


def parse_number_list(value: str, cast_type: type[int] | type[float]) -> list[int] | list[float]:
    parsed = []
    for token in value.split(","):
        stripped = token.strip()
        if stripped:
            parsed.append(cast_type(stripped))
    if not parsed:
        raise ValueError("At least one value is required")
    return parsed


def format_float(value: float | None, digits: int = 6) -> str:
    if value is None or not math.isfinite(value):
        return ""
    return f"{value:.{digits}f}"


def residue_key(atom: PDBAtom) -> tuple[str, str, int | None, str]:
    return atom.chain_id, atom.residue_name, atom.residue_number, atom.insertion_code


def residue_keys_for_chain(atoms: list[PDBAtom]) -> list[tuple[str, int | None, str]]:
    seen: OrderedDict[tuple[str, int | None, str], None] = OrderedDict()
    for atom in atoms:
        seen.setdefault((atom.residue_name, atom.residue_number, atom.insertion_code), None)
    return list(seen)


def group_by_chain(atoms: list[PDBAtom] | tuple[PDBAtom, ...]) -> OrderedDict[str, list[PDBAtom]]:
    grouped: OrderedDict[str, list[PDBAtom]] = OrderedDict()
    for atom in atoms:
        grouped.setdefault(atom.chain_id, []).append(atom)
    return grouped


def build_unit_lookup(atoms: list[PDBAtom], unit_count: int, path: Path) -> dict[tuple[str, str, int | None, str], int]:
    lookup: dict[tuple[str, str, int | None, str], int] = {}
    for chain_id, atoms_for_chain in group_by_chain(atoms).items():
        residues = residue_keys_for_chain(atoms_for_chain)
        if len(residues) != unit_count * 2:
            raise ValueError(
                f"{path} chain {chain_id} has {len(residues)} residues, expected {unit_count * 2} for {unit_count} units"
            )
        for residue_offset, (residue_name, residue_number, insertion_code) in enumerate(residues):
            lookup[(chain_id, residue_name, residue_number, insertion_code)] = residue_offset // 2 + 1
    return lookup


def load_variant(path: Path, unit_count: int, variant_id: str, truncation_rule: str) -> Variant:
    if not path.exists():
        raise FileNotFoundError(f"Mini-hexaplex PDB not found: {path}")
    raw_atoms = load_pdb_atoms(path)
    if not raw_atoms:
        raise ValueError(f"No ATOM/HETATM records found in {path}")
    atoms = dedupe_exact_atoms(raw_atoms)
    warnings: list[str] = []
    if len(atoms) != len(raw_atoms):
        warnings.append(f"removed {len(raw_atoms) - len(atoms)} exact duplicate atoms")
    chains = chain_ids(atoms)
    if len(chains) != 6:
        raise ValueError(f"{path} has {len(chains)} chains ({','.join(chains)}), expected six")
    unit_lookup = build_unit_lookup(atoms, unit_count, path)
    return Variant(
        variant_id=variant_id,
        unit_count=unit_count,
        truncation_rule=truncation_rule,
        path=path,
        atoms=tuple(atoms),
        unit_by_residue=unit_lookup,
        warnings=tuple(warnings),
    )


def residue_class(atom: PDBAtom) -> str:
    return "CYP_MEP" if atom.residue_name in BASE_RESIDUES else atom.residue_name


def contact_category(atom_a: PDBAtom, atom_b: PDBAtom) -> str:
    class_a = residue_class(atom_a)
    class_b = residue_class(atom_b)
    if class_a == "CYP_MEP" and class_b == "CYP_MEP":
        return "CYP_MEP_vs_CYP_MEP"
    if {class_a, class_b} == {"CYP_MEP", "GLU"}:
        return "CYP_MEP_vs_GLU"
    if class_a == "GLU" and class_b == "GLU":
        return "GLU_vs_GLU"
    return "other"


def is_backbone_like_contact(atom_a: PDBAtom, atom_b: PDBAtom) -> bool:
    return atom_a.atom_name.strip().upper() in BACKBONE_LIKE_ATOMS and atom_b.atom_name.strip().upper() in BACKBONE_LIKE_ATOMS


def find_contacts(atoms: tuple[PDBAtom, ...] | list[PDBAtom], unit_lookup: dict[tuple[str, str, int | None, str], int], cutoff: float) -> list[Contact]:
    indexed_heavy = [(index, atom) for index, atom in enumerate(atoms) if not is_hydrogen(atom)]
    contacts: list[Contact] = []
    seen_pairs: set[tuple[int, int]] = set()
    for left_position, (index_a, atom_a) in enumerate(indexed_heavy):
        for index_b, atom_b in indexed_heavy[left_position + 1 :]:
            if atom_a.chain_id == atom_b.chain_id:
                continue
            pair_key = (min(index_a, index_b), max(index_a, index_b))
            if pair_key in seen_pairs or distance(atom_a, atom_b) > cutoff:
                continue
            seen_pairs.add(pair_key)
            chain_pair = tuple(sorted((atom_a.chain_id, atom_b.chain_id)))
            unit_a = unit_lookup.get(residue_key(atom_a))
            unit_b = unit_lookup.get(residue_key(atom_b))
            if unit_a is None or unit_b is None:
                continue
            left_node = (atom_a.chain_id, unit_a)
            right_node = (atom_b.chain_id, unit_b)
            unit_pair = tuple(sorted((left_node, right_node)))
            contacts.append(
                Contact(
                    atom_index_a=pair_key[0],
                    atom_index_b=pair_key[1],
                    atom_a=atom_a,
                    atom_b=atom_b,
                    chain_pair=chain_pair,  # type: ignore[arg-type]
                    unit_pair=unit_pair,  # type: ignore[arg-type]
                    category=contact_category(atom_a, atom_b),
                    backbone_like=is_backbone_like_contact(atom_a, atom_b),
                )
            )
    return contacts


def connected_components(nodes: set[str] | set[tuple[str, int]], edges: set[tuple[object, object]]) -> list[set[object]]:
    adjacency: dict[object, set[object]] = {node: set() for node in nodes}
    for left, right in edges:
        adjacency.setdefault(left, set()).add(right)
        adjacency.setdefault(right, set()).add(left)
    components: list[set[object]] = []
    remaining = set(adjacency)
    while remaining:
        start = remaining.pop()
        stack = [start]
        component = {start}
        while stack:
            node = stack.pop()
            for neighbor in adjacency.get(node, set()):
                if neighbor in component:
                    continue
                component.add(neighbor)
                remaining.discard(neighbor)
                stack.append(neighbor)
        components.append(component)
    return components


def edge_category_counts(contacts: list[Contact]) -> dict[str, int]:
    counts = {
        "CYP_MEP_vs_CYP_MEP_count": 0,
        "CYP_MEP_vs_GLU_count": 0,
        "GLU_vs_GLU_count": 0,
        "backbone_like_count": 0,
        "non_backbone_count": 0,
    }
    for contact in contacts:
        if contact.category == "CYP_MEP_vs_CYP_MEP":
            counts["CYP_MEP_vs_CYP_MEP_count"] += 1
        elif contact.category == "CYP_MEP_vs_GLU":
            counts["CYP_MEP_vs_GLU_count"] += 1
        elif contact.category == "GLU_vs_GLU":
            counts["GLU_vs_GLU_count"] += 1
        if contact.backbone_like:
            counts["backbone_like_count"] += 1
        else:
            counts["non_backbone_count"] += 1
    return counts


def axial_contact_span(contacts: list[Contact], atoms: tuple[PDBAtom, ...]) -> float | None:
    if not contacts:
        return None
    values: list[float] = []
    for contact in contacts:
        values.append(contact.atom_a.z)
        values.append(contact.atom_b.z)
    return max(values) - min(values) if values else None


def summarize_network(
    variant: Variant,
    contacts: list[Contact],
    cutoff: float,
    max_contacts_per_unit: float,
    reference_axial_span: float | None,
    perturbation: dict[str, float | int | None] | None,
) -> tuple[dict[str, str], list[dict[str, str]]]:
    chain_nodes = set(chain_ids(variant.atoms))
    all_chain_pairs = {(left, right) for index, left in enumerate(sorted(chain_nodes)) for right in sorted(chain_nodes)[index + 1 :]}
    chain_edge_contacts: dict[tuple[str, str], list[Contact]] = defaultdict(list)
    unit_edge_contacts: dict[tuple[tuple[str, int], tuple[str, int]], list[Contact]] = defaultdict(list)
    for contact in contacts:
        chain_edge_contacts[contact.chain_pair].append(contact)
        unit_edge_contacts[contact.unit_pair].append(contact)

    chain_edges = set(chain_edge_contacts)
    chain_components = connected_components(chain_nodes, chain_edges)  # type: ignore[arg-type]
    degrees = {chain_id: 0 for chain_id in chain_nodes}
    for left, right in chain_edges:
        degrees[left] += 1
        degrees[right] += 1
    chain_pair_counts = [len(chain_edge_contacts.get(pair, [])) for pair in sorted(all_chain_pairs)]

    unit_nodes = {(chain_id, unit_index) for chain_id in chain_nodes for unit_index in range(1, variant.unit_count + 1)}
    unit_edges = set(unit_edge_contacts)
    unit_components = connected_components(unit_nodes, unit_edges)  # type: ignore[arg-type]
    largest_unit_component = max((len(component) for component in unit_components), default=0)

    category_counts = edge_category_counts(contacts)
    cyp_mep_count = category_counts["CYP_MEP_vs_CYP_MEP_count"] + category_counts["CYP_MEP_vs_GLU_count"]
    glu_count = category_counts["CYP_MEP_vs_GLU_count"] + category_counts["GLU_vs_GLU_count"]
    total_contacts = len(contacts)
    contacts_per_unit = total_contacts / variant.unit_count if variant.unit_count else 0.0
    axial_span = axial_contact_span(contacts, variant.atoms)
    contact_redundancy = min(1.0, contacts_per_unit / max_contacts_per_unit) if max_contacts_per_unit > 0 else None
    axial_score = (
        min(1.0, axial_span / reference_axial_span)
        if axial_span is not None and reference_axial_span is not None and reference_axial_span > 0
        else None
    )
    components = [
        1.0 if len(chain_nodes) == 6 and len(chain_components) == 1 else 0.0,
        max((len(component) for component in chain_components), default=0) / 6.0,
        (sum(degrees.values()) / len(degrees) / 5.0) if degrees else None,
        contact_redundancy,
        largest_unit_component / len(unit_nodes) if unit_nodes else None,
        min(1.0, cyp_mep_count / max(1.0, max_contacts_per_unit * variant.unit_count)),
        axial_score,
    ]
    available_components = [value for value in components if value is not None and math.isfinite(value)]
    network_score = sum(available_components) / len(available_components) if available_components else None

    perturbation = perturbation or {}
    summary = {
        "variant_id": variant.variant_id,
        "units_per_chain": str(variant.unit_count),
        "truncation_rule": variant.truncation_rule,
        "contact_cutoff_A": format_float(cutoff, 3),
        "chain_count": str(len(chain_nodes)),
        "total_atom_count": str(len(variant.atoms)),
        "total_interchain_contacts": str(total_contacts),
        "contacts_per_unit": format_float(contacts_per_unit),
        "contacts_per_chain_pair_mean": format_float(sum(chain_pair_counts) / len(chain_pair_counts) if chain_pair_counts else None),
        "contacts_per_chain_pair_min": str(min(chain_pair_counts) if chain_pair_counts else 0),
        "contacts_per_chain_pair_max": str(max(chain_pair_counts) if chain_pair_counts else 0),
        "number_of_chain_pair_edges": str(len(chain_edges)),
        "chain_graph_connected": str(len(chain_nodes) == 6 and len(chain_components) == 1),
        "chain_graph_average_degree": format_float(sum(degrees.values()) / len(degrees) if degrees else None),
        "chain_graph_min_degree": str(min(degrees.values()) if degrees else 0),
        "chain_graph_edge_density": format_float(len(chain_edges) / len(all_chain_pairs) if all_chain_pairs else None),
        "chain_graph_weighted_density": format_float(total_contacts / len(all_chain_pairs) if all_chain_pairs else None),
        "number_of_connected_components": str(len(chain_components)),
        "largest_component_chain_count": str(max((len(component) for component in chain_components), default=0)),
        "CYP_MEP_contact_count": str(cyp_mep_count),
        "CYP_MEP_contact_fraction": format_float(cyp_mep_count / total_contacts if total_contacts else None),
        "GLU_contact_count": str(glu_count),
        "backbone_like_contact_count": str(category_counts["backbone_like_count"]),
        "unit_graph_node_count": str(len(unit_nodes)),
        "unit_graph_edge_count": str(len(unit_edges)),
        "unit_graph_connected_component_count": str(len(unit_components)),
        "unit_graph_largest_component_fraction": format_float(largest_unit_component / len(unit_nodes) if unit_nodes else None),
        "axial_contact_span_A": format_float(axial_span),
        "contact_redundancy_score": format_float(contact_redundancy),
        "nucleation_network_score": format_float(network_score),
        "perturbation_sample_count": str(int(perturbation.get("sample_count") or 0)),
        "perturbation_contact_fraction_vs_reference_mean": format_float(perturbation.get("contact_fraction_mean")),  # type: ignore[arg-type]
        "perturbation_contact_fraction_vs_reference_std": format_float(perturbation.get("contact_fraction_std")),  # type: ignore[arg-type]
        "perturbation_chain_graph_connected_probability": format_float(perturbation.get("connected_probability")),  # type: ignore[arg-type]
        "perturbation_all_six_one_component_probability": format_float(perturbation.get("all_six_probability")),  # type: ignore[arg-type]
        "notes": "candidate arm/core contact classes unavailable; contacts are heavy-atom interchain pairs in static coordinates",
        "warnings": "; ".join(variant.warnings),
    }

    edge_rows: list[dict[str, str]] = []
    for (left, right), edge_contacts in sorted(chain_edge_contacts.items()):
        counts = edge_category_counts(edge_contacts)
        edge_rows.append(
            {
                "variant_id": variant.variant_id,
                "units_per_chain": str(variant.unit_count),
                "truncation_rule": variant.truncation_rule,
                "contact_cutoff_A": format_float(cutoff, 3),
                "edge_scope": "chain",
                "node_a": left,
                "node_b": right,
                "chain_a": left,
                "chain_b": right,
                "unit_a": "",
                "unit_b": "",
                "contact_count": str(len(edge_contacts)),
                **{key: str(value) for key, value in counts.items()},
            }
        )
    for (left, right), edge_contacts in sorted(unit_edge_contacts.items()):
        counts = edge_category_counts(edge_contacts)
        edge_rows.append(
            {
                "variant_id": variant.variant_id,
                "units_per_chain": str(variant.unit_count),
                "truncation_rule": variant.truncation_rule,
                "contact_cutoff_A": format_float(cutoff, 3),
                "edge_scope": "unit",
                "node_a": f"{left[0]}:{left[1]}",
                "node_b": f"{right[0]}:{right[1]}",
                "chain_a": left[0],
                "chain_b": right[0],
                "unit_a": str(left[1]),
                "unit_b": str(right[1]),
                "contact_count": str(len(edge_contacts)),
                **{key: str(value) for key, value in counts.items()},
            }
        )
    return summary, edge_rows


def load_formed_perturbed_samples(ensemble_dir: Path, unit_count: int) -> list[Path]:
    if not ensemble_dir.exists():
        return []
    pattern = re.compile(rf"central{unit_count}_formed_perturbed_\d+\.pdb$")
    return sorted(path for path in ensemble_dir.glob(f"central{unit_count}_formed_perturbed_*.pdb") if pattern.match(path.name))


def perturbation_robustness(
    variant: Variant,
    reference_contacts: list[Contact],
    cutoff: float,
    ensemble_dir: Path,
) -> dict[str, float | int | None]:
    paths = load_formed_perturbed_samples(ensemble_dir, variant.unit_count)
    if not paths:
        return {"sample_count": 0}
    reference_pairs = {(contact.atom_index_a, contact.atom_index_b) for contact in reference_contacts}
    fractions: list[float] = []
    connected_values: list[float] = []
    all_six_values: list[float] = []
    for path in paths:
        atoms = tuple(dedupe_exact_atoms(load_pdb_atoms(path)))
        if len(atoms) != len(variant.atoms):
            continue
        contacts = find_contacts(atoms, variant.unit_by_residue, cutoff)
        pairs = {(contact.atom_index_a, contact.atom_index_b) for contact in contacts}
        fractions.append(len(reference_pairs.intersection(pairs)) / len(reference_pairs) if reference_pairs else 0.0)
        chain_edges = {contact.chain_pair for contact in contacts}
        components = connected_components(set(chain_ids(atoms)), chain_edges)  # type: ignore[arg-type]
        connected = len(components) == 1
        connected_values.append(1.0 if connected else 0.0)
        all_six_values.append(1.0 if connected and max((len(component) for component in components), default=0) == 6 else 0.0)
    if not fractions:
        return {"sample_count": 0}
    mean_fraction = sum(fractions) / len(fractions)
    variance = sum((value - mean_fraction) ** 2 for value in fractions) / len(fractions)
    return {
        "sample_count": len(fractions),
        "contact_fraction_mean": mean_fraction,
        "contact_fraction_std": math.sqrt(variance),
        "connected_probability": sum(connected_values) / len(connected_values),
        "all_six_probability": sum(all_six_values) / len(all_six_values),
    }


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def safe_float(value: object) -> float | None:
    try:
        parsed = float(str(value))
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def primary_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return [row for row in rows if row["truncation_rule"] == "central" and abs(float(row["contact_cutoff_A"]) - 4.5) < 1e-6]


def write_plots(summary_rows: list[dict[str, str]], edge_rows: list[dict[str, str]], plot_dir: Path) -> list[Path]:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover - optional dependency
        print(f"WARNING: matplotlib unavailable; skipping plots: {exc}", file=sys.stderr)
        return []

    plot_dir.mkdir(parents=True, exist_ok=True)
    plot_paths: list[Path] = []
    rows = primary_rows(summary_rows)

    def line_plot(metric: str, ylabel: str, filename: str) -> None:
        points = sorted(
            (int(row["units_per_chain"]), safe_float(row.get(metric, "")))
            for row in rows
            if safe_float(row.get(metric, "")) is not None
        )
        fig, ax = plt.subplots(figsize=(7.0, 4.5))
        ax.plot([point[0] for point in points], [point[1] for point in points], marker="o")
        ax.set_xlabel("Units per chain")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.25)
        fig.tight_layout()
        path = plot_dir / filename
        fig.savefig(path, dpi=160)
        plt.close(fig)
        plot_paths.append(path)

    line_plot("total_interchain_contacts", "Total interchain contacts", "seed_contacts_units_vs_total_interchain_contacts.png")
    line_plot("contacts_per_unit", "Contacts per unit", "seed_contacts_units_vs_contacts_per_unit.png")
    line_plot("chain_graph_average_degree", "Chain graph average degree", "seed_contacts_units_vs_chain_graph_average_degree.png")
    line_plot(
        "unit_graph_largest_component_fraction",
        "Unit graph largest component fraction",
        "seed_contacts_units_vs_unit_graph_largest_component_fraction.png",
    )
    line_plot("CYP_MEP_contact_count", "CYP/MEP-involving contacts", "seed_contacts_units_vs_CYP_MEP_contact_count.png")
    line_plot("nucleation_network_score", "Exploratory nucleation network score", "seed_contacts_units_vs_nucleation_network_score.png")
    if any(row.get("perturbation_sample_count") not in {"", "0"} for row in rows):
        line_plot(
            "perturbation_chain_graph_connected_probability",
            "Perturbation connected probability",
            "seed_contacts_units_vs_perturbation_connected_probability.png",
        )

    category_by_unit: dict[int, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for row in edge_rows:
        if row["edge_scope"] != "chain" or row["truncation_rule"] != "central" or abs(float(row["contact_cutoff_A"]) - 4.5) > 1e-6:
            continue
        unit = int(row["units_per_chain"])
        for key in ["CYP_MEP_vs_CYP_MEP_count", "CYP_MEP_vs_GLU_count", "GLU_vs_GLU_count", "non_backbone_count"]:
            category_by_unit[unit][key] += int(row[key])
    if category_by_unit:
        categories = ["CYP_MEP_vs_CYP_MEP_count", "CYP_MEP_vs_GLU_count", "GLU_vs_GLU_count", "non_backbone_count"]
        labels = ["CYP/MEP-CYP/MEP", "CYP/MEP-GLU", "GLU-GLU", "Other/non-backbone"]
        units = sorted(category_by_unit)
        bottoms = [0] * len(units)
        fig, ax = plt.subplots(figsize=(8.0, 4.8))
        for key, label in zip(categories, labels):
            values = [category_by_unit[unit][key] for unit in units]
            ax.bar(units, values, bottom=bottoms, label=label)
            bottoms = [bottom + value for bottom, value in zip(bottoms, values)]
        ax.set_xlabel("Units per chain")
        ax.set_ylabel("Contact count")
        ax.legend(fontsize=8)
        fig.tight_layout()
        path = plot_dir / "seed_contacts_category_stacked_by_unit_count.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        plot_paths.append(path)
    return plot_paths


def mean_rows_by_unit(rows: list[dict[str, str]]) -> dict[int, dict[str, str]]:
    return {int(row["units_per_chain"]): row for row in primary_rows(rows)}


def write_report(summary_rows: list[dict[str, str]], plot_paths: list[Path], report_path: Path) -> None:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    rows_by_unit = mean_rows_by_unit(summary_rows)
    units = sorted(rows_by_unit)
    lines = [
        "# Seed Contact-Network Analysis",
        "",
        "## Purpose",
        "",
        "This workflow analyzes the formed mini-hexaplex seed targets as interchain contact networks. It asks whether short 4-8 unit structures already contain a connected and redundant six-chain contact graph, and whether increasing length adds a plausible nucleation threshold signal.",
        "",
        "## Contact Definition",
        "",
        "Contacts are heavy-atom, interchain atom pairs within the requested cutoff. The primary report uses 4.5 A. Duplicate atom-index pairs are excluded. Contacts are classified as CYP/MEP vs CYP/MEP, CYP/MEP vs GLU, GLU vs GLU, backbone-like, or non-backbone. Candidate arm/core contact classes are not assigned because no stable atom-selection map is available in this workflow.",
        "",
        "## Graph Definitions",
        "",
        "Chain-level graph: nodes are chains A-F. An edge exists when any interchain contact is present between a chain pair, and edge weight is the number of heavy-atom contacts.",
        "",
        "Unit-level graph: nodes are chain plus repeat-unit index. An edge exists when residues/units from different chains contact, and edge weight is the number of heavy-atom contacts.",
        "",
        "## Exploratory Score",
        "",
        "`nucleation_network_score` is a transparent average of chain connectivity, largest chain component size, normalized chain average degree, normalized contacts per unit, unit-graph largest component fraction, normalized CYP/MEP contact count, and coordinate-z contact span normalized to the 8-unit central reference. It is exploratory, not a validated free-energy or kinetic coordinate.",
        "",
        "## Summary Table",
        "",
        "| units | chain connected | contacts | contacts/unit | avg degree | unit largest frac | CYP/MEP contacts | axial span A | network score | perturb connected |",
        "|---:|:---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for unit in units:
        row = rows_by_unit[unit]
        lines.append(
            "| "
            + " | ".join(
                [
                    str(unit),
                    row["chain_graph_connected"],
                    row["total_interchain_contacts"],
                    format_float(safe_float(row["contacts_per_unit"]), 2),
                    format_float(safe_float(row["chain_graph_average_degree"]), 2),
                    format_float(safe_float(row["unit_graph_largest_component_fraction"]), 3),
                    row["CYP_MEP_contact_count"],
                    format_float(safe_float(row["axial_contact_span_A"]), 2),
                    format_float(safe_float(row["nucleation_network_score"]), 3),
                    format_float(safe_float(row["perturbation_chain_graph_connected_probability"]), 3),
                ]
            )
            + " |"
        )
    lines.extend(["", "## Plots", ""])
    for path in plot_paths:
        lines.append(f"- `{path}`")
    lines.extend(
        [
            "",
            "## Conservative Interpretation",
            "",
            "A connected six-chain contact network in the 4-unit target means the coordinate-truncated formed seed is already topologically closed by this contact definition. Additional units should therefore be interpreted mainly through redundancy, contact density, axial span, unit-graph component growth, and perturbation robustness rather than binary chain connectivity alone.",
            "",
            "If the 4-8 unit curves grow smoothly, the contact-network evidence supports gradual redundancy and axial-span growth rather than a sharp length threshold. A 7-unit coherent geometry flag is supported only if the 7-unit point shows a distinct increase in redundancy, unit-graph consolidation, axial span, or perturbation robustness compared with 4-6 units.",
            "",
            "Metrics most relevant for later Schrödinger bridge modeling are contact counts by category, chain-pair weighted edges, unit-level largest component fraction, axial contact span, and perturbation contact retention because these can define progress coordinates for closure and register.",
            "",
            "## Limitations",
            "",
            "- Static coordinate-truncated structures.",
            "- Contact cutoff sensitivity.",
            "- No energy model.",
            "- No solvent or counterions.",
            "- No dynamics.",
            "- Does not prove spontaneous formation or stability.",
            "",
        ]
    )
    report_path.write_text("\n".join(lines), encoding="utf-8")


def build_variants(structures_dir: Path, unit_counts: list[int], include_lower_end: bool) -> list[Variant]:
    variants: list[Variant] = []
    for unit_count in unit_counts:
        variants.append(
            load_variant(
                structures_dir / f"mini_hexaplex_central{unit_count}_units.pdb",
                unit_count,
                f"central{unit_count}_units",
                "central",
            )
        )
        if include_lower_end:
            lower_path = structures_dir / f"mini_hexaplex_lower_end_first{unit_count}_units.pdb"
            if lower_path.exists():
                variants.append(load_variant(lower_path, unit_count, f"lower_end_first{unit_count}_units", "lower_end_first"))
    return variants


def run(args: argparse.Namespace) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    unit_counts = parse_number_list(args.unit_counts, int)
    cutoffs = parse_number_list(args.cutoffs, float)
    variants = build_variants(args.structures_dir, unit_counts, args.include_lower_end)
    contacts_by_variant_cutoff: dict[tuple[str, float], list[Contact]] = {}
    contacts_per_unit_values: list[float] = []
    axial_spans: dict[tuple[str, float], float | None] = {}
    for variant in variants:
        for cutoff in cutoffs:
            contacts = find_contacts(variant.atoms, variant.unit_by_residue, cutoff)
            contacts_by_variant_cutoff[(variant.variant_id, cutoff)] = contacts
            contacts_per_unit_values.append(len(contacts) / variant.unit_count)
            axial_spans[(variant.variant_id, cutoff)] = axial_contact_span(contacts, variant.atoms)

    max_contacts_per_unit = max(contacts_per_unit_values) if contacts_per_unit_values else 1.0
    central_reference_spans = {
        cutoff: axial_spans.get((f"central{max(unit_counts)}_units", cutoff))
        for cutoff in cutoffs
    }
    summary_rows: list[dict[str, str]] = []
    edge_rows: list[dict[str, str]] = []
    for variant in variants:
        for cutoff in cutoffs:
            contacts = contacts_by_variant_cutoff[(variant.variant_id, cutoff)]
            perturbation = None
            if not args.skip_perturbation_robustness and variant.truncation_rule == "central":
                perturbation = perturbation_robustness(variant, contacts, cutoff, args.ensemble_dir)
            summary, edges = summarize_network(
                variant,
                contacts,
                cutoff,
                max_contacts_per_unit,
                central_reference_spans.get(cutoff),
                perturbation,
            )
            summary_rows.append(summary)
            edge_rows.extend(edges)
    write_csv(args.summary_csv, summary_rows, SUMMARY_COLUMNS)
    write_csv(args.edges_csv, edge_rows, EDGE_COLUMNS)
    plot_paths = write_plots(summary_rows, edge_rows, args.plot_dir)
    write_report(summary_rows, plot_paths, args.out_report)
    return summary_rows, edge_rows


def main() -> None:
    try:
        summary_rows, edge_rows = run(parse_args())
    except Exception as exc:
        raise SystemExit(f"ERROR: {exc}") from exc
    print(f"Wrote {len(summary_rows)} seed contact-network summary rows and {len(edge_rows)} edge rows")


if __name__ == "__main__":
    main()
