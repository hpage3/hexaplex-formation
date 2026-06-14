#!/usr/bin/env python3
"""Build a combined evidence table across ladder geometry, contacts, and Debye windows."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


CONTACT_CATEGORIES = {
    "scaffold_within_block": "scaffold_within_block",
    "scaffold_between_blocks": "scaffold_between_block",
    "scaffold_hexad_or_other": "scaffold_hexad_or_other",
    "hexad_or_other_internal": "hexad_or_other_internal",
}

DIFFRACTION_COLUMNS = [
    "d_5p5_6p0_fraction",
    "d_4p5_fraction",
    "d_4p1_fraction",
    "d_3p4_fraction",
    "d_3p0_fraction",
]

OUTPUT_COLUMNS = [
    "model_name",
    "included_blocks",
    "block_count",
    "includes_hexads",
    "atom_count",
    "residue_count",
    "angular_coverage_rad_zaxis",
    "angular_coverage_rad_fitted",
    "approximate_turns",
    "approximate_pitch_per_turn",
    "axial_span",
    "z_axis_angle_degrees",
    "contact_count_4p5A",
    "scaffold_within_block_contacts",
    "scaffold_between_block_contacts",
    "scaffold_hexad_or_other_contacts",
    "hexad_or_other_internal_contacts",
    "motif_GLU_GLU",
    "motif_GLU_any",
    "scaffold_within_block_GLU_GLU",
    "scaffold_between_block_GLU_GLU",
    "scaffold_within_block_GLU_involved",
    "scaffold_between_block_GLU_involved",
    "scaffold_hexad_or_other_GLU_involved",
    "top_block_pair_by_contacts",
    "top_block_pair_contact_count",
    "top_block_pair_GLU_GLU_count",
    "d_5p5_6p0_fraction",
    "d_4p5_fraction",
    "d_4p1_fraction",
    "d_3p4_fraction",
    "d_3p0_fraction",
    "d4p5_per_GLU_GLU",
    "d4p5_per_contact",
    "scaffold_between_fraction",
    "GLU_GLU_between_fraction",
    "missing_block_contact_summary",
    "missing_block_pair_summary",
]

REPORT_COLUMNS = [
    "model_name",
    "included_blocks",
    "includes_hexads",
    "angular_coverage_rad_fitted",
    "approximate_turns",
    "contact_count_4p5A",
    "motif_GLU_GLU",
    "scaffold_between_block_contacts",
    "scaffold_hexad_or_other_contacts",
    "d_4p5_fraction",
    "d_3p4_fraction",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--ladder-comparison",
        type=Path,
        default=Path("outputs/metrics/intermediate_ladder/intermediate_ladder_comparison.csv"),
    )
    parser.add_argument(
        "--fitted-helical-comparison",
        type=Path,
        default=Path("outputs/metrics/fitted_helical_order/fitted_helical_comparison.csv"),
    )
    parser.add_argument("--block-contact-dir", type=Path, default=Path("outputs/metrics/block_contacts"))
    parser.add_argument(
        "--diffraction-comparison",
        type=Path,
        default=Path("outputs/metrics/ladder_diffraction/ladder_diffraction_comparison.csv"),
    )
    parser.add_argument("--out-csv", type=Path, default=Path("outputs/metrics/combined_evidence_table.csv"))
    parser.add_argument("--out-md", type=Path, default=Path("outputs/reports/combined_evidence_report.md"))
    return parser.parse_args()


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def normalize_model_name(value: str) -> str:
    name = Path(value).stem if "/" in value else value
    suffixes = [
        ".pdb",
        "_block_contact_summary",
        "_block_pair_summary",
        "_block_contact_decomposition",
        "_debye_profile",
        "_window_scores",
    ]
    for suffix in suffixes:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    return name


def block_count_from_included_blocks(value: str) -> str:
    parts = [part.strip() for part in value.split(",") if part.strip()]
    return str(len(parts)) if parts else ""


def safe_float(value: str | None) -> float | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def safe_int(value: str | None) -> int:
    parsed = safe_float(value)
    return int(parsed) if parsed is not None else 0


def format_float(value: float | None, digits: int = 6) -> str:
    if value is None or not math.isfinite(value):
        return ""
    return f"{value:.{digits}f}"


def divide_or_blank(numerator: str | None, denominator: str | None) -> str:
    num = safe_float(numerator)
    den = safe_float(denominator)
    if num is None or den is None or den == 0.0:
        return ""
    return format_float(num / den)


def pearson(xs: list[float], ys: list[float]) -> float | None:
    if len(xs) < 3 or len(xs) != len(ys):
        return None
    x_mean = sum(xs) / len(xs)
    y_mean = sum(ys) / len(ys)
    numerator = sum((x - x_mean) * (y - y_mean) for x, y in zip(xs, ys))
    x_var = sum((x - x_mean) ** 2 for x in xs)
    y_var = sum((y - y_mean) ** 2 for y in ys)
    if x_var == 0.0 or y_var == 0.0:
        return None
    return numerator / math.sqrt(x_var * y_var)


def index_by_model(rows: list[dict[str, str]], model_column: str = "model_name") -> dict[str, dict[str, str]]:
    indexed: dict[str, dict[str, str]] = {}
    for row in rows:
        name = row.get(model_column) or row.get("model") or row.get("pdb_file") or ""
        if name:
            indexed[normalize_model_name(name)] = row
    return indexed


def aggregate_contact_summary(path: Path) -> tuple[dict[str, str], bool]:
    fields: dict[str, str] = {}
    for prefix in CONTACT_CATEGORIES.values():
        fields[f"{prefix}_contacts"] = "0"
    for field in [
        "scaffold_within_block_GLU_involved",
        "scaffold_between_block_GLU_involved",
        "scaffold_hexad_or_other_GLU_involved",
        "scaffold_within_block_GLU_GLU",
        "scaffold_between_block_GLU_GLU",
    ]:
        fields[field] = "0"
    if not path.exists():
        return fields, True

    for row in read_csv_rows(path):
        category = CONTACT_CATEGORIES.get(row.get("contact_category", ""))
        if not category:
            continue
        fields[f"{category}_contacts"] = str(safe_int(row.get("contact_count")))
        if category in {"scaffold_within_block", "scaffold_between_block", "scaffold_hexad_or_other"}:
            fields[f"{category}_GLU_involved"] = str(safe_int(row.get("GLU_involved_count")))
        if category in {"scaffold_within_block", "scaffold_between_block"}:
            fields[f"{category}_GLU_GLU"] = str(safe_int(row.get("GLU_GLU_count")))
    return fields, False


def aggregate_pair_summary(path: Path) -> tuple[dict[str, str], bool]:
    fields = {
        "top_block_pair_by_contacts": "",
        "top_block_pair_contact_count": "",
        "top_block_pair_GLU_GLU_count": "",
    }
    rows = read_csv_rows(path)
    if not rows:
        return fields, True
    top = max(rows, key=lambda row: safe_int(row.get("contact_count")))
    fields["top_block_pair_by_contacts"] = top.get("block_pair", "")
    fields["top_block_pair_contact_count"] = str(safe_int(top.get("contact_count")))
    fields["top_block_pair_GLU_GLU_count"] = str(safe_int(top.get("GLU_GLU_count")))
    return fields, False


def add_derived_fields(row: dict[str, str]) -> None:
    row["d4p5_per_GLU_GLU"] = divide_or_blank(row.get("d_4p5_fraction"), row.get("motif_GLU_GLU"))
    row["d4p5_per_contact"] = divide_or_blank(row.get("d_4p5_fraction"), row.get("contact_count_4p5A"))

    within = safe_float(row.get("scaffold_within_block_contacts"))
    between = safe_float(row.get("scaffold_between_block_contacts"))
    if within is not None and between is not None and within + between > 0:
        row["scaffold_between_fraction"] = format_float(between / (within + between))
    else:
        row["scaffold_between_fraction"] = ""

    within_glu = safe_float(row.get("scaffold_within_block_GLU_GLU"))
    between_glu = safe_float(row.get("scaffold_between_block_GLU_GLU"))
    if within_glu is not None and between_glu is not None and within_glu + between_glu > 0:
        row["GLU_GLU_between_fraction"] = format_float(between_glu / (within_glu + between_glu))
    else:
        row["GLU_GLU_between_fraction"] = ""


def build_combined_rows(
    ladder_rows: list[dict[str, str]],
    fitted_rows: list[dict[str, str]],
    block_contact_dir: Path,
    diffraction_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    fitted_by_model = index_by_model(fitted_rows, "model")
    diffraction_by_model = index_by_model(diffraction_rows)
    combined: list[dict[str, str]] = []

    for ladder in ladder_rows:
        if ladder.get("atom_mode") and ladder.get("atom_mode") != "heavy_deduped":
            continue
        model_name = normalize_model_name(ladder.get("model_name", ""))
        if not model_name:
            continue

        row = {column: "" for column in OUTPUT_COLUMNS}
        row.update(
            {
                "model_name": model_name,
                "included_blocks": ladder.get("included_blocks", ""),
                "block_count": block_count_from_included_blocks(ladder.get("included_blocks", "")),
                "includes_hexads": ladder.get("includes_hexads", ""),
                "atom_count": ladder.get("atom_count", ""),
                "residue_count": ladder.get("residue_count", ""),
                "angular_coverage_rad_zaxis": ladder.get("angular_coverage_rad", ""),
                "contact_count_4p5A": ladder.get("contact_count_4p5A", ""),
                "motif_GLU_GLU": ladder.get("motif_GLU_GLU", ""),
                "motif_GLU_any": ladder.get("motif_GLU_any", ""),
            }
        )

        fitted = fitted_by_model.get(model_name, {})
        row["angular_coverage_rad_fitted"] = fitted.get("angular_coverage_rad", "")
        row["approximate_turns"] = fitted.get("approximate_turns", "")
        row["approximate_pitch_per_turn"] = fitted.get("approximate_pitch_per_turn", "")
        row["axial_span"] = fitted.get("axial_span", "")
        row["z_axis_angle_degrees"] = fitted.get("z_axis_angle_degrees", "")

        contact_fields, missing_contact = aggregate_contact_summary(
            block_contact_dir / f"{model_name}_block_contact_summary.csv"
        )
        row.update(contact_fields)
        row["missing_block_contact_summary"] = "yes" if missing_contact else "no"

        pair_fields, missing_pair = aggregate_pair_summary(block_contact_dir / f"{model_name}_block_pair_summary.csv")
        row.update(pair_fields)
        row["missing_block_pair_summary"] = "yes" if missing_pair else "no"

        diffraction = diffraction_by_model.get(model_name, {})
        for column in DIFFRACTION_COLUMNS:
            row[column] = diffraction.get(column, "")

        add_derived_fields(row)
        combined.append(row)

    return combined


def write_csv(rows: list[dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def markdown_table(rows: list[dict[str, str]], columns: list[str]) -> str:
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(row.get(column, "") for column in columns) + " |")
    return "\n".join(lines)


def block_count(row: dict[str, str]) -> int:
    return safe_int(row.get("block_count"))


def model_sort_key(row: dict[str, str]) -> tuple[int, int, str]:
    is_reference = 1 if row["model_name"].startswith("reference_") else 0
    return (is_reference, block_count(row), row["model_name"])


def scaffold_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return sorted(
        [
            row
            for row in rows
            if row["model_name"].startswith("scaffold_blocks_") and row.get("includes_hexads") == "no"
        ],
        key=block_count,
    )


def hexad_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return [row for row in rows if row.get("includes_hexads") == "yes"]


def generated_observations(rows: list[dict[str, str]], diffraction_available: bool) -> list[str]:
    notes: list[str] = []
    scaffolds = scaffold_rows(rows)
    block1 = next((row for row in scaffolds if block_count(row) == 1), None)
    full_scaffold = scaffolds[-1] if scaffolds else None

    if block1:
        coverage = safe_float(block1.get("angular_coverage_rad_fitted"))
        if coverage is not None:
            label = "near-full" if coverage >= 5.5 else "partial"
            notes.append(
                f"Scaffold block 1 has {label} fitted angular coverage ({coverage:.3f} rad), "
                "consistent with an individually folded/twisted path descriptor."
            )
        d4 = safe_float(block1.get("d_4p5_fraction"))
        if d4 is not None:
            notes.append(
                f"The d_4p5_fraction is already present in scaffold block 1 ({d4:.6f}); this is a "
                "reciprocal-space-window score, not a literal contact assignment."
            )
    if block1 and full_scaffold:
        first = safe_float(block1.get("d_4p5_fraction"))
        last = safe_float(full_scaffold.get("d_4p5_fraction"))
        if first is not None and last is not None:
            direction = "increases" if last > first else "does not increase"
            notes.append(f"The d_4p5_fraction {direction} from scaffold block 1 to the full scaffold row.")

        first_motif = safe_float(block1.get("motif_GLU_GLU"))
        last_motif = safe_float(full_scaffold.get("motif_GLU_GLU"))
        if first_motif is not None and last_motif is not None:
            direction = "increases" if last_motif > first_motif else "does not increase"
            notes.append(f"motif_GLU_GLU {direction} with scaffold block count in the scaffold-only ladder.")

        first_between = safe_float(block1.get("scaffold_between_block_contacts"))
        last_between = safe_float(full_scaffold.get("scaffold_between_block_contacts"))
        if first_between is not None and last_between is not None:
            direction = "increase" if last_between > first_between else "do not increase"
            notes.append(f"Between-block scaffold contacts {direction} from scaffold block 1 to the full scaffold.")

    if diffraction_available:
        hexad_values = [safe_float(row.get("d_3p4_fraction")) for row in rows if row.get("includes_hexads") == "yes"]
        scaffold_values = [safe_float(row.get("d_3p4_fraction")) for row in rows if row.get("includes_hexads") == "no"]
        hexad_numbers = [value for value in hexad_values if value is not None]
        scaffold_numbers = [value for value in scaffold_values if value is not None]
        if hexad_numbers and scaffold_numbers:
            hexad_mean = sum(hexad_numbers) / len(hexad_numbers)
            scaffold_mean = sum(scaffold_numbers) / len(scaffold_numbers)
            relation = "higher" if hexad_mean > scaffold_mean else "not higher"
            notes.append(
                f"The mean d_3p4_fraction is {relation} in hexad-containing rows "
                f"({hexad_mean:.6f}) than scaffold-only rows ({scaffold_mean:.6f})."
            )
    else:
        notes.append("Ladder diffraction-window scores were unavailable, so diffraction observations are omitted.")

    hexad_contact_rows = [
        row
        for row in hexad_rows(rows)
        if safe_float(row.get("scaffold_hexad_or_other_contacts")) is not None
        and safe_float(row.get("scaffold_hexad_or_other_contacts")) > 0
    ]
    if hexad_contact_rows:
        notes.append("Scaffold-hexad/other contacts appear in full Hexaplex and hexads-plus-scaffold models.")
    return notes


def correlation_pairs(rows: list[dict[str, str]], x_column: str, y_column: str) -> tuple[list[float], list[float]]:
    xs: list[float] = []
    ys: list[float] = []
    for row in rows:
        x_value = safe_float(row.get(x_column))
        if y_column == "includes_hexads_flag":
            y_value = 1.0 if row.get("includes_hexads") == "yes" else 0.0
        else:
            y_value = safe_float(row.get(y_column))
        if x_value is None or y_value is None:
            continue
        xs.append(x_value)
        ys.append(y_value)
    return xs, ys


def correlation_notes(rows: list[dict[str, str]]) -> list[str]:
    specs = [
        ("d_4p5_fraction", "motif_GLU_GLU", "d_4p5_fraction vs motif_GLU_GLU"),
        ("d_4p5_fraction", "motif_GLU_any", "d_4p5_fraction vs motif_GLU_any"),
        ("d_4p5_fraction", "contact_count_4p5A", "d_4p5_fraction vs contact_count_4p5A"),
        (
            "d_4p5_fraction",
            "scaffold_between_block_contacts",
            "d_4p5_fraction vs scaffold_between_block_contacts",
        ),
        (
            "d_4p5_fraction",
            "angular_coverage_rad_fitted",
            "d_4p5_fraction vs angular_coverage_rad_fitted",
        ),
        ("d_3p4_fraction", "includes_hexads_flag", "d_3p4_fraction vs includes_hexads flag"),
    ]
    notes: list[str] = []
    for x_column, y_column, label in specs:
        xs, ys = correlation_pairs(rows, x_column, y_column)
        corr = pearson(xs, ys)
        if corr is None:
            notes.append(f"{label}: unavailable due to insufficient complete numeric values.")
        else:
            notes.append(
                f"{label}: Pearson r = {corr:.3f} across {len(xs)} complete row(s). "
                "In this simplified comparative dataset, this is consistent with an association but does not establish causality."
            )
    return notes


def write_markdown_report(rows: list[dict[str, str]], path: Path, diffraction_available: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    sorted_rows = sorted(rows, key=model_sort_key)
    lines = [
        "# Combined Evidence Report: Hexaplex Formation Metrics",
        "",
        "## Scientific cautions",
        "",
        "- The d ~= 4.5 A feature is a reciprocal-space feature, not a literal atom-distance or atom-contact assignment.",
        "- Contact maps, block contact decomposition, and GLU motif counts are real-space structural summaries.",
        "- Debye scores are simplified comparative approximations, not replacements for full fiber-diffraction simulations.",
        "- Component intensities are comparative controls, not additive decompositions of full intensity.",
        "- The candidate block map remains unvalidated against PyMOL colored strand paths.",
        "- Candidate ladder models do not prove temporal assembly order.",
        "- Hydrogen-bond candidates are rough geometric candidates only.",
        "",
        "## Main evidence table",
        "",
        markdown_table(sorted_rows, REPORT_COLUMNS),
        "",
        "## Generated observations",
        "",
    ]
    lines.extend(f"- {note}" for note in generated_observations(sorted_rows, diffraction_available))
    lines.extend(
        [
            "",
            "## Correlations",
            "",
            "In this simplified comparative dataset, correlations are descriptive screening summaries. They are consistent with possible associations but do not establish causality, temporal order, or physical sufficiency.",
            "",
        ]
    )
    lines.extend(f"- {note}" for note in correlation_notes(sorted_rows))
    lines.extend(
        [
            "",
            "## Mechanistic interpretation",
            "",
            "Current evidence supports a cautious model in which individual scaffold paths are already folded/twisted and helical-like; the 4.5 A-window score is already detectable in a single scaffold path; multi-block scaffold assembly strengthens GLU-rich motif recurrence and the 4.5 A-window score; hexads contribute strongly to the 3.4 A-window score and add scaffold-hexad/other contacts; and final Hexaplex stabilization is likely cooperative.",
            "",
            "## Next steps",
            "",
            "- Validate the candidate block map against PyMOL colored strand paths.",
            "- Run full fiber-diffraction simulations on selected ladder models.",
            "- Add fitted-axis/block-specific helical metrics.",
            "- Refine hydrogen-bond candidate analysis with angle, protonation, ion, and water context.",
            "- Build a contact-state pathway model before Schrodinger bridge modeling.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    ladder_rows = read_csv_rows(args.ladder_comparison)
    fitted_rows = read_csv_rows(args.fitted_helical_comparison)
    diffraction_rows = read_csv_rows(args.diffraction_comparison)
    rows = build_combined_rows(ladder_rows, fitted_rows, args.block_contact_dir, diffraction_rows)
    write_csv(rows, args.out_csv)
    write_markdown_report(rows, args.out_md, diffraction_available=bool(diffraction_rows))
    print(f"Wrote {args.out_csv} with {len(rows)} model row(s)")
    print(f"Wrote {args.out_md}")
    if not diffraction_rows:
        print(f"Diffraction comparison unavailable: {args.diffraction_comparison}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
