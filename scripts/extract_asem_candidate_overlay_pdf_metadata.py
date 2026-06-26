"""Extract candidate labels and Asem metrics from the overlay PDF."""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import Counter
from pathlib import Path
from typing import Iterable


DEFAULT_PDF = Path("inputs/asem_sidechains_20260625/raw/pNAB_per_candidate_overlays.pdf")
DEFAULT_INVENTORY = Path("outputs/asem_sidechains_20260625/asem_sidechain_candidate_manifest.csv")
DEFAULT_MANIFEST = Path("outputs/asem_sidechains_20260625/asem_candidate_overlay_pdf_manifest.csv")
DEFAULT_SUMMARY = Path("outputs/asem_sidechains_20260625/asem_candidate_overlay_pdf_summary.md")

MANIFEST_FIELDS = [
    "angle_deg",
    "candidate_id",
    "candidate_path",
    "pdf_page",
    "r",
    "rwp",
    "energy_kcal_mol",
    "source_pdf",
    "notes",
]

ENTRY_RE = re.compile(
    r"(?P<candidate_path>(?P<angle>\d+)/(?P<candidate>cand\d+)/initial\.pdb)"
    r"\s+r=(?P<r>[+-]?(?:\d+(?:\.\d*)?|\.\d+|nan))"
    r"\s+Rwp=(?P<rwp>[+-]?(?:\d+(?:\.\d*)?|\.\d+|nan))"
    r"\s+E=(?P<energy>[+-]?(?:\d+(?:\.\d*)?|\.\d+|nan))\s+kcal/mol",
    re.IGNORECASE,
)


def parse_overlay_text(text: str, pdf_page: int, source_pdf: Path) -> list[dict[str, str]]:
    """Parse candidate labels and reported metrics from extracted PDF text."""
    rows: list[dict[str, str]] = []
    source_pdf_text = source_pdf.as_posix()
    for match in ENTRY_RE.finditer(text):
        rows.append(
            {
                "angle_deg": match.group("angle"),
                "candidate_id": match.group("candidate"),
                "candidate_path": match.group("candidate_path"),
                "pdf_page": str(pdf_page),
                "r": match.group("r").lower(),
                "rwp": match.group("rwp").lower(),
                "energy_kcal_mol": match.group("energy").lower(),
                "source_pdf": source_pdf_text,
                "notes": (
                    "Asem PDF visual overlay; raw numeric profile data not present in this PDF"
                ),
            }
        )
    return rows


def extract_pdf_rows(pdf_path: Path) -> list[dict[str, str]]:
    """Extract rows from the PDF text layer without OCR."""
    try:
        from pypdf import PdfReader
    except ImportError as exc:  # pragma: no cover - exercised by environment, not logic.
        raise RuntimeError("pypdf is required to extract metadata from the PDF") from exc

    reader = PdfReader(str(pdf_path))
    rows: list[dict[str, str]] = []
    for page_number, page in enumerate(reader.pages, start=1):
        text = page.extract_text() or ""
        rows.extend(parse_overlay_text(text, pdf_page=page_number, source_pdf=pdf_path))
    return rows


def write_manifest(path: Path, rows: Iterable[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def _candidate_path(row: dict[str, str]) -> str:
    value = row.get("candidate_path") or row.get("structure_file") or ""
    return value.replace("\\", "/")


def load_inventory_paths(path: Path) -> set[str]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return {_candidate_path(row) for row in reader if _candidate_path(row)}


def _metric_value(value: str, *, nan_fallback: float) -> float:
    try:
        number = float(value)
    except ValueError:
        return nan_fallback
    if math.isnan(number):
        return nan_fallback
    return number


def _best_rows(
    rows: list[dict[str, str]],
    field: str,
    *,
    reverse: bool,
    nan_fallback: float,
    limit: int = 10,
) -> list[dict[str, str]]:
    return sorted(
        rows,
        key=lambda row: _metric_value(row[field], nan_fallback=nan_fallback),
        reverse=reverse,
    )[:limit]


def write_summary(path: Path, rows: list[dict[str, str]], inventory_paths: set[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pdf_paths = {_candidate_path(row) for row in rows}
    missing_imported = sorted(inventory_paths - pdf_paths)
    extra_pdf = sorted(pdf_paths - inventory_paths)
    per_angle = Counter(row["angle_deg"] for row in rows)

    best_by_rwp = _best_rows(rows, "rwp", reverse=False, nan_fallback=float("inf"))
    best_by_r = _best_rows(rows, "r", reverse=True, nan_fallback=float("-inf"))

    lines = [
        "# Asem Candidate Overlay PDF Summary",
        "",
        "Asem's PDF provides per-candidate visual Debye powder overlays against the experimental trace. "
        "It maps candidate paths to Asem-reported metrics, but it does not provide raw numeric profile data.",
        "",
        "Asem's Debye profile ranking is not the same as this repo's corrected-diffraction radial scoring.",
        "",
        f"- Candidate plots found in PDF: {len(rows)}",
        f"- Matched to imported candidate manifest: {len(pdf_paths & inventory_paths)}",
        f"- Missing imported candidates: {len(missing_imported)}",
        f"- PDF candidates not found in imported manifest: {len(extra_pdf)}",
        "",
        "## Per-Angle Counts",
        "",
        "| angle_deg | pdf_candidate_count |",
        "| --- | ---: |",
    ]
    for angle in sorted(per_angle, key=lambda value: int(value)):
        lines.append(f"| {angle} | {per_angle[angle]} |")

    lines.extend(
        [
            "",
            "## Best Candidates By Asem Rwp",
            "",
            "Lower Rwp is listed first. These are Asem-reported PDF overlay metrics, not corrected-diffraction scores.",
            "",
            "| candidate_path | r | Rwp | energy_kcal_mol | pdf_page |",
            "| --- | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in best_by_rwp:
        lines.append(
            f"| {row['candidate_path']} | {row['r']} | {row['rwp']} | "
            f"{row['energy_kcal_mol']} | {row['pdf_page']} |"
        )

    lines.extend(
        [
            "",
            "## Best Candidates By Asem r",
            "",
            "Higher r is listed first. These are Asem-reported PDF overlay metrics, not corrected-diffraction scores.",
            "",
            "| candidate_path | r | Rwp | energy_kcal_mol | pdf_page |",
            "| --- | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in best_by_r:
        lines.append(
            f"| {row['candidate_path']} | {row['r']} | {row['rwp']} | "
            f"{row['energy_kcal_mol']} | {row['pdf_page']} |"
        )

    lines.extend(["", "## Missing Imported Candidates", ""])
    lines.extend(f"- {path_value}" for path_value in missing_imported[:50])
    if len(missing_imported) > 50:
        lines.append(f"- ... {len(missing_imported) - 50} more")
    if not missing_imported:
        lines.append("None.")

    lines.extend(["", "## PDF Candidates Not Found In Imported Manifest", ""])
    lines.extend(f"- {path_value}" for path_value in extra_pdf[:50])
    if len(extra_pdf) > 50:
        lines.append(f"- ... {len(extra_pdf) - 50} more")
    if not extra_pdf:
        lines.append("None.")

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def validate_unique(rows: list[dict[str, str]]) -> None:
    seen: set[str] = set()
    duplicates: list[str] = []
    for row in rows:
        candidate_path = _candidate_path(row)
        if candidate_path in seen:
            duplicates.append(candidate_path)
        seen.add(candidate_path)
    if duplicates:
        duplicate_text = ", ".join(sorted(set(duplicates)))
        raise ValueError(f"duplicate candidate paths in PDF extraction: {duplicate_text}")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract per-candidate labels and Asem metrics from the overlay PDF."
    )
    parser.add_argument("--pdf", type=Path, default=DEFAULT_PDF)
    parser.add_argument("--inventory", type=Path, default=DEFAULT_INVENTORY)
    parser.add_argument("--output", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    rows = extract_pdf_rows(args.pdf)
    validate_unique(rows)
    inventory_paths = load_inventory_paths(args.inventory)
    write_manifest(args.output, rows)
    write_summary(args.summary, rows, inventory_paths)

    pdf_paths = {_candidate_path(row) for row in rows}
    print(f"extracted_candidates={len(rows)}")
    print(f"matched_imported_candidates={len(pdf_paths & inventory_paths)}")
    print(f"missing_imported_candidates={len(inventory_paths - pdf_paths)}")
    print(f"pdf_candidates_not_in_inventory={len(pdf_paths - inventory_paths)}")
    print(f"manifest={args.output}")
    print(f"summary={args.summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
