"""Import John Basca/Emory corrected powder profile and compare to Nick import."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from statistics import median


NUMERIC_ROW = re.compile(
    r"^\s*([+-]?[0-9]+(?:\.[0-9]+)?)\s+([+-]?[0-9]+(?:\.[0-9]+)?)\s*$"
)


def parse_profile(path: Path) -> list[tuple[float, float]]:
    rows = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = NUMERIC_ROW.match(line)
        if match:
            rows.append((float(match.group(1)), float(match.group(2))))
    if not rows:
        raise ValueError(f"No numeric two-column rows found in {path}")
    return rows


def write_profile_csv(path: Path, rows: list[tuple[float, float]]) -> None:
    max_intensity = max(intensity for _d, intensity in rows)
    if max_intensity <= 0.0:
        raise ValueError("Maximum intensity must be positive for normalization")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["d_A", "intensity_raw", "intensity_normalized"])
        for d_spacing, intensity in rows:
            writer.writerow(
                [
                    f"{d_spacing:.15g}",
                    f"{intensity:.15g}",
                    f"{intensity / max_intensity:.15g}",
                ]
            )


def read_profile_csv(path: Path) -> list[tuple[float, float, float]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return [
            (
                float(row["d_A"]),
                float(row["intensity_raw"]),
                float(row["intensity_normalized"]),
            )
            for row in csv.DictReader(handle)
        ]


def local_maximum_closest_to(rows: list[tuple[float, float]], target_d: float) -> tuple[float, float]:
    maxima = []
    for index in range(1, len(rows) - 1):
        d_spacing, intensity = rows[index]
        if intensity >= rows[index - 1][1] and intensity >= rows[index + 1][1]:
            maxima.append((d_spacing, intensity))
    candidates = maxima or rows
    return min(candidates, key=lambda item: abs(item[0] - target_d))


def write_note(
    path: Path,
    old_profile_path: Path,
    corrected_profile_path: Path,
    old_rows: list[tuple[float, float, float]],
    corrected_rows: list[tuple[float, float]],
) -> dict[str, object]:
    old_d = [row[0] for row in old_rows]
    corrected_d = [row[0] for row in corrected_rows]
    stats: dict[str, object] = {
        "old_rows": len(old_rows),
        "corrected_rows": len(corrected_rows),
        "old_min_d": min(old_d),
        "old_max_d": max(old_d),
        "corrected_min_d": min(corrected_d),
        "corrected_max_d": max(corrected_d),
    }
    if len(old_rows) == len(corrected_rows):
        shifts = [new - old for old, new in zip(old_d, corrected_d)]
        stats.update(
            {
                "mean_shift": sum(shifts) / len(shifts),
                "median_shift": median(shifts),
                "min_shift": min(shifts),
                "max_shift": max(shifts),
            }
        )
    peak_d, peak_intensity = local_maximum_closest_to(corrected_rows, 3.4)
    stats["closest_3p4_local_max_d"] = peak_d
    stats["closest_3p4_local_max_intensity"] = peak_intensity

    mostly_0p1 = "not assessed"
    if "mean_shift" in stats:
        mostly_0p1 = (
            "yes, row-index shifts are about +0.1 A"
            if 0.07 <= float(stats["mean_shift"]) <= 0.13
            else "no, row-index shifts are not centered near +0.1 A"
        )

    lines = [
        "# Experimental Profile Correction Note",
        "",
        "John Basca at Emory supplied corrected experimental powder diffraction values via Nick's note about a correction factor.",
        "",
        "## Inputs",
        "",
        f"- Old profile path: `{old_profile_path}`",
        f"- Corrected profile path: `{corrected_profile_path}`",
        "",
        "## Profile Statistics",
        "",
        f"- Old row count: {stats['old_rows']}",
        f"- Old d-spacing range: {stats['old_min_d']:.12g} A to {stats['old_max_d']:.12g} A",
        f"- Corrected row count: {stats['corrected_rows']}",
        f"- Corrected d-spacing range: {stats['corrected_min_d']:.12g} A to {stats['corrected_max_d']:.12g} A",
    ]
    if "mean_shift" in stats:
        lines.extend(
            [
                "",
                "## Row-Index D-Spacing Shift",
                "",
                f"- Mean shift: {stats['mean_shift']:.12g} A",
                f"- Median shift: {stats['median_shift']:.12g} A",
                f"- Minimum shift: {stats['min_shift']:.12g} A",
                f"- Maximum shift: {stats['max_shift']:.12g} A",
                f"- Mostly shifted by about 0.1 A: {mostly_0p1}",
            ]
        )
    else:
        lines.extend(
            [
                "",
                "## Row-Index D-Spacing Shift",
                "",
                "- Row counts differ, so row-index shift statistics were not computed.",
            ]
        )
    lines.extend(
        [
            "",
            "## 3.4 A Region",
            "",
            f"- Corrected experimental local maximum closest to 3.4 A: {peak_d:.12g} A",
            f"- Raw intensity at that local maximum: {peak_intensity:.12g}",
            "",
            "The corrected values move peak positions modestly relative to the original Nick archive-derived profile. This note documents the correction without overclaiming structural interpretation.",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return stats


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("reference/emory_corrected_experimental_profile/corrected_experimental_values.txt"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("inputs/experimental/nick_powder_profile_corrected_emory.csv"),
    )
    parser.add_argument(
        "--old-profile",
        type=Path,
        default=Path("inputs/experimental/nick_powder_profile_original_from_nick_archive.csv"),
    )
    parser.add_argument(
        "--note",
        type=Path,
        default=Path("outputs/reports/experimental_profile_correction_note.md"),
    )
    args = parser.parse_args()

    corrected_rows = parse_profile(args.input)
    write_profile_csv(args.output, corrected_rows)
    old_rows = read_profile_csv(args.old_profile)
    stats = write_note(args.note, args.old_profile, args.output, old_rows, corrected_rows)
    print(
        f"Wrote {len(corrected_rows)} corrected rows to {args.output}; "
        f"d_A range {stats['corrected_min_d']:.6g}-{stats['corrected_max_d']:.6g}; "
        f"closest 3.4 A local max {stats['closest_3p4_local_max_d']:.6g}"
    )


if __name__ == "__main__":
    main()
