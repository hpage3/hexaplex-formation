"""Extract target-window peaks from twist/rise radial diffraction profiles.

This script converts profile data into the long-form observed peak table used by
the twist/rise scoring pipeline.

Input profile CSV must contain:

    model_id,d_A,intensity

Input target CSV must contain:

    target_label,target_d_A,target_group,notes

For each model_id and each target_label, the script finds the highest intensity
point whose d_A lies within:

    target_d_A +/- window_half_width

Output CSV:

    model_id,target_label,observed_d_A,observed_intensity,target_d_A,
    delta_d_A,window_min_d_A,window_max_d_A,peak_status
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


REQUIRED_PROFILE_COLUMNS = {"model_id", "d_A", "intensity"}
REQUIRED_TARGET_COLUMNS = {"target_label", "target_d_A"}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def require_columns(rows: list[dict[str, str]], required: set[str], path: Path) -> None:
    if not rows:
        raise ValueError(f"No rows found in {path}")

    missing = required - set(rows[0])
    if missing:
        raise ValueError(f"{path} is missing required columns: {sorted(missing)}")


def load_targets(path: Path) -> list[dict[str, str]]:
    rows = read_csv(path)
    require_columns(rows, REQUIRED_TARGET_COLUMNS, path)

    targets: list[dict[str, str]] = []
    seen: set[str] = set()

    for index, row in enumerate(rows, start=2):
        label = row["target_label"].strip()
        if not label:
            raise ValueError(f"{path} row {index} has empty target_label")
        if label in seen:
            raise ValueError(f"Duplicate target_label {label!r} in {path}")
        seen.add(label)

        target_d = float(row["target_d_A"])
        targets.append(
            {
                **row,
                "target_label": label,
                "target_d_A": f"{target_d:.9f}",
            }
        )

    return targets


def load_profiles(path: Path) -> dict[str, list[dict[str, str]]]:
    rows = read_csv(path)
    require_columns(rows, REQUIRED_PROFILE_COLUMNS, path)

    by_model: dict[str, list[dict[str, str]]] = defaultdict(list)

    for index, row in enumerate(rows, start=2):
        model_id = row["model_id"].strip()
        if not model_id:
            raise ValueError(f"{path} row {index} has empty model_id")

        d_A = float(row["d_A"])
        intensity = float(row["intensity"])

        by_model[model_id].append(
            {
                **row,
                "model_id": model_id,
                "d_A": f"{d_A:.9f}",
                "intensity": f"{intensity:.9f}",
            }
        )

    return dict(by_model)


def choose_peak(
    profile_rows: list[dict[str, str]],
    target_d_A: float,
    window_half_width: float,
) -> dict[str, str] | None:
    window_min = target_d_A - window_half_width
    window_max = target_d_A + window_half_width

    candidates = [
        row
        for row in profile_rows
        if window_min <= float(row["d_A"]) <= window_max
    ]

    if not candidates:
        return None

    return max(candidates, key=lambda row: float(row["intensity"]))


def extract_peaks(
    profiles_by_model: dict[str, list[dict[str, str]]],
    targets: list[dict[str, str]],
    window_half_width: float,
    include_missing: bool = False,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []

    for model_id in sorted(profiles_by_model):
        profile_rows = profiles_by_model[model_id]

        for target in targets:
            target_label = target["target_label"]
            target_d = float(target["target_d_A"])
            window_min = target_d - window_half_width
            window_max = target_d + window_half_width
            peak = choose_peak(profile_rows, target_d, window_half_width)

            if peak is None:
                if include_missing:
                    rows.append(
                        {
                            "model_id": model_id,
                            "target_label": target_label,
                            "observed_d_A": "",
                            "observed_intensity": "",
                            "target_d_A": f"{target_d:.9f}",
                            "delta_d_A": "",
                            "window_min_d_A": f"{window_min:.9f}",
                            "window_max_d_A": f"{window_max:.9f}",
                            "peak_status": "missing",
                        }
                    )
                continue

            observed_d = float(peak["d_A"])
            observed_i = float(peak["intensity"])
            rows.append(
                {
                    "model_id": model_id,
                    "target_label": target_label,
                    "observed_d_A": f"{observed_d:.9f}",
                    "observed_intensity": f"{observed_i:.9f}",
                    "target_d_A": f"{target_d:.9f}",
                    "delta_d_A": f"{observed_d - target_d:.9f}",
                    "window_min_d_A": f"{window_min:.9f}",
                    "window_max_d_A": f"{window_max:.9f}",
                    "peak_status": "found",
                }
            )

    return rows


def output_fields() -> list[str]:
    return [
        "model_id",
        "target_label",
        "observed_d_A",
        "observed_intensity",
        "target_d_A",
        "delta_d_A",
        "window_min_d_A",
        "window_max_d_A",
        "peak_status",
    ]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--profiles",
        type=Path,
        required=True,
        help="Profile CSV with model_id,d_A,intensity columns.",
    )
    parser.add_argument(
        "--targets",
        type=Path,
        default=Path("inputs/experimental_peak_windows/hxc590_twist_rise_scan_targets.csv"),
        help="Target peak CSV.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("outputs/twist_rise_scan/observed_peaks_long.csv"),
        help="Long-form observed peak output CSV.",
    )
    parser.add_argument(
        "--window-half-width",
        type=float,
        default=0.20,
        help="Half-width in Angstrom around each target d-spacing.",
    )
    parser.add_argument(
        "--include-missing",
        action="store_true",
        help="Write missing rows when no profile point is found inside a target window.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()

    if args.window_half_width <= 0:
        raise ValueError("--window-half-width must be positive")

    targets = load_targets(args.targets)
    profiles_by_model = load_profiles(args.profiles)
    rows = extract_peaks(
        profiles_by_model=profiles_by_model,
        targets=targets,
        window_half_width=args.window_half_width,
        include_missing=args.include_missing,
    )
    write_csv(args.output, rows, output_fields())
    print(f"Wrote {len(rows)} extracted peak rows to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
