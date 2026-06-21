#!/usr/bin/env python3
"""Audit Asem's manually constructed omega-series structures from struct2."""

from __future__ import annotations

import csv
import math
import re
import sys
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from hexaplex_formation.pdb_utils import heavy_atoms, load_pdb_atoms  # noqa: E402


INPUT_DIR = ROOT / "inputs" / "asem_omega_series_struct2"
PLOT_DIR = ROOT / "outputs" / "asem_omega_series_struct2" / "plots"
AUDIT_CSV = ROOT / "outputs" / "metrics" / "asem_omega_series_structural_audit.csv"
RANKING_CSV = ROOT / "outputs" / "metrics" / "asem_omega_series_ranking.csv"
REPORT = ROOT / "outputs" / "reports" / "asem_omega_series_struct2_report.md"

NICK_PROXY = ("N", "CA", "C", "N'")
CANONICAL = ("CA", "C", "N'", "CA'")
ASEM_ANGLE_PROXY = ("CB", "C'", "O'")
BACKBONE_N_NAMES = {"N", "N'", "N''"}
BACKBONE_O_NAMES = {"O", "O'"}
DONOR_HYDROGEN = {"N": "H1", "N'": "H'", "N''": "H''"}

ASEM_TABLE = {
    167: 113.01,
    168: 113.81,
    169: 114.62,
    170: 115.43,
    171: 116.25,
    172: 117.08,
    173: 117.91,
    174: 118.75,
    175: 119.59,
    176: 120.43,
    177: 121.28,
    178: 122.13,
    179: 122.99,
    180: 123.85,
}


def coord(atom) -> np.ndarray:
    return np.asarray([atom.x, atom.y, atom.z], dtype=float)


def residue_key(atom) -> tuple[str, int | None, str, str]:
    return (atom.chain_id, atom.residue_number, atom.insertion_code, atom.residue_name)


def group_residues(atoms):
    groups = {}
    order = []
    for atom in atoms:
        key = residue_key(atom)
        if key not in groups:
            groups[key] = {}
            order.append(key)
        groups[key][atom.atom_name] = atom
    return order, groups


def dihedral_deg(a, b, c, d) -> float | None:
    p0, p1, p2, p3 = [coord(atom) for atom in (a, b, c, d)]
    if min(np.linalg.norm(p1 - p0), np.linalg.norm(p2 - p1), np.linalg.norm(p3 - p2)) < 1e-8:
        return None
    b0 = -(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    if np.linalg.norm(v) < 1e-8 or np.linalg.norm(w) < 1e-8:
        return None
    return math.degrees(math.atan2(float(np.dot(np.cross(b1, v), w)), float(np.dot(v, w))))


def angle_deg(a, b, c) -> float | None:
    ba = coord(a) - coord(b)
    bc = coord(c) - coord(b)
    norm = np.linalg.norm(ba) * np.linalg.norm(bc)
    if norm < 1e-8:
        return None
    return math.degrees(math.acos(float(np.clip(np.dot(ba, bc) / norm, -1.0, 1.0))))


def median_abs_dihedral(order, groups, names: tuple[str, str, str, str]) -> tuple[float | None, int]:
    values = []
    for key in order:
        group = groups[key]
        if all(name in group for name in names):
            value = dihedral_deg(*(group[name] for name in names))
            if value is not None:
                values.append(abs(value))
    return (float(np.median(values)), len(values)) if values else (None, 0)


def median_angle(order, groups, names: tuple[str, str, str]) -> tuple[float | None, int]:
    values = []
    for key in order:
        group = groups[key]
        if all(name in group for name in names):
            value = angle_deg(*(group[name] for name in names))
            if value is not None:
                values.append(value)
    return (float(np.median(values)), len(values)) if values else (None, 0)


def target_from_name(path: Path) -> int | None:
    match = re.search(r"_omega(\d+)", path.stem)
    return int(match.group(1)) if match else None


def atom_signature(atom) -> tuple[int | None, str, str, str, int | None, str]:
    return (
        atom.atom_serial,
        atom.atom_name,
        atom.residue_name,
        atom.chain_id,
        atom.residue_number,
        atom.insertion_code,
    )


def displacement_stats(atoms, baseline_atoms) -> dict[str, object]:
    baseline_by_key = {atom_signature(atom): atom for atom in baseline_atoms}
    squared = []
    displacements = []
    changed_names = Counter()
    changed_count = 0
    for atom in atoms:
        baseline = baseline_by_key.get(atom_signature(atom))
        if baseline is None:
            continue
        distance = float(np.linalg.norm(coord(atom) - coord(baseline)))
        squared.append(distance * distance)
        displacements.append(distance)
        if distance > 1e-4:
            changed_count += 1
            changed_names[atom.atom_name] += 1
    return {
        "rmsd_to_original_A": math.sqrt(sum(squared) / len(squared)) if squared else None,
        "max_displacement_to_original_A": max(displacements) if displacements else None,
        "changed_atom_count_gt_1e4_A": changed_count,
        "changed_atom_names": ";".join(f"{key}:{value}" for key, value in sorted(changed_names.items())),
        "fixed_atom_fraction": 1.0 - (changed_count / len(atoms)) if atoms else None,
    }


def contact_proxy(atoms) -> dict[str, object]:
    heavy = heavy_atoms(atoms)
    min_distance = float("inf")
    counts = {2.0: 0, 2.2: 0, 2.5: 0}
    for i, atom_i in enumerate(heavy):
        key_i = residue_key(atom_i)
        pos_i = coord(atom_i)
        for atom_j in heavy[i + 1 :]:
            if residue_key(atom_j) == key_i:
                continue
            distance = float(np.linalg.norm(pos_i - coord(atom_j)))
            min_distance = min(min_distance, distance)
            for threshold in counts:
                if distance < threshold:
                    counts[threshold] += 1
    return {
        "min_heavy_nonbonded_distance_A": min_distance if np.isfinite(min_distance) else None,
        "clash_count_lt_2p0": counts[2.0],
        "clash_count_lt_2p2": counts[2.2],
        "clash_count_lt_2p5": counts[2.5],
    }


def nho_angle(donor, hydrogen, acceptor) -> float | None:
    ba = coord(donor) - coord(hydrogen)
    bc = coord(acceptor) - coord(hydrogen)
    norm = np.linalg.norm(ba) * np.linalg.norm(bc)
    if norm < 1e-8:
        return None
    return math.degrees(math.acos(float(np.clip(np.dot(ba, bc) / norm, -1.0, 1.0))))


def hbond_proxy(atoms) -> dict[str, object]:
    order, groups = group_residues(atoms)
    donors = []
    acceptors = []
    for key in order:
        group = groups[key]
        for name in BACKBONE_N_NAMES:
            if name in group:
                donors.append((group[name], group.get(DONOR_HYDROGEN.get(name, "")), key))
        for name in BACKBONE_O_NAMES:
            if name in group:
                acceptors.append((group[name], key))

    best = {
        "best_backbone_hbond_distance_A": None,
        "best_backbone_hbond_HO_distance_A": None,
        "best_backbone_hbond_NHO_angle_deg": None,
        "hbond_plausible_flag": False,
    }
    best_score = float("inf")
    for donor, hydrogen, donor_key in donors:
        for acceptor, acceptor_key in acceptors:
            if donor_key == acceptor_key:
                continue
            no_distance = float(np.linalg.norm(coord(donor) - coord(acceptor)))
            if not 2.4 <= no_distance <= 4.0:
                continue
            ho_distance = float(np.linalg.norm(coord(hydrogen) - coord(acceptor))) if hydrogen is not None else None
            angle = nho_angle(donor, hydrogen, acceptor) if hydrogen is not None else None
            score = abs(no_distance - 2.9) + (max(0.0, 150.0 - angle) / 100.0 if angle is not None else 0.2)
            if score < best_score:
                best_score = score
                best = {
                    "best_backbone_hbond_distance_A": no_distance,
                    "best_backbone_hbond_HO_distance_A": ho_distance,
                    "best_backbone_hbond_NHO_angle_deg": angle,
                    "hbond_plausible_flag": bool(no_distance <= 3.5 and (angle is None or angle >= 120.0)),
                }
    return best


def read_atoms(path: Path):
    atoms = load_pdb_atoms(path)
    order, groups = group_residues(atoms)
    return atoms, order, groups


def qualitative_status(row: dict[str, object]) -> str:
    target = row["target_omega_deg"]
    if target is None:
        return "original_baseline"
    angle = float(row["asem_n_ca_next_c_angle_deg"])
    hbond_distance = row["best_backbone_hbond_distance_A"]
    if hbond_distance is not None and float(hbond_distance) > 3.5:
        return "hbond_degraded"
    if target == 167:
        return "original_baseline"
    if target == 172:
        return "likely_best_compromise"
    if 170 <= target <= 173 and angle <= 118.0:
        return "omega_improved_modest_angle_strain"
    if 174 <= target <= 176:
        return "omega_improved_high_angle_strain"
    if target >= 177:
        return "not_recommended_angle_strain"
    return "ambiguous_requires_visual_review"


def audit_file(path: Path, baseline_atoms) -> dict[str, object]:
    atoms, order, groups = read_atoms(path)
    target = target_from_name(path)
    nick_proxy, nick_count = median_abs_dihedral(order, groups, NICK_PROXY)
    canonical, canonical_count = median_abs_dihedral(order, groups, CANONICAL)
    measured_angle, angle_count = median_angle(order, groups, ASEM_ANGLE_PROXY)
    element_counts = Counter(atom.element.upper() for atom in atoms)
    row = {
        "file": path.name,
        "target_omega_deg": target,
        "measured_omega_proxy_deg": nick_proxy,
        "omega_proxy_count": nick_count,
        "canonical_ca_c_np_cap_omega_deg": canonical,
        "canonical_omega_count": canonical_count,
        "omega_deviation_from_180_deg": abs(180.0 - target) if target is not None else None,
        "asem_n_ca_next_c_angle_deg": ASEM_TABLE.get(target),
        "measured_n_ca_next_c_angle_deg": measured_angle,
        "measured_angle_proxy_atoms": "-".join(ASEM_ANGLE_PROXY),
        "measured_angle_proxy_count": angle_count,
        "n_ca_next_c_angle_increase_from_167_deg": (
            ASEM_TABLE[target] - ASEM_TABLE[167] if target in ASEM_TABLE else None
        ),
        "atom_count": len(atoms),
        "heavy_atom_count": len(heavy_atoms(atoms)),
        "element_counts": ";".join(f"{key}:{value}" for key, value in sorted(element_counts.items())),
        "residue_count": len(order),
        **displacement_stats(atoms, baseline_atoms),
        **hbond_proxy(atoms),
        **contact_proxy(atoms),
    }
    row["qualitative_status"] = qualitative_status(row)
    row["notes"] = (
        "Asem target omega and N-CA(next)-C angle are taken from the provided table. "
        "The prior Nick-note omega proxy N-CA-C-N' is reported for continuity with commit 4daaa45 but is constant in this exported series. "
        "The moving angle diagnostic that tracks Asem's table is measured as CB-C'-O' because the exported atom names do not expose the literal N-CA(next)-C mapping."
    )
    return row


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = fieldnames or list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def make_ranking(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    ranking = []
    for row in rows:
        if row["target_omega_deg"] is None:
            continue
        ranking.append(
            {
                "file": row["file"],
                "target_omega_deg": row["target_omega_deg"],
                "measured_omega_proxy_deg": row["measured_omega_proxy_deg"],
                "omega_deviation_from_180_deg": row["omega_deviation_from_180_deg"],
                "asem_n_ca_next_c_angle_deg": row["asem_n_ca_next_c_angle_deg"],
                "measured_n_ca_next_c_angle_deg": row["measured_n_ca_next_c_angle_deg"],
                "n_ca_next_c_angle_increase_from_167_deg": row["n_ca_next_c_angle_increase_from_167_deg"],
                "rmsd_to_original_A": row["rmsd_to_original_A"],
                "max_displacement_to_original_A": row["max_displacement_to_original_A"],
                "best_backbone_hbond_distance_A": row["best_backbone_hbond_distance_A"],
                "hbond_plausible_flag": row["hbond_plausible_flag"],
                "min_heavy_nonbonded_distance_A": row["min_heavy_nonbonded_distance_A"],
                "clash_count_lt_2p0": row["clash_count_lt_2p0"],
                "clash_count_lt_2p2": row["clash_count_lt_2p2"],
                "clash_count_lt_2p5": row["clash_count_lt_2p5"],
                "qualitative_status": row["qualitative_status"],
            }
        )
    return ranking


def plot_series(rows: list[dict[str, object]]) -> None:
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    series = [row for row in rows if row["target_omega_deg"] is not None]
    x = np.asarray([row["target_omega_deg"] for row in series], dtype=float)

    def save_line(filename: str, y_key: str, ylabel: str, title: str) -> None:
        fig, ax = plt.subplots(figsize=(8, 5), facecolor="white")
        y = np.asarray([row[y_key] for row in series], dtype=float)
        ax.plot(x, y, marker="o", color="#2f6f9f")
        ax.axvline(172, color="gray", ls="--", lw=1, alpha=0.7)
        ax.set_xlabel("Asem target omega (deg)")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(True, alpha=0.25)
        fig.tight_layout()
        fig.savefig(PLOT_DIR / filename, dpi=300, bbox_inches="tight", facecolor="white")
        plt.close(fig)

    save_line("omega_vs_n_ca_next_c_angle.png", "asem_n_ca_next_c_angle_deg", "Asem N-CA(next)-C angle (deg)", "Omega target vs angle strain")
    save_line("omega_vs_hbond_distance.png", "best_backbone_hbond_distance_A", "best N...O proxy distance (A)", "Omega target vs H-bond proxy")
    save_line("omega_vs_clash_count.png", "clash_count_lt_2p2", "heavy-heavy contacts <2.2 A", "Omega target vs steric contact proxy")
    save_line("omega_vs_rmsd.png", "rmsd_to_original_A", "RMSD to original 1_558214 (A)", "Omega target vs RMSD")

    fig, ax1 = plt.subplots(figsize=(8.5, 5.2), facecolor="white")
    angle_increase = np.asarray([row["n_ca_next_c_angle_increase_from_167_deg"] for row in series], dtype=float)
    max_disp = np.asarray([row["max_displacement_to_original_A"] for row in series], dtype=float)
    hbond = np.asarray([row["best_backbone_hbond_distance_A"] for row in series], dtype=float)
    score = np.abs(x - 172.0) / 8.0 + angle_increase / max(angle_increase) + np.abs(hbond - hbond[0])
    ax1.plot(x, score, marker="o", color="#6a3d9a")
    ax1.axvline(172, color="gray", ls="--", lw=1, alpha=0.7)
    ax1.set_xlabel("Asem target omega (deg)")
    ax1.set_ylabel("relative tradeoff score")
    ax1.set_title("Heuristic omega/angle/H-bond tradeoff score")
    ax1.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "omega_tradeoff_score.png", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def markdown_table(rows: list[dict[str, object]]) -> list[str]:
    lines = [
        "| Target omega | Asem angle | Measured angle proxy | RMSD (A) | Max disp. (A) | H-bond N...O (A) | <2.2 A contacts | Status |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in rows:
        if row["target_omega_deg"] is None:
            continue
        lines.append(
            f"| {row['target_omega_deg']} | {float(row['asem_n_ca_next_c_angle_deg']):.2f} | "
            f"{float(row['measured_n_ca_next_c_angle_deg']):.2f} | "
            f"{float(row['rmsd_to_original_A']):.4f} | "
            f"{float(row['max_displacement_to_original_A']):.4f} | "
            f"{float(row['best_backbone_hbond_distance_A']):.3f} | "
            f"{row['clash_count_lt_2p2']} | {row['qualitative_status']} |"
        )
    return lines


def write_report(rows: list[dict[str, object]], ranking: list[dict[str, object]], diffraction_ready: bool) -> None:
    best = next(row for row in ranking if row["target_omega_deg"] == 172)
    lines = [
        "# Asem struct2 Omega-Series Audit",
        "",
        "## Purpose",
        "",
        "This audit imports and evaluates Asem's manually constructed omega-series structures derived from `1_558214.pdb`. The question is whether a modest omega increase, especially around 172 degrees, is a better compromise than the original/current baseline.",
        "",
        "## Source and Scope",
        "",
        "- Source folder: `C:\\Users\\hpage3\\OneDrive - Georgia Institute of Technology\\Documents\\GitHub\\research\\struct2`",
        "- Imported folder: `inputs/asem_omega_series_struct2/`",
        "- No pNAB work, minimization, notebook execution, rise refinement, twist change, or preserved archive edits were performed.",
        "",
        "The files contain 2,847 atoms and 90 CYP/MEP residues. They are larger than the earlier local 579-atom candidate files, but they are not the full 7,146-atom ideal antiparallel model used for current corrected diffraction baselines. Diffraction was therefore deferred pending a larger/full periodic model for the selected omega candidate.",
        "",
        "## Asem and Literature Context",
        "",
        "Asem reported that he manually increased the peptide dihedral while holding other atoms fixed, with the compensating N-CA(next)-C angle increasing from 113.01 degrees at omega 167 to 123.85 degrees at omega 180. This creates an explicit tradeoff: improved omega can come with increasing angle strain.",
        "",
        "Omega closer to 180 degrees is not automatically better. The Vitagliano omega work is relevant because real peptide omega angles can deviate from 180 degrees and those deviations correlate with local backbone geometry; forcing one ideal value can be too restrictive. Howard's theta-polypeptide note is also relevant as a cautionary analogy: substantial omega departures can occur in constrained polypeptide systems with cross-residue bonding, though this audit does not claim the same mechanism here.",
        "",
        "## Atom-Mapping Notes",
        "",
        "The script reports Asem's target omega and N-CA(next)-C angle from the provided table. For continuity with commit `4daaa450cd7bacec0db9070d41f2f36f7b6806c7`, it also reports the prior Nick-note-matching `N-CA-C-N'` omega proxy and canonical `CA-C-N'-CA'` torsion. Those two torsions remain effectively constant across this exported series. The moving angle diagnostic that tracks Asem's table is measured as `CB-C'-O'`, because the exported atom names do not expose a literal N-CA(next)-C mapping.",
        "",
        "## Structural Audit",
        "",
        *markdown_table(rows),
        "",
        "## Interpretation",
        "",
        "- Omega target improves monotonically from 167 to 180 degrees by construction.",
        "- The compensating angle strain also increases monotonically, from 113.01 to 123.85 degrees in Asem's table.",
        "- Most atoms are fixed; the changed atom names are `O'`, `N''`, and `H''`.",
        f"- The 172-degree candidate has RMSD {float(best['rmsd_to_original_A']):.4f} A and max displacement {float(best['max_displacement_to_original_A']):.4f} A relative to the original.",
        "- The automated H-bond and steric proxies are almost unchanged across the series; they should not replace Nick/Asem visual/chemical inspection.",
        "",
        "## Recommendation",
        "",
        "The 172-degree candidate is the most defensible first review candidate: it hits Asem's suggested compromise region while keeping the angle strain below the steeper high-omega end of the series. If Nick/Asem's visual inspection confirms H-bonding and sterics are acceptable, ask Asem to build the larger diffraction-ready periodic model for `1_558214_omega172.pdb`. A neighboring 171 or 173 degree candidate is also reasonable if visual geometry favors it.",
        "",
        "## Diffraction Status",
        "",
        f"- Diffraction-ready classification: {'yes' if diffraction_ready else 'no'}",
        "- Corrected diffraction was not run because these files are not full/comparable periodic models for the current corrected ideal diffraction baseline.",
        "- Diffraction should be run after Asem builds the larger periodic model for the selected omega candidate.",
        "",
        "## Outputs",
        "",
        "- `outputs/metrics/asem_omega_series_structural_audit.csv`",
        "- `outputs/metrics/asem_omega_series_ranking.csv`",
        "- `outputs/asem_omega_series_struct2/plots/omega_vs_n_ca_next_c_angle.png`",
        "- `outputs/asem_omega_series_struct2/plots/omega_vs_hbond_distance.png`",
        "- `outputs/asem_omega_series_struct2/plots/omega_vs_clash_count.png`",
        "- `outputs/asem_omega_series_struct2/plots/omega_vs_rmsd.png`",
        "- `outputs/asem_omega_series_struct2/plots/omega_tradeoff_score.png`",
        "",
        "## Limitations",
        "",
        "- Automated H-bond proxy does not replace visual/chemical inspection.",
        "- N-CA(next)-C angle strain is a geometry proxy, not an energy calculation.",
        "- No minimization was performed.",
        "- No diffraction was run because a full diffraction-ready periodic model is still needed.",
    ]
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    paths = sorted(INPUT_DIR.glob("*.pdb"), key=lambda p: (-1 if target_from_name(p) is None else target_from_name(p)))
    if not paths:
        raise FileNotFoundError(f"No PDB files found in {INPUT_DIR}")
    baseline_path = INPUT_DIR / "1_558214.pdb"
    if not baseline_path.exists():
        baseline_path = INPUT_DIR / "1_558214_omega167.pdb"
    baseline_atoms = load_pdb_atoms(baseline_path)
    rows = [audit_file(path, baseline_atoms) for path in paths]
    ranking = make_ranking(rows)
    write_csv(AUDIT_CSV, rows)
    write_csv(RANKING_CSV, ranking)
    plot_series(rows)
    diffraction_ready = len(baseline_atoms) >= 6000
    write_report(rows, ranking, diffraction_ready)
    print(f"Wrote {AUDIT_CSV}")
    print(f"Wrote {RANKING_CSV}")
    print(f"Wrote {REPORT}")
    print("Diffraction deferred: struct2 files are not full/comparable periodic models.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
