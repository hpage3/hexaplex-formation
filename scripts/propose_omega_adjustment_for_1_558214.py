#!/usr/bin/env python3
"""Assess whether this repo can generate a modest omega adjustment for 1_558214."""

from __future__ import annotations

import csv
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT_CSV = ROOT / "outputs" / "metrics" / "asem_peptide_dihedral_candidate_structural_audit.csv"
FEASIBILITY_CSV = ROOT / "outputs" / "metrics" / "asem_1_558214_omega_adjustment_feasibility.csv"
REPORT = ROOT / "outputs" / "reports" / "asem_peptide_dihedral_candidate_screen_report.md"


def read_audit() -> list[dict[str, str]]:
    if not AUDIT_CSV.exists():
        raise FileNotFoundError(f"Run analyze_asem_peptide_dihedral_candidates.py first: {AUDIT_CSV}")
    with AUDIT_CSV.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def float_or_none(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def markdown_table(rows: list[dict[str, str]]) -> list[str]:
    lines = [
        "| Candidate | Nick-note omega proxy (deg) | Canonical amide omega (deg) | Min heavy contact (A) | <2.2 A contacts | Best N...O (A) | H-bond plausible | Status |",
        "| --- | ---: | ---: | ---: | ---: | ---: | --- | --- |",
    ]
    for row in rows:
        lines.append(
            "| {candidate_file} | {omega} | {canonical} | {contact} | {clash} | {hbond} | {flag} | {status} |".format(
                candidate_file=row["candidate_file"],
                omega=f"{float(row['nick_note_omega_proxy_N_CA_C_Np_deg']):.2f}",
                canonical=f"{float(row['canonical_amide_omega_CA_C_Np_CAp_deg']):.2f}",
                contact=f"{float(row['min_heavy_nonbonded_distance_A']):.2f}",
                clash=row["clash_count_lt_2p2"],
                hbond=f"{float(row['best_backbone_hbond_distance_A']):.2f}" if row["best_backbone_hbond_distance_A"] else "",
                flag=row["hbond_plausible_flag"],
                status=row["qualitative_status"],
            )
        )
    return lines


def write_report(rows: list[dict[str, str]], feasibility_row: dict[str, object]) -> None:
    baseline = next(row for row in rows if row["candidate_file"] == "1_558214.pdb")
    sorted_rows = sorted(rows, key=lambda row: row["candidate_file"])
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    content = [
        "# Asem Peptide-Dihedral Candidate Screen",
        "",
        "## Purpose",
        "",
        "This note audits Asem's peptide-dihedral candidate structures in light of Nick's latest guidance. The immediate question is not whether omega can be forced to 180 degrees, but whether the original H-bond-preserving candidate `1_558214.pdb` can be modestly moved by about +5 degrees toward roughly 172-173 degrees without losing the structural features that made it attractive.",
        "",
        "No corrected diffraction was run in this step. Nick's instruction was to stop after structural and feasibility analysis unless a clean adjusted PDB could be generated and pass geometry checks.",
        "",
        "## Nick and Literature Context",
        "",
        "- Nick's latest note says only `1_558214.pdb` appears to preserve good backbone-backbone H-bonding.",
        "- Nick does not require exact 180 degree omega; a modest shift toward 172-173 degrees would be easier to defend.",
        "- The Vitagliano omega paper is relevant because peptide omega angles in real protein structures deviate from 180 degrees, those deviations correlate with local backbone geometry, and beta-sheet regions can shift omega below 180 degrees.",
        "- Therefore, a roughly 167.5 degree omega-related geometry is strained but not automatically disqualifying.",
        "- Helical rise refinement near 3.38 A is a separate follow-up and was not mixed into this task.",
        "",
        "## Atom-Naming Limitation",
        "",
        "The candidate PDB atom names expose multiple backbone-like torsions. The canonical amide `CA-C-N'-CA'` torsion is near 180 degrees in all candidate files and therefore does not reproduce Nick's cited differences. The audit reports `N-CA-C-N'` as a Nick-note omega proxy because it is the local torsion that separates `1_558214.pdb` from the higher-omega candidates in this file set. This should be confirmed with Asem before using the value as a chemically exact omega.",
        "",
        "## Candidate Audit",
        "",
        *markdown_table(sorted_rows),
        "",
        "## Feasibility of a +5 Degree Adjustment",
        "",
        f"- Baseline file: `{feasibility_row['candidate_file']}`",
        f"- Current Nick-note omega proxy: {float(feasibility_row['current_omega_proxy_deg']):.2f} degrees",
        f"- Target range: {feasibility_row['target_omega_range_deg']}",
        f"- Pipeline decision: `{feasibility_row['pipeline_decision']}`",
        "",
        str(feasibility_row["reason"]),
        "",
        "Because this repo does not have an explicit local bonding/constraint model for coordinated phi/psi-like rotations, and because the task says not to perturb base/hexad stack geometry or helical rise/twist, the script did not generate adjusted PDB variants. A naive coordinate rotation would risk inventing chemistry and changing the wrong degrees of freedom.",
        "",
        "## Recommendation",
        "",
        "`1_558214.pdb` should remain the structural baseline if it is the only candidate that preserves the desired backbone-backbone H-bond. The most defensible next step is to ask Asem to attempt a constrained phi/psi-only adjustment or minimization in his modeling workflow, targeting a modest +5 degree omega shift while explicitly preserving the backbone-backbone H-bond and periodic geometry.",
        "",
        "## Outputs",
        "",
        f"- Structural audit CSV: `{AUDIT_CSV.relative_to(ROOT)}`",
        f"- Feasibility CSV: `{FEASIBILITY_CSV.relative_to(ROOT)}`",
        "",
        "## Checks and Scope",
        "",
        "- Candidate PDBs were copied into `inputs/asem_peptide_dihedral_candidates/`.",
        "- No notebooks were executed.",
        "- Preserved source archives were not modified.",
        "- pNAB was not used or edited.",
        "- Corrected diffraction was not rerun.",
    ]
    REPORT.write_text("\n".join(content) + "\n", encoding="utf-8")


def main() -> int:
    rows = read_audit()
    baseline = next((row for row in rows if row["candidate_file"] == "1_558214.pdb"), None)
    if baseline is None:
        raise RuntimeError("1_558214.pdb was not found in the structural audit")

    omega = float_or_none(baseline["nick_note_omega_proxy_N_CA_C_Np_deg"])
    if omega is None:
        raise RuntimeError("Could not read Nick-note omega proxy for 1_558214.pdb")

    feasibility_row = {
        "candidate_file": "1_558214.pdb",
        "current_omega_proxy_deg": omega,
        "target_omega_range_deg": "172-173",
        "target_delta_deg": 172.5 - omega,
        "variants_generated": 0,
        "rmsd_retained_atoms_A": "",
        "best_backbone_hbond_distance_A": baseline["best_backbone_hbond_distance_A"],
        "min_heavy_nonbonded_distance_A": baseline["min_heavy_nonbonded_distance_A"],
        "pipeline_decision": "no_adjusted_pdb_generated",
        "reason": (
            "The candidate atom naming supports structural auditing, but the coordinated phi/psi-like torsion controls "
            "needed to move only the peptide omega-related geometry while preserving base geometry, rise, twist, and "
            "periodicity are not defined in this repo. A constrained modeling workflow from Asem is required."
        ),
    }
    write_csv(FEASIBILITY_CSV, [feasibility_row])
    write_report(rows, feasibility_row)
    print(f"Wrote {FEASIBILITY_CSV}")
    print(f"Wrote {REPORT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
