# HXC590 S1 Twist/Rise Grid Workflow Plan

## Purpose

The purpose of this plan is to define a future twist/rise grid workflow that moves from "best available screened candidate" toward a stronger test of what the HXC590 S1 powder peak list supports.

The workflow should distinguish among three levels of interpretation:

- support for the broad stacked-hexad / hexaflex-like topology;
- support for a close 30-degree-class antiparallel geometry;
- support for the nominal 30-degree idealized model specifically.

This plan does not generate coordinate variants, add model-generation code, run pNAB, or change existing scored conclusions.

## Current Evidence Summary

Within the current screened candidate set, `central8_units_30deg` is the best available screened candidate. It matches the diagnostic HXC590 S1 powder windows around approximately 3.35 A, 4.33 A, and the 3.7-3.9 A region.

Available wrong-geometry and partial-structure controls fail the same survival criteria used for the stacked-hexad / hexaflex-like candidates.

The current result remains a falsification-style screen against available alternatives, not a definitive phase assignment. Adjacent full-length twist/rise variants remain untested.

## Model-Generation Boundary

Model generation belongs outside this branch.

pNAB remains a separate external utility repository. This HXC590 diffraction/falsification branch should not own the ideal-model generator, add pNAB as a branch-local requirement, vendor pNAB, add pNAB execution code, or generate pNAB-derived twist/rise variants.

This branch should only accept externally generated candidate coordinate files after provenance and baseline-reproduction checks have been completed in a separate workflow. The role of this branch is to score candidate profiles, summarize survival/failure behavior, and report conservative interpretations.

Generic public pNAB examples such as `Hexad.yaml` and `Hexad_Antiparallel.yaml` should not be treated as HXC590 production inputs unless a separate provenance record connects them to the current HXC590/hexaplex baseline.

## Required Baseline Reproduction Check

The first gate is baseline reproduction. An external model-generation workflow should:

1. Regenerate the nominal 30-degree idealized model outside this branch.
2. Convert it to the same coordinate and scoring form as the existing baseline.
3. Generate its radial profile using the established diffraction workflow.
4. Compare the regenerated 30-degree radial profile against the existing full-length 30-degree baseline radial profile.
5. Proceed to nearby twist/rise variants only if the regenerated 30-degree model reproduces the existing baseline profile closely enough for the intended screening use.

The pass/fail criteria should be defined before scoring variants. Conceptually, the regenerated baseline should preserve the diagnostic peak-window behavior, major radial-profile features, atom/element accounting, chain/strand conventions, and coordinate-frame assumptions used by the current scoring workflow.

If the regenerated 30-degree model does not reproduce the current baseline profile, stop and resolve model provenance, coordinate conversion, or diffraction-profile generation before any grid variants are scored.

## Twist/Rise Grid Design

The proposed first-pass grid is:

- Twist: 24, 26, 28, 30, 32, 34, and 36 degrees.
- Rise: 3.20, 3.25, 3.30, 3.35, 3.40, 3.45, 3.50, 3.55, and 3.60 A.

The exact grid can be refined, but this range is sufficient to test whether the current powder data prefer the nominal 30-degree model specifically, a nearby 30-degree-class neighborhood, or the broader stacked-hexad family.

The nominal 30-degree baseline should remain the reference condition until the grid shows whether a nearby twist/rise combination scores better.

## Candidate Naming Convention

Candidate names should encode the model family, geometry class, twist, and rise:

- `hxc590_ideal_antiparallel_twist_24_rise_3p40`
- `hxc590_ideal_antiparallel_twist_30_rise_3p40`
- `hxc590_ideal_antiparallel_twist_36_rise_3p40`

For rise values, use `p` instead of `.` so filenames remain portable, for example `3p35` for 3.35 A.

## Handoff Requirements

Each externally generated candidate should be handed into the diffraction/falsification repo with a manifest row containing:

- coordinate file path;
- generator name, version, or commit;
- source input file;
- twist value;
- rise value;
- atom count;
- element counts;
- coordinate-frame convention;
- whether hydrogens are included or excluded;
- conversion and deduplication steps;
- provenance note;
- checksum, if practical.

The handoff should also state whether the candidate belongs to a baseline-reproduction run, the first twist/rise grid, or a later refined grid.

## Scoring Plan Once Candidates Are Available

Once externally generated candidates pass the baseline-reproduction gate and are handed into this repo, the diffraction/falsification workflow can:

- generate radial profiles using the established diffraction workflow;
- compare each profile against the HXC590 S1 powder peak windows;
- score diagnostic-window matches around approximately 3.35 A, 4.33 A, and 3.7-3.9 A;
- penalize missing diagnostic peaks;
- track spurious predicted peaks where the experimental profile shows local minima, if the current experimental profile data support that comparison;
- run narrow, current, and broad tolerance survival audits;
- compare the grid against existing wrong-geometry and partial-structure controls.

The scoring should preserve the current conservative peak-window framing unless better calibrated experimental data justify a more detailed treatment.

## Possible Outcomes

### A. 30-Degree Neighborhood Supported

The nominal 30-degree model is best, nearby 28-degree and 32-degree models are close, and distant twists score worse.

Claim: the powder peak list supports a 30-degree-class antiparallel stacked-hexad topology.

### B. 30 Degrees Preferred Within The Tested Grid

The nominal 30-degree model is clearly best across the tested twist/rise grid.

Claim: the powder peak list preferentially supports the nominal 30-degree idealized model among tested variants.

### C. Broad Family Only

Many twist/rise variants survive under the same criteria.

Claim: the powder peak list supports the stacked-hexad topology but does not resolve twist/rise within the tested grid.

### D. Different Nearby Variant Preferred

A nearby variant, such as 28 degrees or 32 degrees, scores better than the nominal 30-degree baseline.

Claim: the powder peak list supports a close member of the model family, not necessarily the nominal 30-degree baseline.

### E. Grid Fails

Nearby variants and the baseline fail under consistent scoring.

Claim: revisit model generation, coordinate conversion, diffraction-profile generation, or experimental peak interpretation.

## Risks And Safeguards

- Avoid overfitting low-resolution powder data.
- Do not use generic pNAB examples as HXC590 inputs without provenance.
- Avoid treating idealized models as experimentally solved structures.
- Keep model generation and diffraction scoring separate.
- Require baseline reproduction before scoring variants.
- Preserve all limitations in report language.
- Treat 30 degrees as the nominal tested baseline unless the grid supports a more specific interpretation.

## Reviewer-Aligned Claim Language

"The HXC590 S1 powder X-ray diffraction peak list is compatible with the structural topology of the idealized, antiparallel stacked-hexad / hexaflex-like model family. Within the current screened candidate set, central8_units_30deg yields the highest match rate against key diagnostic powder windows, whereas available partial-structure and wrong-geometry controls fail identical survival criteria. Because adjacent full-length conformational variants remain untested, this result represents a successful falsification-style exclusion of simple alternative architectures rather than a unique, refined phase assignment."

## Acceptance Criteria For This Planning Task

- Adds this workflow plan report.
- Adds a compact CSV planning table.
- Does not generate coordinate variants.
- Does not add pNAB execution code.
- Does not add pNAB as a branch-local requirement or environment step.
- Does not change existing scored conclusions.
- Keeps language conservative.
- Tests pass.

