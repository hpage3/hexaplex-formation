# Anti-Parallel Glu Side-Chain Pipeline Update

Date: 2026-06-25

## Revised Screening Plan

Nick's revised guidance is the working plan for the next anti-parallel screen:

- Fix the helical rise at 3.38 A.
- Start with integer twists 29, 30, and 31 degrees.
- Add 28 and 32 degrees only if the first three profiles do not resolve the
  comparison, or if the boundary profiles remain about as competitive as 30
  degrees.
- Consider half-degree values near 30 degrees, especially 30.5 degrees, only
  after reviewing the integer-twist results.
- Use chemically complete Glu side chains in the primary models.
- Reserve no-side-chain structures for later ablation tests.
- Treat powder/radial-profile agreement as the primary comparison.
- Treat corrected 2D diffraction images as diagnostics for locating the
  sources of powder features, not as the primary ranking target.
- Defer parallel-strand models until the anti-parallel side-chain comparison
  is complete.

The existing Asem side-chain Stage 1 and Stage 2 screens use 3.40 A structures.
They validate the import, conversion, corrected-diffraction, radial-profile,
and scoring pipeline, but they are not final results for the 3.38 A question.

## Local Feasibility Audit

### What is already present

`inputs/pnab_asem_30deg/` contains a reproducible copy of Asem's anti-parallel
pNAB workflow:

- `options.yaml` has alternating strand orientations, six enabled strands,
  CYP/MEP side-chain-aligned component files, and a fixed 3.38 A rise.
- `run.py` runs pNAB and writes `fixed.pdb`.
- `initial.in`, CYP/MEP Amber libraries, and CYP/MEP parameter files are
  present for the subsequent AmberTools step.
- The public `pnab` conda environment is installed and imports pNAB.

The pNAB template and the imported 3.40 A side-chain package do not have the
same length. The template uses 15 base/GLU repeat units per strand (180 total
residues after `run.py`), while Asem's current side-chain candidate package
uses 16 units per strand (192 total residues). A future production set should
hold length fixed across twists. To compare directly with the Stage 1/2
package, the preferred request is the same 16-unit setup.

### What is not already present

There is no complete 3.38 A Glu-side-chain set for twists 29, 30, and 31.

The prior 3.38 A pNAB twist-series manifest records a 30-degree and 31-degree
coordinate, but not a successful 29-degree coordinate, and the generated
structures are not currently retained under the manifest's structure paths.
More importantly, the same pNAB intermediate format seen in the retained
pilot has 2,835 atoms and GLU atom names only through `CB`. It lacks the
`CG`, `CD`, `OE1`, and `OE2` atoms required for complete glutamate side
chains.

By comparison, an imported completed Asem side-chain candidate has 3,810
atoms, 192 residues, and GLU atoms including `CG`, `CD`, `OE1`, and `OE2`.
The imported 29/30/31 reference PDBs are at the earlier approximately 3.40 A
setup, not the requested 3.38 A setup.

### Current blocker

The local pNAB environment is available at:

`C:\Users\Public\hexaplex-tools\miniforge3\envs\pnab\python.exe`

AmberTools `tleap` was not found on `PATH`, in the public pNAB environment, or
in the checked public Miniforge installation. Therefore this machine can
generate the pNAB `fixed.pdb` intermediate, but it cannot currently perform
the documented `tleap -f initial.in` completion step without a separately
reviewed AmberTools installation or environment.

No dependency installation or model generation was performed during this
audit.

## Reproducible Partial Local Workflow

The following command would request only the three 3.38 A anti-parallel pNAB
intermediates from the existing 15-unit template. It is recorded for
reproducibility, but was not run as part of this planning task:

```powershell
.\.venv\Scripts\python.exe scripts\generate_pnab_twist_series.py `
  --workflow-dir inputs\pnab_asem_30deg `
  --output-root outputs\pnab_antiparallel_glu_rise3p38_twist29_30_31 `
  --rise 3.38 `
  --twists 29 30 31 `
  --timeout-seconds 1800
```

If a reviewed AmberTools environment becomes available, each successful work
directory would then be completed independently:

```powershell
cd outputs\pnab_antiparallel_glu_rise3p38_twist29_30_31\work\twist_29
tleap -f initial.in
```

Repeat for `twist_30` and `twist_31`. Before diffraction, each `initial.pdb`
must be audited for the expected chain/residue count and for complete GLU
atoms (`CG`, `CD`, `OE1`, and `OE2`). This partial command uses the 15-unit
template and therefore is not a drop-in replacement for Asem's 16-unit
Stage 1/2 candidate set.

## Recommended Request to Asem

Because the locally installed tools cannot complete the Glu side chains, the
recommended next action is to request the finished structures from Asem:

> Please generate a minimal anti-parallel Hexaplex set using the same complete
> Glu-side-chain workflow and the same 16-unit-per-strand length as the
> side-chain package you sent. Fix the helical rise at 3.38 A and generate
> twists 29, 30, and 31 degrees. Please run the AmberTools/tleap completion so
> every Glu contains the full side chain through CG, CD, OE1, and OE2. For each
> twist, please send the final `initial.pdb`, the pNAB `fixed.pdb`,
> `options.yaml`, `results.csv`, `prefix.yaml`, the tleap input/log, and any
> Amber parameter/library files not identical to the previous package. Please
> identify the candidate selected at each twist and note any pNAB or tleap
> warnings. Multiple 30-degree candidates are welcome, but one audited
> candidate per twist is enough for the first screen.

This request changes rise while preserving side-chain chemistry, strand
arrangement, and model length, which makes the comparison to the existing
3.40 A Stage 1/2 screen interpretable.

## Proposed File Layout

After receipt and audit, source/provenance files should live under:

`inputs/asem_sidechains_rise3p38_antiparallel_20260625/`

Suggested subdirectories are `raw/29/`, `raw/30/`, and `raw/31/`, with one
candidate directory per supplied structure. Generated corrected diffraction,
radial profiles, scores, and diagnostic images should live under:

`outputs/asem_sidechains_rise3p38_antiparallel_20260625/`

Heavy generated XYZ, NPY, radial-profile directories, and individual
diffraction images should remain untracked until explicitly reviewed.

## Commit Guidance

Appropriate to commit now:

- this planning note;
- the planned manifest template, which contains no claimed structure paths.

Do not commit now:

- generated pNAB work directories;
- converted XYZ or NPY files;
- radial-profile directories;
- individual diffraction images;
- unreviewed structures;
- existing unrelated smoke and screening outputs.

The first 29/30/31 comparison remains a screening experiment, not a final
ranking.

## Candidate-to-Profile Correlation Audit

The imported package does not contain Asem's per-candidate fast 1D powder
profile files. A recursive file inventory found:

- 185 expanded `initial.pdb` structures;
- one `options.yaml`;
- one `per_angle_overlays.png`;
- one `tleap_structures.zip`;
- no per-candidate CSV, TXT, NPY, log, JSON, or other profile/metadata files.

`tleap_structures.zip` contains exactly 185 entries, all named
`angle/candN/initial.pdb`. It contains no diffraction or profile output.
Archive timestamps show groups of files written within a short interval, but
there is no accompanying calculation log or index that would map those times
to plotted traces.

Visual inspection of `per_angle_overlays.png` shows panels labeled by angle
and candidate count, for example `30 deg (n=12)`. The blue traces have no
legend, candidate label, line annotation, or other identifier. The image
therefore establishes that multiple candidate profiles were plotted at each
angle, but it does not establish which blue trace corresponds to `cand0`,
`cand1`, and so on.

The repository can correlate candidate names to profiles generated by our own
corrected-diffraction pipeline:

- 124 candidates from angles 27 through 32 have Stage 1 radial profiles.
- 30 selected candidates have Stage 2 refined radial profiles.
- Those mappings follow the explicit `angleNN_candN` model IDs and source PDB
  paths in the Stage 1 and Stage 2 run manifests.

The correlation manifest is:

`outputs/asem_sidechains_20260625/asem_candidate_profile_correlation_manifest.csv`

Its `asem_profile_file` values are intentionally blank. Status values
distinguish `no_asem_profile_file_found`,
`our_corrected_profile_available`, and `not_selected_for_stage2`.

Consequently, candidate structures can be traced by name to our corrected
profiles, but not to Asem's original fast 1D profiles from the files provided.

### Request to Asem for the Missing Mapping

> Could you please send the numerical fast 1D powder profile produced for
> every side-chain candidate, preserving the same `angle/candN` identifiers as
> the structure package? A CSV or text file per candidate is ideal, named for
> example `30/cand0/profile.csv`, with the x-axis definition and units
> (d-spacing, q, 2-theta, or detector radius) and intensity column clearly
> labeled. Please also send the script or command and calculation settings
> used to make `per_angle_overlays.png`, plus a manifest that maps each plotted
> trace to its exact `angle/candN/initial.pdb`. If the original calculation
> used an internal candidate order different from `candN`, please include that
> index mapping. A regenerated overlay with a candidate legend, or separate
> labeled panels/files, would also resolve the ambiguity.
