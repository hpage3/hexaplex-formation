# Asem Side-Chain Candidate Import

## Provenance

- Source folder:
  `C:\Users\hpage3\OneDrive - Georgia Institute of Technology\Documents\GitHub\research\sidechains`
- Imported on: 2026-06-25
- Raw copied package: `inputs\asem_sidechains_20260625\raw\`
- Source files were copied; the original source folder was not modified.

The raw import preserves the expanded angle/candidate folders, `options.yaml`,
`per_angle_overlays.png`, and `tleap_structures.zip`. A checked ZIP entry
(`30/cand0/initial.pdb`) was byte-identical to the corresponding expanded file.
The overlay image is retained as Asem's visual provenance/reference.

After Nick asked whether candidate structures could be correlated by name with
their powder profiles, Asem supplied `pNAB_per_candidate_overlays.pdf`. The PDF
provides per-candidate visual overlays, with each panel labeled by candidate
path such as `29/cand0/initial.pdb`. The panels include Asem's reported `r`,
`Rwp`, and energy values. The PDF does not by itself provide raw numeric profile
data; raw profile tables should be treated as absent unless they are supplied
separately later.

## Asem's Note

Asem reported that he generated pNAB samples with added side chains for twist
angles from 0 to 60 degrees in 1-degree increments at 3.4 A rise. He adjusted
the energy filters to generate more candidates, and not all angles produced
candidates. He did not visually inspect every structure. He also calculated a
fast 1D diffraction pattern, but cautioned that it may not be identical to the
2D diffraction workflow. He observed substantial candidate-to-candidate
variation at a single angle and noted that candidates from roughly 27 to 32
degrees appeared to contain peaks near the experimental feature.

## Cautions

- The imported structures have not been visually validated.
- Only angle folders with produced candidates are present in this package.
- Asem's fast 1D diffraction output is not assumed equivalent to this repo's
  corrected 2D Asem/Nick diffraction and radial-average workflow.
- `options.yaml` is preserved as the pNAB configuration Asem supplied. It
  references component files by paths from Asem's environment; those component
  files are not included in this imported package.
- This directory is raw provenance, not a set of accepted structures or final
  scored results.
- Generated diffraction, radial-profile, and scoring outputs belong under
  `outputs\asem_sidechains_20260625\` and should remain untracked unless
  explicitly reviewed and approved.
