# Asem struct2 Omega-Series Audit

## Purpose

This audit imports and evaluates Asem's manually constructed omega-series structures derived from `1_558214.pdb`. The question is whether a modest omega increase, especially around 172 degrees, is a better compromise than the original/current baseline.

## Source and Scope

- Source folder: `C:\Users\hpage3\OneDrive - Georgia Institute of Technology\Documents\GitHub\research\struct2`
- Imported folder: `inputs/asem_omega_series_struct2/`
- No pNAB work, minimization, notebook execution, rise refinement, twist change, or preserved archive edits were performed.

The files contain 2,847 atoms and 90 CYP/MEP residues. They are larger than the earlier local 579-atom candidate files, but they are not the full 7,146-atom ideal antiparallel model used for current corrected diffraction baselines. Diffraction was therefore deferred pending a larger/full periodic model for the selected omega candidate.

## Asem and Literature Context

Asem reported that he manually increased the peptide dihedral while holding other atoms fixed, with the compensating N-CA(next)-C angle increasing from 113.01 degrees at omega 167 to 123.85 degrees at omega 180. This creates an explicit tradeoff: improved omega can come with increasing angle strain.

Omega closer to 180 degrees is not automatically better. The Vitagliano omega work is relevant because real peptide omega angles can deviate from 180 degrees and those deviations correlate with local backbone geometry; forcing one ideal value can be too restrictive. Howard's theta-polypeptide note is also relevant as a cautionary analogy: substantial omega departures can occur in constrained polypeptide systems with cross-residue bonding, though this audit does not claim the same mechanism here.

## Atom-Mapping Notes

The script reports Asem's target omega and N-CA(next)-C angle from the provided table. For continuity with commit `4daaa450cd7bacec0db9070d41f2f36f7b6806c7`, it also reports the prior Nick-note-matching `N-CA-C-N'` omega proxy and canonical `CA-C-N'-CA'` torsion. Those two torsions remain effectively constant across this exported series. The moving angle diagnostic that tracks Asem's table is measured as `CB-C'-O'`, because the exported atom names do not expose a literal N-CA(next)-C mapping.

## Structural Audit

| Target omega | Asem angle | Measured angle proxy | RMSD (A) | Max disp. (A) | H-bond N...O (A) | <2.2 A contacts | Status |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 167 | 113.01 | 113.25 | 0.0396 | 0.2293 | 2.774 | 84 | original_baseline |
| 168 | 113.81 | 113.67 | 0.0304 | 0.1746 | 2.790 | 84 | ambiguous_requires_visual_review |
| 169 | 114.62 | 114.23 | 0.0225 | 0.1199 | 2.808 | 84 | ambiguous_requires_visual_review |
| 170 | 115.43 | 115.11 | 0.0181 | 0.0664 | 2.825 | 126 | omega_improved_modest_angle_strain |
| 171 | 116.25 | 115.93 | 0.0193 | 0.0804 | 2.842 | 126 | omega_improved_modest_angle_strain |
| 172 | 117.08 | 116.81 | 0.0256 | 0.1028 | 2.861 | 126 | likely_best_compromise |
| 173 | 117.91 | 117.66 | 0.0341 | 0.1248 | 2.878 | 126 | omega_improved_modest_angle_strain |
| 174 | 118.75 | 118.51 | 0.0435 | 0.1596 | 2.896 | 126 | omega_improved_high_angle_strain |
| 175 | 119.59 | 119.35 | 0.0534 | 0.2137 | 2.913 | 126 | omega_improved_high_angle_strain |
| 176 | 120.43 | 120.18 | 0.0635 | 0.2678 | 2.931 | 126 | omega_improved_high_angle_strain |
| 177 | 121.28 | 121.04 | 0.0738 | 0.3227 | 3.709 | 126 | hbond_degraded |
| 178 | 122.13 | 121.87 | 0.0841 | 0.3768 | 3.709 | 126 | hbond_degraded |
| 179 | 122.99 | 122.72 | 0.0945 | 0.4309 | 3.709 | 126 | hbond_degraded |
| 180 | 123.85 | 123.56 | 0.1049 | 0.4843 | 3.709 | 126 | hbond_degraded |

## Interpretation

- Omega target improves monotonically from 167 to 180 degrees by construction.
- The compensating angle strain also increases monotonically, from 113.01 to 123.85 degrees in Asem's table.
- Most atoms are fixed; the changed atom names are `O'`, `N''`, and `H''`.
- The 172-degree candidate has RMSD 0.0256 A and max displacement 0.1028 A relative to the original.
- The automated H-bond and steric proxies are almost unchanged across the series; they should not replace Nick/Asem visual/chemical inspection.

## Recommendation

The 172-degree candidate is the most defensible first review candidate: it hits Asem's suggested compromise region while keeping the angle strain below the steeper high-omega end of the series. If Nick/Asem's visual inspection confirms H-bonding and sterics are acceptable, ask Asem to build the larger diffraction-ready periodic model for `1_558214_omega172.pdb`. A neighboring 171 or 173 degree candidate is also reasonable if visual geometry favors it.

## Diffraction Status

- Diffraction-ready classification: no
- Corrected diffraction was not run because these files are not full/comparable periodic models for the current corrected ideal diffraction baseline.
- Diffraction should be run after Asem builds the larger periodic model for the selected omega candidate.

## Outputs

- `outputs/metrics/asem_omega_series_structural_audit.csv`
- `outputs/metrics/asem_omega_series_ranking.csv`
- `outputs/asem_omega_series_struct2/plots/omega_vs_n_ca_next_c_angle.png`
- `outputs/asem_omega_series_struct2/plots/omega_vs_hbond_distance.png`
- `outputs/asem_omega_series_struct2/plots/omega_vs_clash_count.png`
- `outputs/asem_omega_series_struct2/plots/omega_vs_rmsd.png`
- `outputs/asem_omega_series_struct2/plots/omega_tradeoff_score.png`

## Limitations

- Automated H-bond proxy does not replace visual/chemical inspection.
- N-CA(next)-C angle strain is a geometry proxy, not an energy calculation.
- No minimization was performed.
- No diffraction was run because a full diffraction-ready periodic model is still needed.
