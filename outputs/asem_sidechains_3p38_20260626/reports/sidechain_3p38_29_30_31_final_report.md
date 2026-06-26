# Draft Report: Asem 3.38 A Glu Side-Chain 29/30/31 Screen

## 1. Executive Summary

We screened Asem's complete Glu side-chain 3.38 A rise structures for twists 29, 30, and 31 deg using the corrected-diffraction radial-window scoring workflow. This is a screening result, not structural proof.

The current result is:

- 29 deg is disfavored relative to 30/31 under the current 3.38 A radial-window score.
- 30 deg and 31 deg remain tied/ambiguous under the current scoring.
- The best 3.38 A rows for 29/30/31 improve relative to the prior 3.40 A Stage 2 subset under the same peak-window scoring.
- This is a by-twist comparison. Candidate IDs should not be assumed identical across the 3.38 A and 3.40 A sets.
- Broad 28/32 expansion is not the most urgent immediate next step.
- The more informative next move is a focused half-degree test around the 30/31 boundary, especially 30.5 deg, if Asem or the generation pipeline can produce equivalent 3.38 A Glu-side-chain candidates.

## 2. What Changed Since The Prior 3.40 A Screen

The previous side-chain screen used Asem's 3.40 A rise candidate set. Asem then provided a separate 3.38 A rise structure set, which is closer to Nick's base-stacking spacing target. The new source contains complete Glu side-chain structures across angles 18-32, but the initial Nick-guided analysis intentionally used only 29, 30, and 31 deg.

The main change is therefore the input rise and candidate set, not the scoring convention. The comparison below uses the same target-window peak-position scoring framework, while keeping clear that the 3.38 A candidate IDs are not mapped one-to-one to 3.40 A candidate IDs.

## 3. Input And Provenance Summary

The 3.38 A structure import is documented in:

- `inputs/asem_sidechains_3p38_20260626/README.md`
- `outputs/asem_sidechains_3p38_20260626/asem_sidechain_3p38_candidate_manifest.csv`
- `outputs/asem_sidechains_3p38_20260626/asem_sidechain_3p38_source_inventory_summary.md`
- `outputs/asem_sidechains_3p38_20260626/sidechain_3p38_29_30_31_run_manifest.csv`

Source provenance:

- Source folder: `C:\Users\hpage3\OneDrive - Georgia Institute of Technology\Documents\GitHub\research\tleap_structures_rise3.38`
- Imported on: 2026-06-26
- Raw repo path: `inputs/asem_sidechains_3p38_20260626/raw/`
- Asem provided these as the rise 3.38 A structure set.
- The source did not include YAML, CSV, TXT, ZIP, PDF, log, or raw powder-profile metadata.
- The 3.38 A provenance comes from Asem's communication and source folder name, not from an included YAML file.
- Raw PDBs are preserved under `raw/`.

Total candidate count in the 3.38 A import: 488 structures. The Nick-guided initial subset contained:

| twist_deg | candidate_count |
| ---: | ---: |
| 29 | 14 |
| 30 | 10 |
| 31 | 26 |
| **Total** | **50** |

Representative 29/30/31 structures had 3810 PDB atoms, 2310 heavy atoms after hydrogen removal, residues CYP/GLU/MEP, and GLU side-chain atom names CG, CD, OE1, and OE2 present.

## 4. Methods

The workflow used existing repo conventions:

1. Convert complete Glu side-chain PDBs to heavy-atom XYZ files using the audited converter in `reference/asem_corrected_diffraction_engine/pdb_to_xyz.py`.
2. Run corrected diffraction with `scripts/run_asem_corrected_diffraction.py`.
3. Generate radial averages with `reference/asem_corrected_diffraction_engine/radial_average.py`.
4. Extract observed peaks from the radial profiles using `scripts/extract_twist_rise_profile_peaks.py`.
5. Score the observed peaks against `inputs/experimental_peak_windows/hxc590_twist_rise_scan_targets.csv` with `scripts/score_twist_rise_scan.py`.
6. Rank rows by completeness, combined RMSD, and helical RMSD with `scripts/summarize_twist_rise_scored_manifest.py`.

Corrected-diffraction settings:

- Grid size: 65
- Grid limit: 100 mm
- Theta count: 4
- Phi count: 8
- Psi count: 1
- Theta max: 180 deg
- Orientations: 32
- Radial bins: 120

The 2D corrected-diffraction images are diagnostic only. The primary comparison here is powder/radial-profile agreement through the target-window score.

## 5. 3.38 A 29/30/31 Screening Results

The 3.38 A batch produced and validated:

- 50 corrected-diffraction NPY files
- 50 corrected-diffraction PNG files
- 50 corrected-diffraction JSON files
- 50 radial CSV profiles
- 250 long observed-peak rows, corresponding to 50 models x 5 target windows
- 50 wide observed-peak rows
- 50 scored manifest rows
- 50 ranked manifest rows

All scored rows were complete across the five target windows: base, A, B, C, and D.

The main screening read is that 29 deg is weaker than 30/31, while 30 and 31 remain unresolved under this scoring. This should not be interpreted as proof that 30 deg wins, and it should not be interpreted as excluding 31 deg.

## 6. Best Candidates And Best-Per-Twist Table

3.38 A best per twist:

| twist_deg | best_model | combined_rmsd | helical_rmsd | base_rmsd | current_read |
| ---: | --- | ---: | ---: | ---: | --- |
| 29 | angle29_cand6 | 0.088310743 | 0.098397704 | 0.016294279 | weaker than 30/31 |
| 30 | angle30_cand9 | 0.061332690 | 0.068086325 | 0.016294279 | tied with best 31 |
| 31 | angle31_cand19 | 0.061332690 | 0.068086325 | 0.016294279 | tied with best 30 |

Top-ranked 3.38 A rows begin with `angle30_cand9`, followed by several 31 deg candidates with identical current scores. This makes 30/31 ambiguous under the current radial-window metric.

## 7. 3.38 A vs 3.40 A Comparison

The comparison uses the prior 3.40 A Stage 2 top-five-per-angle subset, not a candidate-identity-matched paired design. Candidate IDs should not be assumed identical across 3.38 and 3.40.

| twist | 3.38 best | 3.38 combined | 3.40 best | 3.40 combined | delta 3.38-3.40 |
| ---: | --- | ---: | --- | ---: | ---: |
| 29 | angle29_cand6 | 0.088310743 | angle29_cand11 | 0.101368623 | -0.013057880 |
| 30 | angle30_cand9 | 0.061332690 | angle30_cand10 | 0.088310743 | -0.026978053 |
| 31 | angle31_cand19 | 0.061332690 | angle31_cand10 | 0.088310743 | -0.026978053 |

Negative deltas indicate lower combined RMSD for the 3.38 A row. Under the same peak-window scoring, the best 3.38 A by-twist rows improve relative to the prior 3.40 A Stage 2 subset for 29, 30, and 31 deg.

This is encouraging for the 3.38 A rise hypothesis, but it remains a screening comparison rather than final structural proof.

## 8. Side-Chain Interpretation And No-Glu Ablation Caution

The complete-Glu side-chain structures are the primary comparison requested here. Prior no-Glu ablation at 3.40 A was a sensitivity analysis, not a replacement scientific model. In that prior ablation, removing the Glu side chains slightly improved mean peak-position RMSD for some rows and worsened others:

- Matched candidates: 30
- Mean delta base RMSD, no-Glu minus with-Glu: 0.000000 A
- Mean delta helical RMSD: -0.005556 A
- Mean delta combined RMSD: -0.004961 A
- Combined RMSD improved/worsened/unchanged after stripping: 11/6/13

That result says the side-chain contribution is not trivially resolved by this peak-position metric alone. It does not establish a chemical mechanism, and it does not mean complete Glu side chains should be discarded. The current 3.38 A analysis should therefore be read as a complete-side-chain screen with a side-chain sensitivity caveat.

## 9. Relationship To Asem's Debye Overlay Provenance

Asem also supplied a per-candidate overlay PDF for the earlier side-chain package. That PDF visually maps candidate names to Asem-reported Debye profile overlays and reports r, Rwp, and energy values. It does not provide raw numeric profile tables by itself.

The repo extracted 185 labeled candidate panels from that PDF and matched all of them to the imported 3.40 A candidate manifest. Asem's Debye overlay ranking is useful provenance, but it is not the same as this repo's corrected-diffraction radial scoring. The current 3.38 A results should be compared using the corrected radial profiles and target-window scores generated in this repo.

## 10. Decision: 28/32 Expansion vs 30.5 Deg Half-Degree Test

Based on the current 29/30/31-only 3.38 A screen:

- 29 deg is disfavored relative to 30/31 under the current radial-window score.
- 30 deg and 31 deg remain tied/ambiguous.
- Broad 28/32 expansion is not the most urgent immediate next step.
- The more informative next scientific move is a focused half-degree test around the 30/31 boundary, especially 30.5 deg, if Asem or the generation pipeline can produce equivalent 3.38 A Glu-side-chain candidates.

A broad 28/32 expansion could still be useful later if plot review suggests edge behavior remains important, but it is not the first thing this screen points to.

## 11. Limitations

- This is a screening workflow, not a final refinement or structural proof.
- The peak-window score uses radial peak positions; it does not capture every aspect of powder-profile intensity or molecular plausibility.
- 2D corrected-diffraction images are diagnostic only.
- Candidate IDs are not assumed to map across the 3.38 A and 3.40 A structure sets.
- The 3.38 A source folder did not include YAML or raw numeric powder-profile metadata.
- The no-Glu ablation suggests side-chain effects are metric-sensitive and should be interpreted cautiously.

## 12. Recommended Next Step

The next scientific step should be a focused 3.38 A half-degree test near the 30/31 boundary, especially 30.5 deg, using complete Glu side-chain candidates generated equivalently to the current Asem set. After that, compare the 30, 30.5, and 31 deg radial profiles directly. Broad 28/32 expansion can wait unless plot review indicates edge behavior remains competitive.

## 13. Suggested Email/Update Wording For Nick

Hi Nick,

I reran the side-chain screen using Asem's 3.38 A rise structures, focusing first on the 29, 30, and 31 deg candidates as discussed. The 3.38 A set looks better than the prior 3.40 A subset under the same radial peak-window scoring, but the result is still a screen rather than a structural proof.

The best 29 deg candidate is weaker than the best 30/31 candidates. The best 30 deg candidate and several 31 deg candidates are tied under the current score, so I would not say 30 deg wins or that 31 deg is excluded. The result mainly narrows the question to the 30/31 boundary.

My cautious recommendation is not to jump immediately to broad 28/32 expansion. A more informative next test would be equivalent 3.38 A Glu-side-chain candidates around 30.5 deg, if Asem or the generation pipeline can produce them. That should help decide whether the apparent 30/31 tie reflects binning/score resolution or a real broad optimum.

## 14. Appendix: Committed File Paths And Commit IDs

Recent commits:

- `be425ec` Add Asem 3.38 A 29-31 side-chain screening results
- `15102b4` Import Asem 3.38 A side-chain candidate structures
- `c173e7a` Add Asem per-candidate overlay provenance
- `afcdc72` Add side-chain no-Glu ablation comparison
- `2baef89` Add side-chain Stage 2 scatter plots
- `b2052fa` Document 3.38 A anti-parallel side-chain plan
- `1c251d8` Add side-chain Stage 2 refined screening summaries

Primary committed inputs and summaries:

- `inputs/asem_sidechains_3p38_20260626/README.md`
- `outputs/asem_sidechains_3p38_20260626/asem_sidechain_3p38_candidate_manifest.csv`
- `outputs/asem_sidechains_3p38_20260626/asem_sidechain_3p38_source_inventory_summary.md`
- `outputs/asem_sidechains_3p38_20260626/sidechain_3p38_29_30_31_run_manifest.csv`
- `outputs/asem_sidechains_3p38_20260626/stage1_29_30_31_glu/sidechain_3p38_stage1_batch_summary.md`
- `outputs/asem_sidechains_3p38_20260626/stage1_29_30_31_glu/sidechain_3p38_decision_summary.md`
- `outputs/asem_sidechains_3p38_20260626/stage1_29_30_31_glu/scoring/sidechain_3p38_scored_manifest_ranked.csv`
- `outputs/asem_sidechains_3p38_20260626/stage1_29_30_31_glu/scoring/sidechain_3p38_summary.md`
- `outputs/asem_sidechains_3p38_20260626/stage1_29_30_31_glu/comparison_3p38_vs_3p40/sidechain_3p38_vs_3p40_best_per_twist.csv`
- `outputs/asem_sidechains_3p38_20260626/stage1_29_30_31_glu/comparison_3p38_vs_3p40/sidechain_3p38_vs_3p40_summary.md`
- `outputs/asem_sidechains_20260625/stage2_top5_refined/scoring/sidechain_stage2_top5_summary.md`
- `outputs/asem_sidechains_20260625/stage2_top5_refined/scoring/sidechain_stage2_top5_scored_manifest_ranked.csv`
- `outputs/asem_sidechains_20260625/stage2_top5_no_glu_ablation/scoring/sidechain_ablation_with_vs_without_summary.md`
- `outputs/asem_sidechains_20260625/asem_candidate_overlay_pdf_summary.md`
