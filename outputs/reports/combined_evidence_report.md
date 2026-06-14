# Combined Evidence Report: Hexaplex Formation Metrics

## Scientific cautions

- The d ~= 4.5 A feature is a reciprocal-space feature, not a literal atom-distance or atom-contact assignment.
- Contact maps, block contact decomposition, and GLU motif counts are real-space structural summaries.
- Debye scores are simplified comparative approximations, not replacements for full fiber-diffraction simulations.
- Component intensities are comparative controls, not additive decompositions of full intensity.
- The candidate block map remains unvalidated against PyMOL colored strand paths.
- Candidate ladder models do not prove temporal assembly order.
- Hydrogen-bond candidates are rough geometric candidates only.

## Main evidence table

| model_name | included_blocks | includes_hexads | angular_coverage_rad_fitted | approximate_turns | contact_count_4p5A | motif_GLU_GLU | scaffold_between_block_contacts | scaffold_hexad_or_other_contacts | d_4p5_fraction | d_3p4_fraction |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| hexads_plus_scaffold_blocks_1_heavy_deduped | 1 | yes | 6.217255 | 2.256204 | 301 | 14 | 0 | 174 | 0.007370 | 0.061783 |
| scaffold_blocks_1_heavy_deduped | 1 | no | 6.233895 | 1.346288 | 57 | 14 | 0 | 0 | 0.014388 | 0.014846 |
| hexads_plus_scaffold_blocks_1_2_heavy_deduped | 1,2 | yes | 6.106985 | 1.096362 | 374 | 43 | 30 | 174 | 0.009785 | 0.059538 |
| scaffold_blocks_1_2_heavy_deduped | 1,2 | no | 6.075752 | 0.070355 | 158 | 43 | 44 | 0 | 0.014381 | 0.012721 |
| hexads_plus_scaffold_blocks_1_2_3_heavy_deduped | 1,2,3 | yes | 6.184021 | 2.569137 | 430 | 57 | 43 | 174 | 0.011554 | 0.061876 |
| scaffold_blocks_1_2_3_heavy_deduped | 1,2,3 | no | 6.237544 | 1.566860 | 228 | 57 | 57 | 0 | 0.013679 | 0.013829 |
| hexads_plus_scaffold_blocks_1_2_3_4_heavy_deduped | 1,2,3,4 | yes | 6.241806 | 1.406890 | 503 | 86 | 73 | 174 | 0.014856 | 0.064457 |
| scaffold_blocks_1_2_3_4_heavy_deduped | 1,2,3,4 | no | 6.214849 | 0.376982 | 329 | 86 | 101 | 0 | 0.016290 | 0.013753 |
| hexads_plus_scaffold_blocks_1_2_3_4_5_heavy_deduped | 1,2,3,4,5 | yes | 6.256888 | 2.879391 | 559 | 100 | 86 | 174 | 0.018306 | 0.072219 |
| scaffold_blocks_1_2_3_4_5_heavy_deduped | 1,2,3,4,5 | no | 6.174892 | 1.861070 | 399 | 100 | 114 | 0 | 0.017841 | 0.015272 |
| hexads_plus_scaffold_blocks_1_2_3_4_5_6_heavy_deduped | 1,2,3,4,5,6 | yes | 6.156728 | 1.758371 | 645 | 129 | 129 | 174 | 0.021626 | 0.077104 |
| scaffold_blocks_1_2_3_4_5_6_heavy_deduped | 1,2,3,4,5,6 | no | 6.071303 | 0.750925 | 513 | 129 | 171 | 0 | 0.019043 | 0.015259 |
| reference_full_hexaplex_heavy_deduped |  | yes | 6.156728 | 0.758371 | 645 | 129 | 129 | 174 | 0.021626 | 0.077104 |
| reference_hexads_only_heavy_deduped |  | yes | 5.767437 | 0.834589 | 258 | 0 | 0 | 174 | 0.004412 | 0.061171 |
| reference_scaffold_only_complement_heavy_deduped | 1,2,3,4,5,6 | no | 6.071303 | 0.750925 | 513 | 129 | 171 | 0 | 0.019043 | 0.015259 |

## Generated observations

- Scaffold block 1 has near-full fitted angular coverage (6.234 rad), consistent with an individually folded/twisted path descriptor.
- The d_4p5_fraction is already present in scaffold block 1 (0.014388); this is a reciprocal-space-window score, not a literal contact assignment.
- The d_4p5_fraction increases from scaffold block 1 to the full scaffold row.
- motif_GLU_GLU increases with scaffold block count in the scaffold-only ladder.
- Between-block scaffold contacts increase from scaffold block 1 to the full scaffold.
- The mean d_3p4_fraction is higher in hexad-containing rows (0.066907) than scaffold-only rows (0.014420).
- Scaffold-hexad/other contacts appear in full Hexaplex and hexads-plus-scaffold models.

## Correlations

In this simplified comparative dataset, correlations are descriptive screening summaries. They are consistent with possible associations but do not establish causality, temporal order, or physical sufficiency.

- d_4p5_fraction vs motif_GLU_GLU: Pearson r = 0.904 across 15 complete row(s). In this simplified comparative dataset, this is consistent with an association but does not establish causality.
- d_4p5_fraction vs motif_GLU_any: Pearson r = 0.906 across 15 complete row(s). In this simplified comparative dataset, this is consistent with an association but does not establish causality.
- d_4p5_fraction vs contact_count_4p5A: Pearson r = 0.609 across 15 complete row(s). In this simplified comparative dataset, this is consistent with an association but does not establish causality.
- d_4p5_fraction vs scaffold_between_block_contacts: Pearson r = 0.846 across 15 complete row(s). In this simplified comparative dataset, this is consistent with an association but does not establish causality.
- d_4p5_fraction vs angular_coverage_rad_fitted: Pearson r = 0.402 across 15 complete row(s). In this simplified comparative dataset, this is consistent with an association but does not establish causality.
- d_3p4_fraction vs includes_hexads flag: Pearson r = 0.982 across 15 complete row(s). In this simplified comparative dataset, this is consistent with an association but does not establish causality.

## Mechanistic interpretation

Current evidence supports a cautious model in which individual scaffold paths are already folded/twisted and helical-like; the 4.5 A-window score is already detectable in a single scaffold path; multi-block scaffold assembly strengthens GLU-rich motif recurrence and the 4.5 A-window score; hexads contribute strongly to the 3.4 A-window score and add scaffold-hexad/other contacts; and final Hexaplex stabilization is likely cooperative.

## Next steps

- Validate the candidate block map against PyMOL colored strand paths.
- Run full fiber-diffraction simulations on selected ladder models.
- Add fitted-axis/block-specific helical metrics.
- Refine hydrogen-bond candidate analysis with angle, protonation, ion, and water context.
- Build a contact-state pathway model before Schrodinger bridge modeling.
