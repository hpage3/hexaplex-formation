# Asem Twist Extended Stack Diffraction Report

## Purpose

Test whether longer synthetic stack lengths improve corrected powder-diffraction discrimination between Asem 29, 30, and 31 degree twist models.

## Inputs

Existing Asem 29/30/31 models and the Emory-corrected experimental profile were used. The synthetic extensions built by `scripts/build_extended_asem_twist_stacks.py` are layer-equivalent repeats, not new chemically searched candidates.

## Build Method

Complete six-residue layer-equivalent units were repeated with the source model's mean screw transform inferred from layer-centroid PCA axis, rise, and twist. New residue IDs were assigned and heavy-atom deduped XYZ files were generated.

## Diffraction Method

Corrected Asem non-accumulating/vectorized diffraction was run with `tilts = [0]`, `rotations = range(0, 181, 5)`, hydrogens excluded, exact heavy-atom deduplication, and the same grid/radial settings as the corrected baseline.

## Ranking By Length

| length | twist | mean abs offset A | max abs offset A | 3.77 offset A | rank | status |
| --- | --- | --- | --- | --- | --- | --- |
| extended_100layer_equiv | 30 | 0.0481597579226 | 0.0867791603027 | -0.0694748196892 | 1 | survives_near_best_filter |
| extended_100layer_equiv | 29 | 0.0493908632873 | 0.107064464 | -0.0178329726963 | 2 | survives_near_best_filter |
| extended_100layer_equiv | 31 | 0.0574287819563 | 0.157490722046 | -0.000283605539219 | 3 | disfavored_by_peak_offsets |
| extended_32layer_equiv | 29 | 0.0270550946209 | 0.0582217172089 | -0.0178329726963 | 1 | survives_near_best_filter |
| extended_32layer_equiv | 30 | 0.0314350431524 | 0.0582217172089 | -0.000283605539219 | 2 | survives_near_best_filter |
| extended_32layer_equiv | 31 | 0.0407607423262 | 0.0863612112229 | -0.0863612112229 | 3 | disfavored_by_peak_offsets |
| extended_64layer_equiv | 30 | 0.037831388524 | 0.0867791603027 | -0.0178329726963 | 1 | survives_near_best_filter |
| extended_64layer_equiv | 31 | 0.0467688377316 | 0.0863612112229 | -0.0863612112229 | 2 | disfavored_by_peak_offsets |
| extended_64layer_equiv | 29 | 0.0493908632873 | 0.107064464 | -0.0178329726963 | 3 | disfavored_by_peak_offsets |
| original | 29 | 0.0235452211895 | 0.0582217172089 | -0.000283605539219 | 1 | survives_near_best_filter |
| original | 30 | 0.0235452211895 | 0.0582217172089 | -0.000283605539219 | 2 | survives_near_best_filter |
| original | 31 | 0.0407607423262 | 0.0863612112229 | -0.0863612112229 | 3 | disfavored_by_peak_offsets |

## Key Questions

- extended_32layer_equiv: 29/30 tie persists? yes; 31 status: disfavored_by_peak_offsets.
- extended_64layer_equiv: 29/30 tie persists? no; 31 status: disfavored_by_peak_offsets.
- extended_100layer_equiv: 29/30 tie persists? yes; 31 status: disfavored_by_peak_offsets.
- Does extending stack length separate 29 from 30? Partly. They remain within the near-best tolerance at 32 and 100 layer-equivalent lengths, but 64 layer-equivalent stacks favor 30 over 29 by the current mean-offset metric.
- Does 31 become more disfavored as length increases? Mostly yes relative to the original at 64 and 100 layers, while 32 layers approximately preserves the original 31-degree mean offset.
- Does the 3.77 A window remain the main discriminator? At original and 32/64 layers, yes: 31 carries the large 3.77 A offset. At 100 layers, 31's 3.77 A peak returns near the experimental position and the largest 31-degree error moves to the 7.3 A window, which is a warning that synthetic coherent repetition can change the apparent discriminator.
- Are any new windows sensitive to length? Yes. The 7.3 A and 4.4 A windows become more length-sensitive in the 100-layer synthetic stacks, especially for 31 degrees.
- Are peak positions stable but intensities changing? Not entirely. Some peak positions remain stable, but several extended-stack windows jump between nearby radial bins, so this is not only an intensity effect.
- Do longer synthetic stacks sharpen or shift radial features? They can do both. The 64/100 layer-equivalent results show enough peak-position changes to treat the longest stacks as coherence sensitivity tests, not literal refined models.
- Does the result strengthen near-30 family or exact-30 framing? Overall it still supports a near-30 family. The 64-layer result hints at possible 30-degree preference, but 100 layers returns 29 to the near-best set, so exact 30 degrees is not uniquely established.

## Interpretation

The length-extension test is mixed. It does not produce a clean monotonic separation of 29 from 30. It does keep 31 disfavored at every tested length, but the feature responsible for that disfavoring changes in the 100-layer synthetic stack. That pattern is consistent with finite-length and coherence sensitivity, and it cautions against overclaiming either exact 30-degree selection or a literal long-stack structural prediction from these synthetic repeats.

## Synthetic-Length Artifacts

The inferred centroid screw transform is a layer-equivalent proxy, not a chemically generated periodic model. The 100-layer stacks in particular show feature-window jumps that are likely affected by ideal coherent repetition and radial-bin sampling. These artifacts are useful for sensitivity testing but should not be interpreted as minimized molecular behavior.

## Limitations

- Synthetic repeated stacks are not energy-minimized or independently generated.
- No AmberTools/tleap, pNAB, minimization, notebooks, or new candidate search was run.
- Longer ideal repeats may exaggerate coherent diffraction relative to real disordered fibers.
- This model-length test does not replace experimental structural refinement.

## Recommendation

Use these results as a sensitivity screen. The decisive next chemistry-facing test remains multiple completed full 30 degree candidates from Asem, ideally with carboxylates built, screened with the same workflow.
