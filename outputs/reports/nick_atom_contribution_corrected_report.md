# Nick Atom Contribution Comparison: Asem-Corrected

## Purpose

Redo Nick's first powder diffraction model-comparison task using the Asem-corrected non-accumulating/vectorized diffraction path and corrected rotation sampling.

## Inputs

| model | label | atoms | path |
| --- | --- | --- | --- |
| with_coo | With COO/full model | 4620 | inputs\nick_asem_models\nick_hexaplex_with_coo.xyz |
| only_bases | Bases only | 1034 | inputs\nick_asem_models\nick_hexaplex_only_bases.xyz |
| eight_hexads | 8 hexads | 2064 | inputs\nick_asem_models\nick_hexaplex_8hexads.xyz |

- Experimental profile: `inputs\experimental\nick_powder_profile.csv`
- Experimental rows: 442
- Experimental d-spacing range: 2.2003638253 A to 9.73072852774 A

## Corrected Rotation Handling

- Engine: `reference/asem_corrected_diffraction_engine/`.
- Asem correction: orientation stacks are generated independently, so azimuthal rotations do not accumulate.
- Vectorized path: `make_oriented_coords` plus `generate_fiber_diffraction_series`.
- Nick-style tilts: `[0]`.
- Nick-style rotations: `range(0, 181, 5)`, producing 37 rotations from 0 to 180 degrees.
- 180 degrees was included to reproduce the explicit inclusive 0-180 degree interpretation. A later sensitivity check can compare `range(0, 180, 5)` if endpoint duplication after image spinning/symmetrization matters.
- The ambiguous `range(0, 5, 180)` form was not used.

## Detector and Radial Settings

- Grid size: 129 x 129
- Detector half-width: 100 mm
- Detector distance: 338.4 mm
- Wavelength: 0.7749 A
- Radial bins: 420
- Profile normalization: max mean intensity at q >= 0.15 A^-1.

## Primary Feature-Window Tops

| feature_window | top_simulated_model | top_simulated_peak_norm | experimental_peak_d_A |
| --- | --- | --- | --- |
| 3.4 A | only_bases | 0.0394129577368 | 3.3470998911 |
| 3.7 A | eight_hexads | 0.00721467243265 | 3.75740106219 |
| 4.4 A | eight_hexads | 0.00730262798683 | 4.35706110773 |
| 5.6 A | with_coo | 0.0113868832852 | 5.52149027844 |
| 7.25 A | eight_hexads | 0.00734404567679 | 7.25282376927 |

## Full Feature Summary

| model | feature_window | simulated_peak_d_A | simulated_peak_intensity_norm | experimental_peak_d_A | experimental_peak_intensity_norm | peak_offset_d_A | window_area_simulated |
| --- | --- | --- | --- | --- | --- | --- | --- |
| with_coo | 3.4 A | 3.41583473048 | 0.0202216355713 | 3.3470998911 | 0.939593860844 | 0.0687348393797 | 0.00132454084051 |
| with_coo | 3.7 A | 3.79471639446 | 0.00510549339756 | 3.75740106219 | 0.553666157405 | 0.0373153322714 | 0.000678521309208 |
| with_coo | 4.4 A | 4.34277828279 | 0.00657505813318 | 4.35706110773 | 1 | -0.0142828249415 | 0.000912180538346 |
| with_coo | 5.6 A | 5.62272176729 | 0.0113868832852 | 5.52149027844 | 0.506916241991 | 0.101231488852 | 0.00138653295493 |
| with_coo | 7.25 A | 7.41213881028 | 0.00663877663594 | 7.25282376927 | 0.409674037895 | 0.159315041006 | 0.00256161416565 |
| with_coo | 3.0 A | 3.07477148552 | 0.00287346151064 | 3.09467544303 | 0.0797339540278 | -0.0199039575014 | 0.000429020644949 |
| with_coo | 4.1 A | 4.01945598315 | 0.00363785156029 | 4.19767195491 | 0.499087626757 | -0.178215971755 | 0.000527922403055 |
| with_coo | 4.5-5.0 A | 4.53836680174 | 0.00592307233353 | 4.50213684803 | 0.405246245447 | 0.0362299537087 | 0.00151783080465 |
| with_coo | 5.5 A | 5.46824763197 | 0.0110588962294 | 5.52149027844 | 0.506916241991 | -0.0532426464659 | 0.00131414169613 |
| with_coo | 7.0 A | 7.01663454308 | 0.00522007771543 | 7.18442637631 | 0.396703619703 | -0.167791833232 | 0.00134195069451 |
| with_coo | 8.4 A | 8.4458409689 | 0.00595638180683 | 8.37530708339 | 0.104851823753 | 0.0705338855085 | 0.00133828080284 |
| only_bases | 3.4 A | 3.41583473048 | 0.0394129577368 | 3.3470998911 | 0.939593860844 | 0.0687348393797 | 0.00213404793703 |
| only_bases | 3.7 A | 3.70863878878 | 0.00428620412705 | 3.75740106219 | 0.553666157405 | -0.0487622734122 | 0.000384250497499 |
| only_bases | 4.4 A | 4.34277828279 | 0.00291619186106 | 4.35706110773 | 1 | -0.0142828249415 | 0.000335039908103 |
| only_bases | 5.6 A | 5.62272176729 | 0.000995872714961 | 5.52149027844 | 0.506916241991 | 0.101231488852 | 0.00014329750325 |
| only_bases | 7.25 A | 7.01663454308 | 0.000700234825029 | 7.25282376927 | 0.409674037895 | -0.236189226187 | 0.000243371484561 |
| only_bases | 3.0 A | 3.06352544086 | 0.0018275971786 | 3.09467544303 | 0.0797339540278 | -0.0311500021709 | 0.000230733110784 |
| only_bases | 4.1 A | 4.05964776927 | 0.00180691444472 | 4.19767195491 | 0.499087626757 | -0.138024185639 | 0.00020824087953 |
| only_bases | 4.5-5.0 A | 4.69745715956 | 0.00330197265519 | 4.50213684803 | 0.405246245447 | 0.195320311531 | 0.00105510750166 |
| only_bases | 5.5 A | 5.46824763197 | 0.00105232291215 | 5.52149027844 | 0.506916241991 | -0.0532426464659 | 0.000148967089032 |
| only_bases | 7.0 A | 7.01663454308 | 0.000700234825029 | 7.18442637631 | 0.396703619703 | -0.167791833232 | 0.000171192740036 |
| only_bases | 8.4 A | 8.4458409689 | 0.00611080958574 | 8.37530708339 | 0.104851823753 | 0.0705338855085 | 0.00123826159532 |
| eight_hexads | 3.4 A | 3.488064464 | 0.012587427099 | 3.3470998911 | 0.939593860844 | 0.140964572897 | 0.00105728773408 |
| eight_hexads | 3.7 A | 3.62656224876 | 0.00721467243265 | 3.75740106219 | 0.553666157405 | -0.130838813429 | 0.000835853795759 |
| eight_hexads | 4.4 A | 4.46291994427 | 0.00730262798683 | 4.35706110773 | 1 | 0.105858836537 | 0.000988199197028 |
| eight_hexads | 5.6 A | 5.62272176729 | 0.0111330730924 | 5.52149027844 | 0.506916241991 | 0.101231488852 | 0.00136774057855 |
| eight_hexads | 7.25 A | 7.41213881028 | 0.00734404567679 | 7.25282376927 | 0.409674037895 | 0.159315041006 | 0.00293217357255 |
| eight_hexads | 3.0 A | 3.06352544086 | 0.00331858867685 | 3.09467544303 | 0.0797339540278 | -0.0311500021709 | 0.000515405007464 |
| eight_hexads | 4.1 A | 4.14259743655 | 0.00429520806565 | 4.19767195491 | 0.499087626757 | -0.0550745183602 | 0.00062489719917 |
| eight_hexads | 4.5-5.0 A | 4.53836680174 | 0.00630693459223 | 4.50213684803 | 0.405246245447 | 0.0362299537087 | 0.00182913050045 |
| eight_hexads | 5.5 A | 5.46824763197 | 0.012184302949 | 5.52149027844 | 0.506916241991 | -0.0532426464659 | 0.00137027524433 |
| eight_hexads | 7.0 A | 7.01663454308 | 0.00623645613302 | 7.18442637631 | 0.396703619703 | -0.167791833232 | 0.00172148054816 |
| eight_hexads | 8.4 A | 8.26823201483 | 0.00691661152145 | 8.37530708339 | 0.104851823753 | -0.10707506856 | 0.0018018904914 |

## Difference Summary

| comparison | feature_window | difference_peak_d_A | difference_peak_intensity_norm | window_area_difference |
| --- | --- | --- | --- | --- |
| with_coo_minus_bases_only | 3.4 A | 3.43426028921 | 0.000203369801404 | -0.000857100189565 |
| with_coo_minus_bases_only | 3.7 A | 3.7977753059 | 0.00266475760009 | 0.000306865375589 |
| with_coo_minus_bases_only | 4.4 A | 4.46562847608 | 0.0041278909294 | 0.000654359117132 |
| with_coo_minus_bases_only | 5.6 A | 5.62380422692 | 0.0102498016403 | 0.0015174592954 |
| with_coo_minus_bases_only | 7.25 A | 7.40756395996 | 0.00621917670689 | 0.00243209402787 |
| with_coo_minus_8hexads | 3.4 A | 3.41735261402 | 0.0067538269351 | 0.000268179451084 |
| with_coo_minus_8hexads | 3.7 A | 3.7977753059 | 4.28477183141e-05 | -0.000174557695404 |
| with_coo_minus_8hexads | 4.4 A | 4.4149054505 | -0.000248231008523 | -8.75253401321e-05 |
| with_coo_minus_8hexads | 5.6 A | 5.69988876529 | 0.0013575961453 | 3.61283129277e-05 |
| with_coo_minus_8hexads | 7.25 A | 7.34838709677 | -0.000217425722955 | -0.000394063831099 |
| 8hexads_minus_bases_only | 3.4 A | 3.43426028921 | 0.000833262578648 | -0.00112527964065 |
| 8hexads_minus_bases_only | 3.7 A | 3.62869855395 | 0.0048427299083 | 0.000481423070993 |
| 8hexads_minus_bases_only | 4.4 A | 4.46562847608 | 0.00497670956349 | 0.000741884457264 |
| 8hexads_minus_bases_only | 5.6 A | 5.62380422692 | 0.0100146086752 | 0.00148133098247 |
| 8hexads_minus_bases_only | 7.25 A | 7.41601779755 | 0.00691017251865 | 0.00282615785896 |

## Preliminary Interpretation

- 3.4 A: strongest simulated peak is in `only_bases` for this corrected reproduction. Treat this as support for stacked-base contribution only if bases-only remains competitive in area and peak intensity.
- 7.25 A: strongest simulated peak is in `eight_hexads`; compare the full-minus-bases and full-minus-8hexads differences before assigning it to backbone/full-model structure.
- 3.7 A, 4.4 A, and 5.6 A: these windows should be treated as mixed until atom-group controls separate backbone and Glu/COO contributions more directly.

## Limitations

- Nick's source archive predates Asem's rotation fix.
- This is a corrected reproduction/sensitivity analysis, not a final structural fit.
- Compare q/d-spacing positions and relative feature trends, not raw 2D detector image shape alone.
- This is not yet the length-convergence or twist-refinement analysis.
