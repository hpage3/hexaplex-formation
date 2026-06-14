# Hexad Arm Atom Geometry Inspection

## Scope and cautions

- This report inspects atom-level geometry only; it does not assign biological meaning to CYP or MEP.
- The central axis is an approximate fitted axis used for radial summaries.
- Proposed fixed/transformable atoms are operational candidates for a later base-length scaling workflow, not chemistry claims.
- The diffraction sweep is not implemented here.

## Input and inferred axis

- PDB: outputs/intermediates/ai_candidate_inputs/full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain.pdb
- Atom count: 3573
- Residue count: 180
- Axis source: backbone_like_heavy_atoms
- Axis origin: (-0.000037, -0.000007, 23.799858)
- Axis direction: (0.001055, 0.001598, 0.999998)

## Inferred chain structure

| chain_id | residue_count | first_residue | last_residue |
| --- | --- | --- | --- |
| A | 30 | CYP1 | GLU30 |
| B | 30 | MEP31 | GLU60 |
| C | 30 | CYP61 | GLU90 |
| D | 30 | MEP91 | GLU120 |
| E | 30 | CYP121 | GLU150 |
| F | 30 | MEP151 | GLU180 |

## Atom-name inventory

| residue_name | atom_name | classification | atom_count | residue_occurrence_count | chains |
| --- | --- | --- | --- | --- | --- |
| CYP | C | backbone_like | 45 | 45 | A;C;E |
| CYP | C2 | candidate_sidechain_or_arm | 45 | 45 | A;C;E |
| CYP | C4 | candidate_sidechain_or_arm | 45 | 45 | A;C;E |
| CYP | C6 | candidate_sidechain_or_arm | 45 | 45 | A;C;E |
| CYP | CA | backbone_like | 45 | 45 | A;C;E |
| CYP | CB | candidate_sidechain_or_arm | 45 | 45 | A;C;E |
| CYP | CC | candidate_sidechain_or_arm | 45 | 45 | A;C;E |
| CYP | H | candidate_sidechain_or_arm | 42 | 42 | A;C;E |
| CYP | H1 | candidate_sidechain_or_arm | 3 | 3 | A;C;E |
| CYP | H2 | candidate_sidechain_or_arm | 3 | 3 | A;C;E |
| CYP | H3 | candidate_sidechain_or_arm | 3 | 3 | A;C;E |
| CYP | HA | candidate_sidechain_or_arm | 45 | 45 | A;C;E |
| CYP | HB2 | candidate_sidechain_or_arm | 45 | 45 | A;C;E |
| CYP | HB3 | candidate_sidechain_or_arm | 45 | 45 | A;C;E |
| CYP | HC2 | candidate_sidechain_or_arm | 45 | 45 | A;C;E |
| CYP | HC3 | candidate_sidechain_or_arm | 45 | 45 | A;C;E |
| CYP | HN3 | candidate_sidechain_or_arm | 45 | 45 | A;C;E |
| CYP | HN5 | candidate_sidechain_or_arm | 45 | 45 | A;C;E |
| CYP | N | backbone_like | 45 | 45 | A;C;E |
| CYP | N1 | candidate_sidechain_or_arm | 45 | 45 | A;C;E |
| CYP | N3 | candidate_sidechain_or_arm | 45 | 45 | A;C;E |
| CYP | N5 | candidate_sidechain_or_arm | 45 | 45 | A;C;E |
| CYP | O | backbone_like | 45 | 45 | A;C;E |
| CYP | OC2 | candidate_sidechain_or_arm | 45 | 45 | A;C;E |
| CYP | OC4 | candidate_sidechain_or_arm | 45 | 45 | A;C;E |
| CYP | OC6 | candidate_sidechain_or_arm | 45 | 45 | A;C;E |
| GLU | C | backbone_like | 90 | 90 | A;B;C;D;E;F |
| GLU | CA | backbone_like | 90 | 90 | A;B;C;D;E;F |
| GLU | CB | glu_scaffold | 90 | 90 | A;B;C;D;E;F |
| GLU | CD | glu_scaffold | 90 | 90 | A;B;C;D;E;F |
| GLU | CG | glu_scaffold | 90 | 90 | A;B;C;D;E;F |
| GLU | H | glu_scaffold | 90 | 90 | A;B;C;D;E;F |
| GLU | HA | glu_scaffold | 90 | 90 | A;B;C;D;E;F |
| GLU | HB2 | glu_scaffold | 90 | 90 | A;B;C;D;E;F |
| GLU | HB3 | glu_scaffold | 90 | 90 | A;B;C;D;E;F |
| GLU | HG2 | glu_scaffold | 90 | 90 | A;B;C;D;E;F |
| GLU | HG3 | glu_scaffold | 90 | 90 | A;B;C;D;E;F |
| GLU | N | backbone_like | 90 | 90 | A;B;C;D;E;F |
| GLU | O | backbone_like | 90 | 90 | A;B;C;D;E;F |
| GLU | OE1 | glu_scaffold | 90 | 90 | A;B;C;D;E;F |
| GLU | OE2 | glu_scaffold | 90 | 90 | A;B;C;D;E;F |
| GLU | OXT | backbone_like | 6 | 6 | A;B;C;D;E;F |
| MEP | C | backbone_like | 45 | 45 | B;D;F |
| MEP | C1 | candidate_sidechain_or_arm | 45 | 45 | B;D;F |
| MEP | C3 | candidate_sidechain_or_arm | 45 | 45 | B;D;F |
| MEP | C5 | candidate_sidechain_or_arm | 45 | 45 | B;D;F |
| MEP | CA | backbone_like | 45 | 45 | B;D;F |
| MEP | CB | candidate_sidechain_or_arm | 45 | 45 | B;D;F |
| MEP | CC | candidate_sidechain_or_arm | 45 | 45 | B;D;F |
| MEP | H | candidate_sidechain_or_arm | 42 | 42 | B;D;F |
| MEP | H1 | candidate_sidechain_or_arm | 3 | 3 | B;D;F |
| MEP | H2 | candidate_sidechain_or_arm | 3 | 3 | B;D;F |
| MEP | H3 | candidate_sidechain_or_arm | 3 | 3 | B;D;F |
| MEP | HA | candidate_sidechain_or_arm | 45 | 45 | B;D;F |
| MEP | HB2 | candidate_sidechain_or_arm | 45 | 45 | B;D;F |
| MEP | HB3 | candidate_sidechain_or_arm | 45 | 45 | B;D;F |
| MEP | HC2 | candidate_sidechain_or_arm | 45 | 45 | B;D;F |
| MEP | HC3 | candidate_sidechain_or_arm | 45 | 45 | B;D;F |
| MEP | HX11 | candidate_sidechain_or_arm | 45 | 45 | B;D;F |
| MEP | HX21 | candidate_sidechain_or_arm | 45 | 45 | B;D;F |
| MEP | HX22 | candidate_sidechain_or_arm | 45 | 45 | B;D;F |
| MEP | HX31 | candidate_sidechain_or_arm | 45 | 45 | B;D;F |
| MEP | HX32 | candidate_sidechain_or_arm | 45 | 45 | B;D;F |
| MEP | N | backbone_like | 45 | 45 | B;D;F |
| MEP | N2 | candidate_sidechain_or_arm | 45 | 45 | B;D;F |
| MEP | N4 | candidate_sidechain_or_arm | 45 | 45 | B;D;F |
| MEP | N6 | candidate_sidechain_or_arm | 45 | 45 | B;D;F |
| MEP | NX1 | candidate_sidechain_or_arm | 45 | 45 | B;D;F |
| MEP | NX2 | candidate_sidechain_or_arm | 45 | 45 | B;D;F |
| MEP | NX3 | candidate_sidechain_or_arm | 45 | 45 | B;D;F |
| MEP | O | backbone_like | 45 | 45 | B;D;F |

## Proposed fixed atoms

| residue_name | atom_name | classification | mean_radial_distance_A | selection_recommendation |
| --- | --- | --- | --- | --- |
| CYP | C | backbone_like | 11.285181 | fixed_candidate |
| CYP | CA | backbone_like | 10.638375 | fixed_candidate |
| CYP | N | backbone_like | 10.993850 | fixed_candidate |
| CYP | O | backbone_like | 11.736500 | fixed_candidate |
| GLU | C | backbone_like | 11.804717 | fixed_candidate |
| GLU | CA | backbone_like | 12.166915 | fixed_candidate |
| GLU | CB | glu_scaffold | 13.665352 | fixed_candidate |
| GLU | CD | glu_scaffold | 15.903977 | fixed_candidate |
| GLU | CG | glu_scaffold | 14.405072 | fixed_candidate |
| GLU | H | glu_scaffold | 11.170382 | fixed_candidate |
| GLU | HA | glu_scaffold | 12.005182 | fixed_candidate |
| GLU | HB2 | glu_scaffold | 14.061176 | fixed_candidate |
| GLU | HB3 | glu_scaffold | 13.867354 | fixed_candidate |
| GLU | HG2 | glu_scaffold | 14.134704 | fixed_candidate |
| GLU | HG3 | glu_scaffold | 14.199542 | fixed_candidate |
| GLU | N | backbone_like | 11.455098 | fixed_candidate |
| GLU | O | backbone_like | 12.354881 | fixed_candidate |
| GLU | OE1 | glu_scaffold | 16.435608 | fixed_candidate |
| GLU | OE2 | glu_scaffold | 16.562854 | fixed_candidate |
| GLU | OXT | backbone_like | 11.066456 | fixed_candidate |
| MEP | C | backbone_like | 11.285113 | fixed_candidate |
| MEP | CA | backbone_like | 10.638351 | fixed_candidate |
| MEP | N | backbone_like | 10.993893 | fixed_candidate |
| MEP | O | backbone_like | 11.736600 | fixed_candidate |

## Proposed transformable atoms

| residue_name | atom_name | mean_radial_distance_A | mean_delta_from_residue_backbone_A | axis_facing_candidate | outward_candidate | selection_recommendation |
| --- | --- | --- | --- | --- | --- | --- |
| CYP | C2 | 6.397742 | -4.765735 | yes | no | transformable_candidate |
| CYP | C4 | 4.127800 | -7.035677 | yes | no | transformable_candidate |
| CYP | C6 | 6.372402 | -4.791075 | yes | no | transformable_candidate |
| CYP | CB | 9.136541 | -2.026936 | yes | no | transformable_candidate |
| CYP | CC | 8.496655 | -2.666822 | yes | no | transformable_candidate |
| CYP | HB2 | 8.881532 | -2.281944 | yes | no | transformable_candidate |
| CYP | HB3 | 8.882876 | -2.280601 | yes | no | transformable_candidate |
| CYP | HC2 | 8.871218 | -2.292259 | yes | no | transformable_candidate |
| CYP | HC3 | 8.870069 | -2.293408 | yes | no | transformable_candidate |
| CYP | HN3 | 4.846345 | -6.317131 | yes | no | transformable_candidate |
| CYP | HN5 | 4.786595 | -6.376881 | yes | no | transformable_candidate |
| CYP | N1 | 6.976040 | -4.187437 | yes | no | transformable_candidate |
| CYP | N3 | 5.011507 | -6.151970 | yes | no | transformable_candidate |
| CYP | N5 | 4.978217 | -6.185260 | yes | no | transformable_candidate |
| CYP | OC2 | 7.308830 | -3.854647 | yes | no | transformable_candidate |
| CYP | OC4 | 2.853081 | -8.310396 | yes | no | transformable_candidate |
| CYP | OC6 | 7.266010 | -3.897467 | yes | no | transformable_candidate |
| MEP | C1 | 6.412087 | -4.751402 | yes | no | transformable_candidate |
| MEP | C3 | 6.407323 | -4.756166 | yes | no | transformable_candidate |
| MEP | C5 | 4.152934 | -7.010556 | yes | no | transformable_candidate |
| MEP | CB | 9.150969 | -2.012521 | yes | no | transformable_candidate |
| MEP | CC | 8.945009 | -2.218481 | yes | no | transformable_candidate |
| MEP | HB2 | 8.740132 | -2.423357 | yes | no | transformable_candidate |
| MEP | HB3 | 8.740196 | -2.423293 | yes | no | transformable_candidate |
| MEP | HC2 | 9.448975 | -1.714515 | yes | no | transformable_candidate |
| MEP | HC3 | 9.447838 | -1.715651 | yes | no | transformable_candidate |
| MEP | HX11 | 7.389594 | -3.773895 | yes | no | transformable_candidate |
| MEP | HX21 | 7.378739 | -3.784751 | yes | no | transformable_candidate |
| MEP | HX22 | 8.507600 | -2.655890 | yes | no | transformable_candidate |
| MEP | HX31 | 2.288693 | -8.874797 | yes | no | transformable_candidate |
| MEP | HX32 | 2.299018 | -8.864472 | yes | no | transformable_candidate |
| MEP | N2 | 7.001677 | -4.161813 | yes | no | transformable_candidate |
| MEP | N4 | 5.016445 | -6.147044 | yes | no | transformable_candidate |
| MEP | N6 | 5.022769 | -6.140721 | yes | no | transformable_candidate |
| MEP | NX1 | 7.501051 | -3.662438 | yes | no | transformable_candidate |
| MEP | NX2 | 7.492897 | -3.670593 | yes | no | transformable_candidate |
| MEP | NX3 | 2.633146 | -8.530343 | yes | no | transformable_candidate |

## Ambiguous/review atoms

| residue_name | atom_name | classification | mean_radial_distance_A | mean_delta_from_residue_backbone_A |
| --- | --- | --- | --- | --- |
| CYP | H | candidate_sidechain_or_arm | 10.563605 | -0.599872 |
| CYP | H1 | candidate_sidechain_or_arm | 11.889023 | 0.725546 |
| CYP | H2 | candidate_sidechain_or_arm | 10.278712 | -0.884765 |
| CYP | H3 | candidate_sidechain_or_arm | 11.324047 | 0.160571 |
| CYP | HA | candidate_sidechain_or_arm | 11.070345 | -0.093132 |
| MEP | H | candidate_sidechain_or_arm | 10.563692 | -0.599798 |
| MEP | H1 | candidate_sidechain_or_arm | 11.889195 | 0.725705 |
| MEP | H2 | candidate_sidechain_or_arm | 11.315881 | 0.152392 |
| MEP | H3 | candidate_sidechain_or_arm | 10.281988 | -0.881501 |
| MEP | HA | candidate_sidechain_or_arm | 11.070333 | -0.093156 |

## Recommended operational base-length definition

For the later base-length sweep, keep all GLU atoms and all N/CA/C/O/OXT backbone-like atoms fixed. For CYP and MEP, treat non-backbone atoms flagged as `transformable_candidate` as candidate arm atoms whose radial component from the fitted central axis may be scaled. In the current six-chain geometry these candidates are axis-facing rather than outward relative to the CYP/MEP backbone-like radial mean. Leave atoms marked `review` fixed until a visual or chemistry-aware inspection justifies including them.

The base length should be measured as the radial distance from the fitted central axis to the selected transformable CYP/MEP atom group, summarized per residue and chain. Scaling should preserve residue identity and should be reported as a geometric perturbation only.

## Output files

- `outputs/metrics/hexad_arm_atom_inventory.csv`
- `outputs/metrics/hexad_arm_geometry_summary.csv`
