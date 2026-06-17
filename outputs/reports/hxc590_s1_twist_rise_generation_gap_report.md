# HXC590 S1 Twist/Rise Generation Gap Audit

## Purpose

This audit checks whether the repository can safely generate the missing HXC590 S1 twist and rise sensitivity candidates before rerunning the powder falsification screen.

This remains a falsification-style screen, not a definitive phase assignment.

Synthetic twist/rise variants are controls for diffraction sensitivity, not chemically optimized structures.

## Generator inspection

- Existing pNAB twist helper: `scripts/generate_pnab_twist_variants.py`.
- pNAB Python import status: not available.
- Current-model pNAB baseline YAML candidates found: 0.
- Existing length/twist sensitivity workflow keeps non-30-degree full-length twist rows as pending placeholders when the official builder inputs are absent.
- No audited rise-generation workflow or stack-axis rise-transform workflow was found for the requested rise scan.

## Twist candidate status

| candidate_id | parameter_value | status | reason |
|---|---|---|---|
| full_length_twist_24 | 24 deg | unavailable | No coordinate/profile files were found for this twist. The repo contains a pNAB twist helper, but no current-model baseline YAML with HelicalParameters.h_twist and CYP content was found. Also, python cannot import pnab. |
| full_length_twist_26 | 26 deg | unavailable | No coordinate/profile files were found for this twist. The repo contains a pNAB twist helper, but no current-model baseline YAML with HelicalParameters.h_twist and CYP content was found. Also, python cannot import pnab. |
| full_length_twist_28 | 28 deg | unavailable | No coordinate/profile files were found for this twist. The repo contains a pNAB twist helper, but no current-model baseline YAML with HelicalParameters.h_twist and CYP content was found. Also, python cannot import pnab. |
| full_length_twist_30 | 30 deg | located | Coordinate and radial-profile files are present. |
| full_length_twist_32 | 32 deg | unavailable | No coordinate/profile files were found for this twist. The repo contains a pNAB twist helper, but no current-model baseline YAML with HelicalParameters.h_twist and CYP content was found. Also, python cannot import pnab. |
| full_length_twist_34 | 34 deg | unavailable | No coordinate/profile files were found for this twist. The repo contains a pNAB twist helper, but no current-model baseline YAML with HelicalParameters.h_twist and CYP content was found. Also, python cannot import pnab. |
| full_length_twist_36 | 36 deg | unavailable | No coordinate/profile files were found for this twist. The repo contains a pNAB twist helper, but no current-model baseline YAML with HelicalParameters.h_twist and CYP content was found. Also, python cannot import pnab. |

## Rise candidate status

| candidate_id | parameter_value | status | reason |
|---|---|---|---|
| rise_3p20_synthetic_control | 3.20 A | unavailable | No coordinate/profile files were found for this rise. No safe existing rise-generation workflow or audited stack-axis transform was found in this repo. |
| rise_3p30_synthetic_control | 3.30 A | unavailable | No coordinate/profile files were found for this rise. No safe existing rise-generation workflow or audited stack-axis transform was found in this repo. |
| rise_3p35_synthetic_control | 3.35 A | unavailable | No coordinate/profile files were found for this rise. No safe existing rise-generation workflow or audited stack-axis transform was found in this repo. |
| rise_3p40_synthetic_control | 3.40 A | unavailable | No coordinate/profile files were found for this rise. No safe existing rise-generation workflow or audited stack-axis transform was found in this repo. |
| rise_3p50_synthetic_control | 3.50 A | unavailable | No coordinate/profile files were found for this rise. No safe existing rise-generation workflow or audited stack-axis transform was found in this repo. |
| rise_3p60_synthetic_control | 3.60 A | unavailable | No coordinate/profile files were found for this rise. No safe existing rise-generation workflow or audited stack-axis transform was found in this repo. |

## Interpretation

The full-length 30-degree baseline is available, but the missing non-30-degree twist variants remain unavailable in this checkout because the current-model pNAB baseline YAML and matching runtime inputs are absent.

Rise variants remain unavailable because the repository does not contain a safe, audited rise-generation path for the current candidate model.

If nearby twist or rise variants survive under current tolerances, the current powder peak list supports the conformation family but does not uniquely determine those parameters.

Because nearby twist and rise variants are not currently generated or profiled, they cannot be used to refine or falsify the twist/rise parameters in this pass.

## Outputs

- `outputs\metrics\hxc590_s1_twist_rise_sensitivity.csv`
- `outputs\metrics\hxc590_s1_twist_rise_generation_gap.csv`
- `outputs\reports\hxc590_s1_twist_rise_generation_gap_report.md`
