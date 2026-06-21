# Nick Ideal 16-mer Atom-Contribution Comparison

## Purpose

Reproduce Nick's four-trace atom-contribution figure using the corrected Asem diffraction path and the Emory-corrected experimental profile.

## Inputs and Derived Model

- Corrected experimental profile: `inputs\experimental\nick_powder_profile_corrected_emory.csv`
- Full ideal 16-mer radial profile: `outputs\nick_ideal_16mer_corrected_emory_profile\radial_profiles\ideal_16mer_antiparallel_30deg_radial_profile.csv`
- Derived no-CH2-COOH PDB: `inputs\nick_ideal_models\derived\Hexaplex_AntiParallel_30deg_Ideal_no_CH2_COOH.pdb`
- Derived no-CH2-COOH XYZ: `inputs\nick_ideal_models\derived\Hexaplex_AntiParallel_30deg_Ideal_no_CH2_COOH.xyz`
- Side-chain subtraction profile: `outputs\nick_ideal_16mer_sidechain_subtraction_corrected_emory_profile\radial_profiles\sidechain_subtraction_profile.csv`
- Final four-trace plot: `outputs\nick_ideal_16mer_corrected_emory_profile\plots\nick_ideal_16mer_four_trace_nick_style.png`

## Atom-Selection Rule

Remove atoms named CG, CD, OE1, OE2, HG2, and HG3 from every GLU residue. This removes the terminal GLU CH2-COOH group while preserving backbone atoms and CB. Hydrogens are excluded from the diffraction XYZ.

- Source PDB atoms: 7146
- Derived PDB atoms: 6066
- Removed atoms: 1080
- Derived heavy deduped XYZ atoms: 1806
- Removed element counts: `{"C": 360, "H": 360, "O": 360}`
- Remaining derived PDB element counts: `{"C": 1800, "H": 2454, "N": 1170, "O": 642}`

## Corrected Diffraction Settings

- Engine: `reference/asem_corrected_diffraction_engine/`.
- Asem correction: azimuthal rotations are applied to independent tilted coordinate stacks rather than accumulated through a rotation loop.
- Tilts: `[0]`.
- Rotations: `range(0, 181, 5)`, 37 rotations from 0 to 180 degrees.
- Hydrogens excluded from XYZ; exact heavy-atom deduplication applied.
- Grid: 129 x 129; detector half-width 100 mm; detector distance 338.4 mm; wavelength 0.7749 A; radial bins 420.

## Quantitative Feature Attribution

| feature_window | experimental_peak_d_A | full_ideal_peak_d_A | no_ch2_cooh_peak_d_A | subtraction_peak_d_A | attribution_note |
| --- | --- | --- | --- | --- | --- |
| 7.25 A | 7.32457786116 | 7.34333958724 | 7.34333958724 | 7.48405253283 | likely backbone-dominated |
| 5.6 A | 5.57661038149 | 5.698561601 | 5.50469043152 | 5.698561601 | mixed backbone/side-chain |
| 4.4 A | 4.40400250156 | 4.34459036898 | 4.46341463415 | 4.4133833646 | mixed backbone/side-chain |
| 3.7 A | 3.794246404 | 3.794246404 | 3.794246404 | 3.794246404 | mixed backbone/side-chain |
| 3.4 A | 3.38148843027 | 3.41588492808 | 3.41588492808 | 3.4440275172 | primarily bases |

## Interpretation

- The 3.4 A window remains most consistent with stacked-base structure because it is strong in the full and no-CH2-COOH traces and is not uniquely isolated by the subtraction trace.
- The 7.25 A window remains likely backbone-dominated because removing the terminal GLU CH2-COOH group does not create a uniquely dominant subtraction feature there.
- The 3.7 A, 4.4 A, and 5.6 A windows are best treated as mixed contribution windows in this corrected comparison.

## Visualization Notes

The four-trace figure uses d-space Gaussian smoothing with sigma 0.10 A for visual comparison only. The subtraction trace is full ideal minus no-CH2-COOH after interpolation onto a common d-spacing grid. Smoothing and stacked offsets do not replace the quantitative peak-position table.
