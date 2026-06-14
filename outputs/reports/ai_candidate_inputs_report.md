# AlphaFold/ESM Candidate Inputs

## Scientific cautions

- CYP and MEP are nonstandard residues.
- The proxy FASTA uses `X` for nonstandard residues and should be treated as an exploratory placeholder only.
- These inputs preserve residue identity in CSV/PDB form, but they do not make the structure a biologically standard AlphaFold/ESM sequence.
- This workflow validates chain inference and residue-pattern extraction; it does not prove the true assembly pathway or biological fold.

## Input summary

- Source PDB: inputs/structures/full_hexaplex_anti_parallel_30deg_ideal.pdb
- Deduplicated inferred-chain PDB: outputs/intermediates/ai_candidate_inputs/full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain.pdb
- Residues after exact deduplication: 180
- Nonstandard residues retained: 90
- Inferred chains: 6

## Chain patterns

| chain_id | residue_count | nonstandard_residue_count | first_residue_label | last_residue_label | residue_pattern | proxy_sequence |
| --- | --- | --- | --- | --- | --- | --- |
| A | 30 | 15 | A:CYP1 | A:GLU30 | CYP-GLU-CYP-GLU-CYP-GLU-CYP-GLU-CYP-GLU-CYP-GLU-CYP-GLU-CYP-GLU-CYP-GLU-CYP-GLU-CYP-GLU-CYP-GLU-CYP-GLU-CYP-GLU-CYP-GLU | XEXEXEXEXEXEXEXEXEXEXEXEXEXEXE |
| B | 30 | 15 | B:MEP31 | B:GLU60 | MEP-GLU-MEP-GLU-MEP-GLU-MEP-GLU-MEP-GLU-MEP-GLU-MEP-GLU-MEP-GLU-MEP-GLU-MEP-GLU-MEP-GLU-MEP-GLU-MEP-GLU-MEP-GLU-MEP-GLU | XEXEXEXEXEXEXEXEXEXEXEXEXEXEXE |
| C | 30 | 15 | C:CYP61 | C:GLU90 | CYP-GLU-CYP-GLU-CYP-GLU-CYP-GLU-CYP-GLU-CYP-GLU-CYP-GLU-CYP-GLU-CYP-GLU-CYP-GLU-CYP-GLU-CYP-GLU-CYP-GLU-CYP-GLU-CYP-GLU | XEXEXEXEXEXEXEXEXEXEXEXEXEXEXE |
| D | 30 | 15 | D:MEP91 | D:GLU120 | MEP-GLU-MEP-GLU-MEP-GLU-MEP-GLU-MEP-GLU-MEP-GLU-MEP-GLU-MEP-GLU-MEP-GLU-MEP-GLU-MEP-GLU-MEP-GLU-MEP-GLU-MEP-GLU-MEP-GLU | XEXEXEXEXEXEXEXEXEXEXEXEXEXEXE |
| E | 30 | 15 | E:CYP121 | E:GLU150 | CYP-GLU-CYP-GLU-CYP-GLU-CYP-GLU-CYP-GLU-CYP-GLU-CYP-GLU-CYP-GLU-CYP-GLU-CYP-GLU-CYP-GLU-CYP-GLU-CYP-GLU-CYP-GLU-CYP-GLU | XEXEXEXEXEXEXEXEXEXEXEXEXEXEXE |
| F | 30 | 15 | F:MEP151 | F:GLU180 | MEP-GLU-MEP-GLU-MEP-GLU-MEP-GLU-MEP-GLU-MEP-GLU-MEP-GLU-MEP-GLU-MEP-GLU-MEP-GLU-MEP-GLU-MEP-GLU-MEP-GLU-MEP-GLU-MEP-GLU | XEXEXEXEXEXEXEXEXEXEXEXEXEXEXE |

## Usage warning

Use these outputs only as candidate-only exploratory use for AlphaFold/ESM-oriented exploration. The nonstandard residues remain explicit in the CSV and PDB outputs, and the proxy FASTA is intentionally lossy. This does not make the structure a biologically standard AlphaFold/ESM sequence.
