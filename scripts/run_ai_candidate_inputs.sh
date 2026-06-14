#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

python3 scripts/build_ai_candidate_inputs.py

echo "Candidate PDB: outputs/intermediates/ai_candidate_inputs/full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain.pdb"
echo "Residue table: outputs/intermediates/ai_candidate_inputs/full_hexaplex_anti_parallel_30deg_ideal_chain_residue_table.csv"
echo "Chain summary: outputs/intermediates/ai_candidate_inputs/full_hexaplex_anti_parallel_30deg_ideal_chain_pattern_summary.csv"
echo "Proxy FASTA: outputs/intermediates/ai_candidate_inputs/full_hexaplex_anti_parallel_30deg_ideal_alphafold_esm_proxy.fasta"
echo "Report: outputs/reports/ai_candidate_inputs_report.md"
