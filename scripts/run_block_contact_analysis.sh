#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

python3 scripts/normalize_structures.py
python3 scripts/build_strand_map_candidate.py

if [ ! -d outputs/intermediates/ladder_structures ]; then
  python3 scripts/build_intermediate_ladder.py
fi

out_dir="outputs/metrics/block_contacts"
mkdir -p "$out_dir"

analyze_pdb() {
  local pdb_path="$1"
  local base_name
  base_name="$(basename "$pdb_path" .pdb)"
  local decomposition="${out_dir}/${base_name}_block_contact_decomposition.csv"
  local summary="${out_dir}/${base_name}_block_contact_summary.csv"

  python3 scripts/block_contact_decomposition.py \
    --pdb "$pdb_path" \
    --out "$decomposition" \
    --cutoff 4.5
  python3 scripts/summarize_block_contacts.py \
    --decomposition-csv "$decomposition" \
    --out "$summary"
}

analyze_pdb outputs/intermediates/normalized_structures/hexaplex_scaffold_only_complement_heavy_deduped.pdb
analyze_pdb outputs/intermediates/normalized_structures/full_hexaplex_anti_parallel_30deg_ideal_heavy_deduped.pdb

shopt -s nullglob
for pdb_path in outputs/intermediates/ladder_structures/*_heavy_deduped.pdb; do
  analyze_pdb "$pdb_path"
done

python3 scripts/report_block_contact_analysis.py
