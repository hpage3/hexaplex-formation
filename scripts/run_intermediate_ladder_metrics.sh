#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

python3 scripts/normalize_structures.py
python3 scripts/build_strand_map_candidate.py
python3 scripts/build_intermediate_ladder.py

metrics_dir="outputs/metrics/intermediate_ladder"
mkdir -p "$metrics_dir"

shopt -s nullglob
for pdb_path in outputs/intermediates/ladder_structures/*_heavy_deduped.pdb; do
  base_name="$(basename "$pdb_path" .pdb)"
  contact_out="${metrics_dir}/${base_name}_contacts_4p5A.csv"
  python3 scripts/contact_map.py \
    --pdb "$pdb_path" \
    --out "$contact_out" \
    --cutoff 4.5 \
    --heavy-only
  python3 scripts/glu_motif_summary.py \
    --contact-map "$contact_out" \
    --out "${metrics_dir}/${base_name}_glu_motifs.csv"
  python3 scripts/helical_order.py \
    --pdb "$pdb_path" \
    --out "${metrics_dir}/${base_name}_helical_order.csv" \
    --heavy-only
done

for pdb_path in outputs/intermediates/ladder_structures/*_allatom_deduped.pdb; do
  base_name="$(basename "$pdb_path" .pdb)"
  python3 scripts/hbond_candidates.py \
    --pdb "$pdb_path" \
    --out "${metrics_dir}/${base_name}_hbond_candidates.csv"
done

python3 scripts/compare_intermediate_ladder.py
