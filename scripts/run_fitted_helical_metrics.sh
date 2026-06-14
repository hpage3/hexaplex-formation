#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

python3 scripts/normalize_structures.py

if [ ! -d outputs/intermediates/ladder_structures ]; then
  python3 scripts/build_strand_map_candidate.py
  python3 scripts/build_intermediate_ladder.py
fi

out_dir="outputs/metrics/fitted_helical_order"
mkdir -p "$out_dir"

shopt -s nullglob
for pdb_path in outputs/intermediates/normalized_structures/*_heavy_deduped.pdb; do
  base_name="$(basename "$pdb_path" .pdb)"
  python3 scripts/helical_order_fitted_axis.py \
    --pdb "$pdb_path" \
    --out "${out_dir}/${base_name}_fitted_helical_order.csv" \
    --heavy-only
done

for pdb_path in outputs/intermediates/ladder_structures/*_heavy_deduped.pdb; do
  base_name="$(basename "$pdb_path" .pdb)"
  python3 scripts/helical_order_fitted_axis.py \
    --pdb "$pdb_path" \
    --out "${out_dir}/${base_name}_fitted_helical_order.csv" \
    --heavy-only
done

python3 scripts/compare_fitted_helical_metrics.py
