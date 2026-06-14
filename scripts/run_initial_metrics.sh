#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

mkdir -p outputs/metrics
mkdir -p outputs/intermediates/normalized_structures

python3 scripts/normalize_structures.py
python3 scripts/inventory_structures.py \
  --input-dir inputs/structures \
  --out outputs/metrics/structure_inventory.csv
python3 scripts/inventory_structures.py \
  --input-dir outputs/intermediates/normalized_structures \
  --out outputs/metrics/normalized_structure_inventory.csv

shopt -s nullglob
for pdb_path in outputs/intermediates/normalized_structures/*_heavy_deduped.pdb; do
  base_name="$(basename "$pdb_path" .pdb)"
  python3 scripts/pair_distance_histogram.py \
    --pdb "$pdb_path" \
    --out "outputs/metrics/${base_name}_pair_distances.csv"
  python3 scripts/pair_distance_histogram.py \
    --pdb "$pdb_path" \
    --out "outputs/metrics/${base_name}_GLU_pair_distances.csv" \
    --residue-filter GLU
done

for pdb_path in outputs/intermediates/normalized_structures/*_allatom_deduped.pdb; do
  if grep -q " GLU " "$pdb_path"; then
    base_name="$(basename "$pdb_path" .pdb)"
    python3 scripts/pair_distance_histogram.py \
      --pdb "$pdb_path" \
      --out "outputs/metrics/${base_name}_GLU_allatom_pair_distances.csv" \
      --all-atoms \
      --residue-filter GLU
  fi
done
