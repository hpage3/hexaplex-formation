#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

if [ ! -d outputs/intermediates/ladder_structures ]; then
  bash scripts/run_intermediate_ladder_metrics.sh
fi

SAMPLE_ATOMS="${SAMPLE_ATOMS:-}"
MAX_ATOMS="${MAX_ATOMS:-}"

profile_dir="outputs/metrics/ladder_diffraction/profiles"
window_dir="outputs/metrics/ladder_diffraction/window_scores"
mkdir -p "$profile_dir" "$window_dir"

shopt -s nullglob
for pdb_path in outputs/intermediates/ladder_structures/*_heavy_deduped.pdb; do
  base_name="$(basename "$pdb_path" .pdb)"
  profile_out="${profile_dir}/${base_name}_debye_profile.csv"
  score_out="${window_dir}/${base_name}_window_scores.csv"
  args=(
    --pdb "$pdb_path"
    --out "$profile_out"
    --heavy-only
    --method histogram
  )
  if [ -n "$SAMPLE_ATOMS" ]; then
    args+=(--sample-atoms "$SAMPLE_ATOMS")
  elif [ -n "$MAX_ATOMS" ]; then
    echo "WARNING: MAX_ATOMS is a legacy first-N truncation fallback and may bias structures; prefer SAMPLE_ATOMS." >&2
    args+=(--max-atoms "$MAX_ATOMS")
  fi
  python3 scripts/debye_radial_profile.py "${args[@]}"
  python3 scripts/score_radial_windows.py \
    --profile "$profile_out" \
    --out "$score_out"
done

python3 scripts/compare_ladder_diffraction.py
