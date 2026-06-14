#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

fiber_dir="${FIBER_DIFFRACTION_DIR:-../fiber-diffraction}"
python_bin="${PYTHON_BIN:-python3}"
baseline_pdb="outputs/intermediates/ai_candidate_inputs/full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain.pdb"
variants=""
unit_counts="8"
analysis_mode="${MINI_HEXAPLEX_ANALYSIS_MODE:-auto}"

usage() {
  cat <<'EOF'
Usage: bash scripts/run_mini_hexaplex_analysis.sh [--units N|--unit-counts CSV] [--variants CSV] [--analysis-mode auto|fiber|debye]

Generates N-unit six-strand mini-hexaplex truncation variants, runs radial
diffraction/scattering profiles, scores reference d-spacing windows, and writes
comparison plots/reports.

Environment overrides:
  FIBER_DIFFRACTION_DIR  sibling diffraction tool directory, default ../fiber-diffraction
  MINI_HEXAPLEX_ANALYSIS_MODE auto, fiber, or debye; default auto
  PYTHON_BIN             Python interpreter for repo and sibling scripts, default python3
  GRID_SIZE              detector grid size for fiber mode, default 31
  GRID_LIMIT             detector half-width in mm for fiber mode, default 100.0
  THETA_COUNT            powder theta samples for fiber mode, default 2
  PHI_COUNT              powder phi samples for fiber mode, default 6
  PSI_COUNT              powder psi samples for fiber mode, default 2
  THETA_MAX              powder theta max in degrees for fiber mode, default 180.0
  WORKERS                powder workers for fiber mode, default 1
  RADIAL_BINS            radial bins for fiber mode, default 240
  DEBYE_Q_MIN            Debye fallback q min, default 0.2
  DEBYE_Q_MAX            Debye fallback q max, default 2.5
  DEBYE_Q_STEP           Debye fallback q step, default 0.01
EOF
}

while [ "$#" -gt 0 ]; do
  case "$1" in
    --units)
      unit_counts="${2:-}"
      shift 2
      ;;
    --unit-counts)
      unit_counts="${2:-}"
      shift 2
      ;;
    --variants)
      variants="${2:-}"
      shift 2
      ;;
    --analysis-mode)
      analysis_mode="${2:-}"
      shift 2
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [ ! -f "$baseline_pdb" ]; then
  echo "Missing cleaned six-chain baseline PDB: $baseline_pdb" >&2
  echo "Run bash scripts/run_ai_candidate_inputs.sh first." >&2
  exit 1
fi

case "$analysis_mode" in
  auto|fiber|debye) ;;
  *)
    echo "--analysis-mode must be one of: auto, fiber, debye" >&2
    exit 2
    ;;
esac

variant_unit_count() {
  local variant_id="$1"
  "$python_bin" -c '
import re
import sys
value = sys.argv[1]
match = re.search(r"(?:first|central)([0-9]+)_units$", value)
print(match.group(1) if match else "")
' "$variant_id"
}

append_unit_count() {
  local count="$1"
  if [ -z "$count" ]; then
    return
  fi
  case ",$unit_counts," in
    *,"$count",*) ;;
    *) unit_counts="${unit_counts},${count}" ;;
  esac
}

if [ -n "$variants" ]; then
  IFS=',' read -r -a requested_variant_array <<< "$variants"
  for raw_variant in "${requested_variant_array[@]}"; do
    requested_variant="$(echo "$raw_variant" | xargs)"
    append_unit_count "$(variant_unit_count "$requested_variant")"
  done
fi

if [ -z "$unit_counts" ]; then
  unit_counts="4,5,6,7,8,12"
fi

default_variants_for_counts() {
  local counts_csv="$1"
  local out=""
  IFS=',' read -r -a counts_array <<< "$counts_csv"
  for raw_count in "${counts_array[@]}"; do
    count="$(echo "$raw_count" | xargs)"
    if [ -z "$count" ]; then
      continue
    fi
    for variant in "lower_end_first${count}_units" "central${count}_units"; do
      if [ -z "$out" ]; then
        out="$variant"
      else
        out="${out},${variant}"
      fi
    done
  done
  echo "$out"
}

if [ -z "$variants" ]; then
  variants="$(default_variants_for_counts "$unit_counts")"
fi

"$python_bin" scripts/generate_mini_hexaplex_variants.py --pdb "$baseline_pdb" --unit-counts "$unit_counts"

diffraction_dir="outputs/mini_hexaplex/diffraction"
radial_dir="outputs/mini_hexaplex/radial_profiles"
plot_dir="outputs/mini_hexaplex/plots"
report_dir="outputs/mini_hexaplex/reports"
mkdir -p "$diffraction_dir" "$radial_dir" "$plot_dir" "$report_dir" outputs/metrics outputs/reports
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/hexaplex_matplotlib}"
mkdir -p "$MPLCONFIGDIR"

mode="$analysis_mode"
if [ "$mode" = "auto" ]; then
  if [ -f "$fiber_dir/run_powder_benchmark.py" ] && [ -f "$fiber_dir/radial_average.py" ] && [ -f "$fiber_dir/pdb_to_xyz.py" ] && "$python_bin" -c "import numpy" >/dev/null 2>&1; then
    mode="fiber"
  else
    mode="debye"
  fi
fi

if [ "$mode" = "fiber" ] && ! "$python_bin" -c "import numpy" >/dev/null 2>&1; then
  echo "Fiber mode requires numpy in the Python environment used by ../fiber-diffraction." >&2
  echo "Use --analysis-mode debye, or install the sibling workflow dependencies." >&2
  exit 1
fi

echo "Mini-hexaplex analysis"
echo "Baseline: $baseline_pdb"
echo "Unit counts: $unit_counts"
echo "Variants: $variants"
echo "Analysis mode: $mode"

run_fiber_profile() {
  local variant_id="$1"
  local pdb_path="$2"
  local xyz_path="${diffraction_dir}/${variant_id}.xyz"
  local output_prefix="${diffraction_dir}/${variant_id}"
  local radial_prefix="${radial_dir}/${variant_id}_radial"
  local grid_size="${GRID_SIZE:-31}"
  local grid_limit="${GRID_LIMIT:-100.0}"
  local theta_count="${THETA_COUNT:-2}"
  local phi_count="${PHI_COUNT:-6}"
  local psi_count="${PSI_COUNT:-2}"
  local theta_max="${THETA_MAX:-180.0}"
  local workers="${WORKERS:-1}"
  local radial_bins="${RADIAL_BINS:-240}"

  "$python_bin" "$fiber_dir/pdb_to_xyz.py" \
    --input-pdb "$pdb_path" \
    --output-xyz "$xyz_path" \
    --include-hetatm \
    --keep-hydrogens

  "$python_bin" "$fiber_dir/run_powder_benchmark.py" \
    --coordinate-file "$xyz_path" \
    --grid-size "$grid_size" \
    --grid-limit "$grid_limit" \
    --theta-count "$theta_count" \
    --phi-count "$phi_count" \
    --psi-count "$psi_count" \
    --theta-max "$theta_max" \
    --workers "$workers" \
    --output-prefix "$output_prefix"

  "$python_bin" "$fiber_dir/radial_average.py" \
    --input-npy "${output_prefix}.npy" \
    --grid-limit "$grid_limit" \
    --bins "$radial_bins" \
    --output-prefix "$radial_prefix" \
    --plot-space both \
    --normalize-plot
}

run_debye_profile() {
  local variant_id="$1"
  local pdb_path="$2"
  "$python_bin" scripts/debye_radial_profile.py \
    --pdb "$pdb_path" \
    --out "${radial_dir}/${variant_id}_radial.csv" \
    --q-min "${DEBYE_Q_MIN:-0.2}" \
    --q-max "${DEBYE_Q_MAX:-2.5}" \
    --q-step "${DEBYE_Q_STEP:-0.01}" \
    --method histogram
}

run_profile() {
  local variant_id="$1"
  local pdb_path="$2"
  if [ "$mode" = "fiber" ]; then
    run_fiber_profile "$variant_id" "$pdb_path"
  else
    run_debye_profile "$variant_id" "$pdb_path"
  fi
}

run_profile "full_length_baseline" "$baseline_pdb"

IFS=',' read -r -a variant_array <<< "$variants"
for raw_variant in "${variant_array[@]}"; do
  variant_id="$(echo "$raw_variant" | xargs)"
  if [ -z "$variant_id" ]; then
    continue
  fi
  pdb_path="outputs/mini_hexaplex/structures/mini_hexaplex_${variant_id}.pdb"
  if [ ! -f "$pdb_path" ]; then
    echo "Missing mini variant PDB: $pdb_path" >&2
    exit 1
  fi
  echo "Running profile for $variant_id"
  run_profile "$variant_id" "$pdb_path"
done

"$python_bin" scripts/analyze_mini_hexaplex_features.py \
  --profile-dir "$radial_dir" \
  --manifest outputs/mini_hexaplex/mini_hexaplex_variant_manifest.csv \
  --geometry outputs/metrics/mini_hexaplex_geometry_summary.csv \
  --plot-dir "$plot_dir" \
  --baseline-pdb "$baseline_pdb" \
  --variants "$variants" \
  --analysis-mode "$mode"

cp outputs/reports/mini_hexaplex_length_response_report.md "${report_dir}/mini_hexaplex_length_response_report.md"
cp outputs/reports/mini_hexaplex_length_response_report.md outputs/reports/mini_hexaplex_analysis_report.md

echo "Mini-hexaplex outputs:"
echo "- Structures: outputs/mini_hexaplex/structures"
echo "- Detector/fiber intermediates: $diffraction_dir"
echo "- Radial profiles: $radial_dir"
echo "- Comparison plots: $plot_dir"
echo "- Feature summary: outputs/metrics/mini_hexaplex_feature_summary.csv"
echo "- Report: outputs/reports/mini_hexaplex_length_response_report.md"
