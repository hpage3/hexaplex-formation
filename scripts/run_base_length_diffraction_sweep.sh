#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

fiber_dir="${FIBER_DIFFRACTION_DIR:-../fiber-diffraction}"
default_scales="0.85,0.90,0.95,1.00,1.05,1.10,1.15"
scales="$default_scales"
include_stress_test=0

usage() {
  cat <<'EOF'
Usage: bash scripts/run_base_length_diffraction_sweep.sh [--scales CSV] [--include-stress-test]

Runs a powder-style detector simulation and radial averaging workflow for generated
base-length variant PDBs. Scale 1.20 is excluded by default because geometry
sanity checks found suspicious heavy-atom overlaps; --include-stress-test appends
1.20 for flagged stress-test inspection.

Environment overrides:
  FIBER_DIFFRACTION_DIR  sibling diffraction tool directory, default ../fiber-diffraction
  GRID_SIZE              detector grid size, default 31
  GRID_LIMIT             detector half-width in mm, default 100.0
  THETA_COUNT            powder theta samples, default 2
  PHI_COUNT              powder phi samples, default 6
  PSI_COUNT              powder psi samples, default 2
  THETA_MAX              powder theta max in degrees, default 180.0
  WORKERS                powder workers, default 1
  RADIAL_BINS            radial bins, default 240
EOF
}

while [ "$#" -gt 0 ]; do
  case "$1" in
    --scales)
      scales="${2:-}"
      shift 2
      ;;
    --include-stress-test)
      include_stress_test=1
      shift
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

if [ -z "$scales" ]; then
  echo "--scales must not be empty" >&2
  exit 2
fi

if [ "$include_stress_test" -eq 1 ]; then
  case ",$scales," in
    *,1.20,*) ;;
    *) scales="${scales},1.20" ;;
  esac
fi

if [ ! -f "$fiber_dir/run_powder_benchmark.py" ] || [ ! -f "$fiber_dir/radial_average.py" ]; then
  echo "Could not find sibling diffraction tools in $fiber_dir" >&2
  echo "Set FIBER_DIFFRACTION_DIR to the directory containing run_powder_benchmark.py and radial_average.py." >&2
  exit 1
fi

manifest="outputs/base_length_sweep/structures/base_length_variant_manifest.csv"
geometry="outputs/metrics/base_length_variant_geometry.csv"
baseline_pdb="outputs/intermediates/ai_candidate_inputs/full_hexaplex_anti_parallel_30deg_ideal_deduped_6chain.pdb"

if [ ! -f "$manifest" ] || [ ! -f "$geometry" ]; then
  python3 scripts/generate_base_length_variants.py
fi

diffraction_dir="outputs/base_length_sweep/diffraction"
radial_dir="outputs/base_length_sweep/radial_profiles"
plot_dir="outputs/base_length_sweep/plots"
report_dir="outputs/base_length_sweep/reports"
mkdir -p "$diffraction_dir" "$radial_dir" "$plot_dir" "$report_dir" outputs/metrics outputs/reports
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/hexaplex_matplotlib}"
mkdir -p "$MPLCONFIGDIR"

grid_size="${GRID_SIZE:-31}"
grid_limit="${GRID_LIMIT:-100.0}"
theta_count="${THETA_COUNT:-2}"
phi_count="${PHI_COUNT:-6}"
psi_count="${PSI_COUNT:-2}"
theta_max="${THETA_MAX:-180.0}"
workers="${WORKERS:-1}"
radial_bins="${RADIAL_BINS:-240}"

echo "Base-length diffraction sweep"
echo "Scales: $scales"
echo "Diffraction tools: $fiber_dir"
echo "Detector grid: ${grid_size}x${grid_size}; grid limit ${grid_limit} mm"
echo "Powder samples: theta=${theta_count}, phi=${phi_count}, psi=${psi_count}, theta_max=${theta_max}, workers=${workers}"

IFS=',' read -r -a scale_array <<< "$scales"
for raw_scale in "${scale_array[@]}"; do
  scale="$(echo "$raw_scale" | xargs)"
  if [ -z "$scale" ]; then
    continue
  fi
  scale_token="$(python3 -c 'import sys; print(f"{float(sys.argv[1]):.2f}".replace(".", "p"))' "$scale")"
  variant_id="hexaplex_base_length_scale_${scale_token}"
  pdb_path="outputs/base_length_sweep/structures/${variant_id}.pdb"
  xyz_path="${diffraction_dir}/${variant_id}.xyz"
  output_prefix="${diffraction_dir}/${variant_id}"
  radial_prefix="${radial_dir}/${variant_id}_radial"

  if [ ! -f "$pdb_path" ]; then
    echo "Missing variant PDB: $pdb_path" >&2
    exit 1
  fi

  if [ "$scale_token" = "1p20" ]; then
    echo "WARNING: including scale 1.20 stress-test variant with known suspicious overlap warning." >&2
  fi

  echo "Converting $pdb_path to $xyz_path"
  python3 "$fiber_dir/pdb_to_xyz.py" \
    --input-pdb "$pdb_path" \
    --output-xyz "$xyz_path" \
    --include-hetatm \
    --keep-hydrogens

  echo "Running powder detector simulation for $variant_id"
  python3 "$fiber_dir/run_powder_benchmark.py" \
    --coordinate-file "$xyz_path" \
    --grid-size "$grid_size" \
    --grid-limit "$grid_limit" \
    --theta-count "$theta_count" \
    --phi-count "$phi_count" \
    --psi-count "$psi_count" \
    --theta-max "$theta_max" \
    --workers "$workers" \
    --output-prefix "$output_prefix"

  echo "Radial averaging $variant_id"
  python3 "$fiber_dir/radial_average.py" \
    --input-npy "${output_prefix}.npy" \
    --grid-limit "$grid_limit" \
    --bins "$radial_bins" \
    --output-prefix "$radial_prefix" \
    --plot-space both \
    --normalize-plot
done

python3 scripts/analyze_base_length_sweep_features.py \
  --profile-dir "$radial_dir" \
  --manifest "$manifest" \
  --geometry "$geometry" \
  --plot-dir "$plot_dir" \
  --scales "$scales" \
  --baseline-pdb "$baseline_pdb"

echo "Sweep outputs:"
echo "- Detector .npy files: $diffraction_dir"
echo "- Radial profile CSVs: $radial_dir"
echo "- Comparison plots: $plot_dir"
echo "- Feature summary: outputs/metrics/base_length_sweep_feature_summary.csv"
echo "- Report: outputs/reports/base_length_sweep_report.md"
