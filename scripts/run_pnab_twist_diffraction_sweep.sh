#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

python_bin="${PYTHON_BIN:-python3}"
fiber_dir="${FIBER_DIFFRACTION_DIR:-../fiber-diffraction}"
baseline_yaml=""
twists="24,26,28,30,32,34,36"
dry_run=0

usage() {
  cat <<'EOF'
Usage: bash scripts/run_pnab_twist_diffraction_sweep.sh --baseline-yaml PATH [--twists CSV] [--dry-run]

Generates pNAB h_twist variants from an explicit baseline YAML. If pNAB and
the sibling diffraction tools are available, selected pNAB structures can then
be passed through the existing detector simulation and radial averaging path.

This script intentionally requires --baseline-yaml; it does not infer chemistry
or fabricate pNAB inputs from the current PDB.
EOF
}

while [ "$#" -gt 0 ]; do
  case "$1" in
    --baseline-yaml)
      baseline_yaml="${2:-}"
      shift 2
      ;;
    --twists)
      twists="${2:-}"
      shift 2
      ;;
    --dry-run)
      dry_run=1
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

if [ -z "$baseline_yaml" ]; then
  echo "--baseline-yaml is required." >&2
  usage >&2
  exit 2
fi

generator_args=(scripts/generate_pnab_twist_variants.py --baseline-yaml "$baseline_yaml" --twists "$twists")
if [ "$dry_run" -eq 1 ]; then
  generator_args+=(--dry-run)
fi
"$python_bin" "${generator_args[@]}"

if [ "$dry_run" -eq 1 ]; then
  "$python_bin" scripts/report_pnab_twist_integration.py --baseline-yaml "$baseline_yaml"
  exit 0
fi

if [ ! -f "$fiber_dir/run_powder_benchmark.py" ] || [ ! -f "$fiber_dir/radial_average.py" ] || [ ! -f "$fiber_dir/pdb_to_xyz.py" ]; then
  echo "pNAB structures were generated or attempted, but sibling diffraction tools were not found in $fiber_dir." >&2
  echo "Set FIBER_DIFFRACTION_DIR to the directory containing run_powder_benchmark.py and radial_average.py." >&2
  "$python_bin" scripts/report_pnab_twist_integration.py --baseline-yaml "$baseline_yaml"
  exit 1
fi

diffraction_dir="outputs/pnab_twist_sweep/diffraction"
radial_dir="outputs/pnab_twist_sweep/radial_profiles"
mkdir -p "$diffraction_dir" "$radial_dir"

grid_size="${GRID_SIZE:-31}"
grid_limit="${GRID_LIMIT:-100.0}"
theta_count="${THETA_COUNT:-2}"
phi_count="${PHI_COUNT:-6}"
psi_count="${PSI_COUNT:-2}"
theta_max="${THETA_MAX:-180.0}"
workers="${WORKERS:-1}"
radial_bins="${RADIAL_BINS:-240}"

"$python_bin" - <<'PY' | while IFS=, read -r model_id structure_path; do
import csv
from pathlib import Path
manifest = Path("outputs/metrics/pnab_twist_model_manifest.csv")
if not manifest.exists():
    raise SystemExit(0)
for row in csv.DictReader(manifest.open()):
    if row.get("model_status") == "selected_structure" and row.get("output_structure"):
        print(row["model_id"] + "," + row["output_structure"])
PY
  xyz_path="${diffraction_dir}/${model_id}.xyz"
  output_prefix="${diffraction_dir}/${model_id}"
  radial_prefix="${radial_dir}/${model_id}_radial"
  "$python_bin" "$fiber_dir/pdb_to_xyz.py" --input-pdb "$structure_path" --output-xyz "$xyz_path" --include-hetatm --keep-hydrogens
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
done

"$python_bin" scripts/report_pnab_twist_integration.py --baseline-yaml "$baseline_yaml"
