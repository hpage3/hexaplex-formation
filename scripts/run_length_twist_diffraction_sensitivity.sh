#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

python_bin="${PYTHON_BIN:-python3}"
regenerate_length=0
analysis_mode="${MINI_HEXAPLEX_ANALYSIS_MODE:-auto}"
unit_counts="4,7,8,12"

usage() {
  cat <<'EOF'
Usage: bash scripts/run_length_twist_diffraction_sensitivity.sh [--regenerate-length] [--analysis-mode auto|fiber|debye]

Builds the controlled length/twist sensitivity summary from the current radial
profiles. With --regenerate-length, first runs the existing mini-hexaplex length
workflow for 4,7,8,12 units plus the full-length baseline.

This script does not generate non-30-degree twist structures unless the official
Proto-Nucleic Acids Building / PNAB workflow is provided separately.
EOF
}

while [ "$#" -gt 0 ]; do
  case "$1" in
    --regenerate-length)
      regenerate_length=1
      shift
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

case "$analysis_mode" in
  auto|fiber|debye) ;;
  *)
    echo "--analysis-mode must be one of: auto, fiber, debye" >&2
    exit 2
    ;;
esac

if [ "$regenerate_length" -eq 1 ]; then
  MINI_HEXAPLEX_ANALYSIS_MODE="$analysis_mode" "$python_bin" -V >/dev/null
  MINI_HEXAPLEX_ANALYSIS_MODE="$analysis_mode" bash scripts/run_mini_hexaplex_analysis.sh \
    --unit-counts "$unit_counts" \
    --variants central4_units,central7_units,central8_units,central12_units \
    --analysis-mode "$analysis_mode"
fi

"$python_bin" scripts/analyze_length_twist_diffraction_sensitivity.py

echo "Length/twist sensitivity outputs:"
echo "- outputs/metrics/length_twist_model_manifest.csv"
echo "- outputs/metrics/length_twist_feature_summary.csv"
echo "- outputs/metrics/length_twist_peak_widths.csv"
echo "- outputs/reports/length_twist_diffraction_sensitivity_report.md"
echo "- scripts/run_pnab_twist_builder_template.sh"
