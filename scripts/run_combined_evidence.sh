#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

if [ ! -f outputs/metrics/intermediate_ladder/intermediate_ladder_comparison.csv ]; then
  bash scripts/run_intermediate_ladder_metrics.sh
fi

if [ ! -f outputs/metrics/fitted_helical_order/fitted_helical_comparison.csv ]; then
  bash scripts/run_fitted_helical_metrics.sh
fi

if ! compgen -G "outputs/metrics/block_contacts/*_block_contact_summary.csv" > /dev/null; then
  bash scripts/run_block_contact_analysis.sh
fi

if [ ! -f outputs/metrics/ladder_diffraction/ladder_diffraction_comparison.csv ]; then
  bash scripts/run_ladder_diffraction.sh
fi

python3 scripts/combined_evidence_table.py

echo "Combined evidence CSV: outputs/metrics/combined_evidence_table.csv"
echo "Combined evidence report: outputs/reports/combined_evidence_report.md"
