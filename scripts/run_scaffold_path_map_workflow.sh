#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

if [ ! -f inputs/metadata/strand_map_candidate.csv ]; then
  python3 scripts/build_strand_map_candidate.py
fi

python3 scripts/convert_strand_map_to_scaffold_path_map.py
python3 scripts/generate_manual_scaffold_path_template.py
python3 scripts/generate_pymol_strand_map_helper.py
python3 scripts/validate_scaffold_path_map.py
python3 scripts/compare_scaffold_path_maps.py

echo "Candidate scaffold path map: inputs/metadata/scaffold_path_map_candidate.csv"
echo "Manual map template: inputs/metadata/scaffold_path_map_manual_template.csv"
echo "PyMOL helper script: outputs/reports/pymol_scaffold_path_map_helper.pml"
echo "Validation report: outputs/reports/scaffold_path_map_validation.md"
echo "Comparison report: outputs/reports/scaffold_path_map_comparison.md"
