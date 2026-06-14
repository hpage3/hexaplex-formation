#!/usr/bin/env bash
set -euo pipefail

cat >&2 <<'EOF'
The official Proto-Nucleic Acids Building / PNAB program was not found locally.

Provide the official builder and baseline input parameters, then replace the
PNAB_BUILDER placeholder below with the actual executable and option names.

Required inputs:
- current 30 degree six-strand baseline parameter file from the original builder
- sequence/residue specification matching full_hexaplex_anti_parallel_30deg_ideal
- strand count fixed at six
- helical twist values, for example: 24,26,28,30,32,34,36
- output PDB path for each generated twist variant

Proposed command template:
  PNAB_BUILDER --input baseline_builder_parameters.json \
    --strand-count 6 \
    --helical-twist-deg ${TWIST_DEG} \
    --output-pdb outputs/length_twist_diffraction/structures/full_length_twist_${TWIST_DEG}.pdb

This stub intentionally exits without generating coordinates.
EOF
exit 1
