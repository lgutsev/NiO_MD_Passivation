#!/bin/bash

set -euo pipefail

cd "$(dirname "$0")/.."

primary="prepared/me-4pacz-alone/held-300K.data"
if [[ ! -s "${primary}" ]]; then
    echo "ERROR: ${primary} is missing; finish the Me-4PACz 300 K hold first." >&2
    exit 1
fi

build_job="$(sbatch --parsable scripts/build_sequential_stage2.sbatch)"
build_job="${build_job%%;*}"

deposit_job="$(sbatch --parsable --dependency="afterok:${build_job}" scripts/run_sequential_deposition_array.sbatch)"
deposit_job="${deposit_job%%;*}"

hold_job="$(sbatch --parsable --dependency="afterok:${deposit_job}" scripts/run_sequential_hold_array.sbatch)"
hold_job="${hold_job%%;*}"

echo "Sequential build job: ${build_job}"
echo "Sequential deposition array: ${deposit_job} (afterok:${build_job})"
echo "Sequential 300 K hold array: ${hold_job} (afterok:${deposit_job})"
