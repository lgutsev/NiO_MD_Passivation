#!/bin/bash
# Submit the hydroxylation-series pilots (3 seeds each via the array in run_pilot_array.sbatch).
# Build first:  STUDIES="$(grep -v '^#' studies/ito/hydroxylation-series.txt | tr '\n' ' ')" sbatch scripts/ito/build_pilots.sbatch
# Then:         bash scripts/ito/submit_hydroxylation_series.sh            # everything (see pilot-plan.md for cost)
#               FILTER=me-4pacz bash scripts/ito/submit_hydroxylation_series.sh   # a subset by substring
# Nothing is submitted when DRY_RUN=1.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
for study in $(grep -v '^#' studies/ito/hydroxylation-series.txt); do
    [[ -n "${FILTER:-}" && "${study}" != *"${FILTER}"* ]] && continue
    if [[ "${DRY_RUN:-0}" == 1 ]]; then echo "would submit ${study}"; continue; fi
    STUDY="${study}" sbatch --export=ALL,STUDY="${study}" -J "ito.${study}" scripts/ito/run_pilot_array.sbatch
done
