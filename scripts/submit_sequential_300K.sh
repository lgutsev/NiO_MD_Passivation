#!/bin/bash

set -euo pipefail

cd "$(dirname "$0")/.."

primary="prepared/me-4pacz-alone/held-300K.data"
if [[ ! -s "${primary}" ]]; then
    echo "ERROR: ${primary} is missing; finish the Me-4PACz 300 K hold first." >&2
    exit 1
fi

submit_job() {
    local output job_id

    if ! output="$(sbatch --parsable "$@" 2>&1)"; then
        printf '%s\n' "${output}" >&2
        return 1
    fi

    # QBD's Lua job-submit plugin writes allocation messages alongside the
    # parsable result. Preserve those messages for the user, but return only
    # the numeric Slurm job ID to the dependency chain.
    printf '%s\n' "${output}" >&2
    job_id="$(
        printf '%s\n' "${output}" |
            awk '
                /^[0-9]+([.;][^[:space:]]+)?$/ {
                    id = $0
                    sub(/[.;].*$/, "", id)
                }
                /Submitted job [0-9]+/ {
                    for (i = 1; i <= NF; i++) {
                        if ($i == "job" && $(i + 1) ~ /^[0-9]+$/) {
                            id = $(i + 1)
                        }
                    }
                }
                END {
                    if (id != "") print id
                }
            '
    )"

    if [[ ! "${job_id}" =~ ^[0-9]+$ ]]; then
        echo "ERROR: could not extract a numeric Slurm job ID from sbatch output." >&2
        return 1
    fi

    printf '%s\n' "${job_id}"
}

build_job="$(submit_job scripts/build_sequential_stage2.sbatch)"

deposit_job="$(submit_job --dependency="afterok:${build_job}" scripts/run_sequential_deposition_array.sbatch)"

hold_job="$(submit_job --dependency="afterok:${deposit_job}" scripts/run_sequential_hold_array.sbatch)"

echo "Sequential build job: ${build_job}"
echo "Sequential deposition array: ${deposit_job} (afterok:${build_job})"
echo "Sequential 300 K hold array: ${hold_job} (afterok:${deposit_job})"
