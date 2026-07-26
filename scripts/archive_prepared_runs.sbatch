#!/bin/bash
#SBATCH -p single
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 04:00:00
#SBATCH -A loni_perovsk27
#SBATCH -J nio.archive
#SBATCH -o nio.archive.%j.out

set -euo pipefail

module purge
module load python/3.11.5-anaconda

cd "${SLURM_SUBMIT_DIR}"
source .venv/bin/activate

if ! command -v nio-md-prep >/dev/null 2>&1; then
    echo "ERROR: nio-md-prep is unavailable after activating .venv." >&2
    exit 1
fi

if [[ -n "${ARCHIVE_ROOTS:-}" ]]; then
    read -r -a roots <<< "${ARCHIVE_ROOTS}"
else
    roots=()
    for candidate in prepared*; do
        if [[ -d "${candidate}" ]]; then
            roots+=("${candidate}")
        fi
    done
fi

if (( ${#roots[@]} == 0 )); then
    echo "ERROR: no prepared/work roots were found." >&2
    echo "Set ARCHIVE_ROOTS to a space-separated list of directories." >&2
    exit 1
fi

archive_output="${ARCHIVE_OUTPUT:-archives/nio-md-rerun-files-$(date -u +%Y%m%dT%H%M%SZ).zip}"

echo "Archiving roots:"
printf '  %s\n' "${roots[@]}"
echo "Output: ${archive_output}"

nio-md-prep archive-runs "${roots[@]}" --output "${archive_output}"

echo "Archive complete:"
ls -lh "${archive_output}"
