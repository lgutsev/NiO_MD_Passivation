#!/bin/bash
#SBATCH -p single
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 04:00:00
#SBATCH -A loni_perovsk27
#SBATCH -J nio.access
#SBATCH -o nio.access.%A_%a.out
#
# Post-processing only: reads existing hold/relax trajectories, launches no MD.
#
# One array task per run directory listed (one per line, relative to the
# repository root) in RUN_LIST, e.g. produced by:
#   find prepared prepared-rerun prepared-rerun2 prepared-lego prepared-lego2 \
#     -mindepth 1 -maxdepth 1 -type d -exec test -s {}/relax-300K.lammpstrj \; -print \
#     | sort > accessibility_runs.txt
#   sbatch --array=1-$(wc -l < accessibility_runs.txt) \
#     --export=ALL,RUN_LIST=accessibility_runs.txt scripts/analyze_accessibility_array.sbatch
#
# After all tasks finish, build the cross-run comparison on a login node:
#   python scripts/collect_accessibility_manifest.py --output accessibility_manifest.csv
#   nio-md-prep compare-routes --workbooks . \
#     --accessibility-manifest accessibility_manifest.csv --output route-comparison

set -euo pipefail

module purge
module load python/3.11.5-anaconda
cd "${SLURM_SUBMIT_DIR}"
source .venv/bin/activate

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MPLCONFIGDIR="${SLURM_TMPDIR:-/tmp}/nio-access-mpl-${SLURM_JOB_ID}"
mkdir -p "${MPLCONFIGDIR}"

run_list="${RUN_LIST:?set RUN_LIST to a file of run directories}"
run_directory="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${run_list}")"
if [[ -z "${run_directory}" || ! -d "${run_directory}" ]]; then
    echo "ERROR: no run directory for task ${SLURM_ARRAY_TASK_ID} in ${run_list}" >&2
    exit 1
fi

# Default: every frame of each 500 ps hold (dump every 2000 steps = 1 ps), so
# persistence and correlation times are measured over the whole window.
last_frames="${LAST_FRAMES:-0}"
stages="${STAGES:-hold-300K relax-300K hold-400K relax-400K}"

echo "Commit: $(git rev-parse --short HEAD)"
echo "Run: ${run_directory}  LAST_FRAMES=${last_frames}  STAGES=${stages}"

failed=0
for stage in ${stages}; do
    trajectory="${run_directory}/${stage}.lammpstrj"
    if [[ ! -s "${trajectory}" ]]; then
        echo "SKIP: ${trajectory} missing"
        continue
    fi
    SECONDS=0
    if nio-md-prep analyze-accessibility "${run_directory}" \
        --trajectory "${stage}.lammpstrj" \
        --output "${run_directory}/accessibility-analysis-${stage}" \
        --last-frames "${last_frames}" \
        --timestep-fs 0.5; then
        echo "COMPLETED: ${run_directory}/${stage} in ${SECONDS} s"
    else
        echo "FAILED: ${run_directory}/${stage}" >&2
        failed=$((failed + 1))
    fi
done
exit "${failed}"
