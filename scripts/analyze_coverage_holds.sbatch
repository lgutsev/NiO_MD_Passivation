#!/bin/bash
#SBATCH -p single
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 12:00:00
#SBATCH -A loni_perovsk27
#SBATCH -J nio.coverage
#SBATCH -o nio.coverage.%j.out

set -euo pipefail

module purge
module load python/3.11.5-anaconda

cd "${SLURM_SUBMIT_DIR}"
source .venv/bin/activate

# Keep NumPy, SciPy, and plotting libraries inside the single allocated CPU.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

export MPLCONFIGDIR="${SLURM_TMPDIR:-/tmp}/nio-coverage-mpl-${SLURM_JOB_ID}"
mkdir -p "${MPLCONFIGDIR}"

if ! command -v nio-md-prep >/dev/null 2>&1; then
    echo "ERROR: nio-md-prep is unavailable after activating .venv." >&2
    exit 1
fi

if ! python -c "import numpy, scipy, matplotlib, openpyxl" >/dev/null 2>&1; then
    echo "ERROR: coverage-analysis dependencies are missing." >&2
    echo "Install once from the repository root with:" >&2
    echo "  python -m pip install -e '.[analysis]'" >&2
    exit 1
fi

prepared_root="${SLURM_SUBMIT_DIR}/prepared"
if [[ ! -d "${prepared_root}" ]]; then
    echo "ERROR: ${prepared_root} does not exist." >&2
    exit 1
fi

# These can be overridden at submission, for example:
# sbatch --export=ALL,LAST_FRAMES=0 scripts/analyze_coverage_holds.sbatch
last_frames="${LAST_FRAMES:-100}"
grid_spacing="${GRID_SPACING:-0.20}"
blocks="${BLOCKS:-5}"

if [[ ! "${last_frames}" =~ ^[0-9]+$ ]]; then
    echo "ERROR: LAST_FRAMES must be a non-negative integer." >&2
    exit 1
fi

echo "Repository: ${SLURM_SUBMIT_DIR}"
echo "Commit: $(git rev-parse --short HEAD)"
echo "Python: $(python --version 2>&1)"
echo "LAST_FRAMES: ${last_frames} (0 means all frames)"
echo "GRID_SPACING: ${grid_spacing} A"
echo "BLOCKS: ${blocks}"

found=0
completed=0
failed=0
total_seconds=0

while IFS= read -r -d '' trajectory; do
    found=$((found + 1))
    run_directory="$(dirname "${trajectory}")"
    trajectory_name="$(basename "${trajectory}")"
    stage="${trajectory_name%.lammpstrj}"
    system_name="$(basename "${run_directory}")"
    output_directory="${run_directory}/coverage-analysis-${stage}"

    echo
    echo "ANALYZING: ${system_name}/${trajectory_name}"
    echo "OUTPUT: ${output_directory}"

    if [[ ! -s "${run_directory}/topology_output.lmp" ]]; then
        echo "ERROR: ${run_directory}/topology_output.lmp is missing." >&2
        failed=$((failed + 1))
        continue
    fi

    if [[ ! -s "${run_directory}/assembly_manifest.json" ]]; then
        echo "ERROR: ${run_directory}/assembly_manifest.json is missing." >&2
        failed=$((failed + 1))
        continue
    fi

    SECONDS=0
    if nio-md-prep analyze-coverage "${run_directory}" \
        --trajectory "${trajectory_name}" \
        --output "${output_directory}" \
        --last-frames "${last_frames}" \
        --grid-spacing "${grid_spacing}" \
        --blocks "${blocks}" \
        --timestep-fs 0.5; then
        elapsed="${SECONDS}"
        total_seconds=$((total_seconds + elapsed))
        completed=$((completed + 1))
        echo "COMPLETED: ${system_name}/${stage} in ${elapsed} seconds"
    else
        elapsed="${SECONDS}"
        total_seconds=$((total_seconds + elapsed))
        failed=$((failed + 1))
        echo "FAILED: ${system_name}/${stage} after ${elapsed} seconds" >&2
    fi
done < <(
    find "${prepared_root}" -mindepth 2 -maxdepth 2 -type f \
        \( -name 'hold-300K.lammpstrj' -o -name 'hold-400K.lammpstrj' \) \
        -size +0c -print0 | sort -z
)

if (( found == 0 )); then
    echo "ERROR: no non-empty hold-300K.lammpstrj or hold-400K.lammpstrj files were found." >&2
    exit 1
fi

workbook="${prepared_root}/coverage_summary.xlsx"
echo
echo "Creating consolidated workbook: ${workbook}"
if nio-md-prep summarize-coverage "${prepared_root}" --output "${workbook}"; then
    echo "WORKBOOK: ${workbook}"
else
    failed=$((failed + 1))
    echo "FAILED: consolidated coverage workbook" >&2
fi

echo
echo "Coverage batch finished."
echo "Found: ${found}"
echo "Completed: ${completed}"
echo "Failed: ${failed}"
echo "Total analyzer time: ${total_seconds} seconds"

if (( failed > 0 )); then
    exit 1
fi
