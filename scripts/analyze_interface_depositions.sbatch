#!/bin/bash
#SBATCH -p single
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 24:00:00
#SBATCH -A loni_perovsk27
#SBATCH -J nio.dep.analysis
#SBATCH -o nio.dep.analysis.%j.out

set -euo pipefail

module purge
module load python/3.11.5-anaconda

cd "${SLURM_SUBMIT_DIR}"
source .venv/bin/activate

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

if ! command -v nio-md-prep >/dev/null 2>&1; then
    echo "ERROR: nio-md-prep is unavailable after activating .venv." >&2
    exit 1
fi

if ! python -c "import numpy, scipy, openpyxl" >/dev/null 2>&1; then
    echo "ERROR: interfacial-analysis dependencies are missing." >&2
    echo "Install once from the repository root with:" >&2
    echo "  python -m pip install -e '.[analysis]'" >&2
    exit 1
fi

prepared_root="${PREPARED_ROOT:-${SLURM_SUBMIT_DIR}/prepared}"
if [[ ! -d "${prepared_root}" ]]; then
    echo "ERROR: ${prepared_root} does not exist." >&2
    exit 1
fi

# Full deposition trajectories are required for a comparable exchange rate.
# Submission-time example:
# sbatch --export=ALL,PREPARED_ROOT="${PWD}/prepared-rerun",EXPERIMENTAL_CSV="/path/to/wbg_pscs_jv_results.csv" \
#   scripts/analyze_interface_depositions.sbatch
last_frames="${LAST_FRAMES:-0}"
blocks="${BLOCKS:-5}"
contact_cutoff="${CONTACT_CUTOFF:-3.25}"
sensitivity_cutoffs="${SENSITIVITY_CUTOFFS:-3.0,3.25,3.5}"
surface_coordination_cutoff="${SURFACE_COORDINATION_CUTOFF:-2.8}"
persistence_threshold="${PERSISTENCE_THRESHOLD:-0.50}"
rdf_rmax="${RDF_RMAX:-15.0}"
exchange_min_dwell_frames="${EXCHANGE_MIN_DWELL_FRAMES:-2}"
exchange_max_vacancy_gap_frames="${EXCHANGE_MAX_VACANCY_GAP_FRAMES:-3}"

if [[ ! "${last_frames}" =~ ^[0-9]+$ ]]; then
    echo "ERROR: LAST_FRAMES must be a non-negative integer." >&2
    exit 1
fi
if [[ ! "${exchange_min_dwell_frames}" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: EXCHANGE_MIN_DWELL_FRAMES must be a positive integer." >&2
    exit 1
fi
if [[ ! "${exchange_max_vacancy_gap_frames}" =~ ^[0-9]+$ ]]; then
    echo "ERROR: EXCHANGE_MAX_VACANCY_GAP_FRAMES must be a non-negative integer." >&2
    exit 1
fi

echo "Repository: ${SLURM_SUBMIT_DIR}"
echo "Commit: $(git rev-parse --short HEAD)"
echo "Python: $(python --version 2>&1)"
echo "PREPARED_ROOT: ${prepared_root}"
echo "LAST_FRAMES: ${last_frames} (0 means all frames)"
echo "BLOCKS: ${blocks}"
echo "CONTACT_CUTOFF: ${contact_cutoff} A"
echo "EXCHANGE_MIN_DWELL_FRAMES: ${exchange_min_dwell_frames}"
echo "EXCHANGE_MAX_VACANCY_GAP_FRAMES: ${exchange_max_vacancy_gap_frames}"

found=0
completed=0
failed=0
total_seconds=0

while IFS= read -r -d '' trajectory; do
    found=$((found + 1))
    run_directory="$(dirname "${trajectory}")"
    system_name="${run_directory#${prepared_root}/}"
    output_directory="${run_directory}/interface-analysis-deposition"

    echo
    echo "ANALYZING: ${system_name}/deposition.lammpstrj"
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
    if nio-md-prep analyze-interface "${run_directory}" \
        --trajectory deposition.lammpstrj \
        --output "${output_directory}" \
        --last-frames "${last_frames}" \
        --blocks "${blocks}" \
        --contact-cutoff "${contact_cutoff}" \
        --sensitivity-cutoffs "${sensitivity_cutoffs}" \
        --surface-coordination-cutoff "${surface_coordination_cutoff}" \
        --persistence-threshold "${persistence_threshold}" \
        --exchange-min-dwell-frames "${exchange_min_dwell_frames}" \
        --exchange-max-vacancy-gap-frames "${exchange_max_vacancy_gap_frames}" \
        --rdf-rmax "${rdf_rmax}" \
        --timestep-fs 0.5; then
        elapsed="${SECONDS}"
        total_seconds=$((total_seconds + elapsed))
        completed=$((completed + 1))
        echo "COMPLETED: ${system_name} in ${elapsed} seconds"
    else
        elapsed="${SECONDS}"
        total_seconds=$((total_seconds + elapsed))
        failed=$((failed + 1))
        echo "FAILED: ${system_name} after ${elapsed} seconds" >&2
    fi
done < <(
    find "${prepared_root}" -mindepth 2 -type f \
        -name 'deposition.lammpstrj' -size +0c -print0 | sort -z
)

if (( found == 0 )); then
    echo "ERROR: no non-empty deposition trajectories were found." >&2
    exit 1
fi

interface_workbook="${prepared_root}/interface_structure_summary.xlsx"
echo
echo "Creating consolidated workbook: ${interface_workbook}"
if nio-md-prep summarize-interface "${prepared_root}" \
    --output "${interface_workbook}"; then
    echo "WORKBOOK: ${interface_workbook}"
else
    failed=$((failed + 1))
    echo "FAILED: consolidated interfacial workbook" >&2
fi

publication_workbook="${prepared_root}/publication_summary.xlsx"
publication_args=()
if [[ -n "${EXPERIMENTAL_CSV:-}" ]]; then
    if [[ ! -s "${EXPERIMENTAL_CSV}" ]]; then
        echo "ERROR: EXPERIMENTAL_CSV=${EXPERIMENTAL_CSV} is missing or empty." >&2
        failed=$((failed + 1))
    else
        publication_args+=(--experimental "${EXPERIMENTAL_CSV}")
    fi
fi
echo
echo "Creating publication-ready workbook: ${publication_workbook}"
if nio-md-prep summarize-publication "${prepared_root}" \
    --output "${publication_workbook}" "${publication_args[@]}"; then
    echo "PUBLICATION WORKBOOK: ${publication_workbook}"
else
    failed=$((failed + 1))
    echo "FAILED: publication-ready workbook" >&2
    echo "Run the hold coverage and interface analyses first." >&2
fi

echo
echo "Deposition interfacial-analysis batch finished."
echo "Found: ${found}"
echo "Completed: ${completed}"
echo "Failed: ${failed}"
echo "Total analyzer time: ${total_seconds} seconds"

if (( failed > 0 )); then
    exit 1
fi
