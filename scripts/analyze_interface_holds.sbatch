#!/bin/bash
#SBATCH -p single
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 24:00:00
#SBATCH -A loni_perovsk27
#SBATCH -J nio.interface
#SBATCH -o nio.interface.%j.out

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

prepared_root="${SLURM_SUBMIT_DIR}/prepared"
if [[ ! -d "${prepared_root}" ]]; then
    echo "ERROR: ${prepared_root} does not exist." >&2
    exit 1
fi

# Submission-time overrides:
# sbatch --export=ALL,LAST_FRAMES=0,CONTACT_CUTOFF=3.25 \
#   scripts/analyze_interface_holds.sbatch
last_frames="${LAST_FRAMES:-100}"
blocks="${BLOCKS:-5}"
contact_cutoff="${CONTACT_CUTOFF:-3.25}"
sensitivity_cutoffs="${SENSITIVITY_CUTOFFS:-3.0,3.25,3.5}"
surface_coordination_cutoff="${SURFACE_COORDINATION_CUTOFF:-2.8}"
persistence_threshold="${PERSISTENCE_THRESHOLD:-0.50}"
rdf_rmax="${RDF_RMAX:-15.0}"

if [[ ! "${last_frames}" =~ ^[0-9]+$ ]]; then
    echo "ERROR: LAST_FRAMES must be a non-negative integer." >&2
    exit 1
fi

echo "Repository: ${SLURM_SUBMIT_DIR}"
echo "Commit: $(git rev-parse --short HEAD)"
echo "Python: $(python --version 2>&1)"
echo "LAST_FRAMES: ${last_frames} (0 means all frames)"
echo "BLOCKS: ${blocks}"
echo "CONTACT_CUTOFF: ${contact_cutoff} A"
echo "SENSITIVITY_CUTOFFS: ${sensitivity_cutoffs} A"
echo "SURFACE_COORDINATION_CUTOFF: ${surface_coordination_cutoff} A"
echo "PERSISTENCE_THRESHOLD: ${persistence_threshold}"
echo "RDF_RMAX: ${rdf_rmax} A"

found=0
completed=0
failed=0
total_seconds=0

while IFS= read -r -d '' trajectory; do
    found=$((found + 1))
    run_directory="$(dirname "${trajectory}")"
    trajectory_name="$(basename "${trajectory}")"
    stage="${trajectory_name%.lammpstrj}"
    system_name="${run_directory#${prepared_root}/}"
    output_directory="${run_directory}/interface-analysis-${stage}"

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
    if nio-md-prep analyze-interface "${run_directory}" \
        --trajectory "${trajectory_name}" \
        --output "${output_directory}" \
        --last-frames "${last_frames}" \
        --blocks "${blocks}" \
        --contact-cutoff "${contact_cutoff}" \
        --sensitivity-cutoffs "${sensitivity_cutoffs}" \
        --surface-coordination-cutoff "${surface_coordination_cutoff}" \
        --persistence-threshold "${persistence_threshold}" \
        --rdf-rmax "${rdf_rmax}" \
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
    find "${prepared_root}" -mindepth 2 -type f \
        \( -name 'hold-300K.lammpstrj' -o -name 'hold-400K.lammpstrj' \
           -o -name 'relax-300K.lammpstrj' -o -name 'relax-400K.lammpstrj' \) \
        -size +0c -print0 | sort -z
)

if (( found == 0 )); then
    echo "ERROR: no non-empty compressed-hold or relaxed-hold trajectories were found." >&2
    exit 1
fi

workbook="${prepared_root}/interface_structure_summary.xlsx"
echo
echo "Creating consolidated workbook: ${workbook}"
if nio-md-prep summarize-interface "${prepared_root}" --output "${workbook}"; then
    echo "WORKBOOK: ${workbook}"
else
    failed=$((failed + 1))
    echo "FAILED: consolidated interfacial workbook" >&2
fi

publication_workbook="${prepared_root}/publication_summary.xlsx"
echo
echo "Creating publication-ready workbook: ${publication_workbook}"
if nio-md-prep summarize-publication "${prepared_root}" \
    --output "${publication_workbook}"; then
    echo "PUBLICATION WORKBOOK: ${publication_workbook}"
else
    failed=$((failed + 1))
    echo "FAILED: publication-ready workbook" >&2
    echo "Run scripts/analyze_coverage_holds.sbatch first so both analyses are available." >&2
fi

echo
echo "Interfacial-analysis batch finished."
echo "Found: ${found}"
echo "Completed: ${completed}"
echo "Failed: ${failed}"
echo "Total analyzer time: ${total_seconds} seconds"

if (( failed > 0 )); then
    exit 1
fi
