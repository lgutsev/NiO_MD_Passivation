#!/bin/bash
#SBATCH -p single
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 12:00:00
#SBATCH -A loni_perovsk27
#SBATCH -J nio.tier0
#SBATCH -o nio.tier0.%j.out

# Tier-0 electronic decision gate (docs/electronic-proxy.md). Analysis only:
# reads relaxed-hold trajectories, launches no MD or DFT.
#
# Submission-time overrides:
# sbatch --export=ALL,SECONDARY=dcz-4p,REFERENCE_DIPOLES=$PWD/reference_dipoles.json \
#   scripts/analyze_electronic_tier0.sbatch

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

if ! python -c "import numpy, scipy" >/dev/null 2>&1; then
    echo "ERROR: analysis dependencies are missing." >&2
    echo "Install once from the repository root with:" >&2
    echo "  python -m pip install -e '.[analysis]'" >&2
    exit 1
fi

secondary="${SECONDARY:-dcz-4p}"
stage="${STAGE:-relax-300K}"
last_frames="${LAST_FRAMES:-100}"
export_conformers="${EXPORT_CONFORMERS:-50}"
include_reference="${INCLUDE_ME4PACZ_ALONE:-1}"
output_root="${OUTPUT_ROOT:-${SLURM_SUBMIT_DIR}/electronic-tier0/${secondary}}"
roots=(${CONVENTIONAL_ROOTS:-prepared prepared-rerun prepared-rerun2})

extra_args=()
if [[ -n "${REFERENCE_DIPOLES:-}" ]]; then
    if [[ ! -s "${REFERENCE_DIPOLES}" ]]; then
        echo "ERROR: REFERENCE_DIPOLES=${REFERENCE_DIPOLES} is missing or empty." >&2
        exit 1
    fi
    extra_args+=(--reference-dipoles "${REFERENCE_DIPOLES}")
fi

# label|run directory
runs=()
for root in "${roots[@]}"; do
    runs+=("cosam|${root}/me-4pacz-${secondary}-cosam")
    runs+=("sequential|${root}/me-4pacz-then-${secondary}")
    if [[ "${include_reference}" == "1" ]]; then
        runs+=("me-4pacz-alone|${root}/me-4pacz-alone")
    fi
done
runs+=("gap-seeded-i|prepared-lego/me-4pacz-then-${secondary}-lego-seeded")
runs+=("gap-seeded-ii|prepared-lego2/me-4pacz-then-${secondary}-lego2-seeded")

echo "Repository: ${SLURM_SUBMIT_DIR}"
echo "Commit: $(git rev-parse --short HEAD)"
echo "SECONDARY: ${secondary}"
echo "STAGE: ${stage}"
echo "LAST_FRAMES: ${last_frames}"
echo "EXPORT_CONFORMERS: ${export_conformers}"
echo "REFERENCE_DIPOLES: ${REFERENCE_DIPOLES:-<none; charge dipoles only>}"
echo "OUTPUT_ROOT: ${output_root}"

mkdir -p "${output_root}"
summaries=()
completed=0
missing=0
failed=0
for entry in "${runs[@]}"; do
    label="${entry%%|*}"
    run_directory="${SLURM_SUBMIT_DIR}/${entry#*|}"
    trajectory="${run_directory}/${stage}.lammpstrj"
    output_directory="${run_directory}/electronic-proxy-${stage}"
    echo
    echo "RUN: ${label} ${run_directory}"
    if [[ ! -s "${trajectory}" ]]; then
        echo "MISSING: ${trajectory}"
        missing=$((missing + 1))
        continue
    fi
    SECONDS=0
    if nio-md-prep analyze-electronic-proxy "${run_directory}" \
        --trajectory "${stage}.lammpstrj" \
        --output "${output_directory}" \
        --last-frames "${last_frames}" \
        --component-alias stage1_primary=me-4pacz \
        --export-conformers "${export_conformers}" \
        "${extra_args[@]}"; then
        completed=$((completed + 1))
        summaries+=("${label}=${output_directory}")
        echo "COMPLETED in ${SECONDS} seconds"
    else
        failed=$((failed + 1))
        echo "FAILED after ${SECONDS} seconds" >&2
    fi
done

if (( completed > 0 )); then
    echo
    nio-md-prep summarize-electronic-proxy "${summaries[@]}" \
        --compare cosam,sequential --output "${output_root}" || failed=$((failed + 1))
    nio-md-prep select-patch-slots "${summaries[@]}" \
        --compare cosam,sequential --reference-slots 2 \
        --output "${output_root}" || failed=$((failed + 1))
fi

echo
echo "Tier-0 batch finished. Completed: ${completed}; missing: ${missing}; failed: ${failed}"
echo "Gate: ${output_root}/electronic_gate.json"
echo "Slots: ${output_root}/patch_slots.csv"
if (( failed > 0 || completed == 0 )); then
    exit 1
fi
