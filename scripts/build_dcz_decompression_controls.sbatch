#!/bin/bash
#SBATCH -p single
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 02:00:00
#SBATCH -A loni_perovsk27
#SBATCH -J nio.dcz.ctrl.build
#SBATCH --array=0-3
#SBATCH -o nio.dcz.ctrl.build.%A_%a.out

set -euo pipefail

module purge
module load python/3.11.5-anaconda

cd "${SLURM_SUBMIT_DIR}"
source .venv/bin/activate

if ! command -v python >/dev/null 2>&1; then
    echo "ERROR: python is unavailable after module and virtual-environment setup." >&2
    exit 1
fi

roots=(
    prepared
    prepared-rerun
    prepared-rerun2
    prepared-lego
)
suffixes=(
    ""
    ""
    ""
    "-lego-seeded"
)

prepared_root="${roots[${SLURM_ARRAY_TASK_ID}]}"
system_suffix="${suffixes[${SLURM_ARRAY_TASK_ID}]}"
name="me-4pacz-then-dcz-4p${system_suffix}"
output="${prepared_root}/${name}"

if [[ ! -s "${output}/topology_output.lmp" ]]; then
    echo "ERROR: ${output}/topology_output.lmp is missing." >&2
    exit 1
fi

echo "Preparing DCZ-4P decompression controls in ${output}"
PYTHONPATH="${SLURM_SUBMIT_DIR}/src${PYTHONPATH:+:${PYTHONPATH}}" \
    python -m nio_md_prep refresh-inputs \
    studies/me-4pacz-then-dcz-4p.toml \
    --output "${output}"

for temperature in 300K 400K; do
    for mode in nonsticky nowall; do
        input="${output}/decompress-control-${mode}-${temperature}.in"
        if [[ ! -s "${input}" ]]; then
            echo "ERROR: ${input} was not produced." >&2
            exit 1
        fi
    done
done

echo "READY: ${output}"
echo "Existing data, restarts, trajectories, and logs were not modified."
