#!/bin/bash
#SBATCH -p single
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 02:00:00
#SBATCH -A loni_perovsk27
#SBATCH -J nio.decomp.build
#SBATCH --array=0-7
#SBATCH -o nio.decomp.build.%A_%a.out

set -euo pipefail

module purge
module load python/3.11.5-anaconda

cd "${SLURM_SUBMIT_DIR}"
source .venv/bin/activate

if ! command -v python >/dev/null 2>&1; then
    echo "ERROR: python is unavailable after module and virtual-environment setup." >&2
    exit 1
fi

systems=(
    me-4pacz-alone
    me-4pacz-meo-2pacz-cosam
    me-4pacz-meo-4padbc-cosam
    me-4pacz-dcz-4p-cosam
    me-4pacz-high-dose
    me-4pacz-then-dcz-4p
    me-4pacz-then-meo-2pacz
    me-4pacz-then-meo-4padbc
)

base_name="${systems[${SLURM_ARRAY_TASK_ID}]}"
name="${base_name}${SYSTEM_SUFFIX:-}"
prepared_root="${PREPARED_ROOT:-prepared}"
output="${prepared_root}/${name}"

if [[ ! -s "${output}/topology_output.lmp" ]]; then
    echo "ERROR: ${output}/topology_output.lmp is missing; this system has not been built." >&2
    exit 1
fi

echo "REFRESHING DECOMPRESSION INPUTS: ${name} (array index ${SLURM_ARRAY_TASK_ID})"
PYTHONPATH="${SLURM_SUBMIT_DIR}/src${PYTHONPATH:+:${PYTHONPATH}}" \
    python -m nio_md_prep refresh-inputs "studies/${base_name}.toml" \
    --output "${output}"

for required in decompress-300K.in decompress-400K.in; do
    if [[ ! -s "${output}/${required}" ]]; then
        echo "ERROR: ${output}/${required} was not produced." >&2
        exit 1
    fi
done

echo "READY: ${name}"
echo "Existing data, restarts, trajectories, and logs were not modified."
