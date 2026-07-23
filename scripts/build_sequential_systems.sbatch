#!/bin/bash
#SBATCH -p single
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 02:00:00
#SBATCH -A loni_perovsk27
#SBATCH -J nio.seq.build
#SBATCH --array=0-2
#SBATCH -o nio.seq.build.%A_%a.out

set -euo pipefail

module purge
module load python/3.11.5-anaconda

cd "${SLURM_SUBMIT_DIR}"
source .venv/bin/activate

for executable in python packmol nio-md-prep; do
    if ! command -v "${executable}" >/dev/null 2>&1; then
        echo "ERROR: ${executable} is not available after module and virtual-environment setup." >&2
        exit 1
    fi
done

primary="prepared/me-4pacz-alone/held-300K.data"
if [[ ! -s "${primary}" ]]; then
    echo "ERROR: ${primary} is missing; finish the Me-4PACz 300 K hold first." >&2
    exit 1
fi

systems=(
    me-4pacz-then-dcz-4p
    me-4pacz-then-meo-2pacz
    me-4pacz-then-meo-4padbc
)

name="${systems[${SLURM_ARRAY_TASK_ID}]}"
output="prepared/${name}"

mkdir -p prepared

if [[ -e "${output}" ]]; then
    echo "ERROR: ${output} already exists; remove or archive it before rebuilding." >&2
    exit 1
fi

echo "BUILDING SEQUENTIAL STAGE 2: ${name} (array index ${SLURM_ARRAY_TASK_ID})"
nio-md-prep prepare-sequential-stage2 "studies/${name}.toml" \
    --primary-final "${primary}" \
    --output "${output}"

for required in topology_output.lmp force_field_settings.lmp deposition.in hold-300K.in; do
    if [[ ! -s "${output}/${required}" ]]; then
        echo "ERROR: ${output}/${required} was not produced." >&2
        exit 1
    fi
done

echo "BUILT: ${name}"
echo "No deposition or hold job was submitted."
