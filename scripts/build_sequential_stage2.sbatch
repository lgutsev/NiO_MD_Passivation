#!/bin/bash
#SBATCH -p workq
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 02:00:00
#SBATCH -A loni_perovsk27
#SBATCH -J nio.seq.build
#SBATCH -o nio.seq.build.%j.out

set -euo pipefail

module purge
module load python/3.11.5-anaconda
module load lammps/22Jul2025-r9.0-cpu

cd "${SLURM_SUBMIT_DIR}"
source .venv/bin/activate

for executable in python packmol nio-md-prep lmp; do
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

mkdir -p prepared

for name in "${systems[@]}"; do
    output="prepared/${name}"
    if [[ -e "${output}/topology_output.lmp" ]]; then
        echo "ERROR: ${output}/topology_output.lmp already exists; refusing to overwrite a prepared stage." >&2
        exit 1
    fi

    echo "BUILDING SEQUENTIAL STAGE 2: ${name}"
    nio-md-prep prepare-sequential-stage2 "studies/${name}.toml" \
        --primary-final "${primary}" \
        --output "${output}"

    echo "VALIDATING: ${name}"
    nio-md-prep validate "${output}"

    if [[ ! -s "${output}/topology_output.lmp" ]]; then
        echo "ERROR: ${output}/topology_output.lmp was not produced." >&2
        exit 1
    fi
done

echo "All three sequential stage-2 systems were built and deposition-smoke validated."
