#!/bin/bash
#SBATCH -p single
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 02:00:00
#SBATCH -A loni_perovsk27
#SBATCH -J nio.build
#SBATCH -o nio.build.%j.out

set -euo pipefail

module purge
module load python/3.11.5-anaconda
module load lammps/22Jul2025-r9.7-cpu

cd "${SLURM_SUBMIT_DIR}"
source .venv/bin/activate

for executable in python packmol nio-md-prep lmp_mpi; do
    if ! command -v "${executable}" >/dev/null 2>&1; then
        echo "ERROR: ${executable} is not available after module and virtual-environment setup." >&2
        exit 1
    fi
done

echo "Repository: $(pwd)"
echo "Commit: $(git rev-parse --short HEAD)"
echo "Python: $(python --version 2>&1)"
echo "Packmol: $(command -v packmol)"
echo "LAMMPS: $(command -v lmp_mpi)"

default_systems=(
    me-4pacz-alone
    me-4pacz-meo-2pacz-cosam
    me-4pacz-meo-4padbc-cosam
    me-4pacz-dcz-4p-cosam
    me-4pacz-high-dose
)

if [[ -n "${SYSTEMS:-}" ]]; then
    read -r -a systems <<< "${SYSTEMS}"
else
    systems=("${default_systems[@]}")
fi
prepared_root="${PREPARED_ROOT:-prepared}"
mkdir -p "${prepared_root}"

for name in "${systems[@]}"; do
    echo "BUILDING: ${name}"
    build_args=(
        "studies/${name}.toml"
        --output "${prepared_root}/${name}"
    )
    if [[ -n "${PACKMOL_SEED:-}" ]]; then
        build_args+=(--packmol-seed "${PACKMOL_SEED}")
    fi
    if [[ -n "${VELOCITY_SEED:-}" ]]; then
        build_args+=(--velocity-seed "${VELOCITY_SEED}")
    fi
    nio-md-prep build "${build_args[@]}"

    echo "VALIDATING: ${name}"
    nio-md-prep validate "${prepared_root}/${name}"

    if [[ ! -s "${prepared_root}/${name}/topology_output.lmp" ]]; then
        echo "ERROR: ${prepared_root}/${name}/topology_output.lmp was not produced." >&2
        exit 1
    fi
done

echo "All requested primary-deposition systems were built and validated."
