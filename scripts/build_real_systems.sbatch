#!/bin/bash
#SBATCH -p workq
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 02:00:00
#SBATCH -A loni_perovsk26
#SBATCH -J nio.build
#SBATCH -o nio.build.%j.out

set -euo pipefail

module purge
module load python/3.11.5-anaconda
module load lammps/04Jul2026-cpu

cd "${SLURM_SUBMIT_DIR}"
source .venv/bin/activate

for executable in python packmol nio-md-prep lmp; do
    if ! command -v "${executable}" >/dev/null 2>&1; then
        echo "ERROR: ${executable} is not available after module and virtual-environment setup." >&2
        exit 1
    fi
done

echo "Repository: $(pwd)"
echo "Commit: $(git rev-parse --short HEAD)"
echo "Python: $(python --version 2>&1)"
echo "Packmol: $(command -v packmol)"
echo "LAMMPS: $(command -v lmp)"

systems=(
    me-4pacz-alone
    me-4pacz-meo-2pacz-cosam
    me-4pacz-meo-4padbc-cosam
    me-4pacz-dcz-4p-cosam
)

mkdir -p prepared

for name in "${systems[@]}"; do
    echo "BUILDING: ${name}"
    nio-md-prep build "studies/${name}.toml" --output "prepared/${name}"

    echo "VALIDATING: ${name}"
    nio-md-prep validate "prepared/${name}"

    if [[ ! -s "prepared/${name}/topology_output.lmp" ]]; then
        echo "ERROR: prepared/${name}/topology_output.lmp was not produced." >&2
        exit 1
    fi
done

echo "All four real CoSAM systems were built and deposition-smoke validated."
