#!/bin/bash
#SBATCH -p single
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 04:00:00
#SBATCH -A loni_perovsk27
#SBATCH -J nio.rep.build
#SBATCH --array=0-5
#SBATCH -o nio.rep.build.%A_%a.out

set -euo pipefail

module purge
module load python/3.11.5-anaconda
module load lammps/22Jul2025-r9.0-cpu

cd "${SLURM_SUBMIT_DIR}"
source .venv/bin/activate

for executable in python packmol nio-md-prep lmp_mpi; do
    if ! command -v "${executable}" >/dev/null 2>&1; then
        echo "ERROR: ${executable} is unavailable after environment setup." >&2
        exit 1
    fi
done

systems=(
    me-4pacz-alone
    me-4pacz-dcz-4p-cosam
    me-4pacz-meo-2pacz-cosam
)
seed_labels=(02 03)
packmol_seeds=(202607242 202607243)
velocity_seeds=(31415927 27182818)

replicate_index=$((SLURM_ARRAY_TASK_ID / ${#systems[@]}))
system_index=$((SLURM_ARRAY_TASK_ID % ${#systems[@]}))

seed_label="${seed_labels[${replicate_index}]}"
packmol_seed="${packmol_seeds[${replicate_index}]}"
velocity_seed="${velocity_seeds[${replicate_index}]}"
name="${systems[${system_index}]}"
output="prepared/replicates/seed-${seed_label}/${name}"

if [[ -e "${output}" ]]; then
    if [[ -s "${output}/topology_output.lmp" ]]; then
        echo "RECOVERING: ${output} is already assembled."
        echo "Rerunning validation through a one-rank Slurm job step."
        nio-md-prep validate "${output}"
        echo "VALIDATED: ${output}"
        exit 0
    fi
    echo "ERROR: ${output} exists but has no completed topology." >&2
    echo "Inspect it and move it aside before rebuilding; it will not be overwritten." >&2
    exit 1
fi

mkdir -p "$(dirname "${output}")"
echo "Building ${name}, independent seed ${seed_label}"
echo "Packmol seed: ${packmol_seed}"
echo "LAMMPS velocity seed: ${velocity_seed}"

nio-md-prep build "studies/${name}.toml" \
    --output "${output}" \
    --packmol-seed "${packmol_seed}" \
    --velocity-seed "${velocity_seed}"

nio-md-prep validate "${output}"

if [[ ! -s "${output}/topology_output.lmp" ]]; then
    echo "ERROR: ${output}/topology_output.lmp was not produced." >&2
    exit 1
fi

echo "BUILT: ${output}"
echo "No deposition or hold job was submitted."
