#!/bin/bash
#SBATCH -p single
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 04:00:00
#SBATCH -A loni_perovsk27
#SBATCH -J nio.rep.seqb
#SBATCH --array=0-3
#SBATCH -o nio.rep.seqb.%A_%a.out

set -euo pipefail

module purge
module load python/3.11.5-anaconda

cd "${SLURM_SUBMIT_DIR}"
source .venv/bin/activate

for executable in python packmol nio-md-prep; do
    if ! command -v "${executable}" >/dev/null 2>&1; then
        echo "ERROR: ${executable} is unavailable after environment setup." >&2
        exit 1
    fi
done

systems=(
    me-4pacz-then-dcz-4p
    me-4pacz-then-meo-2pacz
)
seed_labels=(02 03)
packmol_seeds=(202608242 202608243)
velocity_seeds=(31416927 27183818)

replicate_index=$((SLURM_ARRAY_TASK_ID / ${#systems[@]}))
system_index=$((SLURM_ARRAY_TASK_ID % ${#systems[@]}))

seed_label="${seed_labels[${replicate_index}]}"
packmol_seed="${packmol_seeds[${replicate_index}]}"
velocity_seed="${velocity_seeds[${replicate_index}]}"
name="${systems[${system_index}]}"
replicate_root="prepared/replicates/seed-${seed_label}"
primary="${replicate_root}/me-4pacz-alone/held-300K.data"
output="${replicate_root}/${name}"

if [[ ! -s "${primary}" ]]; then
    echo "ERROR: ${primary} is missing; complete that seed's primary hold first." >&2
    exit 1
fi
if [[ -e "${output}" ]]; then
    echo "ERROR: ${output} already exists; refusing to overwrite it." >&2
    exit 1
fi

echo "Building ${name}, independent seed ${seed_label}"
echo "Primary substrate: ${primary}"
echo "Stage-2 Packmol seed: ${packmol_seed}"
echo "Stage-2 LAMMPS velocity seed: ${velocity_seed}"

nio-md-prep prepare-sequential-stage2 "studies/${name}.toml" \
    --primary-final "${primary}" \
    --output "${output}" \
    --packmol-seed "${packmol_seed}" \
    --velocity-seed "${velocity_seed}"

for required in topology_output.lmp deposition.in hold-300K.in; do
    if [[ ! -s "${output}/${required}" ]]; then
        echo "ERROR: ${output}/${required} was not produced." >&2
        exit 1
    fi
done

echo "BUILT: ${output}"
echo "No deposition or hold job was submitted."
