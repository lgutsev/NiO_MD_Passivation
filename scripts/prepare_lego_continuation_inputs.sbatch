#!/bin/bash
#SBATCH -p single
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 00:30:00
#SBATCH -A loni_perovsk27
#SBATCH -J nio.lego.cont.prepare
#SBATCH --array=0-2
#SBATCH -o nio.lego.cont.prepare.%A_%a.out

set -euo pipefail

module purge
module load python/3.11.5-anaconda

cd "${SLURM_SUBMIT_DIR}"
source .venv/bin/activate

systems=(
    me-4pacz-then-dcz-4p
    me-4pacz-then-meo-2pacz
    me-4pacz-then-meo-4padbc
)

base_name="${systems[${SLURM_ARRAY_TASK_ID}]}"
prepared_root="${PREPARED_ROOT:-prepared-lego}"
system_suffix="${SYSTEM_SUFFIX:--lego-seeded}"
run_directory="${prepared_root}/${base_name}${system_suffix}"

if [[ ! -s "${run_directory}/topology_output.lmp" ]]; then
    echo "ERROR: ${run_directory}/topology_output.lmp is missing." >&2
    exit 1
fi

nio-md-prep prepare-lego-continuation \
    "studies/${base_name}.toml" \
    --output "${run_directory}" \
    --additional-drop "${LEGO_ADDITIONAL_DROP:-15.0}" \
    --steps "${LEGO_CONTINUATION_STEPS:-300000}"

for required in continue-deposition.in hold-300K.in hold-400K.in decompress-300K.in decompress-400K.in; do
    if [[ ! -s "${run_directory}/${required}" ]]; then
        echo "ERROR: ${run_directory}/${required} was not produced." >&2
        exit 1
    fi
done

echo "PREPARED: ${run_directory}"
echo "No LAMMPS job was submitted."
