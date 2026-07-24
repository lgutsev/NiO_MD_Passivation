#!/bin/bash
#SBATCH -p workq
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -c 1
#SBATCH -t 72:00:00
#SBATCH -A loni_perovsk27
#SBATCH -J nio.rep.sdep
#SBATCH --array=0-3
#SBATCH -o nio.rep.sdep.%A_%a.out

set -euo pipefail

module purge
module load lammps/22Jul2025-r9.0-cpu

export OMP_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1

if ! LAMMPS_EXE="$(command -v lmp_mpi)"; then
    echo "ERROR: lmp_mpi is unavailable after loading LAMMPS." >&2
    exit 1
fi

systems=(
    me-4pacz-then-dcz-4p
    me-4pacz-then-meo-2pacz
)
seed_labels=(02 03)

replicate_index=$((SLURM_ARRAY_TASK_ID / ${#systems[@]}))
system_index=$((SLURM_ARRAY_TASK_ID % ${#systems[@]}))
seed_label="${seed_labels[${replicate_index}]}"
name="${systems[${system_index}]}"
run_directory="${SLURM_SUBMIT_DIR}/prepared/replicates/seed-${seed_label}/${name}"

if [[ ! -s "${run_directory}/topology_output.lmp" ]]; then
    echo "ERROR: ${run_directory}/topology_output.lmp is missing; build first." >&2
    exit 1
fi
if [[ -e "${run_directory}/deposited.data" ]]; then
    echo "ERROR: deposited.data already exists in ${run_directory}." >&2
    exit 1
fi

cd "${run_directory}"
echo "Running ${name}, seed ${seed_label}, sequential deposition"
SECONDS=0
srun -n "${SLURM_NTASKS}" "${LAMMPS_EXE}" -in deposition.in \
    -log "log.deposition.${SLURM_JOB_ID}.lammps"
echo "Completed in ${SECONDS} seconds."
