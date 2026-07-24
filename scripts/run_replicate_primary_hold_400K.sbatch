#!/bin/bash
#SBATCH -p workq
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -c 1
#SBATCH -t 48:00:00
#SBATCH -A loni_perovsk27
#SBATCH -J nio.rep.h400
#SBATCH --array=0-5
#SBATCH -o nio.rep.h400.%A_%a.out

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
    me-4pacz-alone
    me-4pacz-dcz-4p-cosam
    me-4pacz-meo-2pacz-cosam
)
seed_labels=(02 03)

replicate_index=$((SLURM_ARRAY_TASK_ID / ${#systems[@]}))
system_index=$((SLURM_ARRAY_TASK_ID % ${#systems[@]}))
seed_label="${seed_labels[${replicate_index}]}"
name="${systems[${system_index}]}"
run_directory="${SLURM_SUBMIT_DIR}/prepared/replicates/seed-${seed_label}/${name}"

if [[ ! -s "${run_directory}/deposited.data" ]]; then
    echo "ERROR: ${run_directory}/deposited.data is missing; deposit first." >&2
    exit 1
fi
if [[ ! -s "${run_directory}/hold-400K.in" ]]; then
    echo "ERROR: ${run_directory}/hold-400K.in is missing; rebuild or refresh inputs." >&2
    exit 1
fi
if [[ -e "${run_directory}/held-400K.data" ]]; then
    echo "ERROR: held-400K.data already exists in ${run_directory}." >&2
    exit 1
fi

cd "${run_directory}"
echo "Running ${name}, seed ${seed_label}, gradual 300-to-400 K heating and compressed hold"
SECONDS=0
srun -n "${SLURM_NTASKS}" "${LAMMPS_EXE}" -in hold-400K.in \
    -log "log.hold-400K.${SLURM_JOB_ID}.lammps"
echo "Completed in ${SECONDS} seconds."
