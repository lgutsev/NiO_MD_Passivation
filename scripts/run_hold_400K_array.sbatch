#!/bin/bash
#SBATCH -p workq
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -c 1
#SBATCH -t 48:00:00
#SBATCH -A loni_perovsk27
#SBATCH -J nio.hold400
#SBATCH --array=0-3
#SBATCH -o nio.hold400.%A_%a.out

set -euo pipefail

module purge
module load lammps/04Jul2026-cpu

export OMP_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1

systems=(
    me-4pacz-alone
    me-4pacz-meo-2pacz-cosam
    me-4pacz-meo-4padbc-cosam
    me-4pacz-dcz-4p-cosam
)

name="${systems[${SLURM_ARRAY_TASK_ID}]}"
run_directory="${SLURM_SUBMIT_DIR}/prepared/${name}"

if [[ ! -s "${run_directory}/deposited.data" ]]; then
    echo "ERROR: ${run_directory}/deposited.data is missing; finish deposition first." >&2
    exit 1
fi

if [[ ! -s "${run_directory}/hold-400K.in" ]]; then
    echo "ERROR: ${run_directory}/hold-400K.in is missing; regenerate the prepared inputs." >&2
    exit 1
fi

if [[ -e "${run_directory}/held-400K.data" ]]; then
    echo "ERROR: ${run_directory}/held-400K.data already exists; refusing to overwrite a completed stage." >&2
    exit 1
fi

cd "${run_directory}"
echo "Running ${name} gradual 300-to-400 K heating and compressed-film hold in ${run_directory}"
echo "LAMMPS: $(command -v lmp)"
echo "MPI tasks: ${SLURM_NTASKS}"

SECONDS=0
srun -n "${SLURM_NTASKS}" lmp -in hold-400K.in -log "log.hold-400K.${SLURM_JOB_ID}.lammps"
echo "400 K heating and hold for ${name} took ${SECONDS} seconds."
