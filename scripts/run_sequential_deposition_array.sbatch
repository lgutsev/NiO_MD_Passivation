#!/bin/bash
#SBATCH -p workq
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -c 1
#SBATCH -t 72:00:00
#SBATCH -A loni_perovsk27
#SBATCH -J nio.seq.deposit
#SBATCH --array=0-2
#SBATCH -o nio.seq.deposit.%A_%a.out

set -euo pipefail

module purge
module load lammps/22Jul2025-r9.3-cpu

export OMP_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1

if ! LAMMPS_EXE="$(command -v lmp_mpi)"; then
    echo "ERROR: lmp_mpi is unavailable after loading the QBD LAMMPS module." >&2
    exit 1
fi

systems=(
    me-4pacz-then-dcz-4p
    me-4pacz-then-meo-2pacz
    me-4pacz-then-meo-4padbc
)

name="${systems[${SLURM_ARRAY_TASK_ID}]}"
prepared_root="${PREPARED_ROOT:-prepared}"
run_directory="${SLURM_SUBMIT_DIR}/${prepared_root}/${name}"

if [[ ! -s "${run_directory}/topology_output.lmp" ]]; then
    echo "ERROR: ${run_directory}/topology_output.lmp is missing; run the sequential build job first." >&2
    exit 1
fi

if [[ ! -s "${run_directory}/deposition.in" ]]; then
    echo "ERROR: ${run_directory}/deposition.in is missing; regenerate the prepared inputs." >&2
    exit 1
fi

if [[ -e "${run_directory}/deposited.data" ]]; then
    echo "ERROR: ${run_directory}/deposited.data already exists; refusing to overwrite a completed stage." >&2
    exit 1
fi

cd "${run_directory}"
echo "Running ${name} sequential stage-2 deposition in ${run_directory}"
echo "LAMMPS: ${LAMMPS_EXE}"
echo "MPI tasks: ${SLURM_NTASKS}"

SECONDS=0
srun -n "${SLURM_NTASKS}" "${LAMMPS_EXE}" -in deposition.in -log "log.deposition.${SLURM_JOB_ID}.lammps"
echo "Sequential deposition for ${name} took ${SECONDS} seconds."
