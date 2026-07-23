#!/bin/bash
#SBATCH -p workq
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -c 1
#SBATCH -t 24:00:00
#SBATCH -A loni_perovsk27
#SBATCH -J nio.seq.hold300
#SBATCH --array=0-2
#SBATCH -o nio.seq.hold300.%A_%a.out

set -euo pipefail

module purge
module load lammps/22Jul2025-r9.0-cpu

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
run_directory="${SLURM_SUBMIT_DIR}/prepared/${name}"

if [[ ! -s "${run_directory}/deposited.data" ]]; then
    echo "ERROR: ${run_directory}/deposited.data is missing; finish sequential deposition first." >&2
    exit 1
fi

if [[ ! -s "${run_directory}/hold-300K.in" ]]; then
    echo "ERROR: ${run_directory}/hold-300K.in is missing; regenerate the prepared inputs." >&2
    exit 1
fi

if [[ -e "${run_directory}/held-300K.data" ]]; then
    echo "ERROR: ${run_directory}/held-300K.data already exists; refusing to overwrite a completed stage." >&2
    exit 1
fi

cd "${run_directory}"
echo "Running ${name} sequential compressed-film hold at 300 K in ${run_directory}"
echo "LAMMPS: ${LAMMPS_EXE}"
echo "MPI tasks: ${SLURM_NTASKS}"

SECONDS=0
srun -n "${SLURM_NTASKS}" "${LAMMPS_EXE}" -in hold-300K.in -log "log.hold-300K.${SLURM_JOB_ID}.lammps"
echo "Sequential 300 K hold for ${name} took ${SECONDS} seconds."
