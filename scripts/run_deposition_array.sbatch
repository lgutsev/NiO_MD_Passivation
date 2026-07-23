#!/bin/bash
#SBATCH -p workq
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -c 1
#SBATCH -t 72:00:00
#SBATCH -A loni_perovsk27
#SBATCH -J nio.deposit
#SBATCH --array=0-4
#SBATCH -o nio.deposit.%A_%a.out

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
    me-4pacz-alone
    me-4pacz-meo-2pacz-cosam
    me-4pacz-meo-4padbc-cosam
    me-4pacz-dcz-4p-cosam
    me-4pacz-high-dose
)

name="${systems[${SLURM_ARRAY_TASK_ID}]}"
run_directory="${SLURM_SUBMIT_DIR}/prepared/${name}"

if [[ ! -s "${run_directory}/topology_output.lmp" ]]; then
    echo "ERROR: ${run_directory}/topology_output.lmp is missing; run the build job first." >&2
    exit 1
fi

if [[ -e "${run_directory}/deposited.data" ]]; then
    echo "ERROR: ${run_directory}/deposited.data already exists; refusing to overwrite a completed stage." >&2
    exit 1
fi

cd "${run_directory}"
echo "Running ${name} deposition in ${run_directory}"
echo "LAMMPS: ${LAMMPS_EXE}"
echo "MPI tasks: ${SLURM_NTASKS}"

SECONDS=0
srun -n "${SLURM_NTASKS}" "${LAMMPS_EXE}" -in deposition.in -log "log.deposition.${SLURM_JOB_ID}.lammps"
echo "Deposition for ${name} took ${SECONDS} seconds."
