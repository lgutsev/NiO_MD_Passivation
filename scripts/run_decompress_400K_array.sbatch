#!/bin/bash
#SBATCH -p workq
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -c 1
#SBATCH -t 72:00:00
#SBATCH -A loni_perovsk27
#SBATCH -J nio.decomp400
#SBATCH --array=0-7
#SBATCH -o nio.decomp400.%A_%a.out

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
    me-4pacz-alone
    me-4pacz-meo-2pacz-cosam
    me-4pacz-meo-4padbc-cosam
    me-4pacz-dcz-4p-cosam
    me-4pacz-high-dose
    me-4pacz-then-dcz-4p
    me-4pacz-then-meo-2pacz
    me-4pacz-then-meo-4padbc
)

name="${systems[${SLURM_ARRAY_TASK_ID}]}"
run_directory="${SLURM_SUBMIT_DIR}/prepared/${name}"

if [[ ! -s "${run_directory}/held-400K.data" ]]; then
    echo "ERROR: ${run_directory}/held-400K.data is missing; finish the compressed 400 K hold first." >&2
    exit 1
fi

if [[ ! -s "${run_directory}/decompress-400K.in" ]]; then
    echo "ERROR: ${run_directory}/decompress-400K.in is missing; run the decompression-input builder first." >&2
    exit 1
fi

if [[ -e "${run_directory}/relaxed-400K.data" ]]; then
    echo "ERROR: ${run_directory}/relaxed-400K.data already exists; refusing to overwrite a completed stage." >&2
    exit 1
fi

cd "${run_directory}"
echo "Running ${name} gradual decompression and relaxed-film hold at 400 K in ${run_directory}"
echo "LAMMPS: ${LAMMPS_EXE}"
echo "MPI tasks: ${SLURM_NTASKS}"

SECONDS=0
srun -n "${SLURM_NTASKS}" "${LAMMPS_EXE}" -in decompress-400K.in \
    -log "log.decompress-400K.${SLURM_JOB_ID}.lammps"
echo "400 K decompression and relaxed hold for ${name} took ${SECONDS} seconds."
