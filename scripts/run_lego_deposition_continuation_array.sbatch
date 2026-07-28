#!/bin/bash
#SBATCH -p workq
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -c 1
#SBATCH -t 72:00:00
#SBATCH -A loni_perovsk27
#SBATCH -J nio.lego.cont
#SBATCH --array=0-2
#SBATCH -o nio.lego.cont.%A_%a.out

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

base_name="${systems[${SLURM_ARRAY_TASK_ID}]}"
prepared_root="${PREPARED_ROOT:-prepared-lego}"
system_suffix="${SYSTEM_SUFFIX:--lego-seeded}"
run_directory="${SLURM_SUBMIT_DIR}/${prepared_root}/${base_name}${system_suffix}"
completion_marker="${run_directory}/deposition-continuation.complete"

if [[ ! -s "${run_directory}/deposited.data" ]]; then
    echo "ERROR: ${run_directory}/deposited.data is missing." >&2
    exit 1
fi
if [[ ! -s "${run_directory}/continue-deposition.in" ]]; then
    echo "ERROR: ${run_directory}/continue-deposition.in is missing; prepare the continuation inputs first." >&2
    exit 1
fi
if [[ -e "${completion_marker}" ]]; then
    echo "ERROR: ${completion_marker} exists; refusing to repeat the final compression." >&2
    exit 1
fi
if [[ -e "${run_directory}/deposited-continued.data" ]]; then
    echo "ERROR: stale deposited-continued.data exists; inspect it before retrying." >&2
    exit 1
fi
for downstream in held-300K.data held-400K.data; do
    if [[ -e "${run_directory}/${downstream}" ]]; then
        echo "ERROR: ${run_directory}/${downstream} already exists; final compression must precede holds." >&2
        exit 1
    fi
done

cd "${run_directory}"
echo "Continuing lego deposition for ${base_name}${system_suffix}"
echo "LAMMPS: ${LAMMPS_EXE}"
echo "MPI tasks: ${SLURM_NTASKS}"

SECONDS=0
srun -n "${SLURM_NTASKS}" "${LAMMPS_EXE}" \
    -in continue-deposition.in \
    -log "log.deposition-continuation.${SLURM_JOB_ID}.lammps"

if [[ ! -s deposited-continued.data ]]; then
    echo "ERROR: LAMMPS returned successfully but deposited-continued.data is missing." >&2
    exit 1
fi

archive="deposited.pre-continuation.${SLURM_JOB_ID}.data"
cp -p -- deposited.data "${archive}"
mv -- deposited-continued.data deposited.data
printf 'job_id=%s\nsource_backup=%s\nelapsed_seconds=%s\n' \
    "${SLURM_JOB_ID}" "${archive}" "${SECONDS}" > deposition-continuation.complete

echo "Final lego compression took ${SECONDS} seconds."
echo "Promoted deposited-continued.data to deposited.data."
echo "Original deposited.data preserved as ${archive}."
