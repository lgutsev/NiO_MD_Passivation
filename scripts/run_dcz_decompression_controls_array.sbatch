#!/bin/bash
#SBATCH -p workq
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -c 1
#SBATCH -t 72:00:00
#SBATCH -A loni_perovsk27
#SBATCH -J nio.dcz.ctrl
#SBATCH --array=0-15
#SBATCH -o nio.dcz.ctrl.%A_%a.out

set -euo pipefail

module purge
module load lammps/22Jul2025-r9.3-cpu

export OMP_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1

if ! LAMMPS_EXE="$(command -v lmp_mpi)"; then
    echo "ERROR: lmp_mpi is unavailable after loading the QBD LAMMPS module." >&2
    exit 1
fi

roots=(
    prepared
    prepared-rerun
    prepared-rerun2
    prepared-lego
)
suffixes=(
    ""
    ""
    ""
    "-lego-seeded"
)
temperatures=(300K 400K)
modes=(nonsticky nowall)

task="${SLURM_ARRAY_TASK_ID}"
mode_index=$((task % 2))
temperature_index=$(((task / 2) % 2))
root_index=$((task / 4))

mode="${modes[${mode_index}]}"
temperature="${temperatures[${temperature_index}]}"
prepared_root="${roots[${root_index}]}"
system_suffix="${suffixes[${root_index}]}"
name="me-4pacz-then-dcz-4p${system_suffix}"
run_directory="${SLURM_SUBMIT_DIR}/${prepared_root}/${name}"
source_data="${run_directory}/held-${temperature}.data"
input="${run_directory}/decompress-control-${mode}-${temperature}.in"
final_data="${run_directory}/relaxed-control-${mode}-${temperature}.data"

if [[ ! -s "${source_data}" ]]; then
    echo "ERROR: ${source_data} is missing; finish the compressed hold first." >&2
    exit 1
fi
if [[ ! -s "${input}" ]]; then
    echo "ERROR: ${input} is missing; run build_dcz_decompression_controls.sbatch first." >&2
    exit 1
fi
if [[ -e "${final_data}" ]]; then
    echo "ERROR: ${final_data} exists; refusing to overwrite a completed control." >&2
    exit 1
fi

cd "${run_directory}"
echo "Running DCZ-4P ${mode} decompression control at ${temperature}"
echo "Prepared tree: ${prepared_root}"
echo "System: ${name}"
echo "LAMMPS: ${LAMMPS_EXE}"
echo "MPI tasks: ${SLURM_NTASKS}"

SECONDS=0
srun -n "${SLURM_NTASKS}" "${LAMMPS_EXE}" \
    -in "decompress-control-${mode}-${temperature}.in" \
    -log "log.decompress-control-${mode}-${temperature}.${SLURM_JOB_ID}.lammps"

if [[ ! -s "relaxed-control-${mode}-${temperature}.data" ]]; then
    echo "ERROR: LAMMPS returned successfully but the final control data is missing." >&2
    exit 1
fi

echo "DCZ-4P ${mode} ${temperature} control took ${SECONDS} seconds."
