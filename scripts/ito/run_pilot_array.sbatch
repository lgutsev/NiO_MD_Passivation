#!/bin/bash
#SBATCH -p workq
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -c 1
#SBATCH -t 24:00:00
#SBATCH -A loni_perovsk27
#SBATCH -J ito.pilot
#SBATCH --array=0-2
#SBATCH -o ito.pilot.%A_%a.out
# Run deposition -> hold -> relax for one (study, seed) per array task, then analyze.
# Default = the pipeline smoke (3 seeds). For the pilots, after review:
#   STUDY=pilot-me-4pacz sbatch scripts/ito/run_pilot_array.sbatch
# Completed stages are never rerun or overwritten.

set -euo pipefail

module purge
module load python/3.11.5-anaconda
module load lammps/22Jul2025-r9.3-cpu
export OMP_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1
LAMMPS_EXE="$(command -v lmp_mpi)" || { echo "ERROR: lmp_mpi unavailable" >&2; exit 1; }

study="${STUDY:-smoke-me-4pacz}"
seeds=(${SEEDS:-11 12 13})
seed="${seeds[${SLURM_ARRAY_TASK_ID}]}"
prepared_root="${PREPARED_ROOT:-prepared/ito}"
run_directory="${SLURM_SUBMIT_DIR}/${prepared_root}/${study}/s${seed}"
[[ -s "${run_directory}/topology_output.lmp" ]] || { echo "ERROR: build ${run_directory} first (scripts/ito/build_pilots.sbatch)" >&2; exit 1; }

cd "${run_directory}"
echo "Study ${study} seed ${seed} in ${run_directory}; commit $(git -C "${SLURM_SUBMIT_DIR}" rev-parse --short HEAD); ${SLURM_NTASKS} MPI tasks"

run_stage () {  # name, input, product
    if [[ -e "$3" ]]; then echo "SKIP $1: $3 exists"; return; fi
    SECONDS=0
    srun -n "${SLURM_NTASKS}" "${LAMMPS_EXE}" -in "$2" -log "log.$1.${SLURM_JOB_ID}.lammps"
    echo "STAGE $1 took ${SECONDS} s"
    grep -h "Performance:" "log.$1.${SLURM_JOB_ID}.lammps" | tail -1 || true
}
run_stage deposition deposition.in deposited.data
run_stage hold hold.in held.data
run_stage relax relax.in relaxed.data

cd "${SLURM_SUBMIT_DIR}"
source .venv/bin/activate
cd src
for traj in hold relax; do
    python -m nio_md_prep.ito analyze "${run_directory}" --trajectory "${run_directory}/${traj}.lammpstrj" --last-frames 20
done
