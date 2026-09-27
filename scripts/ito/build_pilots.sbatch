#!/bin/bash
#SBATCH -p single
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 02:00:00
#SBATCH -A loni_perovsk27
#SBATCH -J ito.build
#SBATCH -o ito.build.%j.out
# Build every ITO pilot study x Packmol seed into prepared/ito/<study>/s<seed>.
# Override with STUDIES="a b" SEEDS="11 12" sbatch scripts/ito/build_pilots.sbatch

set -euo pipefail

module purge
module load python/3.11.5-anaconda
module load lammps/22Jul2025-r9.3-cpu

cd "${SLURM_SUBMIT_DIR}"
source .venv/bin/activate

for executable in python packmol lmp_mpi; do
    command -v "${executable}" >/dev/null 2>&1 || { echo "ERROR: ${executable} unavailable" >&2; exit 1; }
done
echo "Commit: $(git rev-parse --short HEAD)"

studies=(${STUDIES:-smoke-me-4pacz pilot-me-4pacz pilot-meo-2pacz-me-4pacz-1to1 pilot-me-4pacz-monolayer})
seeds=(${SEEDS:-11 12 13})
prepared_root="${PREPARED_ROOT:-prepared/ito}"

for study in "${studies[@]}"; do
    for seed in "${seeds[@]}"; do
        out="${prepared_root}/${study}/s${seed}"
        if [[ -e "${out}/topology_output.lmp" ]]; then
            echo "SKIP ${out}: already built"; continue
        fi
        echo "BUILD ${study} seed ${seed} -> ${out}"
        (cd src && python -m nio_md_prep.ito build-pilot "../studies/ito/${study}.toml" \
            --output "../${out}" --packmol-seed "${seed}" --velocity-seed "$((100 + seed))")
    done
done
