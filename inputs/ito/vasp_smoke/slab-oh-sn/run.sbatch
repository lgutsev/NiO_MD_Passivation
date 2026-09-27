#!/bin/bash
# Adapted from lgutsev/InterfaceForge@4501e34081a63e3735336e693bea1eb40ac35cb1 launch_scripts/runvasp.sh and src/interfaceforge/vasp.py (resolve_potcar_root, assemble_potcar) (MIT)
# ITO/SAM DFT reference point: slab-oh-sn (smoke mode). Submit from this directory: sbatch run.sbatch
#SBATCH -p workq
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -c 1
#SBATCH -t 12:00:00
#SBATCH -A loni_perovsk27
#SBATCH -J ito_slab-oh-sn
#SBATCH -o vasp.cpu.%j.out

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"

# 1. Inputs must match the committed hashes (catches CRLF conversion or local edits).
sha256sum -c inputs.sha256

# 2. Assemble POTCAR from the licensed PAW tree (InterfaceForge search order); never overwrite one.
if [ ! -s POTCAR ]; then
    root=""
    for cand in "${IFACE_POTCAR_ROOT:-}" "${VASP_PP_PATH:+$VASP_PP_PATH/potpaw_PBE}" "${VASP_PP_PATH:-}" "$HOME/pot/potpaw_PBE"; do
        if [ -n "$cand" ] && [ -d "$cand" ]; then root="$cand"; break; fi
    done
    if [ -z "$root" ]; then
        echo "No PBE PAW tree: set IFACE_POTCAR_ROOT or VASP_PP_PATH" >&2; exit 2
    fi
    rm -f POTCAR.tmp
    for ds in $(awk '!/^#/ && NF >= 2 {print $2}' POTCAR.spec); do
        src="$root/$ds/POTCAR"
        if [ ! -s "$src" ]; then echo "Missing licensed POTCAR source: $src" >&2; rm -f POTCAR.tmp; exit 2; fi
        cat "$src" >> POTCAR.tmp
    done
    mv POTCAR.tmp POTCAR
    echo "POTCAR assembled from $root"
fi

# 3. Refuse to run if the POTCAR valences disagree with POTCAR.spec.
got=$(grep -E 'POMASS.*ZVAL' POTCAR | sed -E 's/.*ZVAL *= *([0-9.]+).*/\1/' | awk '{printf "%d ", $1+0.5}')
want=$(awk '!/^#/ && NF >= 4 {printf "%d ", $4}' POTCAR.spec)
if [ "$got" != "$want" ]; then echo "ZVAL mismatch: POTCAR [$got] vs POTCAR.spec [$want]" >&2; exit 3; fi
nel=$(paste <(awk '!/^#/ && NF >= 4 {print $3}' POTCAR.spec) <(echo $got | tr ' ' '\n') | awk '{s += $1*$2} END {print s}')
if [ "$nel" != "1464" ]; then echo "NELECT $nel differs from expected 1464" >&2; exit 3; fi
grep -E 'TITEL|ENMAX' POTCAR | tee potcar_titles.txt

# 4. Run (LONI VASP module convention from InterfaceForge runvasp.sh).
module purge
VASP_MODULE="${IFACE_VASP_MODULE:-vasp6/6.5.1-cpu}"
module load "$VASP_MODULE"
export OMP_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1
{
    echo "case=slab-oh-sn"; echo "mode=smoke"; echo "job=${SLURM_JOB_ID:-none}"; echo "host=$(hostname)"
    echo "module=$VASP_MODULE"; echo "start=$(date -Is)"
    echo "repo_commit=$(git -C "$PWD" rev-parse HEAD 2>/dev/null || echo unknown)"
    sha256sum POSCAR INCAR KPOINTS POTCAR.spec POTCAR
} > run_provenance.txt
SECONDS=0
srun -n "${SLURM_NTASKS:-64}" vasp_gam
echo "end=$(date -Is)" >> run_provenance.txt
echo "took $SECONDS sec." | tee -a run_provenance.txt
