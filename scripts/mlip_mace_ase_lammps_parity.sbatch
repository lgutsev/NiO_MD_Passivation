#!/bin/bash
#SBATCH -p gpu2
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --gres=gpu:1
#SBATCH -t 02:00:00
#SBATCH -A loni_perovsk27
#SBATCH -J nio.mlip-parity
#SBATCH -o nio.mlip-parity.%j.out
#
# MACE ASE <-> LAMMPS/ML-IAP single-point parity (Phase 2 preferred acceptance).
#
# The partition, GPU request and account above are LONI examples; override them on the command
# line, e.g. `sbatch -A <allocation> -p <gpu partition> scripts/mlip_mace_ase_lammps_parity.sbatch`.
#
# Environment (export before sbatch or pass with --export=ALL,VAR=...):
#   MACE_MODEL     required; MACE checkpoint (.model)
#   STRUCTURE      required; structure file (POSCAR, extxyz, ...)
#   ELEMENTS       required; comma list covered by the model, e.g. Ni,O,C,H,N,P
#   PBC            optional per-axis override, e.g. T,T,T (VASP slabs are 3-D periodic)
#   DEVICE         cuda (default) or cpu
#   PRECISION      float64 (default) or float32; the export is created with the same dtype
#   PARITY_ENV     optional conda environment with mace-torch and a LAMMPS Python module built with
#                  ML-IAP + PYTHON + KOKKOS (activated with `conda activate "${PARITY_ENV}"`)
#   PARITY_OUT     output directory (default mlip_parity/<UTC timestamp>); must be new or empty
#   EXTRA_ARGS     optional extra arguments (e.g. "--stress" for fully periodic cells)
#
# The job exports the model with mace_create_lammps_model --format mliap, evaluates it through the
# ASE bridge and the LAMMPS/ML-IAP bridge, and writes parity_report.json plus a CANDIDATE tolerance
# file. Candidates are measurements for human review, never acceptance criteria. The job never
# commits, pushes or publishes; LONI checkouts are pull-only.

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"
if [[ -n "${PARITY_ENV:-}" ]]; then
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${PARITY_ENV}"
elif [[ -f .venv/bin/activate ]]; then
    source .venv/bin/activate
fi

for var in MACE_MODEL STRUCTURE ELEMENTS; do
    if [[ -z "${!var:-}" ]]; then
        echo "ERROR: set ${var}." >&2
        exit 1
    fi
done
if ! python -c "import mace, lammps" >/dev/null 2>&1; then
    echo "ERROR: mace-torch and a LAMMPS Python module are required (set PARITY_ENV)." >&2
    exit 1
fi

out="${PARITY_OUT:-mlip_parity/$(date -u +%Y%m%dT%H%M%SZ)}"
args=(--model "${MACE_MODEL}" --structure "${STRUCTURE}" --elements "${ELEMENTS}"
      --device "${DEVICE:-cuda}" --precision "${PRECISION:-float64}" --export --output "${out}")
[[ -n "${PBC:-}" ]] && args+=(--pbc "${PBC}")
if [[ -n "${EXTRA_ARGS:-}" ]]; then
    read -r -a extra <<< "${EXTRA_ARGS}"
    args+=("${extra[@]}")
fi

echo "git $(git rev-parse --short HEAD 2>/dev/null || echo unknown); python $(python --version 2>&1)"
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv 2>/dev/null || true
python scripts/mlip_mace_ase_lammps_parity.py "${args[@]}"
