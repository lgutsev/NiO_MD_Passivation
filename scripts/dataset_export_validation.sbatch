#!/bin/bash
#SBATCH -p single
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 08:00:00
#SBATCH -A loni_perovsk27
#SBATCH -J nio.dataset
#SBATCH -o nio.dataset.%j.out
#
# Dataset export validation on real VASP trees: scan -> export -> split -> audit -> spot check.
#
# The account above is an example; override it on the command line:
#     sbatch -A <your_allocation> scripts/dataset_export_validation.sbatch
#
# Environment (export before sbatch or pass with --export=ALL,VAR=...):
#   DATASET_ROOT     required; space-separated roots, each ALIAS=PATH or PATH (alias = basename).
#                    The alias is part of every run id, so keep it stable between exports.
#   DATASET_OUT      output directory (default dataset_validation/<UTC timestamp>); must not exist
#                    or be empty (nothing is overwritten; no --force is ever passed here).
#   LINEAGE_POLICY   optional run|parent|depth:N for runs without manifest/inventory lineage.
#   INVENTORY        optional selection/lineage TOML (--inventory).
#   SETTINGS_OVERRIDES optional reviewed [[equivalence]] TOML (--settings-overrides).
#   POOL             optional reference-settings pool id (needed when the scan reports several).
#   EXTRA_ARGS       optional extra scan/export flags (e.g. "--include-stress").
#   SPLIT_SEED       split seed (default 11).
#   SPOTCHECK_N      frames compared with ASE's own VASP reader (default 50).
#
# Nothing is written into the calculation roots. The job never commits, pushes or publishes;
# LONI checkouts are pull-only.

set -euo pipefail

module purge
module load python/3.11.5-anaconda

cd "${SLURM_SUBMIT_DIR}"
source .venv/bin/activate

if ! command -v nio-md-prep >/dev/null 2>&1; then
    echo "ERROR: nio-md-prep is unavailable after activating .venv." >&2
    exit 1
fi
if ! python -c "import ase, numpy" >/dev/null 2>&1; then
    echo "ERROR: ASE/numpy missing; install the extra: pip install -e '.[dataset]'" >&2
    exit 1
fi
if [[ -z "${DATASET_ROOT:-}" ]]; then
    echo "ERROR: set DATASET_ROOT to one or more ALIAS=PATH calculation roots." >&2
    exit 1
fi

out="${DATASET_OUT:-dataset_validation/$(date -u +%Y%m%dT%H%M%SZ)}"
split_seed="${SPLIT_SEED:-11}"
spot_n="${SPOTCHECK_N:-50}"

read -r -a roots <<< "${DATASET_ROOT}"
root_args=()
verify_args=()
for root in "${roots[@]}"; do
    if [[ "${root}" == *=* ]]; then
        alias="${root%%=*}"
        path="${root#*=}"
    else
        path="${root}"
        alias="$(basename "${path}")"
    fi
    if [[ ! -d "${path}" ]]; then
        echo "ERROR: root ${path} is not a directory." >&2
        exit 1
    fi
    root_args+=(--root "${alias}=${path}")
    verify_args+=(--verify-sources "${alias}=${path}")
done

common=("${root_args[@]}")
[[ -n "${LINEAGE_POLICY:-}" ]] && common+=(--lineage-policy "${LINEAGE_POLICY}")
[[ -n "${INVENTORY:-}" ]] && common+=(--inventory "${INVENTORY}")
[[ -n "${SETTINGS_OVERRIDES:-}" ]] && common+=(--settings-overrides "${SETTINGS_OVERRIDES}")
if [[ -n "${EXTRA_ARGS:-}" ]]; then
    read -r -a extra <<< "${EXTRA_ARGS}"
    common+=("${extra[@]}")
fi
pool_args=()
[[ -n "${POOL:-}" ]] && pool_args=(--pool "${POOL}")

mkdir -p "${out}"
echo "nio-md-prep $(nio-md-prep --version 2>/dev/null || true) git $(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
echo "roots: ${root_args[*]}"
echo "output: ${out}"

echo "== 1/5 scan (accounting only)"
nio-md-prep dataset scan "${common[@]}" --output "${out}/scan"

echo "== 2/5 export"
if ! nio-md-prep dataset export "${common[@]}" "${pool_args[@]}" --output "${out}/export"; then
    echo "ERROR: export refused; read ${out}/scan/audit.md and exclusions.csv (several pools -> set POOL," >&2
    echo "       unresolved lineage -> set LINEAGE_POLICY or INVENTORY, no exportable frame -> see counts)." >&2
    exit 2
fi

echo "== 3/5 split (lineage-group train/valid/test)"
nio-md-prep dataset split "${out}/export" --output "${out}/split" --seed "${split_seed}"

echo "== 4/5 audit (hashes, accounting, ASE round trip, leakage, source re-hash)"
nio-md-prep dataset audit "${out}/export" --split "${out}/split" "${verify_args[@]}" --output "${out}/audit"

echo "== 5/5 spot check against the original VASP files with ASE's reader"
python scripts/dataset_spotcheck.py "${out}/export" --n "${spot_n}" --seed "${split_seed}" \
    --json "${out}/spotcheck.json"

echo "Validation complete: ${out}"
ls -l "${out}/export"
