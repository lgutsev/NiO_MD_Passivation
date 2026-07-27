#!/bin/bash
#SBATCH -p single
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 02:00:00
#SBATCH -A loni_perovsk27
#SBATCH -J nio.lego.build
#SBATCH --array=0
#SBATCH -o nio.lego.build.%A_%a.out

set -euo pipefail

module purge
module load python/3.11.5-anaconda

cd "${SLURM_SUBMIT_DIR}"
source .venv/bin/activate

for executable in python packmol nio-md-prep; do
    if ! command -v "${executable}" >/dev/null 2>&1; then
        echo "ERROR: ${executable} is not available after module and virtual-environment setup." >&2
        exit 1
    fi
done

source_root="${SOURCE_ROOT:-prepared}"
prepared_root="${PREPARED_ROOT:-prepared-lego}"
primary="${source_root}/me-4pacz-alone/held-300K.data"
coverage_map="${COVERAGE_MAP:-${source_root}/me-4pacz-alone/coverage-analysis-hold-300K/coverage_probability.npz}"

if [[ ! -s "${primary}" ]]; then
    echo "ERROR: ${primary} is missing; finish the original Me-4PACz 300 K hold first." >&2
    exit 1
fi
if [[ ! -s "${coverage_map}" ]]; then
    echo "ERROR: ${coverage_map} is missing; run coverage analysis on the original Me-4PACz hold first." >&2
    exit 1
fi

systems=(
    me-4pacz-then-dcz-4p
    me-4pacz-then-meo-2pacz
    me-4pacz-then-meo-4padbc
)

base_name="${systems[${SLURM_ARRAY_TASK_ID}]}"
lego_suffix="${LEGO_SUFFIX:--lego-seeded}"
name="${base_name}${lego_suffix}"
output="${prepared_root}/${name}"

mkdir -p "${prepared_root}"

if [[ -e "${output}" ]]; then
    echo "ERROR: ${output} already exists; remove or archive it before rebuilding." >&2
    exit 1
fi

echo "BUILDING GAP-SEEDED SEQUENTIAL STAGE 2: ${name}"
echo "Stage-1 source: ${primary}"
echo "Coverage map: ${coverage_map}"
echo "Lateral confinement: ${LEGO_LATERAL_WALL:-0}"

build_args=(
    "studies/${base_name}.toml"
    --primary-final "${primary}"
    --coverage-map "${coverage_map}"
    --output "${output}"
    --occupancy-threshold "${OCCUPANCY_THRESHOLD:-0.20}"
    --packing-inset "${PACKING_INSET:-4.0}"
    --minimum-gap-fraction "${MINIMUM_GAP_FRACTION:-0.12}"
    --minimum-tunnel-width "${MINIMUM_TUNNEL_WIDTH:-12.0}"
    --deposition-clearance "${LEGO_DEPOSITION_CLEARANCE:-30.0}"
    --wall-strength "${LEGO_WALL_STRENGTH:-0.50}"
    --wall-cutoff "${LEGO_WALL_CUTOFF:-3.0}"
)
if [[ "${LEGO_LATERAL_WALL:-0}" == "1" ]]; then
    build_args+=(--lateral-wall)
fi
if [[ -n "${PACKMOL_SEED:-}" ]]; then
    build_args+=(--packmol-seed "${PACKMOL_SEED}")
fi
if [[ -n "${VELOCITY_SEED:-}" ]]; then
    build_args+=(--velocity-seed "${VELOCITY_SEED}")
fi

nio-md-prep prepare-lego-stage2 "${build_args[@]}"

for required in lego_plan.json lego_tunnel_profile.csv lego-stage1-shifted.data topology_output.lmp force_field_settings.lmp deposition.in hold-300K.in; do
    if [[ ! -s "${output}/${required}" ]]; then
        echo "ERROR: ${output}/${required} was not produced." >&2
        exit 1
    fi
done

echo "BUILT: ${name}"
echo "Review ${output}/lego_plan.json before launching deposition."
echo "No deposition or hold job was submitted."
