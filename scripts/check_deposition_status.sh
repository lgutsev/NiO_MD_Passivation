#!/bin/bash
# Read-only completion check for deposition reruns.

set -u

if (( $# )); then
    roots=("$@")
else
    roots=(prepared-rerun prepared-rerun2)
fi

latest_log() {
    find "$1" -maxdepth 1 -type f -name 'log.deposition.*.lammps' \
        -printf '%T@ %p\n' 2>/dev/null |
        sort -nr | head -n 1 | cut -d' ' -f2-
}

target_step() {
    awk '$1=="run" && $2~/^[0-9]+$/ && $3=="start" && $5=="stop" {
        print $2
        exit
    }' "$1"
}

last_step() {
    awk '$1~/^[0-9]+$/ && $2~/^[-+0-9.eE]+$/ {
        if ($1>last) last=$1
        found=1
    } END {
        if (found) print last
        else print 0
    }' "$1"
}

complete_atom_section() {
    local declared actual
    declared="$(awk '$1~/^[0-9]+$/ && $2=="atoms" {print $1; exit}' "$1")"
    actual="$(awk '
        /^Atoms([[:space:]]|$)/ {inside=1; next}
        inside && /^[[:alpha:]]/ {inside=0}
        inside && $1~/^[0-9]+$/ {count++}
        END {print count+0}
    ' "$1")"
    [[ -n "${declared}" && "${actual}" == "${declared}" ]]
}

stable_file() {
    local first second
    first="$(stat -c '%s' "$1" 2>/dev/null)" || return 1
    sleep 2
    second="$(stat -c '%s' "$1" 2>/dev/null)" || return 1
    [[ "${first}" -gt 0 && "${first}" == "${second}" ]]
}

safe_jobs=()
printf '%-22s %-38s %-23s %-16s %-12s\n' \
    "ROOT" "SYSTEM" "STATUS" "STEP/TARGET" "JOB"

for root in "${roots[@]}"; do
    while IFS= read -r -d '' directory; do
        [[ -s "${directory}/deposition.in" ]] || continue
        name="$(basename "${directory}")"
        target="$(target_step "${directory}/deposition.in")"
        [[ -n "${target}" ]] || target=0
        log="$(latest_log "${directory}")"
        step=0
        job="-"
        state=""
        if [[ -n "${log}" ]]; then
            step="$(last_step "${log}")"
            job="$(basename "${log}")"
            job="${job#log.deposition.}"
            job="${job%.lammps}"
            state="$(squeue -h -j "${job}" -o '%T' 2>/dev/null | head -n 1)"
        fi

        if [[ -s "${directory}/deposited.data" ]] &&
           (( step >= target )) &&
           stable_file "${directory}/deposited.data" &&
           complete_atom_section "${directory}/deposited.data"; then
            status="SAFE TO CANCEL"
            [[ -n "${state}" && "${job}" != "-" ]] && safe_jobs+=("${job}")
        elif [[ -s "${directory}/deposited.data" ]]; then
            status="FINALIZING/CHECK"
        elif [[ -n "${state}" ]]; then
            status="${state}"
        elif [[ -n "${log}" ]] &&
             grep -Eqi 'ERROR|MPI_Abort|CANCELLED|Killed|missing on proc' "${log}"; then
            status="FAILED"
        elif [[ -n "${log}" ]]; then
            status="INCOMPLETE"
        else
            status="NOT STARTED"
        fi

        printf '%-22s %-38s %-23s %7s/%-8s %-12s\n' \
            "${root}" "${name}" "${status}" "${step}" "${target}" "${job}"
    done < <(find "${root}" -mindepth 1 -maxdepth 1 -type d -print0 2>/dev/null | sort -z)
done

echo
echo "SAFE TO CANCEL = target step reached + stable deposited.data + complete atom count."
if (( ${#safe_jobs[@]} )); then
    printf 'Suggested command after reviewing the table: scancel'
    printf ' %q' "${safe_jobs[@]}"
    echo
fi
