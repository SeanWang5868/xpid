#!/usr/bin/env bash

# Run XPID for every X-ray and EM resolution-bin list against both the
# standard PDB mirror and the PDB-REDO mirror. Completed tasks are skipped.

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR" || exit 1

PDB_MIRROR="${PDB_MIRROR:-/vault/xhpi_pdb_July2026/pdb/mmCIF}"
REDO_MIRROR="${REDO_MIRROR:-/vault/xhpi_pdb_July2026/pdb-redo}"
OUT_DIR="${OUT_DIR:-$SCRIPT_DIR/xpid_result}"
LOG_DIR="${LOG_DIR:-$SCRIPT_DIR/xpid_logs}"
SUMMARY_FILE="${SUMMARY_FILE:-$SCRIPT_DIR/xpid_summary.tsv}"
XPID_JOBS="${XPID_JOBS:-64}"

detect_jobs() {
    local detected="$XPID_JOBS"

    if [[ ! "$detected" =~ ^[1-9][0-9]*$ ]]; then
        detected=64
    fi
    printf '%s' "$detected"
}

task_is_complete() {
    local output_name="$1"
    local log_file="$2"
    local csv_file="$OUT_DIR/${output_name}.csv"
    local metadata_file="$OUT_DIR/${output_name}_metadata.json"

    [[ -s "$csv_file" ]] || return 1
    [[ -s "$metadata_file" ]] || return 1
    [[ -s "$log_file" ]] || return 1
    grep -qE '\[INFO\] Progress[[:space:]]+:[[:space:]]+[0-9]+/[0-9]+ \(100\.0%\)' \
        "$log_file" || return 1
    grep -qE '\[OUTPUT\] (Merged result saved to:|No interactions found\.)' \
        "$log_file" || return 1
}

extract_file_count() {
    local metadata_file="$1"
    awk -F: '
        /"file_count"/ {
            value=$2
            gsub(/[^0-9]/, "", value)
            print value
            exit
        }
    ' "$metadata_file"
}

append_statistics() {
    local list_name="$1"
    local dataset="$2"
    local output_name="$3"
    local csv_file="$OUT_DIR/${output_name}.csv"
    local metadata_file="$OUT_DIR/${output_name}_metadata.json"
    local detected_files="NA"
    local result_count="NA"
    local mean_count="NA"
    local line_count

    if [[ -s "$metadata_file" ]]; then
        detected_files="$(extract_file_count "$metadata_file")"
        [[ -n "$detected_files" ]] || detected_files="NA"
    fi

    if [[ -f "$csv_file" ]]; then
        line_count="$(wc -l < "$csv_file")"
        if (( line_count > 0 )); then
            result_count=$((line_count - 1))
        else
            result_count=0
        fi
    fi

    if [[ "$detected_files" =~ ^[0-9]+$ ]] &&
       [[ "$result_count" =~ ^[0-9]+$ ]] &&
       (( detected_files > 0 )); then
        mean_count="$(awk -v total="$result_count" -v count="$detected_files" \
            'BEGIN { printf "%.6f", total / count }')"
    fi

    printf '%s\t%s\t%s\t%s\t%s\n' \
        "$list_name" "$dataset" "$detected_files" "$result_count" "$mean_count" \
        >> "$SUMMARY_FILE"
}

if ! command -v xpid >/dev/null 2>&1; then
    printf '[ERROR] xpid is not available in the current environment.\n' >&2
    printf '        Activate xhpi_env and run this script again.\n' >&2
    exit 127
fi

if [[ ! -d "$PDB_MIRROR" ]]; then
    printf '[ERROR] Standard PDB mirror not found: %s\n' "$PDB_MIRROR" >&2
    exit 1
fi

if [[ ! -d "$REDO_MIRROR" ]]; then
    printf '[ERROR] PDB-REDO mirror not found: %s\n' "$REDO_MIRROR" >&2
    exit 1
fi

shopt -s nullglob
list_files=("$SCRIPT_DIR"/Xray_bins/*.txt "$SCRIPT_DIR"/EM_bins/*.txt)
shopt -u nullglob

if (( ${#list_files[@]} == 0 )); then
    printf '[ERROR] No .txt list files found in Xray_bins or EM_bins.\n' >&2
    exit 1
fi

mkdir -p "$OUT_DIR" "$LOG_DIR"
# Rebuild the summary from all complete outputs on every invocation.
printf 'list\tdataset\tdetected_files\tresult_count\tmean_per_structure\n' > "$SUMMARY_FILE"

jobs="$(detect_jobs)"
total_tasks=$(( ${#list_files[@]} * 2 ))
task_index=0
success_count=0
skipped_count=0
failed_tasks=()
warning_tasks=()

printf '[XPID] Work directory : %s\n' "$SCRIPT_DIR"
printf '[XPID] CPU workers    : %s\n' "$jobs"
printf '[XPID] Resolution lists: %s (%s tasks)\n' "${#list_files[@]}" "$total_tasks"
printf '[XPID] Output directory: %s\n' "$OUT_DIR"
printf '[XPID] Log directory   : %s\n\n' "$LOG_DIR"

for list_file in "${list_files[@]}"; do
    list_name="$(basename "$list_file")"
    stem="${list_name%.txt}"

    for dataset in pdb pdbredo; do
        task_index=$((task_index + 1))
        output_name="${stem}_${dataset}"
        log_file="$LOG_DIR/${output_name}.log"

        if task_is_complete "$output_name" "$log_file"; then
            append_statistics "$list_name" "$dataset" "$output_name"
            success_count=$((success_count + 1))
            skipped_count=$((skipped_count + 1))
            printf '[%s/%s] SKIP  %s (%s): complete result found.\n' \
                "$task_index" "$total_tasks" "$list_name" "$dataset"
            if grep -qE '\[WARNING\] [0-9]+ files failed processing' "$log_file"; then
                warning_tasks+=("${stem}:${dataset}")
            fi
            continue
        fi

        if [[ "$dataset" == "pdb" ]]; then
            mirror_option="--pdb-mirror"
            mirror_path="$PDB_MIRROR"
        else
            mirror_option="--redo-mirror"
            mirror_path="$REDO_MIRROR"
        fi

        printf '\n[%s/%s] RUN   %s (%s)\n' \
            "$task_index" "$total_tasks" "$list_name" "$dataset" | tee "$log_file"
        printf '[%s/%s] Output: %s.csv\n' \
            "$task_index" "$total_tasks" "$output_name" | tee -a "$log_file"

        PYTHONUNBUFFERED=1 \
        OMP_NUM_THREADS=1 \
        OPENBLAS_NUM_THREADS=1 \
        MKL_NUM_THREADS=1 \
        NUMEXPR_NUM_THREADS=1 \
        xpid \
            --pdb-list "$list_file" \
            "$mirror_option" "$mirror_path" \
            --out-dir "$OUT_DIR" \
            --output-name "$output_name" \
            --file-type csv \
            --sasa \
            --provenance \
            --h-mode 4 \
            --jobs "$jobs" \
            --sym-contacts \
            -v \
            --include-coordinates \
            2>&1 | tee -a "$log_file"
        task_rc=${PIPESTATUS[0]}

        if (( task_rc == 0 )) && task_is_complete "$output_name" "$log_file"; then
            append_statistics "$list_name" "$dataset" "$output_name"
            success_count=$((success_count + 1))
            printf '[%s/%s] DONE  %s (%s)\n' \
                "$task_index" "$total_tasks" "$list_name" "$dataset" | tee -a "$log_file"
            if grep -qE '\[WARNING\] [0-9]+ files failed processing' "$log_file"; then
                warning_tasks+=("${stem}:${dataset}")
            fi
        else
            printf '%s\t%s\tNA\tNA\tNA\n' "$list_name" "$dataset" >> "$SUMMARY_FILE"
            failed_tasks+=("${stem}:${dataset}(exit=${task_rc})")
            printf '[%s/%s] FAILED %s (%s), exit=%s; continuing.\n' \
                "$task_index" "$total_tasks" "$list_name" "$dataset" "$task_rc" \
                | tee -a "$log_file" >&2
        fi
    done
done

printf '\n[XPID] All tasks finished: %s/%s completed.\n' "$success_count" "$total_tasks"
printf '[XPID] Previously complete tasks skipped: %s\n' "$skipped_count"
printf '[XPID] Summary: %s\n' "$SUMMARY_FILE"
printf '[XPID] Logs   : %s\n' "$LOG_DIR"

if (( ${#warning_tasks[@]} > 0 )); then
    printf '[XPID] Completed with per-structure warnings: %s\n' "${warning_tasks[*]}"
fi

if (( ${#failed_tasks[@]} > 0 )); then
    printf '[XPID] Failed tasks: %s\n' "${failed_tasks[*]}" >&2
    exit 1
fi

exit 0
