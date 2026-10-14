#!/usr/bin/env bash

# Return the scheduler/job identifier when available, else a local fallback.
get_job_id() {
    local job_id=""
    local array_index=""

    if [ -n "${AWS_BATCH_JOB_ID:-}" ]; then
        job_id="$AWS_BATCH_JOB_ID"
        array_index="${AWS_BATCH_JOB_ARRAY_INDEX:-}"
    elif [ -n "${LSB_JOBID:-}" ]; then
        job_id="$LSB_JOBID"
        array_index="${LSB_JOBINDEX:-}"
    elif [ -n "${SLURM_JOB_ID:-}" ]; then
        job_id="$SLURM_JOB_ID"
        array_index="${SLURM_ARRAY_TASK_ID:-}"
    elif [ -n "${PBS_JOBID:-}" ]; then
        job_id="$PBS_JOBID"
        array_index="${PBS_ARRAY_INDEX:-}"
    elif [ -n "${JOB_ID:-}" ]; then
        job_id="$JOB_ID"
        array_index="${SGE_TASK_ID:-}"
    else
        job_id="local_$$"
    fi

    if [ -n "$array_index" ]; then
        echo "${job_id}_${array_index}"
    else
        echo "$job_id"
    fi
}

# Force line-buffering for every subprocess in the current shell so long
# running tools stream stdout/stderr to the log in real time instead of
# sitting in libc's 8KB buffers or Python's stdout buffer.
#
# LD_PRELOAD is inherited by all descendants, so grandchildren (e.g. eddy
# spawned by trac-preproc spawned by trac-all) get the setting automatically.
# PYTHONUNBUFFERED covers Python's own layer on top of libc.
#
# Call this once after the rule's `exec > log_file` redirection, before any
# long-running tool is invoked.
enable_line_buffering() {
    local libstdbuf=/usr/libexec/coreutils/libstdbuf.so
    if [ -f "$libstdbuf" ]; then
        export STDBUF_O=L STDBUF_E=L
        export LD_PRELOAD="${LD_PRELOAD:+$LD_PRELOAD:}$libstdbuf"
    fi
    export PYTHONUNBUFFERED=1
}

# Log scratch usage (GB), report peak RSS, and remove the scratch directory.
# A 120-second timeout guards against frozen/stale network-attached scratch
# filesystems that would otherwise block the job indefinitely on rm -rf.
scratch_cleanup() {
    local scratch_path="${1:-}"

    # Report peak RSS across all stages and total before scratch cleanup
    report_total_peak_mem

    if [ -n "$scratch_path" ] && [ -d "$scratch_path" ]; then
        echo "Cleaning up scratch directory: $scratch_path"
        local scratch_size_gb
        scratch_size_gb=$(du -sk "$scratch_path" 2>/dev/null | awk '{printf "%.2f", $1/1024/1024}')
        if [ -n "$scratch_size_gb" ]; then
            echo "Scratch usage before cleanup: ${scratch_size_gb} GB"
        fi
    elif [ -n "$scratch_path" ]; then
        echo "Cleaning up scratch directory: $scratch_path (not found or already removed)"
    fi

    if [ -z "$scratch_path" ]; then
        return 0
    fi

    local rm_status=0
    if command -v timeout >/dev/null 2>&1; then
        timeout 600 rm -rf "$scratch_path" || rm_status=$?
    else
        rm -rf "$scratch_path" || rm_status=$?
    fi

    if [ "$rm_status" -eq 0 ]; then
        echo "Scratch cleanup complete."
    else
        echo "WARNING: scratch cleanup timed out or failed for: $scratch_path" >&2
    fi
}

# Verify a real GPU is actually present when -platf 1 (CUDA) is requested.
# A no-op for any other platf value (CPU mode never needs a real GPU).
#
# NiftyReg's own reg_aladin/reg_f3d/reg_resample/seg_GIF already refuse to
# proceed in this situation ("CUDA driver is a stub library"), but that
# error surfaces as a raw C++ exception deep inside the binary, potentially
# after real setup work (T1 prep, healthy-cohort materialization, etc.) has
# already run. This catches the same condition immediately, in plain
# English, before any of that starts. Checks nvidia-smi rather than just
# device-file presence -- our bundled libcuda.so.1 stub (see GIF_200826.sh)
# exists precisely so the binary can still dynamic-link and start on a
# non-GPU node, so its mere presence proves nothing about a real GPU.
require_gpu_if_platf1() {
    local platf="$1"
    if [ "$platf" != "1" ]; then
        return 0
    fi
    if ! command -v nvidia-smi >/dev/null 2>&1 || ! nvidia-smi -L 2>/dev/null | grep -q '^GPU'; then
        echo "ERROR: GPU is not available on this node (requested -platf 1 / --gpu, but nvidia-smi found no real GPU). Stopping." >&2
        exit 1
    fi
}

# Track multi-stage peak memory within a script run.
SCRIPT_PEAK_MEM_FILE="${SCRIPT_PEAK_MEM_FILE:-}"

_record_stage_peak_mem() {
    local label="$1"
    local mb="$2"
    local elapsed="$3"
    local cpu="$4"
    local kb="$5"

    if [ -z "${SCRIPT_PEAK_MEM_FILE:-}" ] || [ ! -f "${SCRIPT_PEAK_MEM_FILE:-}" ]; then
        SCRIPT_PEAK_MEM_FILE=$(mktemp -t leukoquant_peak_mem.XXXXXX 2>/dev/null || mktemp)
        export SCRIPT_PEAK_MEM_FILE
    fi
    echo "${label}|${mb}|${elapsed:-N/A}|${cpu:-N/A}|${kb}" >> "$SCRIPT_PEAK_MEM_FILE"
}

# Report a summary of all executed stages and the overall peak RSS.
# Safe and idempotent: can be called multiple times (e.g. at script end and in EXIT trap)
# without duplicate printing.
report_total_peak_mem() {
    if [ "${_REPORT_TOTAL_PEAK_MEM_DONE:-0}" -eq 1 ]; then
        return 0
    fi
    _REPORT_TOTAL_PEAK_MEM_DONE=1

    local mem_file="${SCRIPT_PEAK_MEM_FILE:-}"
    if [ -z "$mem_file" ] || [ ! -f "$mem_file" ] || [ ! -s "$mem_file" ]; then
        [ -n "$mem_file" ] && rm -f "$mem_file"
        return 0
    fi

    local count=0
    local max_kb=0
    local lines=()
    while IFS='|' read -r r_label r_mb r_elapsed r_cpu r_kb; do
        [ -z "$r_kb" ] && continue
        count=$((count + 1))
        if [ "$r_kb" -gt "$max_kb" ]; then
            max_kb="$r_kb"
        fi
        if [ "$r_elapsed" != "N/A" ]; then
            lines+=("  * ${r_label:-stage}: ${r_mb} MB (Elapsed: ${r_elapsed}, CPU: ${r_cpu})")
        else
            lines+=("  * ${r_label:-stage}: ${r_mb} MB")
        fi
    done < "$mem_file"

    rm -f "$mem_file"
    unset SCRIPT_PEAK_MEM_FILE

    if [ "$count" -eq 0 ] || [ "$max_kb" -le 0 ]; then
        return 0
    fi

    local total_mb=$((max_kb / 1024))
    [ "$total_mb" -eq 0 ] && total_mb="<1"

    if [ "$count" -gt 1 ]; then
        local summary_hdr="[mem] --- Peak RSS Summary (${count} stages) ---"
        echo "$summary_hdr"
        echo "$summary_hdr" >&2
        for line in "${lines[@]}"; do
            echo "$line"
            echo "$line" >&2
        done
    fi

    local total_msg="[mem] Total Peak RSS: ${total_mb} MB"
    echo "$total_msg"
    echo "$total_msg" >&2
}

# ---------------------------------------------------------------------------
# Loud failure reporting.
#
# A tool that runs out of memory under the SGE h_vmem limit (address-space cap)
# fails with an allocation error, or is killed (status 137). If a script hides
# that (output sent to /dev/null, "|| true") the rule still "succeeds" and the
# stale or empty outputs look like a result. These helpers make the failure
# visible as text in the log and error files and as a non-zero exit status:
#
#   "OUT OF MEMORY: <stage> ..."   the step ran out of memory
#   "FAILED: <stage> ..."          any other failure
#   "[oom] <stage> ..."            written by run_with_peak_mem for a stage that
#                                  was killed or hit an allocation failure
#   "[mem] Peak RSS (<stage>): ... (stage FAILED, exit N)"
# ---------------------------------------------------------------------------

# Text that means "out of memory": glibc and C++ allocation failures, Python's
# MemoryError, the kernel OOM killer, and our own "OUT OF MEMORY" line.
_OOM_PATTERN="can't allocate memory|cannot allocate memory|std::bad_alloc|MemoryError|out of memory|oom-kill|Killed"
# Stricter variant for a stage that reported success (status 0): a bare word
# like "Killed" is not enough evidence there.
_OOM_PATTERN_STRICT="can't allocate memory|cannot allocate memory|std::bad_alloc|MemoryError|out of memory|oom-kill"

# First line of stdin that matches the out-of-memory text (exit 0), else exit 1.
# Usage: first_oom_line [strict] < text
first_oom_line() {
    local pattern="$_OOM_PATTERN"
    [ "${1:-}" = "strict" ] && pattern="$_OOM_PATTERN_STRICT"
    local line
    line=$(grep -Eim1 "$pattern") || return 1
    printf '%s\n' "$line"
}

# The file this shell's stderr is redirected to (empty if stderr is not a file).
# $LEUKOQUANT_STDERR_FILE overrides the lookup (needed where there is no /proc).
_stderr_file() {
    if [ -n "${LEUKOQUANT_STDERR_FILE:-}" ] && [ -f "$LEUKOQUANT_STDERR_FILE" ]; then
        printf '%s\n' "$LEUKOQUANT_STDERR_FILE"
        return 0
    fi
    local target
    target=$(readlink "/proc/$$/fd/2" 2>/dev/null) || return 0
    [ -f "$target" ] && printf '%s\n' "$target"
    return 0
}

# Size in bytes of a file (portable: no stat flags).
_file_size() {
    wc -c < "$1" 2>/dev/null | tr -d ' ' || echo 0
}

# Print a classified one-line verdict for a failed step (stdout and stderr, like
# the [mem] lines) and return 0 when it was an out-of-memory failure, 1 when not.
# Usage: report_step_failure "stage" exit_status "captured output text"
report_step_failure() {
    local stage="$1" status="${2:-1}" text="${3:-}"
    local hit="" msg
    hit=$(printf '%s\n' "$text" | first_oom_line) || hit=""
    if [ -n "$hit" ] || [ "$status" = "137" ]; then
        msg="OUT OF MEMORY: ${stage} ran out of memory (exit ${status}${hit:+; matched: ${hit}})"
        echo "$msg"
        echo "$msg" >&2
        return 0
    fi
    msg="FAILED: ${stage} (exit ${status})"
    echo "$msg"
    echo "$msg" >&2
    return 1
}

# Skull-strip with mri_synthstrip. On failure print the tool output, the legacy
# "ERROR: mri_synthstrip failed on <what>" line and a classified verdict, and
# return the tool's status. The caller decides whether to stop (they all do).
# Usage: run_synthstrip "T1" input.nii.gz output.nii.gz
run_synthstrip() {
    local what="$1" input="$2" output="$3"
    local text rc=0
    text=$(mri_synthstrip -i "$input" -o "$output" 2>&1) || rc=$?
    [ "$rc" -eq 0 ] && return 0
    echo "ERROR: mri_synthstrip failed on ${what}. Output:" >&2
    printf '%s\n' "$text" >&2
    report_step_failure "mri_synthstrip (${what})" "$rc" "$text" || true
    return "$rc"
}

# ERR trap handler: say which command failed instead of dying silently.
# Usage: trap 'report_failed_command "$?" "$LINENO" "$BASH_COMMAND"' ERR
report_failed_command() {
    local status="$1" line="$2" cmd="$3"
    report_step_failure "${0##*/} line ${line}: ${cmd}" "$status" "" || true
}

# After a stage has run: on failure of a stage whose stderr (or status 137) points to
# memory print the "[oom]" line; on a stage that exited 0 but wrote an allocation
# message to its stderr, flag the swallowed failure. Text only; the status is not changed.
# Usage: _finish_stage_report stage status peak_mb stderr_file stderr_start_offset
_finish_stage_report() {
    local stage="$1" status="$2" peak_mb="${3:-unknown}" err_file="${4:-}" err_start="${5:-0}"
    local hit="" text="" msg

    if [ "$status" -eq 0 ]; then
        # Success: a swallowed allocation failure still deserves a flag (strict text only).
        if [ -n "$err_file" ] && [ -f "$err_file" ]; then
            hit=$(tail -c +"$((err_start + 1))" "$err_file" 2>/dev/null | first_oom_line strict) || hit=""
        fi
        if [ -n "$hit" ]; then
            msg="[oom] ${stage}: exit 0 but the log contains an out-of-memory message (${hit}); the stage may have hidden a failure"
            echo "$msg"
            echo "$msg" >&2
        fi
        return 0
    fi

    if [ -n "$err_file" ] && [ -f "$err_file" ]; then
        text=$(tail -c +"$((err_start + 1))" "$err_file" 2>/dev/null | tail -c 20000)
    fi
    hit=$(printf '%s\n' "$text" | first_oom_line) || hit=""
    if [ -n "$hit" ] || [ "$status" -eq 137 ]; then
        local why="exit ${status}"
        [ "$status" -eq 137 ] && why="exit 137 (SIGKILL: memory limit, time limit or the OOM killer)"
        [ -n "$hit" ] && why="${why}; matched: ${hit}"
        msg="[oom] ${stage}: out of memory (${why}); peak RSS ${peak_mb} MB; address-space limit $(ulimit -v 2>/dev/null || echo unknown) KB"
        echo "$msg"
        echo "$msg" >&2
    fi
}

# Execute a command and log its kernel-measured peak RSS, elapsed time,
# and CPU % using /usr/bin/time (Linux GNU time or macOS BSD time).
#
# Unlike background polling loops, this has ZERO polling overhead, spawns
# zero background processes, and incurs zero cluster slowdown. The Linux
# kernel tracks peak RSS natively in struct rusage via the wait4 syscall.
#
# Usage:
#   run_with_peak_mem "stage_name" command arg1 arg2 ...
#   run_with_peak_mem "stage_name" bash -c "complex shell pipeline"
run_with_peak_mem() {
    local label="${1:-}"
    shift
    local status=0
    local tmp_time
    tmp_time=$(mktemp)

    if command -v /usr/bin/time >/dev/null 2>&1; then
        if /usr/bin/time -f "%M" true >/dev/null 2>&1; then
            # GNU time (Linux, Singularity containers, SGE cluster nodes)
            # Remember where this stage's stderr starts, to scan only its own text afterwards.
            local err_file err_start=0
            err_file=$(_stderr_file)
            [ -n "$err_file" ] && err_start=$(_file_size "$err_file")
            /usr/bin/time -o "$tmp_time" -f "%M %E %P" "$@" || status=$?
            # Belt and braces: GNU time names a SIGKILL in its output even if its exit status were 0.
            if [ "$status" -eq 0 ] && grep -q "terminated by signal 9" "$tmp_time" 2>/dev/null; then
                status=137
            fi
            local peak_kb elapsed cpu
            # On a signal GNU time first writes "Command terminated by signal N": the figures are on the last line.
            read -r peak_kb elapsed cpu <<< "$(tail -n 1 "$tmp_time" 2>/dev/null)"
            rm -f "$tmp_time"
            local peak_mb="unknown"
            if [[ "$peak_kb" =~ ^[0-9]+$ ]] && [ "$peak_kb" -gt 0 ]; then
                peak_mb=$((peak_kb / 1024))
                [ "$peak_mb" -eq 0 ] && peak_mb="<1"
                local msg="[mem] Peak RSS${label:+ ($label)}: ${peak_mb} MB (Elapsed: ${elapsed:-N/A}, CPU: ${cpu:-N/A})"
                # A failed stage must not read as a normal measurement.
                [ "$status" -ne 0 ] && msg="${msg} (stage FAILED, exit ${status})"
                echo "$msg"
                echo "$msg" >&2
                _record_stage_peak_mem "$label" "$peak_mb" "${elapsed:-N/A}" "${cpu:-N/A}" "$peak_kb"
            elif [ "$status" -ne 0 ]; then
                local msg="[mem] Peak RSS${label:+ ($label)}: unknown (stage FAILED, exit ${status})"
                echo "$msg"
                echo "$msg" >&2
            fi
            _finish_stage_report "${label:-stage}" "$status" "$peak_mb" "$err_file" "$err_start"
            return "$status"
        elif [ "$(uname -s)" = "Darwin" ]; then
            # BSD time (macOS local workstations)
            local err_file err_start=0
            err_file=$(_stderr_file)
            [ -n "$err_file" ] && err_start=$(_file_size "$err_file")
            /usr/bin/time -l -o "$tmp_time" "$@" || status=$?
            local peak_bytes
            peak_bytes=$(awk '/maximum resident set size/ {print $1}' "$tmp_time" 2>/dev/null || true)
            rm -f "$tmp_time"
            local peak_mb="unknown"
            if [[ "$peak_bytes" =~ ^[0-9]+$ ]] && [ "$peak_bytes" -gt 0 ]; then
                local peak_kb=$((peak_bytes / 1024))
                peak_mb=$((peak_bytes / 1024 / 1024))
                [ "$peak_mb" -eq 0 ] && peak_mb="<1"
                local msg="[mem] Peak RSS${label:+ ($label)}: ${peak_mb} MB"
                [ "$status" -ne 0 ] && msg="${msg} (stage FAILED, exit ${status})"
                echo "$msg"
                echo "$msg" >&2
                _record_stage_peak_mem "$label" "$peak_mb" "N/A" "N/A" "$peak_kb"
            elif [ "$status" -ne 0 ]; then
                local msg="[mem] Peak RSS${label:+ ($label)}: unknown (stage FAILED, exit ${status})"
                echo "$msg"
                echo "$msg" >&2
            fi
            _finish_stage_report "${label:-stage}" "$status" "$peak_mb" "$err_file" "$err_start"
            return "$status"
        fi
    fi

    rm -f "$tmp_time"
    "$@"
}

# Backward compatibility no-ops / hooks:
start_peak_mem_monitor() {
    :
}
stop_peak_mem_monitor() {
    report_total_peak_mem
}

# Generate a temporary path by splicing "_tmp" before the file extension.
# Dynamically handles compound extensions (.nii.gz, .tar.gz) and single
# extensions (.mgz, .csv), so callers writing to a temp path first and
# renaming into place (the atomic-write pattern used elsewhere in this
# codebase) don't need to special-case compressed formats themselves.
make_tmp_path() {
    local target="$1"
    local dir filename base ext
    dir="$(dirname "$target")"
    filename="$(basename "$target")"

    if [[ "$filename" =~ \.[a-zA-Z0-9]+\.(gz|bz2|xz|z|Z)$ ]]; then
        ext="${BASH_REMATCH[0]}"
        base="${filename%"$ext"}"
    elif [[ "$filename" == *.* ]]; then
        ext=".${filename##*.}"
        base="${filename%"$ext"}"
    else
        ext=""
        base="$filename"
    fi

    if [ "$dir" = "." ]; then
        echo "${base}_tmp${ext}"
    else
        echo "${dir}/${base}_tmp${ext}"
    fi
}
