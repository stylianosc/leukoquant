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
#                                  hit an allocation failure, or was killed with
#                                  its peak near the address-space limit
#   "[killed] <stage> ..."         written by run_with_peak_mem for a stage that
#                                  was killed (SIGKILL or a "Killed" line) while its
#                                  peak was far below the limit: cause unknown
#   "[mem] Peak RSS (<stage>): ... (stage FAILED, exit N)"
#
# Why a kill alone is not called out of memory: a SIGKILL also comes from a time
# limit, the kernel OOM killer or a parent tool. bedpostx_mgh's monitor, for one,
# runs kill -9 on bedpostx as soon as any of its own error files is non-empty; in
# the 021-s-7055/7062/7092 TRACULA runs that happened at a 2.1 GB peak under a 24 GB
# limit, and the old "[oom]" label pointed the investigation at memory.
# ---------------------------------------------------------------------------

# Text that is direct evidence of running out of memory: glibc and C++ allocation
# failures, Python's MemoryError, the kernel OOM killer, our own "OUT OF MEMORY" line.
_OOM_PATTERN_STRICT="can't allocate memory|cannot allocate memory|std::bad_alloc|MemoryError|out of memory|oom-kill"
# Wider variant that also catches the shell's "Killed" line. A kill is evidence that
# the stage died, not of why; see _finish_stage_report for how it is classified.
_OOM_PATTERN="${_OOM_PATTERN_STRICT}|Killed"
# A killed stage counts as out of memory only if its peak reached this share of the
# address-space limit (ulimit -v, set by SGE from h_vmem). Below it the limit cannot
# have been what killed it.
_OOM_NEAR_CAP_PCT=90

# The address-space limit of this shell in KB ("unlimited" or "unknown" when there
# is none). A function, so tests can stand in a limit: `ulimit -v` can neither be
# set on macOS nor raised inside an SGE job.
_address_space_limit_kb() {
    ulimit -v 2>/dev/null || echo unknown
}

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
    # Only direct allocation-failure text counts; a bare "Killed" line says the step died, not why.
    # Status 137 (SIGKILL) stays out of memory here: these steps run under the SGE h_vmem limit, whose
    # kill is the usual cause, and no peak is measured here to tell otherwise.
    hit=$(printf '%s\n' "$text" | first_oom_line strict) || hit=""
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
    _OOM_REPORTED=0     # read by run_with_peak_mem for the usage record's "killed" column
    _KILL_REPORTED=0    # likewise, for a kill whose cause is unknown

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
    local alloc_hit="" kill_hit=""
    alloc_hit=$(printf '%s\n' "$text" | first_oom_line strict) || alloc_hit=""
    [ -z "$alloc_hit" ] && { kill_hit=$(printf '%s\n' "$text" | first_oom_line) || kill_hit=""; }
    # Not killed and no allocation failure: an ordinary failure, already marked on the [mem] line.
    [ -z "$alloc_hit" ] && [ -z "$kill_hit" ] && [ "$status" -ne 137 ] && return 0

    local limit_kb vm_peak_kb="${VM_PEAK_KB:-0}" vm_txt=""
    limit_kb=$(_address_space_limit_kb)
    [[ "$vm_peak_kb" =~ ^[0-9]+$ ]] || vm_peak_kb=0
    [ "$vm_peak_kb" -gt 0 ] && vm_txt="; peak address space $((vm_peak_kb / 1024)) MB (sampled)"
    local tail_txt="peak RSS ${peak_mb} MB${vm_txt}; address-space limit ${limit_kb} KB"

    local why="exit ${status}"
    [ "$status" -eq 137 ] && why="exit 137 (SIGKILL)"

    # Direct evidence (an allocation failure in the stage's own stderr): out of memory, whatever the peak.
    if [ -n "$alloc_hit" ]; then
        _OOM_REPORTED=1
        msg="[oom] ${stage}: out of memory (${why}; matched: ${alloc_hit}); ${tail_txt}"
        echo "$msg"
        echo "$msg" >&2
        return 0
    fi

    # Killed without an allocation message. Out of memory only if the larger of the exact RSS peak and the sampled
    # address-space peak reached _OOM_NEAR_CAP_PCT % of the limit; otherwise the limit cannot have killed it.
    [ -n "$kill_hit" ] && why="${why}; matched: ${kill_hit}"
    local peak_kb=0
    [[ "$peak_mb" =~ ^[0-9]+$ ]] && peak_kb=$((peak_mb * 1024))
    [ "$vm_peak_kb" -gt "$peak_kb" ] && peak_kb=$vm_peak_kb
    if [[ "$limit_kb" =~ ^[0-9]+$ ]] && [ "$limit_kb" -gt 0 ] \
        && [ $((peak_kb * 100)) -ge $((limit_kb * _OOM_NEAR_CAP_PCT)) ]; then
        _OOM_REPORTED=1
        msg="[oom] ${stage}: out of memory (${why}; peak within $((100 - _OOM_NEAR_CAP_PCT))% of the limit); ${tail_txt}"
    else
        _KILL_REPORTED=1
        msg="[killed] ${stage}: killed, cause unknown (${why}); ${tail_txt}"
        if [[ "$limit_kb" =~ ^[0-9]+$ ]]; then
            msg="${msg}; the peak stayed below ${_OOM_NEAR_CAP_PCT}% of the limit, so the limit did not kill it"
        else
            msg="${msg}; no numeric limit to compare the peak with"
        fi
        msg="${msg} (other causes: a time limit, the kernel OOM killer, a parent tool's kill -9)"
    fi
    echo "$msg"
    echo "$msg" >&2
}

# ---------------------------------------------------------------------------
# Peak address space of a stage.
#
# SGE limits the address space (h_vmem -> ulimit -v), but GNU time only measures resident memory (%M). To see how large
# the address space of a stage really gets, the processes of the stage are sampled from /proc every
# LEUKOQUANT_VM_SAMPLE_S seconds (default 30): the largest VmPeak (peak address space) and VmHWM (peak resident) of any
# one process, the unit the per-process limit acts on. Limits: a process that starts and ends between two samples is
# missed, and its peak is lost with it. /proc is read from LEUKOQUANT_PROC_DIR (default /proc; tests use a fake one).
# ---------------------------------------------------------------------------

# _vm_scan <root_pid> [exclude_pid]
# Prints "<VmPeak kB> <name> <VmHWM kB> <name>": the largest values among root_pid and its descendants. Prints nothing
# when /proc is not there. Always succeeds.
_vm_scan() {
    local root="$1" exclude="${2:-0}" proc="${LEUKOQUANT_PROC_DIR:-/proc}"
    [ -d "$proc" ] || return 0
    awk -v root="$root" -v excl="$exclude" '
        FNR == 1 { pid = ""; }
        # The whole name, with blanks turned into underscores and "?" for an empty one, so that it is always exactly
        # one word: the sampler stores the values as whitespace-separated fields and a name with spaces (or none)
        # would shift every field after it.
        /^Name:/   { n = $0; sub(/^Name:[ \t]*/, "", n); gsub(/[ \t]+/, "_", n); if (n == "") n = "?"; name[FILENAME] = n }
        /^Pid:/    { pid = $2; file[pid] = FILENAME }
        /^PPid:/   { ppid[pid] = $2 }
        /^VmPeak:/ { vmpeak[pid] = $2 }
        /^VmHWM:/  { vmhwm[pid] = $2 }
        /^VmSize:/ { vmsize[pid] = $2 }
        END {
            bp = 0; bh = 0; vs = 0; np = "-"; nh = "-"
            for (p in file) {
                x = p; found = 0
                for (depth = 0; depth < 64 && x != "" && x != 0; depth++) {
                    if (x == root) { found = 1; break }
                    x = ppid[x]
                }
                if (!found || p == excl) continue
                n = name[file[p]]
                vs += vmsize[p] + 0
                if (vmpeak[p] + 0 > bp) { bp = vmpeak[p] + 0; np = n }
                if (vmhwm[p] + 0 > bh) { bh = vmhwm[p] + 0; nh = n }
            }
            if (bp > 0 || bh > 0) print bp, np, bh, nh, vs
        }' "$proc"/[0-9]*/status 2>/dev/null
    return 0
}

# start_vm_sampler <label> [root_pid]   (root_pid defaults to this shell)
start_vm_sampler() {
    VM_SAMPLER_LABEL="${1:-stage}"
    local root="${2:-$$}" proc="${LEUKOQUANT_PROC_DIR:-/proc}" interval="${LEUKOQUANT_VM_SAMPLE_S:-30}"
    VM_SAMPLER_PID=""
    VM_SAMPLER_FILE=""
    VM_SAMPLER_INTERVAL="$interval"
    # Forget the previous stage's peaks, so a stage without samples is not judged by them.
    VM_PEAK_KB=""; VM_HWM_KB=""; VM_SUM_KB=""; VM_SCRATCH_KB=""; VM_MINFREE_KB=""; VM_NSAMPLES=""
    [ -d "$proc" ] || return 0
    VM_SAMPLER_FILE=$(mktemp 2>/dev/null) || { VM_SAMPLER_FILE=""; return 0; }
    (
        while :; do
            sample=$(_vm_scan "$root" "$BASHPID")
            if [ -n "$sample" ]; then
                # The trailing "_" takes the job-wide fields that follow; without it the last named variable would
                # swallow them and the next write would carry extra fields (a syntax error one sample later).
                read -r cur_p cur_pn cur_h cur_hn _ < "$VM_SAMPLER_FILE" 2>/dev/null || true
                read -r new_p new_pn new_h new_hn new_vs <<< "$sample"
                [ "${new_p:-0}" -gt "${cur_p:-0}" ] 2>/dev/null || { new_p="${cur_p:-0}"; new_pn="${cur_pn:--}"; }
                [ "${new_h:-0}" -gt "${cur_h:-0}" ] 2>/dev/null || { new_h="${cur_h:-0}"; new_hn="${cur_hn:--}"; }
                # Job-wide figures for the usage record: largest sum of address space, largest scratch, smallest free space.
                read -r _ _ _ _ cur_vs cur_sc cur_fr cur_n < "$VM_SAMPLER_FILE" 2>/dev/null || true
                [ "${new_vs:-0}" -gt "${cur_vs:-0}" ] 2>/dev/null || new_vs="${cur_vs:-0}"
                new_sc="${cur_sc:-0}"; new_fr="${cur_fr:-0}"
                if [ -n "${LEUKOQUANT_USAGE_SCRATCH:-}" ] && [ -d "$LEUKOQUANT_USAGE_SCRATCH" ]; then
                    sc=$(du -sk "$LEUKOQUANT_USAGE_SCRATCH" 2>/dev/null | awk '{print $1}')
                    fr=$(df -Pk "$LEUKOQUANT_USAGE_SCRATCH" 2>/dev/null | awk 'NR==2 {print $4}')
                    [ "${sc:-0}" -gt "${new_sc:-0}" ] 2>/dev/null && new_sc="$sc"
                    if [ -n "$fr" ] && { [ "${new_fr:-0}" -eq 0 ] || [ "$fr" -lt "$new_fr" ]; } 2>/dev/null; then new_fr="$fr"; fi
                fi
                echo "$new_p $new_pn $new_h $new_hn $new_vs $new_sc $new_fr $(( ${cur_n:-0} + 1 ))" > "$VM_SAMPLER_FILE"
            fi
            sleep "$interval"
        done
    ) &
    VM_SAMPLER_PID=$!
    return 0
}

# stop_vm_sampler: stops the sampler and prints one "[vm] <label>: ..." line. Always succeeds.
stop_vm_sampler() {
    local label="${VM_SAMPLER_LABEL:-stage}" msg peak_kb peak_name hwm_kb hwm_name
    if [ -n "${VM_SAMPLER_PID:-}" ]; then
        kill "$VM_SAMPLER_PID" 2>/dev/null || true
        wait "$VM_SAMPLER_PID" 2>/dev/null || true
    fi
    if [ -z "${VM_SAMPLER_FILE:-}" ] || [ ! -s "${VM_SAMPLER_FILE:-/nonexistent}" ]; then
        msg="[vm] ${label}: unavailable (no /proc, or no sample taken)"
    else
        read -r peak_kb peak_name hwm_kb hwm_name VM_SUM_KB VM_SCRATCH_KB VM_MINFREE_KB VM_NSAMPLES < "$VM_SAMPLER_FILE"
        VM_PEAK_KB="$peak_kb"; VM_HWM_KB="$hwm_kb"
        msg="[vm] ${label}: peak VmPeak $((peak_kb / 1024)) MB (${peak_name}), VmHWM $((hwm_kb / 1024)) MB (${hwm_name});"
        msg="${msg} address-space limit $(ulimit -v 2>/dev/null || echo unknown) KB; sampled every ${VM_SAMPLER_INTERVAL:-30} s,"
        msg="${msg} a process that ends between two samples is missed"
    fi
    echo "$msg"
    echo "$msg" >&2
    [ -n "${VM_SAMPLER_FILE:-}" ] && rm -f "$VM_SAMPLER_FILE"
    VM_SAMPLER_PID=""
    VM_SAMPLER_FILE=""
    return 0
}

# _write_usage_record <stage> <status> <peak_rss_kb> <elapsed_s> <killed 0|1>
# One line per stage into the usage tables named in LEUKOQUANT_USAGE_FILES (space-separated), with the context
# (dataset, subject, requests, ...) from LEUKOQUANT_USAGE_CONTEXT ("key=value ..."). See utils/usage_record.py for the
# schema and the limits of the sampled columns. Does nothing when no table is configured; never fails the stage.
_write_usage_record() {
    [ -n "${LEUKOQUANT_USAGE_FILES:-}" ] || return 0
    local script="${LEUKOQUANT_USAGE_SCRIPT:-$(dirname "${BASH_SOURCE[0]}")/usage_record.py}" f files=()
    for f in $LEUKOQUANT_USAGE_FILES; do files+=(--file "$f"); done
    # shellcheck disable=SC2086
    python3 "$script" append "${files[@]}" $LEUKOQUANT_USAGE_CONTEXT stage="$1" exit_status="$2" \
        peak_rss_mb="$(( ${3:-0} / 1024 ))" elapsed_s="$4" killed="$5" \
        peak_vmpeak_mb="$(( ${VM_PEAK_KB:-0} / 1024 ))" peak_vmhwm_mb="$(( ${VM_HWM_KB:-0} / 1024 ))" \
        peak_vm_sum_mb="$(( ${VM_SUM_KB:-0} / 1024 ))" peak_scratch_mb="$(( ${VM_SCRATCH_KB:-0} / 1024 ))" \
        min_free_scratch_mb="$(( ${VM_MINFREE_KB:-0} / 1024 ))" sample_interval_s="${VM_SAMPLER_INTERVAL:-}" \
        n_samples="${VM_NSAMPLES:-0}" host="$(hostname 2>/dev/null)" job_id="${JOB_ID:-}" \
        task_id="${SGE_TASK_ID:-}" >/dev/null 2>&1 || true
    return 0
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
            local t_start=$SECONDS
            start_vm_sampler "$label"
            /usr/bin/time -o "$tmp_time" -f "%M %E %P %x" "$@" || status=$?
            stop_vm_sampler
            local peak_kb elapsed cpu exit_code
            # On a signal newer GNU time first writes "Command terminated by signal N": the figures are on the last line.
            read -r peak_kb elapsed cpu exit_code <<< "$(tail -n 1 "$tmp_time" 2>/dev/null)"
            # A command killed by a signal must be reported as 128 + signal (137 for SIGKILL), as the shell does.
            # GNU time 1.9 (the containers) already returns that and names the signal in its output. GNU time 1.7
            # (the cluster hosts) returns the bare signal number (9) and writes nothing about it, which looks like a
            # plain "exit 9"; its %x (the command's own exit status) tells them apart: 0 after a signal, 9 after exit 9.
            if grep -q "terminated by signal 9" "$tmp_time" 2>/dev/null; then
                status=137
            elif [ "$status" -gt 0 ] && [ "$status" -lt 128 ] && [ "${exit_code:-}" = "0" ]; then
                status=$((128 + status))
            fi
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
            local killed=0
            { [ "$status" -eq 137 ] || [ "${_OOM_REPORTED:-0}" -eq 1 ] || [ "${_KILL_REPORTED:-0}" -eq 1 ]; } && killed=1
            _write_usage_record "${label:-stage}" "$status" "${peak_kb:-0}" "$((SECONDS - t_start))" "$killed"
            return "$status"
        elif [ "$(uname -s)" = "Darwin" ]; then
            # BSD time (macOS local workstations)
            local err_file err_start=0
            err_file=$(_stderr_file)
            [ -n "$err_file" ] && err_start=$(_file_size "$err_file")
            local t_start=$SECONDS
            start_vm_sampler "$label"
            /usr/bin/time -l -o "$tmp_time" "$@" || status=$?
            stop_vm_sampler
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
            local killed=0
            { [ "$status" -eq 137 ] || [ "${_OOM_REPORTED:-0}" -eq 1 ] || [ "${_KILL_REPORTED:-0}" -eq 1 ]; } && killed=1
            _write_usage_record "${label:-stage}" "$status" "${peak_kb:-0}" "$((SECONDS - t_start))" "$killed"
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


# Copy the logs of a failed bedpostx run (dmri.bedpostX/logs, esp. bedpostx_post.e*) out of node scratch.
#   $1  the dmri.bedpostX directory in scratch
#   $2  destination folder; the copy goes to $2/bedpostx_logs (a previous copy is replaced)
# Never fails: a missing directory or a copy error only prints a note, so it cannot hide the original failure.
rescue_bedpostx_logs() {
    local bedpostx_dir="${1:-}" dest_dir="${2:-}"
    if [ -z "$bedpostx_dir" ] || [ -z "$dest_dir" ] || [ ! -d "$bedpostx_dir/logs" ]; then
        echo "[bedpostx] no logs to rescue from '$bedpostx_dir'"
        return 0
    fi
    rm -rf "$dest_dir/bedpostx_logs"
    if cp -r "$bedpostx_dir/logs" "$dest_dir/bedpostx_logs"; then
        echo "[bedpostx] logs of the failed bedpostx run copied to $dest_dir/bedpostx_logs"
    else
        echo "[bedpostx] could not copy the bedpostx logs from $bedpostx_dir/logs to $dest_dir"
    fi
    return 0
}


# Log the node and its free disk space, to tell a kill for lack of node scratch from a memory kill.
#   $1  label of the moment (e.g. "before bedpostx")
#   $2  the job's scratch folder (optional)
# Prints: host, df of the scratch and tmp locations, the size of the job's scratch folder. Never fails.
log_node_disk_state() {
    local label="${1:-}" scratch_dir="${2:-}"
    echo "[node] ${label}: host: $(hostname 2>/dev/null || echo unknown), $(date '+%F %T')"
    df -h "${scratch_dir:-/scratch0}" "${TMPDIR:-/tmp}" /scratch0 2>&1 | sed 's/^/[node]   df /'
    if [ -n "$scratch_dir" ] && [ -d "$scratch_dir" ]; then
        echo "[node]   scratch folder size: $(du -sh "$scratch_dir" 2>/dev/null | cut -f1) ($scratch_dir)"
    else
        echo "[node]   scratch folder size: n/a ($scratch_dir)"
    fi
    return 0
}

# ---------------------------------------------------------------------------
# Provenance lines: which node and CPU ran a job, and a checksum of every intermediate.
#
# Two runs of the same job with the same code and inputs gave slightly different tract
# volumes on some nodes (EPAD, 2026-11-18) while their logs were identical. These lines
# make two logs diffable at the point where the runs part ways:
#   [host] <hostname> cpu=<model> simd=<widest vector extension>
#   [sha] <label> <relative path> <md5 | missing>
# Text only: neither function changes the exit status of the caller (safe under set -e).
# ---------------------------------------------------------------------------

# md5 of one file, portable between Linux (md5sum) and macOS (md5).
_md5_of() {
    if command -v md5sum >/dev/null 2>&1; then
        md5sum < "$1" | cut -c1-32
    else
        md5 -q "$1"
    fi
}

# Print the host name, the CPU model and its widest SIMD extension. The vector width
# matters because numerical libraries pick code paths by it, which can change the last
# bits of a floating-point result.
log_host_cpu() {
    local model="unknown" simd="unknown"
    if [ -r /proc/cpuinfo ]; then
        model=$(grep -m1 "model name" /proc/cpuinfo | cut -d: -f2- | sed 's/^ *//; s/  */ /g')
        local flags
        flags=$(grep -m1 "^flags" /proc/cpuinfo)
        for ext in avx512f avx2 avx sse4_2; do
            if printf '%s\n' "$flags" | grep -qw "$ext"; then simd="$ext"; break; fi
        done
    elif command -v sysctl >/dev/null 2>&1; then
        model=$(sysctl -n machdep.cpu.brand_string 2>/dev/null || echo unknown)
    fi
    echo "[host] $(hostname) cpu=${model:-unknown} simd=${simd}"
    return 0
}

# Make the numerical libraries behave the same on every node.
#
# numpy and scipy ship an OpenBLAS that picks its matrix kernel at run time from the CPU. A kernel
# with fused multiply-add rounds differently from one without, so the last bits of every matrix
# product depend on the node. In the metrics step that decides which voxel a transformed streamline
# point falls in when it lies close to a voxel boundary: the same job gave different tract masks on
# the abner nodes (AVX only) than on the AVX2 and AVX-512 nodes (EPAD, 2026-11-18; see
# tests/test_blas_kernel_pin.py). Pinning one kernel that every node can run makes the result
# independent of the node. Core2 (SSSE3, no AVX, no FMA) is a kernel every node we know can run: SandyBridge, the first
# choice, needs AVX and killed the first python3 call with SIGILL on arbuckle (Xeon E7-8837, SSE4.2 only; ADNI-3 metrics,
# 2026-11-20). Core2 and Nehalem gave byte-identical density maps on arbuckle, abner (AVX), burns and saunders (AVX2/FMA),
# and identical to the SandyBridge ones, so results made under the SandyBridge pin stay valid.
# The variable is exported, so Python started by this job and by its helper scripts inherits it.
# A value set by the caller wins, so another kernel can be tried; the log line says which applied.
# Never fails the job.
pin_blas_kernel() {
    if [ -n "${OPENBLAS_CORETYPE:-}" ]; then
        echo "[blas] OPENBLAS_CORETYPE=${OPENBLAS_CORETYPE} (already set by the caller, not changed)"
        return 0
    fi
    export OPENBLAS_CORETYPE="Core2"
    echo "[blas] OPENBLAS_CORETYPE=${OPENBLAS_CORETYPE} (pinned: one kernel on every node, so results do not depend on the CPU)"
    return 0
}

# Make the skull-strip network give the same mask on as many nodes as possible.
#
# mri_synthstrip is a PyTorch network that runs on the CPU, and its MKL picks its kernels from the CPU at run
# time. On an AVX2 node (Haswell) the T1 mask differed from the one on AVX-512 (Xeon Gold) and AMD EPYC nodes for the
# same input, the DWI-to-T1 affine moved with it, and so did every DTI and NODDI value sampled in the lesions
# (ADNI-3, 2026-11-21; see tests/test_mkl_kernel_pin.py). MKL_CBWR=COMPATIBLE restricts MKL to kernels that give the same
# bits on every x86 CPU. With it, three nodes of three CPU families gave identical b0 and T1 masks and identical
# reg_aladin affines (3 sessions, 2 repeats). Known exception: the AVX-only abner node still gave a different T1 mask
# for one session (subj-130-s-6611, 2026-11-21), so this is not a full node-independence guarantee. The OpenBLAS pin above does not reach MKL, and ONEDNN_MAX_CPU_ISA and a
# single thread did not change the mask. The run time is unchanged. MKL reads the variable when it is loaded, so it is
# exported before any mri_synthstrip process starts. A value set by the caller wins. Never fails the job.
pin_mkl_kernel() {
    if [ -n "${MKL_CBWR:-}" ]; then
        echo "[mkl] MKL_CBWR=${MKL_CBWR} (already set by the caller, not changed)"
        return 0
    fi
    export MKL_CBWR="COMPATIBLE"
    echo "[mkl] MKL_CBWR=${MKL_CBWR} (pinned: restricts MKL to CPU-independent kernels; identical masks on AVX2, AVX-512 and AMD nodes, not guaranteed on AVX-only nodes)"
    return 0
}

# Usage: log_checksums <label> <file-or-directory>...
# One [sha] line per file; a directory is walked recursively and listed sorted, with
# paths relative to it, so the lines of two runs line up. A missing path is reported
# as "missing" instead of failing.
log_checksums() {
    local label="$1"; shift
    local target rel f
    for target in "$@"; do
        if [ -z "$target" ]; then
            echo "[sha] ${label} (unset)"
        elif [ -f "$target" ]; then
            echo "[sha] ${label} $(basename "$target") $(_md5_of "$target")"
        elif [ -d "$target" ]; then
            while IFS= read -r f; do
                rel="${f#"$target"/}"
                echo "[sha] ${label} ${rel} $(_md5_of "$f")"
            done < <(find "$target" -type f | LC_ALL=C sort)
        else
            echo "[sha] ${label} $(basename "$target") missing"
        fi
    done
    return 0
}
