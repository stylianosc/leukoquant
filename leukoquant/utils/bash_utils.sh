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

# Log scratch usage (GB) and remove the scratch directory.
# A 120-second timeout guards against frozen/stale network-attached scratch
# filesystems that would otherwise block the job indefinitely on rm -rf.
scratch_cleanup() {
    local scratch_path="${1:-}"

    if [ -z "$scratch_path" ]; then
        return 0
    fi

    echo "Cleaning up scratch directory: $scratch_path"
    if [ -d "$scratch_path" ]; then
        local scratch_size_gb
        scratch_size_gb=$(du -sb "$scratch_path" 2>/dev/null | awk '{printf "%.2f", $1/1024/1024/1024}')
        if [ -n "$scratch_size_gb" ]; then
            echo "Scratch usage before cleanup: ${scratch_size_gb} GB"
        fi
    fi

    if timeout 600 rm -rf "$scratch_path"; then
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

# Internal: polls `ps` every 5s for the combined RSS (KB) of every
# descendant of root_pid, keeping a running maximum in peak_file. Run in
# the background by start_peak_mem_monitor; not meant to be called
# directly.
_peak_mem_monitor_loop() {
    local root_pid="$1"
    local peak_file="$2"
    while true; do
        local total_kb
        total_kb=$(ps --no-headers -e -o pid,ppid,rss 2>/dev/null | awk -v root="$root_pid" '
            {
                ppid[$1] = $2
                rss[$1]  = $3
            }
            END {
                children[root] = 1
                changed = 1
                while (changed) {
                    changed = 0
                    for (p in ppid) {
                        if (!(p in children) && (ppid[p] in children)) {
                            children[p] = 1
                            changed = 1
                        }
                    }
                }
                sum = 0
                for (p in children) sum += rss[p]
                print sum
            }')
        local peak_kb
        peak_kb=$(cat "$peak_file" 2>/dev/null || echo 0)
        if [ -n "$total_kb" ] && [ "$total_kb" -gt "$peak_kb" ]; then
            echo "$total_kb" > "$peak_file"
        fi
        sleep 5
    done
}

# Starts a background loop tracking the peak combined RSS of the calling
# script's entire process tree (every descendant of $$), so heavy
# multi-tool pipelines can be sized from real usage instead of guesswork.
# Sets PEAK_MEM_PID/PEAK_MEM_FILE in the caller's shell.
#
# Call stop_peak_mem_monitor from the script's OWN existing EXIT/ERR/INT/
# TERM trap function (not a fresh `trap ... EXIT` here) -- bash traps
# don't chain, so registering a second one would silently replace or be
# replaced by whatever cleanup trap the script already has.
start_peak_mem_monitor() {
    PEAK_MEM_FILE=$(mktemp)
    echo 0 > "$PEAK_MEM_FILE"
    _peak_mem_monitor_loop "$$" "$PEAK_MEM_FILE" &
    PEAK_MEM_PID=$!
}

# Stops the loop started by start_peak_mem_monitor and prints the peak
# combined RSS observed, in MB. No-op if the monitor was never started.
stop_peak_mem_monitor() {
    [ -n "${PEAK_MEM_PID:-}" ] || return 0
    kill "$PEAK_MEM_PID" 2>/dev/null || true
    wait "$PEAK_MEM_PID" 2>/dev/null || true
    local peak_kb
    peak_kb=$(cat "$PEAK_MEM_FILE" 2>/dev/null || echo 0)
    echo "[mem] Peak whole-script RSS: $((peak_kb / 1024)) MB" >&2
    rm -f "$PEAK_MEM_FILE"
    unset PEAK_MEM_PID
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
