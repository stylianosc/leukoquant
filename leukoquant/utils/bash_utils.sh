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
