#!/bin/bash
# Modified by Stylianos Charalampous on 11/11/2025

# Exit on error
set -euo pipefail

# Default values
filename=""
mask=""
output_folder=""
verbose=0
threads=1
db_path="/GIF/db_mideface/db.xml"
gif_bin_path="/leukoquant/leukoquant/external/gif/bin"

# Parse command line arguments
while [ $# -gt 0 ]; do
    case $1 in
        -f | --filename) 
            filename="$2"
            shift 2
            ;;
        -m | --mask)
            mask="$2"
            shift 2
            ;;
        -o | --output-folder)
            output_folder="$2"
            shift 2
            ;;
        -d | --db)
            db_path="$2"
            shift 2
            ;;
        -v | --verbose)
            verbose=1
            shift
            ;;
        -t | --threads)
            threads="$2"
            shift 2
            ;;
        -b | --gif-bin)
            gif_bin_path="$2"
            shift 2
            ;;
        *)
            echo "Unknown parameter: $1" >&2
            exit 1
            ;;
    esac
done

# Validate required parameters
if [ -z "$filename" ] || [ -z "$output_folder" ]; then
    echo "Usage: $0 -f|--filename <input_nifti> -o|--output-folder <output_dir> [options]" >&2
    echo "Options:" >&2
    echo "  -m, --mask <mask_file>      Brain mask file" >&2
    echo "  -d, --db <db_file>          GIF database file (default: /GIF/GIF/db_mideface/db.xml)" >&2
    echo "  -v, --verbose              Enable verbose output" >&2
    echo "  -t, --threads <num_threads> Number of threads" >&2
    exit 1
fi

log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1"
}

# Check if input file exists
if [ ! -f "$filename" ]; then
    log "Error: Input file not found: $filename" >&2
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$output_folder"

# Prepare mask parameter
mask_param=""
if [ -n "$mask" ]; then
    if [ ! -f "$mask" ]; then
        log "Warning: Mask file not found: $mask" >&2
    else
        mask_param="-mask $mask"
    fi
fi

log "Starting GIF processing for $filename"
log "Output directory: $output_folder"
[ -n "$mask" ] && log "Using mask: $mask"

log "========== DIAGNOSTICS BEFORE RUN =========="
log "System Memory:"
free -h || echo "(free command unavailable)"
log "System Limits:"
ulimit -a || true
log "OMP/Threading Environment:"
env | grep -i "omp\|thread" || true
log "Input NIfTI Details:"
ls -lh "$filename" || true
log "Checking GIF DB & Binaries:"
ls -lh "$db_path" || true
ls -lh "$gif_bin_path/seg_GIF" || true
# Removing file command as it's not available in the container
log "============================================"

# Ensure all GIF binaries are executable (git checkouts may strip the bit)
chmod +x "$gif_bin_path"/* 2>/dev/null || true

# Run GIF segmentation.
# GNU time's -o flag opens the target path as a fresh file descriptor, which
# on NFS causes it to write the timing block at a different byte offset than
# the main stdout stream, leaving null-byte gaps in the log file.  Writing
# the timing summary to a local temp file and cat-ing it afterwards guarantees
# all output flows through the same fd 1 sequentially and eliminates the nulls.
_time_out=$(mktemp /tmp/gif_time_XXXXXX)
set +e
/usr/bin/time -v -o "$_time_out" \
"$gif_bin_path/seg_GIF" \
    -in "$filename" \
    -db "$db_path" \
    -v 1 \
    -regNMI \
    -segPT 0.1 \
    -out "$output_folder/" \
    -temper 0.05 \
    $mask_param \
    -omp "$threads" \
    -lncc_ker -4 \
    -regBE 0.001 \
    -regJL 0.00005

exit_code=$?
set -e
cat "$_time_out"
rm -f "$_time_out"

if [ $exit_code -ne 0 ]; then
    log "Error: GIF processing failed with exit code $exit_code" >&2
    if [ $exit_code -eq 139 ]; then
        log "CRITICAL ERROR: Exit Code 139 indicates a SIGSEGV (Signal 11 / Segmentation Fault)." >&2
        log "This means seg_GIF accessed invalid memory. Possible causes:" >&2
        log "  1. Thread stack limit exceeded (try increasing OMP_STACKSIZE further or setting ulimit -s unlimited)." >&2
        log "  2. Corrupt, Float64, or incompatible input NIfTI file ($filename)." >&2
        log "  3. Incomplete, corrupted, or unreadable GIF database at $db_path." >&2
        log "  4. Singularity/OS library incompatibilities (glibc mismatch)." >&2
        log "Attempting to retrieve kernel termination logs:" >&2
        dmesg -T | grep -i "seg_gif\|segfault\|error 4" | tail -n 10 || echo "(dmesg logs unavailable; typically requires root)" >&2
    fi
    exit $exit_code
fi

log "GIF processing completed successfully"
log "Results saved to: $output_folder"

exit 0

