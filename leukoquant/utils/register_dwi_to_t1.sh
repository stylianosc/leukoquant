#!/bin/bash
#
# register_dwi_to_t1.sh
#
# Reusable: register a single DWI (4D) image to a T1 image.
#
# Steps:
#   1. Extract b0 volumes via select_dwi_vols (or first volume if --bval not given)
#   2. Extract first b0 via fslroi
#   3. Skull-strip b0           (disable via --skip-skullstrip-dwi true)
#   4. Skull-strip T1           (disable via --skip-skullstrip-t1 true)
#   5. Affine register b0 -> T1 via reg_aladin
#   6. Compute inverse T1 -> DWI transform via reg_transform
#
# Outputs in --output-dir:
#   b0_4d.nii.gz         extracted b0 volumes
#   b0_first.nii.gz      first b0 volume
#   b0_brain.nii.gz      skull-stripped b0
#   t1_brain.nii.gz      skull-stripped T1
#   diff2t1_affine.txt   DWI -> T1 affine
#   t12diff_affine.txt   T1 -> DWI affine (inverse)
#
# Prints on stdout (always, for caller capture):
#   DIFF2T1: <path>
#   T12DIFF: <path>
#
# Exit code: 0 on success, 1 on error.
#

set -Eeuo pipefail

trap 'echo "ERROR: register_dwi_to_t1.sh failed at line ${LINENO}: ${BASH_COMMAND}" >&2' ERR

# ============================================================================
# Defaults
# ============================================================================

DWI=""
BVAL=""
T1=""
OUTPUT_DIR=""
THREADS=4
SKIP_SKULLSTRIP_T1=false
SKIP_SKULLSTRIP_DWI=false
VERBOSE=false
FORCE=false
DWI_WORK=""
BVAL_WORK=""

# ============================================================================
# Argument parsing
# ============================================================================

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dwi)                   DWI="$2";                  shift 2 ;;
        --bval)                  BVAL="$2";                 shift 2 ;;
        --t1)                    T1="$2";                   shift 2 ;;
        --output-dir)            OUTPUT_DIR="$2";           shift 2 ;;
        --threads)               THREADS="$2";              shift 2 ;;
        --skip-skullstrip-t1)    SKIP_SKULLSTRIP_T1="$2";  shift 2 ;;
        --skip-skullstrip-dwi)   SKIP_SKULLSTRIP_DWI="$2"; shift 2 ;;
        --verbose)               VERBOSE="$2";              shift 2 ;;
        --force)                 FORCE="$2";                shift 2 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

# ============================================================================
# Validation
# ============================================================================

if [[ -z "$DWI" || -z "$T1" || -z "$OUTPUT_DIR" ]]; then
    echo "ERROR: --dwi, --t1, and --output-dir are required." >&2
    echo "Usage: $0 --dwi DWI.nii.gz --t1 T1.nii.gz --output-dir DIR" >&2
    echo "         [--bval bvals] [--threads N]" >&2
    echo "         [--skip-skullstrip-t1 true|false] [--skip-skullstrip-dwi true|false]" >&2
    echo "         [--verbose true|false] [--force true|false]" >&2
    exit 1
fi

#echo "  [register] Inputs received"
#echo "  [register] DWI: $DWI"
#echo "  [register] T1 : $T1"
#echo "  [register] OUT: $OUTPUT_DIR"

if [[ ! -f "$DWI" ]]; then
    echo "ERROR: DWI file not found: $DWI" >&2; exit 1
fi
if [[ ! -f "$T1" ]]; then
    echo "ERROR: T1 file not found: $T1" >&2; exit 1
fi

# ============================================================================
# Helpers
# ============================================================================

log() {
    if [[ "$VERBOSE" == "true" || "$VERBOSE" == "True" || "$VERBOSE" == "1" ]]; then
        echo "$@"
    fi
}

is_true() {
    [[ "$1" == "true" || "$1" == "True" || "$1" == "1" ]]
}

print_image_info() {
    local label="$1"
    local file="$2"
    case "$file" in
        *.nii|*.nii.gz)
            local d1 d2 d3 d4 px py pz
            d1=$(fslval "$file" dim1    2>/dev/null || echo "?")
            d2=$(fslval "$file" dim2    2>/dev/null || echo "?")
            d3=$(fslval "$file" dim3    2>/dev/null || echo "?")
            d4=$(fslval "$file" dim4    2>/dev/null || echo "?")
            px=$(fslval "$file" pixdim1 2>/dev/null || echo "?")
            py=$(fslval "$file" pixdim2 2>/dev/null || echo "?")
            pz=$(fslval "$file" pixdim3 2>/dev/null || echo "?")
            echo "  ${label}: $(basename "$file") | dims=${d1}x${d2}x${d3} vols=${d4} | voxel=${px}x${py}x${pz}mm"
            ;;
        *)
            echo "  ${label}: $(basename "$file") (non-NIfTI, no space info)"
            ;;
    esac
}

prepare_dwi_input() {
    # Normalize DWI input using dwi_utils.py when needed (zip/dicom/other).
    local script_dir
    local dwi_utils
    local prep_dir
    local prep_cmd

    script_dir="$(cd "$(dirname "$0")" && pwd)"
    dwi_utils="$script_dir/dwi_utils.py"

    case "$DWI" in
        *.nii|*.nii.gz)
            # NIfTI input can be used directly; bvals optional.
            #echo "  [register] DWI input is NIfTI; using directly"
            DWI_WORK="$DWI"
            if [[ -n "$BVAL" && -f "$BVAL" ]]; then
                BVAL_WORK="$BVAL"
            else
                BVAL_WORK=""
            fi
            return 0
            ;;
        *)
            #echo "  [register] DWI input is non-NIfTI; preprocessing with dwi_utils.py"
            if [[ ! -f "$dwi_utils" ]]; then
                echo "ERROR: dwi_utils.py not found at $dwi_utils (required for non-NIfTI DWI input)" >&2
                exit 1
            fi

            prep_dir="$OUTPUT_DIR/dwi_prep"
            mkdir -p "$prep_dir"

            prep_cmd=(python3 "$dwi_utils" --dwi "$DWI" --outdir "$prep_dir")
            if [[ -n "$BVAL" && -f "$BVAL" ]]; then
                prep_cmd+=(--bvals "$BVAL")
            fi

            #echo "  [register] Running: ${prep_cmd[*]}"
            "${prep_cmd[@]}"

            DWI_WORK="$prep_dir/data.nii.gz"
            BVAL_WORK="$prep_dir/bvals"

            if [[ ! -f "$DWI_WORK" ]]; then
                echo "ERROR: dwi_utils.py did not produce expected DWI file: $DWI_WORK" >&2
                exit 1
            fi
            if [[ ! -f "$BVAL_WORK" ]]; then
                # Keep working even without bvals; fallback path below uses fslroi on first volume.
                BVAL_WORK=""
            fi
            return 0
            ;;
    esac
}

# ============================================================================
# Output paths
# ============================================================================

mkdir -p "$OUTPUT_DIR"

B0_4D="$OUTPUT_DIR/b0_4d.nii.gz"
B0_FIRST="$OUTPUT_DIR/b0_first.nii.gz"
B0_BRAIN="$OUTPUT_DIR/b0_brain.nii.gz"
T1_BRAIN="$OUTPUT_DIR/t1_brain.nii.gz"
DIFF2T1_AFF="$OUTPUT_DIR/diff2t1_affine.txt"
T12DIFF_AFF="$OUTPUT_DIR/t12diff_affine.txt"
B0_BRAIN_IN_T1="$OUTPUT_DIR/b0_brain_in_t1.nii.gz"

# Normalize DWI input to a NIfTI that downstream tools can consume.
prepare_dwi_input

# ============================================================================
# Skip if already done
# ============================================================================

if [[ -f "$DIFF2T1_AFF" && -f "$T12DIFF_AFF" ]] && ! is_true "$FORCE"; then
    echo "[SKIP] DWI->T1 registration already complete."
    #echo "DIFF2T1: $DIFF2T1_AFF"
    #echo "T12DIFF: $T12DIFF_AFF"
    exit 0
fi

# ============================================================================
# Print input image spaces
# ============================================================================

#log "  [DWI->T1] Starting registration"
#print_image_info "DWI input" "$DWI_WORK"
#print_image_info "T1  input" "$T1"

# ============================================================================
# Step 1: Extract b0 volumes
# ============================================================================

if [[ ! -f "$B0_4D" ]]; then
    if [[ -n "$BVAL_WORK" && -f "$BVAL_WORK" ]] && command -v select_dwi_vols >/dev/null 2>&1; then
        log "  [1/6] Extracting b0 volumes via select_dwi_vols ..."
        if ! select_dwi_vols "$DWI_WORK" "$BVAL_WORK" "$B0_4D" 0 > /dev/null 2>&1; then
            log "  [1/6] select_dwi_vols failed; falling back to fslroi (first volume) ..."
            rm -f "$B0_4D"
            fslroi "$DWI_WORK" "$B0_4D" 0 1
        fi
    else
        log "  [1/6] Using first volume as b0 (no bval or select_dwi_vols unavailable) ..."
        fslroi "$DWI_WORK" "$B0_4D" 0 1
    fi
    if [[ ! -f "$B0_4D" ]]; then
        echo "ERROR: b0 4D extraction failed: $B0_4D" >&2; exit 1
    fi
else
    log "  [1/6] [SKIP] b0 4D already exists."
fi

# ============================================================================
# Step 2: Extract first b0
# ============================================================================

if [[ ! -f "$B0_FIRST" ]]; then
    log "  [2/6] Extracting first b0 volume ..."
    fslroi "$B0_4D" "$B0_FIRST" 0 1
    if [[ ! -f "$B0_FIRST" ]]; then
        echo "ERROR: first b0 extraction failed: $B0_FIRST" >&2; exit 1
    fi
else
    log "  [2/6] [SKIP] First b0 already extracted."
fi

# ============================================================================
# Step 3: Skull-strip b0
# ============================================================================

if [[ ! -f "$B0_BRAIN" ]]; then
    if is_true "$SKIP_SKULLSTRIP_DWI"; then
        log "  [3/6] Skipping b0 skull stripping (--skip-skullstrip-dwi) ..."
        cp "$B0_FIRST" "$B0_BRAIN"
    else
        log "  [3/6] Skull stripping b0 ..."
        if ! _SYNTHSTRIP_OUT=$(mri_synthstrip -i "$B0_FIRST" -o "$B0_BRAIN" 2>&1); then
            echo "ERROR: mri_synthstrip failed on b0. Output:" >&2
            echo "$_SYNTHSTRIP_OUT" >&2
            exit 1
        fi
    fi
else
    log "  [3/6] [SKIP] b0 brain already exists."
fi

# ============================================================================
# Step 4: Skull-strip T1
# ============================================================================

if [[ ! -f "$T1_BRAIN" ]]; then
    if is_true "$SKIP_SKULLSTRIP_T1"; then
        log "  [4/6] Skipping T1 skull stripping (--skip-skullstrip-t1) ..."
        cp "$T1" "$T1_BRAIN"
    else
        log "  [4/6] Skull stripping T1 ..."
        if ! _SYNTHSTRIP_OUT=$(mri_synthstrip -i "$T1" -o "$T1_BRAIN" 2>&1); then
            echo "ERROR: mri_synthstrip failed on T1. Output:" >&2
            echo "$_SYNTHSTRIP_OUT" >&2
            exit 1
        fi
    fi
else
    log "  [4/6] [SKIP] T1 brain already exists."
fi

if [[ ! -f "$B0_BRAIN" ]]; then
    echo "ERROR: b0 brain not found after skull stripping: $B0_BRAIN" >&2; exit 1
fi
if [[ ! -f "$T1_BRAIN" ]]; then
    echo "ERROR: T1 brain not found after skull stripping: $T1_BRAIN" >&2; exit 1
fi

# ============================================================================
# Step 5: Affine register b0 -> T1
# ============================================================================
log "  [5/6] Registering b0 to T1 (affine, threads=$THREADS) ..."

reg_aladin \
    -ref "$T1_BRAIN"    \
    -flo "$B0_BRAIN"    \
    -aff "$DIFF2T1_AFF" \
    -res "$B0_BRAIN_IN_T1" \
    -omp "$THREADS"     \
    -voff > /dev/null 2>&1

if [[ ! -f "$DIFF2T1_AFF" ]]; then
    echo "ERROR: reg_aladin failed; affine not found: $DIFF2T1_AFF" >&2; exit 1
fi

# ============================================================================
# Step 6: Compute inverse transform (T1 -> DWI)
# ============================================================================

#echo "  [6/6] Computing inverse transform (T1 -> DWI) ..."
if ! reg_transform -invAff "$DIFF2T1_AFF" "$T12DIFF_AFF" > /dev/null 2>&1; then
    echo "ERROR: reg_transform -invAff failed" >&2; exit 1
fi

if [[ ! -f "$T12DIFF_AFF" ]]; then
    echo "ERROR: reg_transform failed; inverse affine not found: $T12DIFF_AFF" >&2; exit 1
fi

# ============================================================================
# Summary
# ============================================================================

log "  Registration complete."
#print_image_info "b0 brain" "$B0_BRAIN"
#print_image_info "T1 brain" "$T1_BRAIN"

#echo "DIFF2T1: $DIFF2T1_AFF"
#echo "T12DIFF: $T12DIFF_AFF"

# Delete not needed files
rm -f "$B0_4D" "$B0_FIRST"
rm -f "$T1_BRAIN"
# rm -f "$B0_BRAIN"
