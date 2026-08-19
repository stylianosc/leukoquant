#!/bin/bash
set -euo pipefail

# ------------------------------------------------------------------------------
# Script: z_score_calc.sh
# Description: Full Z-score pipeline for one target subject across ALL metrics:
#   Pre-Step: (if metric-space=dwi) Prepare DWI inputs via dwi_utils.py
#   1. Generate GLM design matrix (once -- same for all metrics).
#   2. Register each healthy T1 to target T1 (once -- CPP reused for all metrics).
#   2.5. (if dwi) Extract b0, register to healthy T1 (diff2t1 affine).
#   3. For each metric:
#      a. Resample healthy metrics (compose transforms if dwi space).
#      b. Merge into a 4D image.
#      c. Run FSL GLM; compute residual std dev and predicted map.
#      d. Compute Z-score: (target - predicted) / std_dev.
#
# Named arguments:
#   --target-t1         brain-extracted T1 of the target subject
#   --target-id         subject ID string of the target
#   --demographics      CSV with subject IDs and covariate columns
#   --healthy-ids       plain-text file with one healthy subject ID per line
#   --output-dir        per-subject output directory
#   --shared-cache-dir  dataset-level cache for healthy-cohort DWI/T1 prep,
#                       shared and reused across every target subject's run
#                       (optional; falls back to the old per-target
#                       location under --output-dir/tmp if omitted)
#   --metrics           comma-separated metric names, e.g. "FA,MD,ICVF"
#   --covariates        comma-separated covariate names (optional)
#   --poly-terms        comma-separated covariates for squared terms (optional)
#   --threads           CPU threads for NiftyReg tools (default: 4)
#   --healthy-t1s       space-separated healthy T1 paths (one per healthy subject)
#   --target-metrics    space-separated target metric paths (one per metric)
#   --healthy-metrics   space-separated healthy metric paths
#                       (flat, grouped by metric then subject:
#                        m1_h1, m1_h2, ..., m2_h1, m2_h2, ...)
#   --metric-space      "t1" or "dwi" (default: "t1")
#   --output-space      "t1" or "dwi" (default: "t1")
#   --dwi-paths         space-separated DWI paths (target + N_HEALTHY, required if metric-space=dwi)
#                       Supports: .nii/.nii.gz, DICOM directory, or .zip of DICOMs
#   --bval-paths        (optional) space-separated bval file paths
#   --verbose           enable progress/debug logging
# ------------------------------------------------------------------------------

TARGET_T1=""
TARGET_ID=""
DEMO_CSV=""
HEALTHY_IDS_FILE=""
OUTPUT_DIR=""
SHARED_CACHE_DIR=""
METRIC_NAMES_STR=""
COVARIATES=""
POLY_TERMS=""
THREADS="4"
METRIC_SPACE="t1"
OUTPUT_SPACE="t1"
DWI_PATHS=()
BVAL_PATHS=()
HEALTHY_T1S=()
TARGET_METRICS=()
HEALTHY_METRICS_ALL=()
QC=0
SKIP_SKULLSTRIP_T1=0
SKIP_SKULLSTRIP_DWI=0
VERBOSE=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --target-t1)       TARGET_T1="$2";        shift 2 ;;
        --target-id)       TARGET_ID="$2";        shift 2 ;;
        --demographics)    DEMO_CSV="$2";         shift 2 ;;
        --healthy-ids)     HEALTHY_IDS_FILE="$2"; shift 2 ;;
        --output-dir)      OUTPUT_DIR="$2";       shift 2 ;;
        --shared-cache-dir) SHARED_CACHE_DIR="$2"; shift 2 ;;
        --metrics)         METRIC_NAMES_STR="$2"; shift 2 ;;
        --covariates)      COVARIATES="$2";       shift 2 ;;
        --poly-terms)      POLY_TERMS="$2";       shift 2 ;;
        --threads)         THREADS="$2";          shift 2 ;;
        --metric-space)    METRIC_SPACE="$2";      shift 2 ;;
        --output-space)    OUTPUT_SPACE="$2";      shift 2 ;;
        --skip-skullstrip-t1)  SKIP_SKULLSTRIP_T1=1; shift ;;
        --skip-skullstrip-dwi) SKIP_SKULLSTRIP_DWI=1; shift ;;
        --verbose)         VERBOSE=1;            shift ;;
        --dwi-paths)       shift; while [[ $# -gt 0 && "$1" != --* ]]; do DWI_PATHS+=("$1"); shift; done ;;
        --bval-paths)      shift; while [[ $# -gt 0 && "$1" != --* ]]; do BVAL_PATHS+=("$1"); shift; done ;;
        --healthy-t1s)     shift; while [[ $# -gt 0 && "$1" != --* ]]; do HEALTHY_T1S+=("$1"); shift; done ;;
        --target-metrics)  shift; while [[ $# -gt 0 && "$1" != --* ]]; do TARGET_METRICS+=("$1"); shift; done ;;
        --healthy-metrics) shift; while [[ $# -gt 0 && "$1" != --* ]]; do HEALTHY_METRICS_ALL+=("$1"); shift; done ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# Validate required arguments
# ---------------------------------------------------------------------------

for _req in TARGET_T1 TARGET_ID HEALTHY_IDS_FILE OUTPUT_DIR METRIC_NAMES_STR; do
    if [[ -z "${!_req}" ]]; then
        echo "Error: --${_req//_/-} is required." >&2; exit 1
    fi
done

mkdir -p "$OUTPUT_DIR"

IFS=',' read -r -a METRIC_NAMES <<< "$METRIC_NAMES_STR"
N_METRICS="${#METRIC_NAMES[@]}"
N_HEALTHY="${#HEALTHY_T1S[@]}"

# ---------------------------------------------------------------------------
# Validate list arguments
# ---------------------------------------------------------------------------

if [[ "$N_HEALTHY" -eq 0 ]]; then
    echo "Error: --healthy-t1s is required and must be non-empty." >&2; exit 1
fi
if [[ "${#TARGET_METRICS[@]}" -ne "$N_METRICS" ]]; then
    echo "Error: ${#TARGET_METRICS[@]} --target-metrics paths but $N_METRICS metric names." >&2; exit 1
fi
EXPECTED=$(( N_METRICS * N_HEALTHY ))
if [[ "${#HEALTHY_METRICS_ALL[@]}" -ne "$EXPECTED" ]]; then
    echo "Error: expected $EXPECTED --healthy-metrics paths, got ${#HEALTHY_METRICS_ALL[@]}." >&2; exit 1
fi
if [[ "$METRIC_SPACE" != "t1" && "$METRIC_SPACE" != "dwi" ]]; then
    echo "Error: --metric-space must be 't1' or 'dwi'." >&2; exit 1
fi
if [[ "$OUTPUT_SPACE" != "t1" && "$OUTPUT_SPACE" != "dwi" ]]; then
    echo "Error: --output-space must be 't1' or 'dwi'." >&2; exit 1
fi

# ---------------------------------------------------------------------------
# Locate Python utility (same directory as this script)
# ---------------------------------------------------------------------------

log() {
    if [[ "$VERBOSE" -eq 1 ]]; then
        echo "$@"
    fi
}

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
UTILS_SCRIPT="$SCRIPT_DIR/z_score_utils.py"
if [[ ! -f "$UTILS_SCRIPT" ]]; then
    echo "Error: z_score_utils.py not found at $UTILS_SCRIPT" >&2; exit 1
fi

DWI_UTILS_SCRIPT="$SCRIPT_DIR/dwi_utils.py"
if [[ ! -f "$DWI_UTILS_SCRIPT" ]]; then
    echo "Error: dwi_utils.py not found at $DWI_UTILS_SCRIPT" >&2; exit 1
fi

echo ""
echo "FSLDIR: $FSLDIR"
source "${FSLDIR}/etc/fslconf/fsl.sh"


# Per-target scratch: target's own T1/DWI prep, all healthy->target
# registration output, design matrix, xfm files -- none of it is reusable
# by any other target's run, so it belongs on ephemeral per-job scratch,
# not the SAN. TMPDIR is set by the caller (z_score_workflow.smk) to a
# per-job /scratch0 path; fall back to the old SAN-based location only for
# standalone CLI use without that wrapper.
TMP_DIR="${TMPDIR:-$OUTPUT_DIR/tmp}"
mkdir -p "$TMP_DIR"

# Shared, dataset-level cache for healthy-cohort DWI/T1 prep (the part that
# genuinely is identical across every target's run -- see prepare_dwi_inputs
# and prepare_t1_inputs below). Falls back to the old per-target location
# if the caller doesn't pass --shared-cache-dir, so standalone CLI use still
# works, just without the cross-run sharing.
SHARED_CACHE_DIR="${SHARED_CACHE_DIR:-$OUTPUT_DIR/tmp}"
mkdir -p "$SHARED_CACHE_DIR"

# Registered here, right after TMP_DIR exists, not at the bottom of the
# script -- a trap only fires for signals received after it's registered,
# so defining this at the end (as it was before, disabled) meant a crash
# anywhere in Steps 1-6 would exit without ever cleaning up. Only $TMP_DIR
# is removed: the shared cache above must never be touched here, since its
# entire purpose is to persist and be reused by every other target's run.
function finish {
    echo "Cleaning up scratch directory: $TMP_DIR"
    rm -rf "$TMP_DIR"
}
trap finish EXIT ERR INT TERM

echo "=== Z-Score Pipeline: $TARGET_ID ==="
echo "  Metrics          : $METRIC_NAMES_STR"
echo "  Metric space     : $METRIC_SPACE"
echo "  Output space     : $OUTPUT_SPACE"
echo "  Healthy count    : $N_HEALTHY"
echo "  Covariates       : ${COVARIATES:-<none>}"
echo "  Poly terms       : ${POLY_TERMS:-<none>}"
echo "  Threads          : $THREADS"

if [[ "$SKIP_SKULLSTRIP_T1" -eq 1 ]]; then
    echo "Skipping T1 Skull Stripping: Enabled"
fi
if [[ "$SKIP_SKULLSTRIP_DWI" -eq 1 ]]; then
    echo "Skipping DWI Skull Stripping: Enabled"
fi
echo "------------------------------------------------------------------"

# ===========================================================================
# Step state and counter
# ===========================================================================

STEP_NUM=0

function increment_step {
    STEP_NUM=$((STEP_NUM + 1))
}

function get_step {
    echo "$STEP_NUM"
}


extract_brain_t1() {
    local input="$1"
    local output="$2"

    #robustfov -i "$input" -r "$input"_cropped.nii.gz
    #bet "$input"_cropped.nii.gz "$output" -B -f 0.1

    # Previously "> /dev/null 2>&1" -- discarded synthstrip's own output
    # entirely, so a killed/failed strip (e.g. OOM'd by node-level memory
    # contention, confirmed 2026-08-15 for OASIS3/EPAD's healthy-cohort T1
    # prep) left a completely empty z_score_error.log with zero diagnostic
    # trail. Letting it inherit stdout/stderr routes it into the same
    # z_score.log/z_score_error.log the calling Snakemake rule already sets up.
    mri_synthstrip -i "$input" -o "$output"
}

extract_brain_dwi() {
    local input="$1"
    local output="$2"

    #bet "$input" "$output" -f 0.3
    mri_synthstrip -i "$input" -o "$output"
}

PREPARED_DWIS=()
PREPARED_BVALS=()
PREPARED_TARGET_METRICS=()
HEALTHY_T1_TO_TARGET_T1_CPPS=()
HEALTHY_DIFF2T1_AFFINES=()
PREPARED_HEALTHY_T1S=()
PREPARED_TARGET_T1=""
PREPARED_TARGET_DWI=""
PREPARED_TARGET_BVAL=""
TARGET_B0_BRAIN=""
TARGET_DIFF2T1_AFFINE=""
TARGET_T12DIFF_AFFINE=""
DESIGN_DIR="$TMP_DIR/design_matrix"
DESIGN_MAT="$DESIGN_DIR/design.mat"
TARGET_VALS="$DESIGN_DIR/target_regressors.txt"
XFM_DIR="$TMP_DIR/xfm"

# Find the first path whose full path contains the subject ID (case-insensitive).
find_subject_path() {
    local subject_id="$1"
    shift
    local id_lc
    local candidate
    local candidate_lc

    id_lc=$(echo "$subject_id" | tr '[:upper:]' '[:lower:]')
    for candidate in "$@"; do
        candidate_lc=$(echo "$candidate" | tr '[:upper:]' '[:lower:]')
        if [[ "$candidate_lc" == *"$id_lc"* ]]; then
            echo "$candidate"
            return 0
        fi
    done

    return 1
}

prepare_dwi_inputs() {
    if [[ "$METRIC_SPACE" == "t1" ]]; then
        return 0
    fi

    increment_step
    echo "Step $(get_step): Preparing healthy DWI inputs (skipping target) ..."

    if [[ ${#DWI_PATHS[@]} -eq 0 ]]; then
        echo "Error: --dwi-paths required when --metric-space=dwi" >&2
        exit 1
    fi

    local DWI_PREP_DIR
    local HEALTHY_ID
    local DWI_INPUT
    local BVAL_INPUT
    local BVAL_MATCH
    local IDX
    local DWI_FINAL_DIR
    local DWI_TMP_DIR
    local PREP_CMD

    # Shared across every target subject's run -- see prepare_dwi_inputs'
    # comment above: this step never references the target at all, so its
    # output is identical no matter which target job produces it.
    DWI_PREP_DIR="$SHARED_CACHE_DIR/dwi_prep"
    mkdir -p "$DWI_PREP_DIR"

    # Process only healthy subjects.
    # Target is intentionally skipped: Healthy DIFF -> Healthy T1 -> Target T1.
    for (( i=0; i<N_HEALTHY; i++ )); do
        HEALTHY_ID=$(sed -n "$((i+1))p" "$HEALTHY_IDS_FILE")
        IDX=$(printf "%04d" "$i")

        if ! DWI_INPUT=$(find_subject_path "$HEALTHY_ID" "${DWI_PATHS[@]}"); then
            echo "Error: could not find DWI path for healthy subject '$HEALTHY_ID' in --dwi-paths" >&2
            exit 1
        fi

        BVAL_INPUT=""
        if [[ ${#BVAL_PATHS[@]} -gt 0 ]]; then
            if BVAL_MATCH=$(find_subject_path "$HEALTHY_ID" "${BVAL_PATHS[@]}"); then
                BVAL_INPUT="$BVAL_MATCH"
            else
                echo "Warning: no bval matched for '$HEALTHY_ID'; dwi_utils.py will try to infer bvals" >&2
            fi
        fi

        # Keyed by the real subject ID, not loop position -- stable and
        # self-documenting regardless of how --healthy-ids happens to be
        # ordered on any given run.
        DWI_FINAL_DIR="$DWI_PREP_DIR/${HEALTHY_ID}"
        PREPARED_DWI="$DWI_FINAL_DIR/data.nii.gz"
        PREPARED_BVAL="$DWI_FINAL_DIR/bvals"

        if [[ -f "$PREPARED_DWI" && -f "$PREPARED_BVAL" ]]; then
            echo "  [$IDX] [SKIP] DWI already prepared for healthy subject: ${HEALTHY_ID} ($((i+1))/${N_HEALTHY})"
            PREPARED_DWIS+=("$PREPARED_DWI")
            PREPARED_BVALS+=("$PREPARED_BVAL")
            continue
        fi

        echo "  [$IDX] Preparing DWI for healthy subject: ${HEALTHY_ID} ($((i+1))/${N_HEALTHY})"

        # Every target subject's job processes the same healthy cohort and
        # can race on this exact cache entry. Prepare into a job-unique temp
        # dir first, then atomically rename into place (same filesystem, so
        # this is a single rename() syscall) -- no other job can ever
        # observe a partially-written file here. If another job's rename
        # already landed first, ours fails cleanly and we just discard our
        # own redundant copy: wasted compute once, never a corrupted read.
        DWI_TMP_DIR="$DWI_PREP_DIR/.tmp_${HEALTHY_ID}_$$_${RANDOM}"
        mkdir -p "$DWI_TMP_DIR"

        PREP_CMD="python3 \"$DWI_UTILS_SCRIPT\" --dwi \"$DWI_INPUT\" --outdir \"$DWI_TMP_DIR\""
        if [[ -n "$BVAL_INPUT" ]]; then
            PREP_CMD="$PREP_CMD --bvals \"$BVAL_INPUT\""
        fi

        #echo "    Running: $PREP_CMD"
        eval "$PREP_CMD"

        if [[ ! -f "$DWI_TMP_DIR/data.nii.gz" || ! -f "$DWI_TMP_DIR/bvals" ]]; then
            echo "ERROR: Prepared DWI or bvals not found in $DWI_TMP_DIR" >&2
            rm -rf "$DWI_TMP_DIR"
            exit 1
        fi

        # HEALTHY_ID may itself contain a "/" (subject/session format, e.g.
        # EPAD), making DWI_FINAL_DIR a nested path -- mkdir -p its parent
        # first, or mv fails with "No such file or directory" and gets
        # misread below as "another job already cached this", silently
        # discarding a successfully-prepared DWI (confirmed root cause of
        # EPAD's healthy-cohort DWI-prep never completing, 2026-08-14).
        mkdir -p "$(dirname "$DWI_FINAL_DIR")"
        if [[ -d "$DWI_FINAL_DIR" ]] || ! mv "$DWI_TMP_DIR" "$DWI_FINAL_DIR" 2>/dev/null; then
            echo "  [$IDX] Another job already cached this healthy subject's DWI -- discarding redundant copy."
            rm -rf "$DWI_TMP_DIR"
        fi

        PREPARED_DWIS+=("$PREPARED_DWI")
        PREPARED_BVALS+=("$PREPARED_BVAL")

        if [[ ! -f "${PREPARED_DWIS[$i]}" ]]; then
            echo "ERROR: Prepared DWI not found: ${PREPARED_DWIS[$i]}" >&2
            exit 1
        fi
        if [[ ! -f "${PREPARED_BVALS[$i]}" ]]; then
            echo "ERROR: Prepared bvals not found: ${PREPARED_BVALS[$i]}" >&2
            exit 1
        fi
    done

    echo "  All healthy DWI inputs prepared successfully."
}

prepare_target_dwi_transforms() {
    if [[ "$METRIC_SPACE" != "dwi" && "$OUTPUT_SPACE" != "dwi" ]]; then
        return 0
    fi

    sub_step=$((sub_step + 1))
    echo "Step $(get_step).$sub_step: Preparing target DWI transforms ..."

    if [[ ${#DWI_PATHS[@]} -eq 0 ]]; then
        echo "Error: --dwi-paths required when metric-space or output-space is dwi" >&2
        exit 1
    fi

    local TARGET_DWI_INPUT
    local TARGET_BVAL_INPUT
    local TARGET_BVAL_MATCH
    local TARGET_DWI_OUT_DIR
    local TARGET_B0_4D
    local TARGET_B0_FIRST
    local PREP_CMD

    if ! TARGET_DWI_INPUT=$(find_subject_path "$TARGET_ID" "${DWI_PATHS[@]}"); then
        echo "Error: could not find DWI path for target subject '$TARGET_ID' in --dwi-paths" >&2
        exit 1
    fi

    TARGET_BVAL_INPUT=""
    if [[ ${#BVAL_PATHS[@]} -gt 0 ]]; then
        if TARGET_BVAL_MATCH=$(find_subject_path "$TARGET_ID" "${BVAL_PATHS[@]}"); then
            TARGET_BVAL_INPUT="$TARGET_BVAL_MATCH"
        else
            echo "Warning: no bval matched for target '$TARGET_ID'; dwi_utils.py will try to infer bvals" >&2
        fi
    fi

    TARGET_DWI_OUT_DIR="$TMP_DIR/dwi_prep/target"
    mkdir -p "$TARGET_DWI_OUT_DIR"
    PREPARED_TARGET_DWI="$TARGET_DWI_OUT_DIR/data.nii.gz"
    PREPARED_TARGET_BVAL="$TARGET_DWI_OUT_DIR/bvals"

    if [[ -f "$PREPARED_TARGET_DWI" && -f "$PREPARED_TARGET_BVAL" ]]; then
        echo "  [SKIP] Target DWI already prepared"
    else
        PREP_CMD="python3 \"$DWI_UTILS_SCRIPT\" --dwi \"$TARGET_DWI_INPUT\" --outdir \"$TARGET_DWI_OUT_DIR\""
        if [[ -n "$TARGET_BVAL_INPUT" ]]; then
            PREP_CMD="$PREP_CMD --bvals \"$TARGET_BVAL_INPUT\""
        fi
        eval "$PREP_CMD"
    fi

    if [[ ! -f "$PREPARED_TARGET_DWI" || ! -f "$PREPARED_TARGET_BVAL" ]]; then
        echo "Error: target DWI preparation failed" >&2
        exit 1
    fi

    TARGET_B0_4D="$TARGET_DWI_OUT_DIR/target_b0_4d.nii.gz"
    TARGET_B0_FIRST="$TARGET_DWI_OUT_DIR/target_b0_first.nii.gz"
    TARGET_B0_BRAIN="$TARGET_DWI_OUT_DIR/target_b0_brain.nii.gz"
    TARGET_DIFF2T1_AFFINE="$TARGET_DWI_OUT_DIR/target_diff2t1_affine.txt"
    TARGET_T12DIFF_AFFINE="$TARGET_DWI_OUT_DIR/target_t12diff_affine.txt"
    TARGET_B0_BRAIN_IN_T1="$TARGET_DWI_OUT_DIR/target_b0_brain_in_t1.nii.gz"

    if [[ -f "$TARGET_DIFF2T1_AFFINE" && -f "$TARGET_T12DIFF_AFFINE" && -f "$TARGET_B0_BRAIN" ]]; then
        echo "  [SKIP] Target DWI<->T1 transforms already exist"
        return 0
    fi

    select_dwi_vols "$PREPARED_TARGET_DWI" "$PREPARED_TARGET_BVAL" "$TARGET_B0_4D" 0 > /dev/null 2>&1
    fslroi "$TARGET_B0_4D" "$TARGET_B0_FIRST" 0 1

    if [[ "$SKIP_SKULLSTRIP_DWI" -eq 1 ]]; then
        cp "$TARGET_B0_FIRST" "$TARGET_B0_BRAIN"
    else
        extract_brain_dwi "$TARGET_B0_FIRST" "$TARGET_B0_BRAIN"
    fi

    reg_aladin \
        -ref "$PREPARED_TARGET_T1" \
        -flo "$TARGET_B0_BRAIN" \
        -aff "$TARGET_DIFF2T1_AFFINE" \
        -res "$TARGET_B0_BRAIN_IN_T1" \
        -omp "$THREADS" \
        -voff > /dev/null 2>&1

    reg_transform -invAff "$TARGET_DIFF2T1_AFFINE" "$TARGET_T12DIFF_AFFINE" > /dev/null 2>&1

    if [[ ! -f "$TARGET_DIFF2T1_AFFINE" || ! -f "$TARGET_T12DIFF_AFFINE" ]]; then
        echo "Error: target DWI/T1 transform generation failed" >&2
        exit 1
    fi

    echo "  Target DWI transforms prepared successfully."
}

prepare_target_metrics_for_output_space() {
    sub_step=$((sub_step + 1))
    echo "Step $(get_step).$sub_step: Preparing target metrics in output space ..."

    local METRIC_NAME
    local TARGET_METRIC_RAW
    local METRIC_TMP
    local TARGET_METRIC_PREP

    PREPARED_TARGET_METRICS=()

    for (( m=0; m<N_METRICS; m++ )); do
        METRIC_NAME="${METRIC_NAMES[$m]}"
        TARGET_METRIC_RAW="${TARGET_METRICS[$m]}"
        METRIC_TMP="$TMP_DIR/$METRIC_NAME"
        mkdir -p "$METRIC_TMP"

        if [[ "$METRIC_SPACE" == "$OUTPUT_SPACE" ]]; then
            PREPARED_TARGET_METRICS+=("$TARGET_METRIC_RAW")
            continue
        fi

        TARGET_METRIC_FOLDER="$METRIC_TMP/target"
        mkdir -p "$TARGET_METRIC_FOLDER"

        TARGET_METRIC_PREP="$TARGET_METRIC_FOLDER/target_in_${OUTPUT_SPACE}_${METRIC_NAME}.nii.gz"

        if [[ -f "$TARGET_METRIC_PREP" ]]; then
            echo "  [SKIP] Target metric already in ${OUTPUT_SPACE} space: $METRIC_NAME"
            PREPARED_TARGET_METRICS+=("$TARGET_METRIC_PREP")
            continue
        fi

        if [[ "$METRIC_SPACE" == "dwi" && "$OUTPUT_SPACE" == "t1" ]]; then
            reg_resample \
                -ref "$PREPARED_TARGET_T1" \
                -flo "$TARGET_METRIC_RAW" \
                -trans "$TARGET_DIFF2T1_AFFINE" \
                -res "$TARGET_METRIC_PREP" \
                -inter 1 \
                -voff > /dev/null 2>&1
        elif [[ "$METRIC_SPACE" == "t1" && "$OUTPUT_SPACE" == "dwi" ]]; then
            reg_resample \
                -ref "$TARGET_B0_BRAIN" \
                -flo "$TARGET_METRIC_RAW" \
                -trans "$TARGET_T12DIFF_AFFINE" \
                -res "$TARGET_METRIC_PREP" \
                -inter 1 \
                -voff > /dev/null 2>&1
        else
            echo "Error: unsupported metric-space/output-space combination: $METRIC_SPACE -> $OUTPUT_SPACE" >&2
            exit 1
        fi

        if [[ ! -f "$TARGET_METRIC_PREP" ]]; then
            echo "Error: target metric conversion failed for $METRIC_NAME" >&2
            exit 1
        fi

        PREPARED_TARGET_METRICS+=("$TARGET_METRIC_PREP")
    done
}

generate_design_matrix() {
    increment_step
    echo "Step $(get_step): Generating design matrix for GLM ..."
    if [[ -z "$DEMO_CSV" ]]; then
        echo "  [SKIP] Demographics CSV not provided. Falling back to simple Z-score."
        return 0
    fi
    mkdir -p "$DESIGN_DIR"
    local DM_CMD=(
        python3 "$UTILS_SCRIPT" generate-design-matrix
        --demo        "$DEMO_CSV"
        --target-id   "$TARGET_ID"
        --healthy-ids "$HEALTHY_IDS_FILE"
        --covariates  "$COVARIATES"
        --poly-terms  "$POLY_TERMS"
        --design-mat  "$DESIGN_MAT"
        --target-vals "$TARGET_VALS"
    )

    if [[ "$VERBOSE" -eq 1 ]]; then
        DM_CMD+=(--verbose)
    fi

    if ! "${DM_CMD[@]}"; then
        echo "  [WARNING] Design matrix generation failed. Falling back to simple Z-score."
        rm -f "$DESIGN_MAT" "$TARGET_VALS"
    fi
}

prepare_t1_inputs() {

    sub_step=$((sub_step + 1))
    echo "Step $(get_step).$sub_step: Preparing T1 inputs (skull stripping) ..."

    local T1_PREP_DIR
    local HEALTHY_T1_PREP_DIR
    local IDX
    local HEALTHY_ID
    local HEALTHY_T1
    local HEALTHY_T1_BRAIN
    local HEALTHY_T1_BRAIN_TMP

    PREPARED_HEALTHY_T1S=()
    T1_PREP_DIR="$TMP_DIR/t1_prep"
    mkdir -p "$T1_PREP_DIR"
    # Shared across every target subject's run, same reasoning as
    # prepare_dwi_inputs above: skull-stripping a healthy subject's own T1
    # doesn't depend on the target at all. Only the target's own T1 (right
    # above) is genuinely per-target and stays in $TMP_DIR.
    HEALTHY_T1_PREP_DIR="$SHARED_CACHE_DIR/t1_prep"
    mkdir -p "$HEALTHY_T1_PREP_DIR"

    if [[ "$SKIP_SKULLSTRIP_T1" -eq 1 ]]; then
        echo "  Skipping T1 skull stripping for target and healthy T1s."
        PREPARED_TARGET_T1="$TARGET_T1"
        PREPARED_HEALTHY_T1S=("${HEALTHY_T1S[@]}")
        return 0
    fi

    PREPARED_TARGET_T1="$T1_PREP_DIR/target_t1_brain.nii.gz"
    if [[ -f "$PREPARED_TARGET_T1" ]]; then
        echo "  [SKIP] Target T1 already skull-stripped: $PREPARED_TARGET_T1"
    else
        echo "  Skull stripping target T1 for subject: $TARGET_ID ..."
        extract_brain_t1 "$TARGET_T1" "$PREPARED_TARGET_T1"
        if [[ ! -f "$PREPARED_TARGET_T1" ]]; then
            echo "ERROR: Skull-stripped target T1 not found: $PREPARED_TARGET_T1" >&2
            exit 1
        fi
    fi

    for (( i=0; i<N_HEALTHY; i++ )); do
        IDX=$(printf "%04d" "$i")
        HEALTHY_ID=$(sed -n "$((i+1))p" "$HEALTHY_IDS_FILE")
        HEALTHY_T1="${HEALTHY_T1S[$i]}"
        # Keyed by real subject ID, same reasoning as prepare_dwi_inputs.
        HEALTHY_T1_BRAIN="$HEALTHY_T1_PREP_DIR/${HEALTHY_ID}_healthy_t1_brain.nii.gz"

        if [[ -f "$HEALTHY_T1_BRAIN" ]]; then
            echo "  [$IDX] [SKIP] Healthy T1 already processed: ${HEALTHY_ID} ($((i+1))/${N_HEALTHY})"
        else
            echo "  [$IDX] Skull stripping healthy T1 for subject: ${HEALTHY_ID} ($((i+1))/${N_HEALTHY})"
            # Same race as prepare_dwi_inputs: strip into a job-unique temp
            # file, atomically rename into the shared cache. A regular-file
            # rename() is just as atomic as a directory one, so a losing
            # job's mv fails cleanly and it discards its own redundant copy.
            HEALTHY_T1_BRAIN_TMP="$HEALTHY_T1_PREP_DIR/.tmp_${HEALTHY_ID}_$$_${RANDOM}.nii.gz"
            # HEALTHY_ID contains a literal "/" for session-nested datasets
            # (e.g. "sub-OAS30005/ses-d1274"), turning the line above into a
            # multi-component path whose intermediate directory was never
            # created -- unlike prepare_dwi_inputs' equivalent temp path
            # (a directory, always mkdir -p'd), this one is a single output
            # FILE, so nobody mkdir -p'd its parent before synthstrip tried
            # to write it. Confirmed 2026-08-15: this crashed every single
            # OASIS3/EPAD attempt with FileNotFoundError on the temp path
            # (ADNI3 never hit it -- its healthy IDs are flat, no "/").
            mkdir -p "$(dirname "$HEALTHY_T1_BRAIN_TMP")"
            extract_brain_t1 "$HEALTHY_T1" "$HEALTHY_T1_BRAIN_TMP"
            if [[ ! -f "$HEALTHY_T1_BRAIN_TMP" ]]; then
                echo "ERROR: Skull-stripped healthy T1 not found: $HEALTHY_T1_BRAIN_TMP" >&2
                exit 1
            fi
            # Same nested-HEALTHY_ID issue as prepare_dwi_inputs -- see the
            # comment there.
            mkdir -p "$(dirname "$HEALTHY_T1_BRAIN")"
            if [[ -f "$HEALTHY_T1_BRAIN" ]] || ! mv "$HEALTHY_T1_BRAIN_TMP" "$HEALTHY_T1_BRAIN" 2>/dev/null; then
                echo "  [$IDX] Another job already cached this healthy subject's T1 -- discarding redundant copy."
                rm -f "$HEALTHY_T1_BRAIN_TMP"
            fi
            if [[ ! -f "$HEALTHY_T1_BRAIN" ]]; then
                echo "ERROR: Skull-stripped healthy T1 not found: $HEALTHY_T1_BRAIN" >&2
                exit 1
            fi
        fi

        PREPARED_HEALTHY_T1S+=("$HEALTHY_T1_BRAIN")
    done
}

register_dwi_to_t1() {
    if [[ "$METRIC_SPACE" != "dwi" ]]; then
        return 0
    fi

    # Registration: healthy DWI b0 space -> healthy T1 space
    sub_step=$((sub_step + 1))
    echo "Step $(get_step).$sub_step: Extracting b0 and registering DWI b0 to healthy T1 ..."

    local DIFF2T1_DIR
    local HEALTHY_ID
    local HEALTHY_DWI
    local HEALTHY_BVAL
    local HEALTHY_T1
    local IDX
    local B0_4D
    local B0_FIRST
    local B0_BRAIN
    local DIFF2T1_AFF
    local HEALTY_DIFF2T1_FOLDER
    local HEALTY_DIFF2T1_FOLDER_TMP

    # Shared across every target subject's run: b0 extraction, b0 skull-
    # stripping, and b0->healthy-T1 affine registration only ever touch
    # healthy-subject inputs (PREPARED_DWIS/PREPARED_BVALS/PREPARED_HEALTHY_T1S)
    # -- no TARGET_* variable appears anywhere in this loop -- so the output
    # is identical no matter which target job produces it. Same reasoning
    # and race-safety pattern as prepare_dwi_inputs/prepare_t1_inputs above.
    DIFF2T1_DIR="$SHARED_CACHE_DIR/diff2t1"
    mkdir -p "$DIFF2T1_DIR"
    HEALTHY_DIFF2T1_AFFINES=()

    for (( i=0; i<N_HEALTHY; i++ )); do
        HEALTHY_ID=$(sed -n "$((i+1))p" "$HEALTHY_IDS_FILE")
        HEALTHY_DWI="${PREPARED_DWIS[$i]}"
        HEALTHY_BVAL="${PREPARED_BVALS[$i]}"
        HEALTHY_T1="${PREPARED_HEALTHY_T1S[$i]}"
        IDX=$(printf "%04d" "$i")

        # Keyed by the real subject ID, not loop position -- same reasoning
        # as prepare_dwi_inputs.
        HEALTY_DIFF2T1_FOLDER="$DIFF2T1_DIR/${HEALTHY_ID}"
        HEALTHY_DIFF2T1_AFF="$HEALTY_DIFF2T1_FOLDER/${IDX}_diff2t1_affine.txt"

        if [[ -f "$HEALTHY_DIFF2T1_AFF" ]]; then
            echo "  [$IDX] [SKIP] DWI b0 to T1 already registered: ${HEALTHY_ID} ($((i+1))/${N_HEALTHY})"
        else
            echo "  [$IDX] Processing DWI for: ${HEALTHY_ID} ($((i+1))/${N_HEALTHY}) ($((i+1))/${N_HEALTHY})"

            # Every target subject's job processes the same healthy cohort
            # and can race on this exact cache entry. Prepare into a
            # job-unique temp dir first, then atomically rename into place
            # (same filesystem, so this is a single rename() syscall) -- no
            # other job can ever observe a partially-written folder here.
            HEALTY_DIFF2T1_FOLDER_TMP="$DIFF2T1_DIR/.tmp_${HEALTHY_ID}_$$_${RANDOM}"
            mkdir -p "$HEALTY_DIFF2T1_FOLDER_TMP"

            HEALTHY_B0_4D="$HEALTY_DIFF2T1_FOLDER_TMP/${IDX}_b0_4d.nii.gz"
            HEALTHY_B0_FIRST="$HEALTY_DIFF2T1_FOLDER_TMP/${IDX}_b0_first.nii.gz"
            HEALTHY_B0_BRAIN="$HEALTY_DIFF2T1_FOLDER_TMP/${IDX}_b0_brain.nii.gz"
            HEALTHY_DIFF2T1_AFF_TMP="$HEALTY_DIFF2T1_FOLDER_TMP/${IDX}_diff2t1_affine.txt"
            HEALTHY_B0_BRAIN_IN_T1="$HEALTY_DIFF2T1_FOLDER_TMP/${IDX}_b0_brain_in_t1_aladin.nii.gz"

            #echo "    Extracting b0 volumes ..."
            select_dwi_vols "$HEALTHY_DWI" "$HEALTHY_BVAL" "$HEALTHY_B0_4D" 0 > /dev/null 2>&1

            #echo "    Extracting first b0 volume ..."
            fslroi "$HEALTHY_B0_4D" "$HEALTHY_B0_FIRST" 0 1

            if [[ "$SKIP_SKULLSTRIP_DWI" -eq 1 ]]; then
                #echo "    Skipping b0 skull stripping (using first b0 directly) ..."
                cp "$HEALTHY_B0_FIRST" "$HEALTHY_B0_BRAIN"
            else
                #echo "    Brain extracting b0 ..."
                extract_brain_dwi "$HEALTHY_B0_FIRST" "$HEALTHY_B0_BRAIN"
            fi

            #echo "    Registering b0 to T1 (affine) ..."
            reg_aladin \
                -ref "$HEALTHY_T1"  \
                -flo "$HEALTHY_B0_BRAIN"    \
                -aff "$HEALTHY_DIFF2T1_AFF_TMP" \
                -res "$HEALTHY_B0_BRAIN_IN_T1" \
                -omp "$THREADS"     \
                -voff > /dev/null 2>&1

            if [[ ! -f "$HEALTHY_DIFF2T1_AFF_TMP" ]]; then
                echo "ERROR: diff2t1 affine not found: $HEALTHY_DIFF2T1_AFF_TMP" >&2
                rm -rf "$HEALTY_DIFF2T1_FOLDER_TMP"
                exit 1
            fi

            if [[ "$QC" -eq 1 ]]; then
                # Store Healthy b0 brain / T1 in T1 space for QC, inside the
                # same job-unique temp dir so they move atomically together.
                reg_resample \
                    -ref "$HEALTHY_T1"  \
                    -flo "$HEALTHY_B0_BRAIN" \
                    -trans "$HEALTHY_DIFF2T1_AFF_TMP" \
                    -res "$HEALTY_DIFF2T1_FOLDER_TMP/${IDX}_b0_brain_in_t1.nii.gz" \
                    -inter 1 \
                    -voff > /dev/null 2>&1
                cp "$HEALTHY_T1" "$HEALTY_DIFF2T1_FOLDER_TMP/${IDX}_t1.nii.gz"
            fi

            # Same nested-HEALTHY_ID issue as prepare_dwi_inputs/prepare_t1_inputs
            # above -- HEALTHY_ID can contain a "/" for session-nested datasets,
            # making HEALTY_DIFF2T1_FOLDER's parent directory not exist yet.
            # Confirmed 2026-08-15: this made OASIS3's sub-OAS30001 die here
            # with "diff2t1 affine not found" right after the T1 fix let it
            # progress past skull-stripping.
            mkdir -p "$(dirname "$HEALTY_DIFF2T1_FOLDER")"
            if [[ -d "$HEALTY_DIFF2T1_FOLDER" ]] || ! mv "$HEALTY_DIFF2T1_FOLDER_TMP" "$HEALTY_DIFF2T1_FOLDER" 2>/dev/null; then
                echo "  [$IDX] Another job already cached this healthy subject's b0->T1 registration -- discarding redundant copy."
                rm -rf "$HEALTY_DIFF2T1_FOLDER_TMP"
            fi

            if [[ ! -f "$HEALTHY_DIFF2T1_AFF" ]]; then
                echo "ERROR: diff2t1 affine not found: $HEALTHY_DIFF2T1_AFF" >&2
                exit 1
            fi
            sync
        fi

        HEALTHY_DIFF2T1_AFFINES+=("$HEALTHY_DIFF2T1_AFF")
    done

    echo "  All b0 registrations complete."
    sync
}

qc_register_nodiff_to_t1_folder() {
    if [[ "$METRIC_SPACE" != "dwi" ]]; then
        return 0
    fi

    # QC step: write explicit nodiff->target-T1 resampled images in a dedicated folder.
    sub_step=$((sub_step + 1))
    echo "Step $(get_step).$sub_step (QC): Resampling nodiff brain images to target T1 space ..."

    local NODIFF_QC_DIR
    local HEALTHY_T1
    local HEALTHY_ID
    local IDX
    local HEALTHY_B0_BRAIN
    local DIFF2T1_AFF
    local HEALTHY_T1_TO_TARGET_T1_CPP
    local COMPOSED
    local NODIFF_IN_TARGET_T1
    local reg_transform_cmd

    NODIFF_QC_DIR="$TMP_DIR/nodiff_to_target_t1_qc"
    mkdir -p "$NODIFF_QC_DIR"

    for (( i=0; i<N_HEALTHY; i++ )); do
        HEALTHY_T1="${PREPARED_HEALTHY_T1S[$i]}"
        HEALTHY_ID=$(sed -n "$((i+1))p" "$HEALTHY_IDS_FILE")
        IDX=$(printf "%04d" "$i")

        # register_dwi_to_t1 now writes this into the shared cache, keyed by
        # real subject ID (see that function's comment for why this step is
        # target-independent).
        HEALTY_DIFF2T1_FOLDER="$SHARED_CACHE_DIR/diff2t1/${HEALTHY_ID}"
        HEALTHY_B0_BRAIN="$HEALTY_DIFF2T1_FOLDER/${IDX}_b0_brain.nii.gz"

        NODIFF_QC_SUBJECT_DIR="$NODIFF_QC_DIR/${IDX}"
        mkdir -p "$NODIFF_QC_SUBJECT_DIR"
        
        COMPOSED="$NODIFF_QC_SUBJECT_DIR/${IDX}_nodiff_to_target_composed.nii.gz"
        NODIFF_IN_TARGET_T1="$NODIFF_QC_SUBJECT_DIR/${IDX}_nodiff_brain_in_target_t1.nii.gz"

        if [[ -f "$NODIFF_IN_TARGET_T1" ]]; then
            echo "  [$IDX] [SKIP] QC nodiff->target T1 already exists ($((i+1))/${N_HEALTHY})"
            continue
        fi
        echo "  [$IDX] QC: Resampling nodiff brain to target T1 space for ${HEALTHY_ID} ($((i+1))/${N_HEALTHY}) ..."

        if [[ ! -f "$HEALTHY_B0_BRAIN" ]]; then
            echo "ERROR: nodiff brain image not found for QC step: $HEALTHY_B0_BRAIN" >&2
            exit 1
        fi
        if [[ ! -f "${HEALTHY_T1_TO_TARGET_T1_CPPS[$i]}" ]]; then
            echo "ERROR: healthy->target CPP missing for QC step: ${HEALTHY_T1_TO_TARGET_T1_CPPS[$i]}" >&2
            exit 1
        fi

        #echo "  [$IDX] Composing nodiff->target transform for QC ($((i+1))/${N_HEALTHY})"
        reg_transform_cmd="reg_transform -comp \"${HEALTHY_DIFF2T1_AFFINES[$i]}\" \"${HEALTHY_T1_TO_TARGET_T1_CPPS[$i]}\" \"$COMPOSED\" -ref \"$HEALTHY_T1\""
        #echo "    Running reg_transform with command: $reg_transform_cmd"
        eval "$reg_transform_cmd" > /dev/null 2>&1

        #echo "  [$IDX] Resampling nodiff brain to target T1 for QC ($((i+1))/${N_HEALTHY})"
        reg_resample \
            -ref "$PREPARED_TARGET_T1" \
            -flo "$HEALTHY_B0_BRAIN" \
            -trans "$COMPOSED" \
            -res "$NODIFF_IN_TARGET_T1" \
            -inter 1 \
            -voff > /dev/null 2>&1

        if [[ ! -f "$NODIFF_IN_TARGET_T1" ]]; then
            echo "ERROR: nodiff QC output missing: $NODIFF_IN_TARGET_T1" >&2
            exit 1
        fi
        
        # Force filesystem sync to prevent stale file handles on network mounts
        sync
    done

    #echo "  QC nodiff->T1 outputs written to: $NODIFF_QC_DIR"
    sync
}

register_t1_to_target_t1() {
    # Registration: healthy T1 space -> target T1 space
    sub_step=$((sub_step + 1))
    echo "Step $(get_step).$sub_step: Registering healthy T1s to target T1 space ..."

    local HEALTHY_ID
    local HEALTHY_T1
    local IDX
    local HEALTHY_T1_TO_TARGET_T1_AFFINE
    local CPP
    local RES_T1

    mkdir -p "$XFM_DIR"
    HEALTHY_T1_TO_TARGET_T1_CPPS=()

    for (( i=0; i<N_HEALTHY; i++ )); do
        HEALTHY_ID=$(sed -n "$((i+1))p" "$HEALTHY_IDS_FILE")
        HEALTHY_T1="${PREPARED_HEALTHY_T1S[$i]}"
        IDX=$(printf "%04d" "$i")

        XFM_SUBJECT_DIR="$XFM_DIR/${IDX}"
        mkdir -p "$XFM_SUBJECT_DIR"
        
        HEALTHY_T1_TO_TARGET_T1_AFFINE="$XFM_SUBJECT_DIR/${IDX}_affine.txt"
        CPP="$XFM_SUBJECT_DIR/${IDX}_cpp.nii.gz"
        RES_T1="$XFM_SUBJECT_DIR/${IDX}_res_t1.nii.gz"
        RES_T1_ALADIN="$XFM_SUBJECT_DIR/${IDX}_res_t1_aladin.nii.gz"

        if [[ -f "$CPP" ]]; then
            echo "  [${IDX}] [SKIP] T1 to target T1 already registered: ${HEALTHY_ID}"
            HEALTHY_T1_TO_TARGET_T1_CPPS+=("$CPP")
            continue
        fi

        echo "  [${IDX}] Registering healthy subject: ${HEALTHY_ID} ($((i+1))/${N_HEALTHY}) ..."

        #echo "    Step $(get_step).$((sub_step))a: Running reg_aladin (affine) ..."
        reg_aladin \
            -ref "$PREPARED_TARGET_T1"  \
            -flo "$HEALTHY_T1" \
            -aff "$HEALTHY_T1_TO_TARGET_T1_AFFINE"     \
            -res "$RES_T1_ALADIN" \
            -omp "$THREADS"    \
            -voff > /dev/null 2>&1

        #echo "    Step $(get_step).$((sub_step))b: Running reg_f3d (non-linear CPP) ..."
        reg_f3d \
            -ref "$PREPARED_TARGET_T1"  \
            -flo "$HEALTHY_T1" \
            -aff "$HEALTHY_T1_TO_TARGET_T1_AFFINE"     \
            -cpp "$CPP"        \
            -res "$RES_T1"     \
            -omp "$THREADS"    \
            -voff > /dev/null 2>&1

        if [[ ! -f "$CPP" ]]; then
            echo "ERROR: CPP file missing after registration: $CPP" >&2
            exit 1
        fi
        sync

        HEALTHY_T1_TO_TARGET_T1_CPPS+=("$CPP")

        if [[ "$QC" -eq 1 ]]; then
            echo "    QC: Checking Healthy T1 to target T1 registration for ${HEALTHY_ID} ..."
            # Copy the resampled T1 to the QC folder for visual inspection.
            QC_SUBJECT_DIR="$TMP_DIR/t1_to_target_t1_qc/${IDX}"
            mkdir -p "$QC_SUBJECT_DIR"
            QC_T1_IN_TARGET="$QC_SUBJECT_DIR/${IDX}_t1_in_target_t1.nii.gz"
            QC_TARGET_T1="$QC_SUBJECT_DIR/${IDX}_target_t1.nii.gz"

            if [[ -f "$QC_T1_IN_TARGET" && -f "$QC_TARGET_T1" ]]; then
                echo "    [SKIP] QC T1->target already exists for ${HEALTHY_ID}"
                continue
            fi

            cp "$RES_T1" "$QC_T1_IN_TARGET"

            # Copy the target T1 to the QC folder for visual inspection.
            cp "$PREPARED_TARGET_T1" "$QC_TARGET_T1"
        fi

    done

    echo "  All healthy T1s registered to target T1 space."
    sync
}

register_metrics() {
    sub_step=$((sub_step + 1))
    echo "Step $(get_step).$sub_step: Registering healthy metrics to target space ..."

    local METRIC_NAME
    local METRIC_TMP
    local OFFSET
    local HEALTHY_METRIC
    local HEALTHY_T1
    local IDX
    local REG_METRIC
    local COMPOSED
    local reg_transform_cmd
    local reg_resample_cmd
    
    for (( m=0; m<N_METRICS; m++ )); do
        METRIC_NAME="${METRIC_NAMES[$m]}"
        METRIC_TMP="$TMP_DIR/$METRIC_NAME"
        mkdir -p "$METRIC_TMP"

        echo "--- Metric: $METRIC_NAME ---"

        OFFSET=$(( m * N_HEALTHY ))

        for (( i=0; i<N_HEALTHY; i++ )); do
            HEALTHY_METRIC="${HEALTHY_METRICS_ALL[$(( OFFSET + i ))]}"
            HEALTHY_T1="${PREPARED_HEALTHY_T1S[$i]}"
            IDX=$(printf "%04d" "$i")

            REG_METRIC_FOLDER="$METRIC_TMP/${IDX}"
            mkdir -p "$REG_METRIC_FOLDER"
            REG_METRIC="$REG_METRIC_FOLDER/${IDX}_in_target.nii.gz"

            if [[ -f "$REG_METRIC" ]]; then
                echo "    [SKIP] Metric already registered for ${IDX}"
                continue
            fi

            if [[ "$METRIC_SPACE" == "dwi" ]]; then
                # Registration: healthy metric in DWI space -> target T1 space
                COMPOSED_FOLDER="$METRIC_TMP/${IDX}"
                mkdir -p "$COMPOSED_FOLDER"
                COMPOSED="$COMPOSED_FOLDER/${IDX}_composed.nii.gz"

                echo "    Composing transforms for ${IDX}: HEALTHY DWI -> HEALTHY T1 -> TARGET T1"

                reg_transform_cmd="reg_transform -comp \"${HEALTHY_DIFF2T1_AFFINES[$i]}\" \"${HEALTHY_T1_TO_TARGET_T1_CPPS[$i]}\" \"$COMPOSED\" -ref \"$HEALTHY_T1\""
                #echo "    Running reg_transform with command: $reg_transform_cmd"
                eval "$reg_transform_cmd" > /dev/null 2>&1

                reg_resample_cmd="reg_resample -ref \"$PREPARED_TARGET_T1\" -flo \"$HEALTHY_METRIC\" -trans \"$COMPOSED\" -res \"$REG_METRIC\" -inter 1 -voff"
                #echo "    Running reg_resample with command: $reg_resample_cmd"
                eval "$reg_resample_cmd" > /dev/null 2>&1
            else
                # Registration: healthy metric in T1 space -> target T1 space
                echo "    Registering healthy metric to target T1 space for ${IDX}: HEALTHY T1 -> TARGET T1"
                reg_resample_cmd="reg_resample -ref \"$PREPARED_TARGET_T1\" -flo \"$HEALTHY_METRIC\" -trans \"${HEALTHY_T1_TO_TARGET_T1_CPPS[$i]}\" -res \"$REG_METRIC\" -inter 1 -voff"
                #echo "    Running reg_resample with command: $reg_resample_cmd"
                eval "$reg_resample_cmd" > /dev/null 2>&1
            fi
        done
    done

    echo "  All metrics registered to target space."
    sync
}

project_registered_metrics_to_output_space() {
    if [[ "$OUTPUT_SPACE" != "dwi" ]]; then
        return 0
    fi

    sub_step=$((sub_step + 1))
    echo "Step $(get_step).$sub_step: Projecting registered healthy metrics to target DWI space ..."

    local METRIC_NAME
    local METRIC_TMP
    local IDX
    local REG_METRIC_T1
    local REG_METRIC_DWI

    for (( m=0; m<N_METRICS; m++ )); do
        METRIC_NAME="${METRIC_NAMES[$m]}"
        METRIC_TMP="$TMP_DIR/$METRIC_NAME"

        for (( i=0; i<N_HEALTHY; i++ )); do
            IDX=$(printf "%04d" "$i")
            REG_METRIC_T1="$METRIC_TMP/${IDX}/${IDX}_in_target.nii.gz"
            REG_METRIC_DWI="$METRIC_TMP/${IDX}/${IDX}_in_target_dwi.nii.gz"

            if [[ -f "$REG_METRIC_DWI" ]]; then
                echo "    [SKIP] Metric already projected to DWI for ${IDX}"
                continue
            fi

            reg_resample \
                -ref "$TARGET_B0_BRAIN" \
                -flo "$REG_METRIC_T1" \
                -trans "$TARGET_T12DIFF_AFFINE" \
                -res "$REG_METRIC_DWI" \
                -inter 1 \
                -voff > /dev/null 2>&1

            if [[ ! -f "$REG_METRIC_DWI" ]]; then
                echo "Error: projection to target DWI failed for metric $METRIC_NAME subject $IDX" >&2
                exit 1
            fi
        done
    done
}

merge_metrics() {
    sub_step=$((sub_step + 1))
    echo "Step $(get_step).$sub_step: Merging registered healthy metrics ..."

    local METRIC_NAME
    local METRIC_TMP
    local OFFSET
    local IDX
    local REGISTERED_METRICS
    local MERGED

    for (( m=0; m<N_METRICS; m++ )); do
        METRIC_NAME="${METRIC_NAMES[$m]}"
        METRIC_TMP="$TMP_DIR/$METRIC_NAME"

        echo ""
        echo "--- Metric: $METRIC_NAME ---"

        REGISTERED_METRICS=()
        for (( i=0; i<N_HEALTHY; i++ )); do
            IDX=$(printf "%04d" "$i")
            REG_METRIC_FOLDER="$METRIC_TMP/${IDX}"
            if [[ "$OUTPUT_SPACE" == "dwi" ]]; then
                REG_METRIC="$REG_METRIC_FOLDER/${IDX}_in_target_dwi.nii.gz"
            else
                REG_METRIC="$REG_METRIC_FOLDER/${IDX}_in_target.nii.gz"
            fi
            REGISTERED_METRICS+=("$REG_METRIC")
        done

        MERGED_FOLDER="$METRIC_TMP/merged"
        mkdir -p "$MERGED_FOLDER"

        MERGED="$MERGED_FOLDER/healthy_merged_${METRIC_NAME}.nii.gz"
        
        if [[ -f "$MERGED" ]]; then
            echo "  [SKIP] Merged file already exists: $MERGED"
        else
            echo "  Merging ${N_HEALTHY} registered metrics ..."
            fslmerge -t "$MERGED" "${REGISTERED_METRICS[@]}"

            if [[ ! -f "$MERGED" ]]; then
                echo "ERROR: Merged file missing: $MERGED" >&2
                exit 1
            fi
        fi
    done
}

compute_z_scores() {
    increment_step
    echo "Step $(get_step): Computing GLM and Z-scores ..."

    local METRIC_NAME
    local TARGET_METRIC
    local OUTPUT_FILE
    local METRIC_TMP
    local MERGED
    local BETAS
    local RESIDUALS
    local STD_DEV
    local PREDICTED
    local N_REGRESSORS
    # Peak-memory report lines for fsl_glm (the confirmed OOM culprit --
    # see the -m mask above), one per metric, printed as a summary at the
    # end of this function so the actual mem_mb budget in
    # z_score_workflow.smk can be tuned from real numbers instead of a
    # guess. Remove once the brain-mask fix's real headroom is confirmed
    # across a representative run and the resource budget is settled.
    local -a GLM_PEAK_MEM_REPORT=()

    # Brain mask, computed once and reused for every metric (previously
    # recomputed per-metric from the same reference, and only ever applied
    # AFTER fsl_glm ran -- meaning fsl_glm was regressing across the full
    # head/background bounding box of every registered healthy volume, not
    # just brain tissue. Passing it into fsl_glm via -m bounds the GLM's own
    # memory to brain voxels only (a real "fsl_glm: std::bad_alloc" OOM was
    # confirmed here on 2026-08-13; background voxels typically make up the
    # majority of a whole-head bounding box), independent of healthy-cohort
    # size -- the actual fix, rather than just requesting more memory.
    if [[ "$OUTPUT_SPACE" == "dwi" ]]; then
        BRAIN_MASK_REF="$TARGET_B0_BRAIN"
    else
        BRAIN_MASK_REF="$PREPARED_TARGET_T1"
    fi
    BRAIN_MASK="$TMP_DIR/brain_mask.nii.gz"
    fslmaths "$BRAIN_MASK_REF" -bin "$BRAIN_MASK"

    for (( m=0; m<N_METRICS; m++ )); do
        METRIC_NAME="${METRIC_NAMES[$m]}"
        TARGET_METRIC="${PREPARED_TARGET_METRICS[$m]}"

        # Build filename suffix with covariates and poly terms
        FILENAME_SUFFIX=""
        if [[ -n "$COVARIATES" ]]; then
            COV_SAFE=$(echo "$COVARIATES" | tr ',' '_')
            FILENAME_SUFFIX="${FILENAME_SUFFIX}_cov_${COV_SAFE}"
        fi
        if [[ -n "$POLY_TERMS" ]]; then
            POLY_SAFE=$(echo "$POLY_TERMS" | tr ',' '_')
            FILENAME_SUFFIX="${FILENAME_SUFFIX}_poly_${POLY_SAFE}"
        fi

        # Flat under $OUTPUT_DIR (z-score/outputs/), matching every other
        # stage's convention (dti/outputs/fa.nii.gz etc) -- a prior nested
        # Z_SCORES/ subfolder here didn't match this rule's own Snakemake
        # `output:` declaration (z_score_workflow.smk) or process_all's
        # z_score_outputs() (process_all_workflow.smk), silently breaking
        # the z-score -> metrics dependency chain those rely on.
        OUTPUT_FILE="$OUTPUT_DIR/${METRIC_NAME}_z_score${FILENAME_SUFFIX}.nii.gz"
        METRIC_TMP="$TMP_DIR/$METRIC_NAME"
        MERGED_FOLDER="$METRIC_TMP/merged"
        GLM_FOLDER="$METRIC_TMP/GLM"
        
        mkdir -p "$MERGED_FOLDER" "$GLM_FOLDER"

        MERGED="$MERGED_FOLDER/healthy_merged_${METRIC_NAME}.nii.gz"

        if [[ ! -f "$MERGED" ]]; then
            echo "ERROR: merged healthy metric missing before GLM: $MERGED" >&2
            exit 1
        fi

        #echo ""
        echo "--- Metric: $METRIC_NAME ---"
        
        # Skip check: if final output exists, skip entire metric
        if [[ -f "$OUTPUT_FILE" ]]; then
            echo "  [SKIP] Z-score already computed: $OUTPUT_FILE"
            continue
        fi
        
        if [[ ! -f "$TARGET_METRIC" ]]; then
            echo "  ERROR: Target metric file does not exist: $TARGET_METRIC" >&2
            exit 1
        fi

        BETAS="$GLM_FOLDER/betas_${METRIC_NAME}${FILENAME_SUFFIX}.nii.gz"
        RESIDUALS="$GLM_FOLDER/residuals_${METRIC_NAME}${FILENAME_SUFFIX}.nii.gz"
        STD_DEV="$GLM_FOLDER/std_dev_${METRIC_NAME}${FILENAME_SUFFIX}.nii.gz"
        PREDICTED="$GLM_FOLDER/predicted_${METRIC_NAME}${FILENAME_SUFFIX}.nii.gz"

        if [[ -f "$DESIGN_MAT" ]]; then
            echo "  Running fsl_glm (masked to brain voxels only) ..."
            FSL_GLM_TIME_LOG="$GLM_FOLDER/fsl_glm_time_${METRIC_NAME}${FILENAME_SUFFIX}.log"
            /usr/bin/time -v fsl_glm \
                -i "$MERGED"      \
                -d "$DESIGN_MAT"  \
                -m "$BRAIN_MASK"  \
                --out="$BETAS"    \
                --out_res="$RESIDUALS" \
                2> "$FSL_GLM_TIME_LOG"
            GLM_PEAK_KB=$(grep "Maximum resident set size" "$FSL_GLM_TIME_LOG" | awk '{print $NF}')
            if [[ -n "$GLM_PEAK_KB" ]]; then
                echo "  fsl_glm peak memory: $((GLM_PEAK_KB / 1024)) MB"
                GLM_PEAK_MEM_REPORT+=("${METRIC_NAME}${FILENAME_SUFFIX}: $((GLM_PEAK_KB / 1024)) MB")
            fi

            if [[ ! -f "$BETAS" ]]; then
                echo "ERROR: Betas file missing: $BETAS" >&2
                cat "$FSL_GLM_TIME_LOG" >&2
                exit 1
            fi
            if [[ ! -f "$RESIDUALS" ]]; then
                echo "ERROR: Residuals file missing: $RESIDUALS" >&2
                exit 1
            fi

            fslmaths "$RESIDUALS" -Tstd "$STD_DEV"

            echo "  Computing predicted map ..."

            log "  DEBUG: Splitting betas: $BETAS"
            fslsplit "$BETAS" "$GLM_FOLDER/beta_${METRIC_NAME}${FILENAME_SUFFIX}_" -t

            log "  DEBUG: Reading target values from: $TARGET_VALS"
            if [[ ! -f "$TARGET_VALS" ]]; then
                echo "ERROR: Target values file missing: $TARGET_VALS" >&2
                exit 1
            fi
            TARGET_VALS_LINE=$(grep -vE '^[[:space:]]*#|^[[:space:]]*$' "$TARGET_VALS" | head -n 1)
            if [[ -z "$TARGET_VALS_LINE" ]]; then
                echo "ERROR: No numeric target regressors found in: $TARGET_VALS" >&2
                exit 1
            fi
            read -r -a VALS <<< "$TARGET_VALS_LINE" || { echo "ERROR: Failed to parse regressors from $TARGET_VALS" >&2; exit 1; }
            log "  DEBUG: Target values: ${VALS[*]}"
            log "  DEBUG: Number of values: ${#VALS[@]}"

            log "  DEBUG: Creating predicted map starting with beta_${METRIC_NAME}${FILENAME_SUFFIX}_0000.nii.gz * ${VALS[0]}"
            fslmaths "$GLM_FOLDER/beta_${METRIC_NAME}${FILENAME_SUFFIX}_0000.nii.gz" -mul "${VALS[0]}" "$PREDICTED"

            N_REGRESSORS="${#VALS[@]}"
            for (( i=1; i<N_REGRESSORS; i++ )); do
                REGRESSOR_IDX=$(printf "%04d" "$i")
                log "  DEBUG: Adding beta_${METRIC_NAME}${FILENAME_SUFFIX}_${REGRESSOR_IDX}.nii.gz * ${VALS[$i]}"
                fslmaths "$GLM_FOLDER/beta_${METRIC_NAME}${FILENAME_SUFFIX}_${REGRESSOR_IDX}.nii.gz" \
                    -mul "${VALS[$i]}" \
                    -add "$PREDICTED"  \
                    "$PREDICTED"
            done
        else
            echo "  Computing simple mean and standard deviation from healthy subjects ..."
            fslmaths "$MERGED" -Tmean "$PREDICTED"
            fslmaths "$MERGED" -Tstd "$STD_DEV"
        fi

        # Floor the std at a small, data-driven epsilon before dividing. Near the
        # edge of the brain mask, only a handful of registered healthy volumes
        # overlap, so -Tstd can come out near-zero there; dividing by that
        # produces enormous, meaningless z-scores at a tiny fraction of voxels.
        # The floor is 1% of the median std among voxels that do have real
        # (non-zero) variance, so it scales automatically with this metric's
        # own units rather than using one fixed constant across metric types.
        STD_DEV_FLOORED="$GLM_FOLDER/std_dev_floored_${METRIC_NAME}${FILENAME_SUFFIX}.nii.gz"
        MEDIAN_STD=$(fslstats "$STD_DEV" -P 50)
        EPSILON=$(echo "${MEDIAN_STD:-0} * 0.01" | bc -l 2>/dev/null)
        if [[ -z "$EPSILON" ]] || (( $(echo "$EPSILON <= 0" | bc -l) )); then
            EPSILON=0.000001
        fi
        echo "  Flooring std at epsilon=$EPSILON (1% of median non-zero std=$MEDIAN_STD)"
        fslmaths "$STD_DEV" -max "$EPSILON" "$STD_DEV_FLOORED"

        echo "  Computing final Z-score ..."

        TARGET_DIMS="$(fslval "$TARGET_METRIC" dim1) $(fslval "$TARGET_METRIC" dim2) $(fslval "$TARGET_METRIC" dim3)"
        PREDICTED_DIMS="$(fslval "$PREDICTED" dim1) $(fslval "$PREDICTED" dim2) $(fslval "$PREDICTED" dim3)"
        STD_DIMS="$(fslval "$STD_DEV_FLOORED" dim1) $(fslval "$STD_DEV_FLOORED" dim2) $(fslval "$STD_DEV_FLOORED" dim3)"
        if [[ "$TARGET_DIMS" != "$PREDICTED_DIMS" || "$TARGET_DIMS" != "$STD_DIMS" ]]; then
            echo "ERROR: space mismatch before final z-score for $METRIC_NAME" >&2
            echo "  target dims   : $TARGET_DIMS" >&2
            echo "  predicted dims: $PREDICTED_DIMS" >&2
            echo "  stddev dims   : $STD_DIMS" >&2
            exit 1
        fi

        echo "  Computing: (target metric - predicted) / std_dev ..."
        Z_SCORE_RAW="$GLM_FOLDER/z_score_raw_${METRIC_NAME}${FILENAME_SUFFIX}.nii.gz"
        fslmaths "$TARGET_METRIC" \
            -sub "$PREDICTED"     \
            -div "$STD_DEV_FLOORED" \
            "$Z_SCORE_RAW"

        # Clamp z-score as a safety net for any remaining extreme-but-finite
        # values -- the std floor above prevents literal blow-ups, but this keeps
        # the output in a statistically sane range for any downstream consumer
        # that isn't expecting outliers.
        z_score_abs_max=15
        echo "  Clamping z-score to [-$z_score_abs_max, $z_score_abs_max] ..."
        Z_SCORE_UNMASKED="$GLM_FOLDER/z_score_unmasked_${METRIC_NAME}${FILENAME_SUFFIX}.nii.gz"
        fslmaths "$Z_SCORE_RAW" -max "-$z_score_abs_max" -min "$z_score_abs_max" "$Z_SCORE_UNMASKED"

        echo "  Applying brain mask to Z-score ..."
        # Reuses the single mask computed once at the top of this function.
        fslmaths "$Z_SCORE_UNMASKED" -mas "$BRAIN_MASK" "$OUTPUT_FILE"

        echo "  Output saved at: $OUTPUT_FILE"
    done

    if [[ ${#GLM_PEAK_MEM_REPORT[@]} -gt 0 ]]; then
        echo ""
        echo "=== fsl_glm peak memory usage (this run) ==="
        for _line in "${GLM_PEAK_MEM_REPORT[@]}"; do
            echo "  $_line"
        done
    fi

    echo ""
    echo "=== All metrics complete for $TARGET_ID ==="
}

function prepare_target_space_step {
    increment_step
    echo "Step $(get_step): Preparing target/structural spaces"

    local sub_step=0
    # T1 -> skull-stripped target/healthy T1s (robust registration input)
    prepare_t1_inputs
    # Target DWI <-> T1 -> transforms for cross-space projection (if needed)
    prepare_target_dwi_transforms
    # Target metric -> output space (t1/dwi) before final z-score math
    prepare_target_metrics_for_output_space
}

function register_healthy_to_target_step {
    increment_step
    echo "Step $(get_step): Registering healthy subjects to target space"

    local sub_step=0
    # Healthy DWI b0 -> healthy T1 (enables DWI metric chaining)
    register_dwi_to_t1
    # Healthy T1 -> target T1 (single transform reused for all metrics)
    register_t1_to_target_t1
    # QC nodiff -> target T1 (optional diagnostics only)
    qc_register_nodiff_to_t1_folder
}

function build_metric_stacks_step {
    increment_step
    echo "Step $(get_step): Building metric stacks in output space"

    local sub_step=0
    # Healthy metrics -> target T1 (common intermediate analysis space)
    register_metrics
    # Target T1 metrics -> target DWI (only when output-space=dwi)
    project_registered_metrics_to_output_space
    # Per-subject metrics -> 4D merged stack for GLM
    merge_metrics
}

# ===========================================================================
# Pipeline Execution Plan
# Reorder function calls below if you want to change pipeline flow.
# ===========================================================================

# Step 1: GLM setup (design matrix)
generate_design_matrix

# Step 2: Data preparation (healthy DWI inputs - optional)
prepare_dwi_inputs

# Step 3: Prepare target and structural spaces
prepare_target_space_step

# Step 4: Register healthy subjects to target
register_healthy_to_target_step

# Step 5: Build metric stacks in final output space
build_metric_stacks_step

# Step 6: GLM + predicted map + Z-score computation
compute_z_scores


