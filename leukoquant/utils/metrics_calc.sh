#!/bin/bash
#
# Metrics calculation script - runs inside Singularity container.
#
# Accepts all metric extraction parameters and handles conditional logic.
# Computes statistics along tracts, lesions, and other ROIs.
#

set -e

# ============================================================================
# Parse Arguments
# ============================================================================

SUBJECT=""
TRACTOGRAPHY_PATH=""
TRACT_MODE="tractography-atlas"
OUTPUT_DIR=""
LESION_PATH=""
T1_PATH=""
METRICS_SPEC=""
VERBOSE=false
OUTPUT_CSV=""
DWI=""
BVAL=""
SKIP_SKULLSTRIP_T1=false
SKIP_SKULLSTRIP_DWI=false
NEED_DWI_TO_T1_REG=false
DWI_FILE=""
AUTO_FA_SKELETON_PATH=""
AUTO_FA_SKELETON_SPACE=""
SCRATCH_DIR=""
QC_REPORT=""
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

while [[ $# -gt 0 ]]; do
    case $1 in
        --subject)
            SUBJECT="$2"
            shift 2
            ;;
        --tractography-path)
            TRACTOGRAPHY_PATH="$2"
            shift 2
            ;;
        --tract-mode)
            TRACT_MODE="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --lesion-path)
            LESION_PATH="$2"
            shift 2
            ;;
        --t1-path)
            T1_PATH="$2"
            shift 2
            ;;
        --metrics)
            METRICS_SPEC="$2"
            shift 2
            ;;
        --verbose)
            VERBOSE=true
            shift 1
            ;;
        --output-csv)
            OUTPUT_CSV="$2"
            shift 2
            ;;
        --dwi)
            DWI="$2"
            shift 2
            ;;
        --bval)
            BVAL="$2"
            shift 2
            ;;
        --skip-skullstrip-t1)
            SKIP_SKULLSTRIP_T1="$2"
            shift 2
            ;;
        --skip-skullstrip-dwi)
            SKIP_SKULLSTRIP_DWI="$2"
            shift 2
            ;;
        --qc-report)
            QC_REPORT="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Prefer the conda base Python (where dipy/nibabel are installed) over any
# tool-specific python3 that may be prepended to PATH (e.g. FSL's python3).
if [ -x "/opt/conda/bin/python3" ]; then
    PYTHON3="/opt/conda/bin/python3"
else
    PYTHON3="python3"
fi

# Set up scratch directory now that $SUBJECT is known.
# Use /scratch0 which is bind-mounted into the container; fall back to /tmp.
_scratch_base="/scratch0"
if [ ! -d "$_scratch_base" ] || [ ! -w "$_scratch_base" ]; then
    _scratch_base="/tmp"
fi
SCRATCH_DIR="${_scratch_base}/metrics_calc_${SUBJECT}"
mkdir -p "$SCRATCH_DIR"

function cleanup {
    rm -rf "$SCRATCH_DIR"
}

trap cleanup EXIT ERR INT TERM

# ============================================================================
# Validation
# ============================================================================

# ============================================================================
# Helper Functions
# ============================================================================

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
            echo "    $(basename "$file") | dims=${d1}x${d2}x${d3} vols=${d4} | voxel=${px}x${py}x${pz}mm"
            ;;
        *)
            #echo "    $(basename "$file") (non-NIfTI)"
            ;;
    esac
}

init_fsl_tools() {
    # Ensure FSL commands are available for header/space inspection.
    if command -v fslval >/dev/null 2>&1; then
        return 0
    fi

    if [ -n "${FSLDIR:-}" ] && [ -f "${FSLDIR}/etc/fslconf/fsl.sh" ]; then
        # shellcheck disable=SC1090
        source "${FSLDIR}/etc/fslconf/fsl.sh" >/dev/null 2>&1 || true
    fi
}

parse_spec() {
    # Parse base:glob[:space] into 3 vars echoed as: base|glob|space
    # Optional second argument sets default space (default: dwi).
    local spec="$1"
    local default_space="${2:-dwi}"
    local base=""
    local glob=""
    local space="$default_space"

    base=$(echo "$spec" | cut -d: -f1)
    glob=$(echo "$spec" | cut -d: -f2)

    # 3rd component is optional space token
    local third
    third=$(echo "$spec" | cut -d: -f3)
    if [ -n "$third" ]; then
        space="$third"
    fi

    echo "${base}|${glob}|${space}"
}

parse_lesion_entry() {
    # Parse [name=]base:glob[:space] -> name|base|glob|space
    # Name is separated from the path spec by '=' (optional; default: 'lesion').
    local entry="$1"
    local default_space="${2:-t1}"
    local name="lesion"
    local spec="$entry"

    if [[ "$entry" == *'='* ]]; then
        name="${entry%%=*}"
        spec="${entry#*=}"
    fi

    local parsed
    parsed=$(parse_spec "$spec" "$default_space")
    local base glob space
    base=$(echo "$parsed" | cut -d'|' -f1)
    glob=$(echo "$parsed" | cut -d'|' -f2)
    space=$(echo "$parsed" | cut -d'|' -f3)

    echo "${name}|${base}|${glob}|${space}"
}

if [ -z "$SUBJECT" ] || [ -z "$TRACTOGRAPHY_PATH" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "ERROR: Missing required parameters"
    echo "Usage: $0 --subject SUBJ --tractography-path PATH --output-dir DIR [--tract-mode MODE] [--lesion-path PATH] [--metrics JSON] [--verbose bool]"
    exit 1
fi

init_fsl_tools

if [ "$VERBOSE" = "True" ] || [ "$VERBOSE" = "true" ]; then
    echo "=========================================="
    echo "Metrics Extraction for $SUBJECT"
    echo "=========================================="
    echo "Tractography: $TRACTOGRAPHY_PATH"
    echo "Tract mode: $TRACT_MODE"
    echo "Output: $OUTPUT_DIR"
    if [ -n "$LESION_PATH" ]; then
        echo "Lesion path: $LESION_PATH"
    fi
    if [ -n "$T1_PATH" ]; then
        echo "T1 path: $T1_PATH"
    fi
    if [ -n "$DWI" ]; then
        echo "DWI: $DWI"
    fi
    if [ -n "$METRICS_SPEC" ]; then
        echo "Metrics: $METRICS_SPEC"
    fi
fi

# ============================================================================
# Main Processing
# ============================================================================

mkdir -p "$OUTPUT_DIR"

T1_FILE=""
INPUT_MANIFEST="$SCRATCH_DIR/t1_space_inputs_manifest.tsv"
: > "$INPUT_MANIFEST"

# Extract tractography base directory (before the colon, if present)
TRACTOGRAPHY_PARSED=$(parse_spec "$TRACTOGRAPHY_PATH" "dwi")
TRACTOGRAPHY_BASE=$(echo "$TRACTOGRAPHY_PARSED" | cut -d'|' -f1)
TRACTOGRAPHY_GLOB=$(echo "$TRACTOGRAPHY_PARSED" | cut -d'|' -f2)
TRACTOGRAPHY_SPACE=$(echo "$TRACTOGRAPHY_PARSED" | cut -d'|' -f3)

# Validate input directories exist
if [ ! -d "$TRACTOGRAPHY_BASE" ]; then
    echo "ERROR: Tractography base path not found: $TRACTOGRAPHY_BASE"
    exit 1
fi

echo ""
# Check tractography files using the glob pattern
TRACT_FILES=$(find "$TRACTOGRAPHY_BASE" -path "*$TRACTOGRAPHY_GLOB" 2>/dev/null | sort)
if [ -z "$TRACT_FILES" ]; then
    echo "WARNING: No tractography files found matching $TRACTOGRAPHY_GLOB in $TRACTOGRAPHY_BASE"
else
    TRACT_FILE_COUNT=$(echo "$TRACT_FILES" | wc -l | tr -d ' ')
    echo "✓ Found $TRACT_FILE_COUNT tractography file(s) [space=${TRACTOGRAPHY_SPACE}]"
    while IFS= read -r _f; do
        #print_image_info "tract" "$_f"
        _tract_name=$(basename "$(dirname "$_f")")
        echo -e "tractography\t${_tract_name}\t${TRACTOGRAPHY_SPACE}\t${_f}" >> "$INPUT_MANIFEST"
    done <<< "$TRACT_FILES"
fi

# Process lesion(s) if provided  (format: [name:]base:glob[:space] comma-separated)
if [ -n "$LESION_PATH" ] && [ "$LESION_PATH" != "None" ]; then
    IFS=',' read -ra _LESION_ENTRIES <<< "$LESION_PATH"
    for _LESION_ENTRY in "${_LESION_ENTRIES[@]}"; do
        _LESION_ENTRY=$(echo "$_LESION_ENTRY" | xargs)  # trim whitespace
        [ -z "$_LESION_ENTRY" ] && continue

        _LESION_EPARSED=$(parse_lesion_entry "$_LESION_ENTRY" "t1")
        _LESION_NAME=$(echo  "$_LESION_EPARSED" | cut -d'|' -f1)
        _LESION_BASE=$(echo  "$_LESION_EPARSED" | cut -d'|' -f2)
        _LESION_GLOB=$(echo  "$_LESION_EPARSED" | cut -d'|' -f3)
        _LESION_SPACE=$(echo "$_LESION_EPARSED" | cut -d'|' -f4)

        if [ "$_LESION_SPACE" = "dwi" ]; then
            NEED_DWI_TO_T1_REG=true
        fi
        if [ ! -d "$_LESION_BASE" ]; then
            echo "ERROR: Lesion base path not found: $_LESION_BASE"
            exit 1
        fi
        _LESION_FILES=$(find "$_LESION_BASE" -path "*$_LESION_GLOB" 2>/dev/null | sort)
        if [ -z "$_LESION_FILES" ]; then
            echo "WARNING: No lesion files found matching $_LESION_GLOB in $_LESION_BASE"
        else
            _LESION_COUNT=$(echo "$_LESION_FILES" | wc -l | tr -d ' ')
            echo "✓ Found $_LESION_COUNT lesion file(s) [name=${_LESION_NAME}, space=${_LESION_SPACE}]"
            while IFS= read -r _f; do
                echo -e "lesion\t${_LESION_NAME}\t${_LESION_SPACE}\t${_f}" >> "$INPUT_MANIFEST"
            done <<< "$_LESION_FILES"
        fi
    done
fi

# Process T1 if provided
if [ -n "$T1_PATH" ] && [ "$T1_PATH" != "None" ]; then
    T1_PARSED=$(parse_spec "$T1_PATH" "t1")
    T1_BASE=$(echo "$T1_PARSED" | cut -d'|' -f1)
    T1_GLOB=$(echo "$T1_PARSED" | cut -d'|' -f2)
    T1_SPACE=$(echo "$T1_PARSED" | cut -d'|' -f3)

    if [ -f "$T1_BASE" ]; then
        # T1 is a direct file path - use it as-is
        T1_FILE="$T1_BASE"
        echo "✓ Found T1 file (direct) [space=${T1_SPACE}]: $T1_FILE"
        echo -e "t1\tmain\t${T1_SPACE}\t${T1_FILE}" >> "$INPUT_MANIFEST"
    elif [ -d "$T1_BASE" ]; then
        # T1 is a directory - search with glob
        T1_FILES_LIST=$(find "$T1_BASE" -path "*$T1_GLOB" 2>/dev/null | sort)
        if [ -z "$T1_FILES_LIST" ]; then
            echo "WARNING: No T1 files found matching $T1_GLOB in $T1_BASE"
        else
            T1_COUNT=$(echo "$T1_FILES_LIST" | wc -l | tr -d ' ')
            T1_FILE=$(echo "$T1_FILES_LIST" | head -n 1)
            echo "✓ Found $T1_COUNT T1 file(s) [space=${T1_SPACE}]"
            while IFS= read -r _f; do
                echo -e "t1\tmain\t${T1_SPACE}\t${_f}" >> "$INPUT_MANIFEST"
            done <<< "$T1_FILES_LIST"
        fi
    else
        echo "ERROR: T1 base path not found: $T1_BASE"
        exit 1
    fi
fi

# Process metrics if provided
if [ -n "$METRICS_SPEC" ] && [ "$METRICS_SPEC" != "None" ] && [ "$METRICS_SPEC" != "{}" ]; then
    echo ""
    echo "Processing metric maps:"
    
    # Parse metrics spec (format: metric1=base:glob:space,metric2=base:glob:space,...)
    # Split by comma to handle multiple metrics
    IFS=',' read -ra METRICS_ARRAY <<< "$METRICS_SPEC"
    
    for METRIC_ENTRY in "${METRICS_ARRAY[@]}"; do
        # Parse metric_name=base:glob:space
        METRIC_NAME=$(echo "$METRIC_ENTRY" | cut -d'=' -f1 | xargs)
        METRIC_SPEC=$(echo "$METRIC_ENTRY" | cut -d'=' -f2- | xargs)
        
        if [ -z "$METRIC_NAME" ] || [ -z "$METRIC_SPEC" ]; then
            continue
        fi
        
        # Extract base (first component before first colon)
        METRIC_PARSED=$(parse_spec "$METRIC_SPEC" "dwi")
        METRIC_BASE=$(echo "$METRIC_PARSED" | cut -d'|' -f1)
        METRIC_GLOB=$(echo "$METRIC_PARSED" | cut -d'|' -f2)
        METRIC_SPACE=$(echo "$METRIC_PARSED" | cut -d'|' -f3)

        if [ "$METRIC_SPACE" = "dwi" ]; then
            NEED_DWI_TO_T1_REG=true
        fi
        
        if [ -z "$METRIC_GLOB" ]; then
            echo "  WARNING: Metric '$METRIC_NAME' has no glob pattern in spec: $METRIC_SPEC"
            continue
        fi
        
        if [ ! -d "$METRIC_BASE" ]; then
            echo "  WARNING: Metric '$METRIC_NAME' base path not found: $METRIC_BASE"
            continue
        fi
        
        # Check for metric files matching pattern
        METRIC_FILES_LIST=$(find "$METRIC_BASE" -path "*$METRIC_GLOB" 2>/dev/null | sort)
        
        if [ -n "$METRIC_FILES_LIST" ]; then
            METRIC_COUNT=$(echo "$METRIC_FILES_LIST" | wc -l | tr -d ' ')
            echo "  ✓ Metric '$METRIC_NAME' [space=${METRIC_SPACE}]: found $METRIC_COUNT file(s)"

            # Auto-detect FA metric by name and keep first matching file for skeleton generation.
            if [ -z "$AUTO_FA_SKELETON_PATH" ]; then
                METRIC_NAME_LOWER=$(echo "$METRIC_NAME" | tr '[:upper:]' '[:lower:]')
                if echo "$METRIC_NAME_LOWER" | grep -q "fa"; then
                    AUTO_FA_SKELETON_PATH=$(echo "$METRIC_FILES_LIST" | head -n 1)
                    AUTO_FA_SKELETON_SPACE="$METRIC_SPACE"
                fi
            fi

            while IFS= read -r _f; do
                #print_image_info "$METRIC_NAME" "$_f"
                echo -e "metric\t${METRIC_NAME}\t${METRIC_SPACE}\t${_f}" >> "$INPUT_MANIFEST"
            done <<< "$METRIC_FILES_LIST"
        fi
    done
fi


# ============================================================================
# Build FA Skeleton 
# ============================================================================

if [ -n "$AUTO_FA_SKELETON_PATH" ] && [ -f "$AUTO_FA_SKELETON_PATH" ]; then
    echo ""
    echo "--- Building FA Skeleton ---"
    SCRIPT_DIR_LOCAL="$(cd "$(dirname "$0")" && pwd)"
    SKELETON_SCRIPT="$SCRIPT_DIR_LOCAL/skeleton.sh"
    SKELETON_OUT_DIR="$SCRATCH_DIR/skeleton"

    if [ ! -f "$SKELETON_SCRIPT" ]; then
        echo "WARNING: skeleton.sh not found at $SKELETON_SCRIPT, skipping skeleton generation."
    else
        SKEL_ARGS=(
            --fa-path "$AUTO_FA_SKELETON_PATH"
            --output-dir "$SKELETON_OUT_DIR"
        )

        skeleton_cmd="bash $SKELETON_SCRIPT ${SKEL_ARGS[*]}"
        #echo "Running FA skeleton script:"
        #echo "  $skeleton_cmd"
        eval "$skeleton_cmd" > /dev/null 2>&1
        
        echo "✓ FA skeleton generated."
    fi
    # Add skeleton to manifest. Separate folder: skeleton, and space the same as the FA metric used to generate it
    skeleton_file="$SKELETON_OUT_DIR/skeleton_mask_direct.nii.gz"
    skeleton_space="$AUTO_FA_SKELETON_SPACE"
    echo -e "skeleton\tskeleton\t${skeleton_space}\t${skeleton_file}" >> "$INPUT_MANIFEST"
fi

# ============================================================================
# Register FA to MGH35_HCP_FA_template to use tract atlas
# ============================================================================
if [ -n "$AUTO_FA_SKELETON_PATH" ] && [ -f "$AUTO_FA_SKELETON_PATH" ]; then
    echo ""
    echo "--- Preparing Tract Atlas ---"

    MGH35_HCP_FA_TEMPLATE=""
    if [ -d "$FREESURFER_HOME" ]; then
        _template_candidate="$FREESURFER_HOME/trctrain/hcp/MGH35_HCP_FA_template.nii.gz"
        if [ -f "$_template_candidate" ]; then
            MGH35_HCP_FA_TEMPLATE="$_template_candidate"
        else
            echo "WARNING: MGH35_HCP_FA_template not found at $FREESURFER_HOME/trctrain/hcp/. Skipping atlas tract processing."
        fi
    else
        echo "WARNING: FREESURFER_HOME not set or invalid. Skipping atlas tract processing."
    fi

    if [ -n "$MGH35_HCP_FA_TEMPLATE" ]; then
        FA_REGISTERED_DIR="$SCRATCH_DIR/MGH35_HCP_FA_SUBJECT_SPACE"
        mkdir -p "$FA_REGISTERED_DIR"

        FA_FILE="$AUTO_FA_SKELETON_PATH"
        MGH2SUBJ_AFF="$FA_REGISTERED_DIR/mgh2subj_affine.txt"
        MGH2SUBJ_CPP="$FA_REGISTERED_DIR/mgh2subj_cpp.nii.gz"
        RES_FA="$FA_REGISTERED_DIR/fa_in_template.nii.gz"
        FA_TO_MGH35_DEF="$FA_REGISTERED_DIR/fa_to_mgh35_deformation.nii.gz"

        # reg_f3d's gradient-descent optimisation runs multi-threaded (-omp), and
        # floating-point reductions across threads are not bit-reproducible run to
        # run; occasionally that variance is enough to tip the nonlinear optimiser
        # into a bad local optimum even on data it usually registers correctly.
        # Confirmed on a subject whose cluster run produced a catastrophically
        # mislocated atlas tract (128mm centroid error): re-running this exact
        # registration offline against the same FA file and template succeeded
        # cleanly (NCC 0.88 vs. the failed run's near-zero/negative NCC). So a
        # failed attempt is retried a bounded number of times before giving up,
        # rather than accepted or abandoned after a single try.
        MAX_ATLAS_REG_ATTEMPTS=4  # 1 initial attempt + 3 retries
        ATLAS_REG_STATUS="FAIL"
        for _atlas_reg_attempt in $(seq 1 $MAX_ATLAS_REG_ATTEMPTS); do
            echo "  FA→MGH35 registration attempt $_atlas_reg_attempt/$MAX_ATLAS_REG_ATTEMPTS"

            # Affine register Subject FA to MGH35 template.
            echo "  Affine registration (reg_aladin): FA → MGH35 template"
            reg_aladin \
            -ref "$FA_FILE"    \
            -flo "$MGH35_HCP_FA_TEMPLATE"    \
            -aff "$MGH2SUBJ_AFF" \
            -res "$FA_REGISTERED_DIR/fa_aladin_result.nii.gz" \
            -omp "$THREADS"     \
            -voff > /dev/null

            # Non-linear CPP refinement.
            echo "  Non-linear registration (reg_f3d): FA → MGH35 template"
            reg_f3d \
            -ref "$FA_FILE"  \
            -flo "$MGH35_HCP_FA_TEMPLATE"  \
            -aff "$MGH2SUBJ_AFF"     \
            -cpp "$MGH2SUBJ_CPP"        \
            -res "$RES_FA"     \
            -omp "$THREADS"    \
            -voff > /dev/null

            # Convert CPP to a dense forward deformation field on the FA grid.
            # Each FA voxel stores the corresponding MGH35 world-mm coordinate
            # (FA_world → MGH35_world).  Used by transform_trk.py to warp atlas
            # streamlines to subject space via fixed-point inversion.
            echo "  Generating forward deformation field (reg_transform)"
            # NiftyReg v2.0.0 syntax: -def <input_transform> <output_field>
            # The CPP file is auto-detected as a recognised transformation type.
            reg_transform \
            -ref "$FA_FILE" \
            -def "$MGH2SUBJ_CPP" "$FA_TO_MGH35_DEF" > /dev/null

            # Validate registration accuracy before trusting it for the atlas fallback.
            # See check_atlas_registration.py's docstring for why this check exists
            # (transform_trk.py's --verify-nifti alone cannot catch a genuinely bad
            # reg_aladin/reg_f3d registration).
            ATLAS_REG_QC=$($PYTHON3 "$SCRIPT_DIR/check_atlas_registration.py" "$FA_FILE" "$RES_FA" --threshold 0.4 2>/dev/null)
            ATLAS_REG_STATUS="${ATLAS_REG_QC%% *}"
            ATLAS_REG_NCC="${ATLAS_REG_QC#* }"

            if [ "$ATLAS_REG_STATUS" = "PASS" ]; then
                echo "  ✓ FA→MGH35 registration NCC=$ATLAS_REG_NCC (passed >= 0.4 threshold, attempt $_atlas_reg_attempt)"
                break
            else
                echo "  ✗ FA→MGH35 registration NCC=${ATLAS_REG_NCC:-unavailable} below threshold (attempt $_atlas_reg_attempt/$MAX_ATLAS_REG_ATTEMPTS)"
            fi
        done

        if [ "$ATLAS_REG_STATUS" != "PASS" ]; then
            echo "WARNING: FA→MGH35 template registration quality too low after $MAX_ATLAS_REG_ATTEMPTS attempts (NCC=${ATLAS_REG_NCC:-unavailable})."
            echo "         Skipping tract atlas generation for this subject/session -- the atlas fallback"
            echo "         would silently place tracts at the wrong anatomical location otherwise."
            MGH35_HCP_FA_TEMPLATE=""  # gates the atlas-tract loop below off via the existing guard
        fi

    if [ -n "$MGH35_HCP_FA_TEMPLATE" ]; then
    for tract_atlas_temp in "$FREESURFER_HOME"/trctrain/hcp/syn/*display.trk; do

        output_tracts_atlas_dir="$SCRATCH_DIR/tracts_atlas"
        mkdir -p "$output_tracts_atlas_dir"

        tract_name_temp=$(basename "$tract_atlas_temp" .display.trk)

        track_MGH35_space="$FA_REGISTERED_DIR/${tract_name_temp}_MGH35_space.trk"
        density_map_MGH35_space="$FA_REGISTERED_DIR/${tract_name_temp}_density_MGH35.nii.gz"
        pd_trk_dest="$output_tracts_atlas_dir/${tract_name_temp}.nii.gz"
        trk_subj_dest="$output_tracts_atlas_dir/${tract_name_temp}_in_subject.trk"
        density_map_subj_from_trk="$FA_REGISTERED_DIR/${tract_name_temp}_density_subject_from_trk.nii.gz"

        convert_trk_script="$SCRIPT_DIR/transform_trk.py"
        resample_trk_to_MGH35_space_cmd="$PYTHON3 $convert_trk_script \
        --input-trk '$tract_atlas_temp' \
        --output-trk '$track_MGH35_space' \
        --in-reference '$MGH35_HCP_FA_TEMPLATE' \
        --convert-affine
        "
        eval "$resample_trk_to_MGH35_space_cmd" > /dev/null

        density_map_cmd="$PYTHON3 $convert_trk_script \
        --input-trk '$track_MGH35_space' \
        --output-vol '$density_map_MGH35_space' \
        "
        eval "$density_map_cmd" > /dev/null

        # Resample density map to subject space using inverted CPP
        reg_resample \
        -ref "$FA_FILE" \
        -flo "$density_map_MGH35_space" \
        -trans "$MGH2SUBJ_CPP" \
        -res "$pd_trk_dest" \
        -omp "$THREADS" \
        -inter 1 \
        -voff > /dev/null

        # Make sure output values are larger than 0
        fslmaths "$pd_trk_dest" -thr 0 "$pd_trk_dest" > /dev/null

        # Warp atlas trk from MGH35 space to subject FA space using the full
        # non-linear deformation field (affine + CPP).  The affine is passed
        # separately as the initialisation for the fixed-point inversion.
        # --verify-nifti compares the resulting trk CoM against the already-
        # registered atlas density map so we can catch a wrong warp direction.
        # --verify-output-vol saves the density map generated from the warped TRK.
        transform_atlas_trk_cmd="$PYTHON3 $convert_trk_script \
        --input-trk '$track_MGH35_space' \
        --output-trk '$trk_subj_dest' \
        --warp '$FA_TO_MGH35_DEF' \
        --affine '$MGH2SUBJ_AFF' \
        --out-reference '$FA_FILE' \
        --verify-nifti '$pd_trk_dest' \
        --verify-output-vol '$density_map_subj_from_trk'
        "
        eval "$transform_atlas_trk_cmd" 2>&1

        tract_atlas_space="$AUTO_FA_SKELETON_SPACE"
        echo -e "tract_atlas\t${tract_name_temp}\t${tract_atlas_space}\t${pd_trk_dest}" >> "$INPUT_MANIFEST"
    done
        echo "✓ Tract atlas registration to subject space completed."
    fi  # [ -n "$MGH35_HCP_FA_TEMPLATE" ] (post-registration-QC gate)
    fi  # [ -n "$MGH35_HCP_FA_TEMPLATE" ]
fi


# ============================================================================
# DWI -> T1 Registration (if DWI and T1 are both available)
# ============================================================================

if [ "$NEED_DWI_TO_T1_REG" = "true" ] && [ -n "$DWI" ] && [ "$DWI" != "None" ] && [ -n "$T1_FILE" ]; then
    echo ""
    echo "--- DWI -> T1 Registration ---"
    REG_DIR="$SCRATCH_DIR/reg_dwi_to_t1"
    REGISTER_SCRIPT="$SCRIPT_DIR/register_dwi_to_t1.sh"
    if [ ! -f "$REGISTER_SCRIPT" ]; then
        echo "WARNING: register_dwi_to_t1.sh not found at $REGISTER_SCRIPT, skipping registration."
    else
        # Resolve --dwi value to an actual file path if provided as base:glob[:space]
        if echo "$DWI" | grep -q ":"; then
            DWI_PARSED=$(parse_spec "$DWI" "dwi")
            DWI_BASE=$(echo "$DWI_PARSED" | cut -d'|' -f1)
            DWI_GLOB=$(echo "$DWI_PARSED" | cut -d'|' -f2)
            DWI_FILES=$(find "$DWI_BASE" -path "*$DWI_GLOB" 2>/dev/null | sort)
            DWI_FILE=$(echo "$DWI_FILES" | head -n 1)
        else
            DWI_FILE="$DWI"
        fi

        if [ -z "$DWI_FILE" ] || [ ! -f "$DWI_FILE" ]; then
            echo "WARNING: DWI-space inputs detected, but resolved DWI file is missing. Skipping DWI->T1 registration."
            echo "  DWI input: $DWI"
            DWI_FILE=""
        fi

        # Run registration at most once for this subject/output directory.
        if [ -n "$DWI_FILE" ] && [ -f "$REG_DIR/diff2t1_affine.txt" ] && [ -f "$REG_DIR/t12diff_affine.txt" ]; then
            echo "[SKIP] DWI -> T1 registration already exists in $REG_DIR"
            #echo "DIFF2T1: $REG_DIR/diff2t1_affine.txt"
            #echo "T12DIFF: $REG_DIR/t12diff_affine.txt"
        elif [ -n "$DWI_FILE" ]; then
        REG_ARGS=(
            --dwi "$DWI_FILE"
            --t1 "$T1_FILE"
            --output-dir "$REG_DIR"
            --skip-skullstrip-t1 "$SKIP_SKULLSTRIP_T1"
            --skip-skullstrip-dwi "$SKIP_SKULLSTRIP_DWI"
            --verbose "$VERBOSE"
        )
        if [ -n "$BVAL" ] && [ "$BVAL" != "None" ]; then
            REG_ARGS+=(--bval "$BVAL")
        fi
        bash "$REGISTER_SCRIPT" "${REG_ARGS[@]}"
        fi

        echo "✓ DWI -> T1 registration completed."
    fi
elif [ "$NEED_DWI_TO_T1_REG" = "true" ]; then
    MISSING=()
    if [ -z "$DWI" ] || [ "$DWI" = "None" ]; then
        MISSING+=("--dwi")
    fi
    if [ -z "$T1_FILE" ]; then
        MISSING+=("subject-specific T1 (from --t1-path)")
    fi
    echo "WARNING: DWI-space inputs detected, missing ${MISSING[*]}. Skipping DWI->T1 registration."
fi

# ============================================================================
# Create T1-space inputs folder for downstream metrics use
# ============================================================================

T1_SPACE_DIR="$SCRATCH_DIR/t1_space_inputs"
REG_DIR="$SCRATCH_DIR/reg_dwi_to_t1"
DIFF2T1_AFF="$REG_DIR/diff2t1_affine.txt"
B0_BRAIN="$REG_DIR/b0_brain.nii.gz"

mkdir -p "$T1_SPACE_DIR/tractography" "$T1_SPACE_DIR/lesion" "$T1_SPACE_DIR/t1" "$T1_SPACE_DIR/metrics" "$T1_SPACE_DIR/skeleton" "$T1_SPACE_DIR/tract_atlas"


if [ -s "$INPUT_MANIFEST" ]; then
    echo ""
    echo "--- Building T1-space input folder ---"

    while IFS=$'\t' read -r category subcat space src; do
        [ -z "$src" ] && continue

        case "$category" in
            tractography)
                # Preserve tract identity by creating one folder per tract label.
                dest_dir="$T1_SPACE_DIR/tractography/$subcat"
                ;;
            lesion)
                # Each lesion type goes into its own named subdirectory.
                dest_dir="$T1_SPACE_DIR/lesion/$subcat"
                ;;
            t1)
                dest_dir="$T1_SPACE_DIR/t1"
                ;;
            metric)
                dest_dir="$T1_SPACE_DIR/metrics/$subcat"
                ;;
            skeleton)
                dest_dir="$T1_SPACE_DIR/skeleton"
                ;;
            tract_atlas)
                dest_dir="$T1_SPACE_DIR/tract_atlas"
                ;;
            *)
                continue
                ;;
        esac

        mkdir -p "$dest_dir"
        base_name=$(basename "$src")
        dest="$dest_dir/$base_name"

        if [ "$space" = "dwi" ]; then
            case "$src" in
                *.nii|*.nii.gz)
                    ;;
                *)
                    echo "WARNING: Skipping non-NIfTI dwi-space input: $src"
                    continue
                    ;;
            esac

            if [ -z "$T1_FILE" ] || [ ! -f "$T1_FILE" ]; then
                echo "WARNING: Cannot project dwi-space input (missing T1 reference): $src"
                continue
            fi
            if [ ! -f "$DIFF2T1_AFF" ]; then
                echo "WARNING: Cannot project dwi-space input (missing transform $DIFF2T1_AFF): $src"
                continue
            fi

            # Skip resampling path.pd.nii.gz directly; we only create tracts.nii.gz downstream
            if [ "$category" != "tractography" ]; then
                if [ ! -f "$dest" ]; then
                    reg_resample \
                        -ref "$T1_FILE" \
                        -flo "$src" \
                        -trans "$DIFF2T1_AFF" \
                        -res "$dest" \
                        -inter 1 \
                        -voff > /dev/null 2>&1
                fi
            fi
        else
            if [ ! -f "$dest" ]; then
                cp "$src" "$dest"
            fi
        fi

        # Copy the warped atlas TRK file if it exists, renaming it for metrics_calc.py
        if [ "$category" = "tract_atlas" ]; then
            src_parent=$(dirname "$src")
            _atlas_trk_src="$src_parent/${subcat}_in_subject.trk"
            _atlas_trk_dest="$dest_dir/${subcat}_in_subject.trk"
            if [ -f "$_atlas_trk_src" ] && [ ! -f "$_atlas_trk_dest" ]; then
                cp "$_atlas_trk_src" "$_atlas_trk_dest"
            fi
        fi

        # Keep only selected tract sidecar files available for downstream metrics.
        if [ "$category" = "tractography" ]; then
            src_parent=$(dirname "$src")
            if [ -d "$src_parent" ]; then
                TRACT_SIDECAR_FILES=(
                    "length.samples.txt"
                )
                for _sidecar_name in "${TRACT_SIDECAR_FILES[@]}"; do
                    _aux="$src_parent/$_sidecar_name"
                    aux_dest="$dest_dir/$_sidecar_name"
                    if [ -f "$_aux" ] && [ ! -f "$aux_dest" ]; then
                        cp "$_aux" "$aux_dest"
                    fi
                done
                if [ -f "$B0_BRAIN" ]; then
                    # Transform "path.pd.trk" files to T1 space for use in metrics calc.
                    _pd_trk="$src_parent/path.pd.trk"
                    if [ -f "$_pd_trk" ]; then
                        pd_trk_dest="$dest_dir/tracts.trk"
                        pd_nii_dest="$dest_dir/tracts.nii.gz"
                        pd_nii_dwi="$SCRATCH_DIR/tractography_dwi_${subcat}.nii.gz"
                        pd_nii_t1_from_trk="$SCRATCH_DIR/tractography_from_trk_${subcat}.nii.gz"
                        convert_trk_script="$SCRIPT_DIR/transform_trk.py"

                        # 1. Generate the B0/DWI space density map in scratch
                        generate_dwi_density_cmd="$PYTHON3 $convert_trk_script \
                        --input-trk '$_pd_trk' \
                        --output-vol '$pd_nii_dwi' \
                        "
                        eval "$generate_dwi_density_cmd" > /dev/null

                        # 2. Resample the B0 density map to T1 space (used directly for metrics)
                        reg_resample \
                        -ref "$T1_FILE" \
                        -flo "$pd_nii_dwi" \
                        -trans "$DIFF2T1_AFF" \
                        -res "$pd_nii_dest" \
                        -inter 1 \
                        -voff > /dev/null

                        # Resample original path.pd.nii.gz to T1 space for comparison
                        pd_nii_t1_tracula_registered="$SCRATCH_DIR/tractography_tracula_pd_t1_${subcat}.nii.gz"
                        if [ -f "$src_parent/path.pd.nii.gz" ]; then
                            reg_resample \
                            -ref "$T1_FILE" \
                            -flo "$src_parent/path.pd.nii.gz" \
                            -trans "$DIFF2T1_AFF" \
                            -res "$pd_nii_t1_tracula_registered" \
                            -inter 1 \
                            -voff > /dev/null
                        fi

                        # 3. Transform TRK and generate density map directly from transformed TRK
                        transform_trk_cmd="$PYTHON3 $convert_trk_script \
                        --input-trk '$_pd_trk' \
                        --output-trk '$pd_trk_dest' \
                        --in-reference '$B0_BRAIN' \
                        --out-reference '$T1_FILE' \
                        --affine '$DIFF2T1_AFF' \
                        --verify-nifti '$pd_nii_dest' \
                        --verify-output-vol '$pd_nii_t1_from_trk' \
                        "

                        #echo "  Transforming tractography sidecar to T1 space: $pd_trk_dest"
                        #echo "    Command: $transform_trk_cmd"
                        eval "$transform_trk_cmd" 2>&1

                    fi
                fi

            fi
        fi
    done < "$INPUT_MANIFEST"

    echo "✓ T1-space inputs prepared."
fi



# Skull-strip T1 for TIV computation (unless explicitly skipped).
T1_BRAIN_FOR_METRICS=""
if [ "$SKIP_SKULLSTRIP_T1" != "true" ] && [ "$SKIP_SKULLSTRIP_T1" != "True" ]; then
    T1_IN_SPACE=$(find "$T1_SPACE_DIR/t1" \( -name "*.nii.gz" -o -name "*.nii" \) 2>/dev/null | sort | head -n 1)
    if [ -n "$T1_IN_SPACE" ]; then
        T1_BRAIN_FOR_METRICS="$T1_SPACE_DIR/t1/t1_brain.nii.gz"
        if [ ! -f "$T1_BRAIN_FOR_METRICS" ]; then
            #echo "  [metrics] Skull stripping T1 for TIV computation ..."
            mkdir -p "$T1_SPACE_DIR/t1"
            if command -v mri_synthstrip >/dev/null 2>&1; then
                mri_synthstrip -i "$T1_IN_SPACE" -o "$T1_BRAIN_FOR_METRICS" >/dev/null 2>&1
            else
                echo "WARNING: mri_synthstrip not found; TIV will be skipped."
                T1_BRAIN_FOR_METRICS=""
            fi
        fi
    fi
fi


# Snapshot exact inputs to metrics_calc.py for debugging and standalone re-runs.
DEBUG_METRICS_DIR="$(dirname "$OUTPUT_DIR")/debug/metrics_inputs"
rm -rf "$DEBUG_METRICS_DIR"
mkdir -p "$DEBUG_METRICS_DIR"
cp -r "$T1_SPACE_DIR/." "$DEBUG_METRICS_DIR/"
if [ -n "$QC_REPORT" ] && [ -f "$QC_REPORT" ]; then
    cp "$QC_REPORT" "$DEBUG_METRICS_DIR/qc_report.json"
fi
echo "✓ Debug inputs saved to debug/metrics_inputs"

echo ""
echo "--- Computing Metrics ---"
# Prepare metrics CSV via metrics_calc.py.
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PY_METRICS_PREP="$SCRIPT_DIR/metrics_calc.py"
if [ -f "$PY_METRICS_PREP" ]; then
    MAPS_PREP_SPEC=""
    for map_dir in "$T1_SPACE_DIR"/metrics/*; do
        [ -d "$map_dir" ] || continue
        map_name="$(basename "$map_dir")"
        entry="${map_name}=${map_dir}"
        if [ -z "$MAPS_PREP_SPEC" ]; then
            MAPS_PREP_SPEC="$entry"
        else
            MAPS_PREP_SPEC="${MAPS_PREP_SPEC},${entry}"
        fi
    done

    PREP_ARGS=(
        --subject "$SUBJECT"
        --t1-path "$T1_SPACE_DIR/t1"
        --output "$OUTPUT_DIR"
        --tract-mode "$TRACT_MODE"
    )

    if [ -n "$T1_BRAIN_FOR_METRICS" ] && [ -f "$T1_BRAIN_FOR_METRICS" ]; then
        PREP_ARGS+=(--t1-brain-path "$T1_BRAIN_FOR_METRICS")
    else
        PREP_ARGS+=(--skip-skullstrip)
    fi

    # Build tracts spec: directory path
    if [ -d "$T1_SPACE_DIR/tractography" ]; then
        PREP_ARGS+=(--tracts-path "$T1_SPACE_DIR/tractography")
    fi

    # Build named lesion spec: name1=dir1,name2=dir2
    _LESION_SPEC=""
    for _ldir in "$T1_SPACE_DIR"/lesion/*/; do
        [ -d "$_ldir" ] || continue
        _lname=$(basename "$_ldir")
        _entry="${_lname}=${_ldir%/}"
        if [ -z "$_LESION_SPEC" ]; then
            _LESION_SPEC="$_entry"
        else
            _LESION_SPEC="${_LESION_SPEC},${_entry}"
        fi
    done
    if [ -n "$_LESION_SPEC" ]; then
        PREP_ARGS+=(--lesion-path "$_LESION_SPEC")
    fi
    if [ -n "$MAPS_PREP_SPEC" ]; then
        PREP_ARGS+=(--maps "$MAPS_PREP_SPEC")
    fi
    # Add skeleton if it exists
    if [ -d "$T1_SPACE_DIR/skeleton" ]; then
        PREP_ARGS+=(--skeleton-path "$T1_SPACE_DIR/skeleton")
    fi

    # Add tract atlas if it exists
    if [ -d "$T1_SPACE_DIR/tract_atlas" ]; then
        PREP_ARGS+=(--tract-atlas-path "$T1_SPACE_DIR/tract_atlas")
    fi

    # Add verbose flag if set
    if [ "$VERBOSE" = "True" ] || [ "$VERBOSE" = "true" ]; then
        PREP_ARGS+=(--verbose)
    fi

    # Pass pre-computed tract QC report when available so metrics_calc.py
    # does not need to re-run QC inline for every tract.
    if [ -n "$QC_REPORT" ] && [ -f "$QC_REPORT" ]; then
        PREP_ARGS+=(--qc-report "$QC_REPORT")
    fi

    metrics_cmd="$PYTHON3 $PY_METRICS_PREP ${PREP_ARGS[*]}"
    #echo "Running metrics preparation script:"
    #echo "  $metrics_cmd"
    eval "$metrics_cmd"

    echo "✓ Metrics CSV prepared."
fi

# Delete the DWI -> T1 registration directory and the input manifest.

# if [ -d "$REG_DIR" ]; then
#     rm -rf "$REG_DIR"
# fi

# if [ -f "$INPUT_MANIFEST" ]; then
#     rm -f "$INPUT_MANIFEST"
# fi

# Remove all image files
# if [ -d "$T1_SPACE_DIR" ]; then
#     rm -rf "$T1_SPACE_DIR"
# fi

echo ""
echo "Metrics extraction completed successfully for $SUBJECT"