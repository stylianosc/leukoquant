#!/bin/bash
#
# Build a study FA skeleton in MNI152 space using FMRIB58_FA_1mm and TBSS tools.
#
# Inputs are an FA image specified full path
#
# Outputs:
#   - skeleton.nii.gz
#

set -e

FA_PATH=""
OUTPUT_DIR=""
VERBOSE="false"
METHOD="2"

parse_args() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --fa-path)
                FA_PATH="$2"
                shift 2
                ;;
            --output-dir)
                OUTPUT_DIR="$2"
                shift 2
                ;;
            --verbose)
                VERBOSE="$2"
                shift 2
                ;;
            --method)
                METHOD="$2"
                shift 2
                ;;
            *)
                echo "Unknown option: $1"
                exit 1
                ;;
        esac
    done
}

check_fsl_env() {
    if [ -z "${FSLDIR:-}" ] || [ ! -f "${FSLDIR}/etc/fslconf/fsl.sh" ]; then
        echo "ERROR: FSLDIR is not set or invalid."
        echo "Please export FSLDIR and source fsl.sh before running."
        exit 1
    fi

    # shellcheck disable=SC1090
    source "${FSLDIR}/etc/fslconf/fsl.sh" >/dev/null 2>&1 || true

    if ! command -v fslmaths >/dev/null 2>&1; then
        echo "ERROR: fslmaths not found in PATH after sourcing FSL."
        exit 1
    fi

    if ! command -v tbss_skeleton >/dev/null 2>&1; then
        echo "ERROR: tbss_skeleton not found in PATH after sourcing FSL."
        exit 1
    fi
}

check_inputs() {
    if [ -z "$FA_PATH" ] || [ -z "$OUTPUT_DIR" ]; then
        echo "ERROR: Missing required parameters."
        echo "Usage: $0 --fa-path /full/path/to/fa.nii.gz --output-dir DIR [--verbose true]"
        exit 1
    fi

    if [ ! -f "$FA_PATH" ]; then
        echo "ERROR: FA image not found: $FA_PATH"
        exit 1
    fi

    mkdir -p "$OUTPUT_DIR"

    if [ "$VERBOSE" = "true" ] || [ "$VERBOSE" = "True" ]; then
        echo "=========================================="
        echo "FA Skeleton Builder"
        echo "=========================================="
        echo "FA image: $FA_PATH"
        echo "Output dir: $OUTPUT_DIR"
    fi
}

build_fmrib58_skeleton() {
    TEMP_DIR=$(mktemp -d)
    trap 'rm -rf "$TEMP_DIR"' EXIT
    FA_TEMP="$TEMP_DIR/FA.nii.gz"
    cp "$FA_PATH" "$FA_TEMP"

    OUTPUT_FILE="$TEMP_DIR/skeleton_mask.nii.gz"

    (
        cd "$TEMP_DIR"
        if [ "$METHOD" = "1" ]; then
            METHOD_STRING="FMRIB58 Template"
        elif [ "$METHOD" = "2" ]; then
            METHOD_STRING="Direct"
        elif [ "$METHOD" = "3" ]; then
            METHOD_STRING="Intersection of FMRIB58 and Direct"
        else
            echo "ERROR: Invalid method specified: $METHOD"
            exit 1
        fi

        echo "Building FA skeleton using TBSS method (Method $METHOD_STRING)..."
        
        # Method 1: Use TBSS to get the skeleton in MNI space from FRMIB58 template, then warp it back to native space
        # 1. Preprocess
        # Creates directory: origdata/
        # Creates directory: FA/
        # Creates file: origdata/FA.nii.gz (copy of original)
        # Creates file: FA/FA_FA.nii.gz (slightly eroded, preprocessed FA)
        # Creates file: FA/FA_FA_mask.nii.gz (brain mask)
        tbss_1_preproc "FA.nii.gz" 

        if [ "$METHOD" = "1" ] || [ "$METHOD" = "3" ]; then
            # 2. Register to FMRIB58 standard space
            # Creates file: FA/FA_FA_to_target_warp.nii.gz (the non-linear warp to MNI space)
            tbss_2_reg -T

            # 3. Post-registration and Skeletonization
            # Creates directory: stats/
            # Creates file: stats/all_FA.nii.gz (FA image mapped into standard MNI space)
            # Creates file: stats/mean_FA.nii.gz (same as above since N=1)
            # Creates file: stats/mean_FA_skeleton.nii.gz (the raw, unthresholded white matter skeleton in standard space)
            tbss_3_postreg -T

            # 4. Threshold the skeleton at 0.2
            # Creates file: stats/mean_FA_skeleton_mask.nii.gz (the cleaned, binary skeleton mask in standard space)
            # Creates file: stats/all_FA_skeletonised.nii.gz (FA values projected onto the skeleton)
            tbss_4_prestats 0.2

            # 5. Create the inverse warp and map the skeleton back to native space
            WARP_FILE=$(ls FA/*_to_target_warp.nii.gz | head -n 1)
            # Create the inverse warp (Standard space -> Subject space)
            invwarp -w "$WARP_FILE" \
                    -o FA/target_to_FA_warp.nii.gz \
                    -r FA/FA_FA.nii.gz

            # Map the MNI skeleton (mean_FA_skeleton) back to Subject Space
            applywarp -i stats/mean_FA_skeleton.nii.gz \
                    -r FA/FA_FA.nii.gz \
                    -w FA/target_to_FA_warp.nii.gz \
                    -o skeleton_MNI_in_Native.nii.gz \
                    --interp=trilinear

            # Threshold the skeleton in native space to create a binary mask. Image was multiplied by 10000 during the TBSS process, so we need to scale the image back before.
            fslmaths skeleton_MNI_in_Native.nii.gz \
                    -mul 0.0001 \
                    -thr 0.2 \
                    -bin \
                    skeleton_mask_template.nii.gz
            
            cp skeleton_mask_template.nii.gz "$OUTPUT_DIR/skeleton_mask_template.nii.gz"
        fi
        if [ "$METHOD" = "2" ] || [ "$METHOD" = "3" ]; then
            # Method 2
            # Run tbss_skeleton directly on the preprocessed FA image to get the skeleton in native space
            tbss_skeleton -i FA/FA_FA.nii.gz \
                          -o skeleton_mask.nii.gz

            # Threshold the skeleton to create a binary mask
            fslmaths skeleton_mask.nii.gz \
                    -thr 0.2 \
                    -bin \
                    skeleton_mask_direct.nii.gz
            
            cp skeleton_mask_direct.nii.gz "$OUTPUT_DIR/skeleton_mask_direct.nii.gz"
        fi
        if [ "$METHOD" = "3" ]; then
            # Method 3: Run both methods and compare
            fslmaths skeleton_mask_template.nii.gz -mul skeleton_mask_direct.nii.gz \
                    -bin \
                    skeleton_mask_intersection.nii.gz
            
            # Union of the two skeletons (for visualization)
            fslmaths skeleton_mask_template.nii.gz -add skeleton_mask_direct.nii.gz \
                    -bin \
                    skeleton_mask_union.nii.gz
            
            cp skeleton_mask_intersection.nii.gz "$OUTPUT_DIR/skeleton_mask_intersection.nii.gz"
            cp skeleton_mask_union.nii.gz "$OUTPUT_DIR/skeleton_mask_union.nii.gz"
        fi

    )

}


main() {
    parse_args "$@"
    check_fsl_env
    check_inputs
    build_fmrib58_skeleton
}

main "$@"
