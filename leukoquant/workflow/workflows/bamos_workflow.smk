"""
Snakemake workflow for BaMoS (Bayesian Model Selection) processing.
Supports one or more subjects via the {subject} wildcard.

Output structure: {output_dir}/{subject}/bamos/outputs/CorrectLesion.nii.gz
                  {output_dir}/{subject}/bamos/intermediate/  (all other working files)
                  {output_dir}/{subject}/bamos/logs/
"""
import glob as _glob
import os
import re
import sys

subjects        = config.get("subjects", [config.get("subject_id", "")])
if isinstance(subjects, str):
    subjects = [subjects]

flair_files     = config.get("flair_files", [config.get("flair_file", "")])
if isinstance(flair_files, str):
    flair_files = [flair_files]

flair_files_sing = config.get("flair_files_singularity", [config.get("flair_file_singularity", "")])
if isinstance(flair_files_sing, str):
    flair_files_sing = [flair_files_sing]

t1_files        = config.get("t1_files", [config.get("t1_file", "")])
if isinstance(t1_files, str):
    t1_files = [t1_files]

t1_files_sing   = config.get("t1_files_singularity", [config.get("t1_file_singularity", "")])
if isinstance(t1_files_sing, str):
    t1_files_sing = [t1_files_sing]

gif_results_paths      = config.get("gif_results_paths", [config.get("gif_results_path", "")])
if isinstance(gif_results_paths, str):
    gif_results_paths = [gif_results_paths]

gif_results_paths_sing = config.get("gif_results_paths_singularity", [config.get("gif_results_path_singularity", "/gif_input")])
if isinstance(gif_results_paths_sing, str):
    gif_results_paths_sing = [gif_results_paths_sing]

OUTPUT_DIR      = config.get("output_dir")
OUTPUT_DIR_SING = config.get("output_dir_singularity", "/bamos_output")
SCRATCH_DIR_SING = config.get("scratch_dir_singularity", "/scratch0")
KEEP_INTERMEDIATE = config.get("keep_intermediate", False)
TOOL_NAME       = "bamos"

# Per-subject lookup dicts
flair_map      = dict(zip(subjects, flair_files))
flair_sing_map = dict(zip(subjects, flair_files_sing))
t1_map         = dict(zip(subjects, t1_files))
t1_sing_map    = dict(zip(subjects, t1_files_sing))
gif_map        = dict(zip(subjects, gif_results_paths))
gif_sing_map   = dict(zip(subjects, gif_results_paths_sing))

def _gif_parcellation_filename(subject):
    """Return the actual parcellation filename produced by GIF for this subject.

    For pre-existing external GIF directories, glob for the real filename.
    For the internal path (GIF will be run as part of this workflow), derive
    the expected filename from the T1 stem - the same stem the GIF workflow uses.
    """
    gif_path = gif_map.get(subject, "")
    if gif_path:
        results = _glob.glob(os.path.join(gif_path, "*_NeuroMorph_Parcellation.nii.gz"))
        if results:
            return results[0].split("/")[-1]
        # process_all_workflow.smk always sets gif_results_paths to the
        # internal {output_dir}/{s}/gif/outputs path for every subject,
        # whether or not GIF has actually run yet -- so a truthy gif_path
        # with no match here is ambiguous: it could be a genuinely external
        # directory with mismatched naming (real error), or process-all's
        # own internal path for a subject whose GIF hasn't run yet, so the
        # directory doesn't exist or is still empty (expected, not an
        # error). Only raise when the directory actually has *some* files in
        # it but none match our naming -- that's the real mismatch case.
        # Directory missing/empty falls through to the same T1-stem
        # prediction used for the genuinely-internal (gif_path empty) case
        # below, exactly as before this check was added.
        if os.path.isdir(gif_path) and os.listdir(gif_path):
            raise FileNotFoundError(
                f"No file matching '*_NeuroMorph_Parcellation.nii.gz' found in "
                f"external GIF results directory for subject '{subject}': {gif_path}\n"
                "Ensure the parcellation file's name ends with "
                "'_NeuroMorph_Parcellation.nii.gz', or check that --gif-results-dir "
                "points to the correct directory."
            )
    # Internal path: GIF names outputs after the primary input file stem (T1).
    image_stem = t1_sing_map[subject].split("/")[-1].split(".")[0]
    return f"{image_stem}_NeuroMorph_Parcellation.nii.gz"

gif_parc_filename_map = {s: _gif_parcellation_filename(s) for s in subjects}

def _gif_parcellation_path(subject):
    """Return the full path to the GIF parcellation file."""
    gif_path = gif_map.get(subject, "")
    if gif_path:
        return f"{gif_path}/{gif_parc_filename_map[subject]}"
    else:
        # Expected output from the gif module
        return f"{OUTPUT_DIR}/{subject}/gif/outputs/{t1_sing_map[subject].split('/')[-1].split('.')[0]}_NeuroMorph_Parcellation.nii.gz"

def _gif_output_dir(subject):
    gif_path = gif_map.get(subject, "")
    if gif_path:
        return gif_path
    else:
        return f"{OUTPUT_DIR}/{subject}/gif/outputs"


SPACE      = int(config.get("space", 1))
JUMP_START = int(config.get("jump_start", 0))
OPT        = config.get("opt", "TA")

SCRIPTS_DIR              = "/leukoquant/leukoquant/external/bamos/scripts"
BAMOS_SCRIPT             = f"{SCRIPTS_DIR}/BaMoS_WMH_080526_local.sh"
BAMOS_CORRECTIONS_SCRIPT = f"{SCRIPTS_DIR}/correction_lesions_111125.py"

LEUKOQUANT_PARENT_DIR = config.get("leukoquant_parent_dir")

# When the user has already run GIF separately, redirect the GIF module to look
# for parcellation files in those external directories instead of re-running GIF.
gif_output_dirs = [
    gif_map[s] if gif_map.get(s) else f"{OUTPUT_DIR}/{s}/gif/outputs"
    for s in subjects
]

# image_ids: stem of the primary input image used by GIF.
gif_image_ids = [gif_parc_filename_map[s].replace("_NeuroMorph_Parcellation.nii.gz", "") for s in subjects]

config_gif = {
    "subjects": subjects,
    "t1_files": t1_files,
    "image_ids": gif_image_ids,
    "t1_files_singularity": t1_files_sing,
    "gif_output_dirs": gif_output_dirs,
    "output_dir": OUTPUT_DIR,
    "output_dir_singularity": OUTPUT_DIR_SING,
    "leukoquant_parent_dir": LEUKOQUANT_PARENT_DIR,
    "keep_intermediate": KEEP_INTERMEDIATE,
    "gif_home_host": config.get("gif_home_host", ""),
}
module gif:
    snakefile: "gif_workflow.smk"
    config: config_gif
# Import all GIF rules under the gif_* namespace so gif_all can be referenced below.
use rule * from gif as gif_*


sys.path.insert(0, LEUKOQUANT_PARENT_DIR)
from leukoquant.utils.container_utils import ensure_container

CONTAINER_SIF = os.path.join(
    LEUKOQUANT_PARENT_DIR,
    "leukoquant/workflow/containers/miniconda_unified_container.sif",
)
ensure_container(CONTAINER_SIF)

container:
    CONTAINER_SIF

# subject→image_id lookup used to resolve the GIF output path for each subject.
subject_to_image_id = dict(zip(subjects, gif_image_ids))

wildcard_constraints:
    subject="|".join(re.escape(s) for s in subjects),
    image_id="|".join(re.escape(iid) for iid in gif_image_ids),


rule all:
    input:
        # GIF parcellation files are a BaMoS prerequisite, but only include the
        # internal GIF targets for subjects that did NOT supply an external
        # --gif-results-dir. Including them unconditionally would force
        # gif_run_gif into the DAG even when parcellation files already exist
        # externally, causing spurious SGE submissions.
        [
            f"{OUTPUT_DIR}/{s}/gif/outputs/{gif_image_ids[i]}_NeuroMorph_Parcellation.nii.gz"
            for i, s in enumerate(subjects) if not gif_map.get(s)
        ],
        # BaMoS final corrected lesion masks
        [f"{OUTPUT_DIR}/{s}/{TOOL_NAME}/outputs/CorrectLesion.nii.gz" for s in subjects],

rule run_bamos:
    input:
        flair=lambda wildcards: flair_map[wildcards.subject],
        t1=lambda wildcards: t1_map[wildcards.subject],
        # When the user supplied --gif-results-dir for this subject, point at
        # the pre-existing external parcellation file directly so gif_run_gif
        # is never added to the DAG (no spurious SGE submission). Otherwise,
        # reference rules.gif_run_gif.output so Snakemake wires the SGE job
        # dependency at DAG-build time (required for --immediate-submit).
        gif_parcellation_file=lambda wildcards: (
            _gif_parcellation_path(wildcards.subject)
            if gif_map.get(wildcards.subject)
            else expand(
                rules.gif_run_gif.output.parcellation,
                subject=wildcards.subject,
                image_id=subject_to_image_id[wildcards.subject],
            )[0]
        ),
    output:
        bamos_lesion_file=f"{OUTPUT_DIR}/{{subject}}/{TOOL_NAME}/intermediate/essential/Correct_WS3WT3WC1Lesion_corr.nii.gz",
        bamos_connectivity_file=f"{OUTPUT_DIR}/{{subject}}/{TOOL_NAME}/intermediate/essential/Connect_WS3WT3WC1Lesion_corr.nii.gz",
        bamos_label_file=f"{OUTPUT_DIR}/{{subject}}/{TOOL_NAME}/intermediate/essential/TxtLesion_WS3WT3WC1Lesion_corr.txt",
    resources:
        mem_mb=12*1024,
        scratch_size=5*1024,
        time="168:00:00" if JUMP_START == 0 else "2:00:00",
        name="BaMoS" if JUMP_START == 0 else f"BaMoSLes_{wildcards.subject}",
        workdir=lambda wildcards: f"{OUTPUT_DIR}/{wildcards.subject}/{TOOL_NAME}",
    params:
        bamos_script=BAMOS_SCRIPT,
        flair_sing=lambda wildcards: flair_sing_map[wildcards.subject],
        t1_sing=lambda wildcards: t1_sing_map[wildcards.subject],
        gif_path_sing=lambda wildcards: gif_sing_map[wildcards.subject] if gif_map.get(wildcards.subject) else f"{OUTPUT_DIR_SING}/{wildcards.subject}/gif/outputs",
        scratch_dir_sing=SCRATCH_DIR_SING,
        bamos_essential_final=lambda wildcards: f"{OUTPUT_DIR_SING}/{wildcards.subject}/{TOOL_NAME}/intermediate/essential",
        bamos_intermediate_final=lambda wildcards: f"{OUTPUT_DIR_SING}/{wildcards.subject}/{TOOL_NAME}/intermediate/other",
        jump_start=JUMP_START,
        opt=OPT,
        space=SPACE,
        keep_intermediate=KEEP_INTERMEDIATE,
        log_file=lambda wildcards: f"{OUTPUT_DIR_SING}/{wildcards.subject}/{TOOL_NAME}/logs/bamos.log",
        error_file=lambda wildcards: f"{OUTPUT_DIR_SING}/{wildcards.subject}/{TOOL_NAME}/logs/bamos_error.log",
    shell:
        """
        set -euo pipefail

        # Verify the log directory exists and is writable BEFORE redirecting.
        # If the bind mount is missing or the path is wrong, fail loudly to
        # SGE's stderr (visible in qacct) rather than silently losing every
        # message from the rest of the job to a broken file descriptor.
        source /leukoquant/leukoquant/utils/bash_utils.sh
        log_dir="$(dirname "{params.log_file}")"
        mkdir -p "$log_dir"
        if [ ! -w "$log_dir" ]; then
            echo "FATAL: bamos log directory not writable: $log_dir" >&2
            exit 1
        fi
        exec > "{params.log_file}"
        exec 2> "{params.error_file}"

        echo "Date: $(date)"
        echo "Date: $(date)" >&2

        enable_line_buffering

        scratch_path=""
        trap 'scratch_cleanup "$scratch_path"' EXIT INT TERM

        chmod +x {params.bamos_script}

        MY_JOB_ID="$(get_job_id)"

        scratch_path="{params.scratch_dir_sing}/$USER/$MY_JOB_ID/BaMoS/{wildcards.subject}"
        bamos_results_dir="$scratch_path/intermediate"
        mkdir -p "$bamos_results_dir"

        echo "Date: $(date)"
        echo "Date: $(date)" >&2
        echo "Running BaMoS for subject {wildcards.subject}"
        echo "Scratch path: $scratch_path"

        # Some raw T1/FLAIR acquisitions are exported as 4D with a duplicate
        # second frame instead of a normal 3D volume; collapse to the first
        # frame before anything else touches the image (no-op if already
        # single-volume -- see image_utils.py:select_first_frame).
        t1_singleframe="$scratch_path/t1_singleframe.nii.gz"
        flair_singleframe="$scratch_path/flair_singleframe.nii.gz"
        echo "Selecting first frame of T1: {params.t1_sing} -> $t1_singleframe"
        python /leukoquant/leukoquant/utils/image_utils.py \
            --input  "{params.t1_sing}" \
            --output "$t1_singleframe" \
            --op first_frame \
            --verbose
        echo "Selecting first frame of FLAIR: {params.flair_sing} -> $flair_singleframe"
        python /leukoquant/leukoquant/utils/image_utils.py \
            --input  "{params.flair_sing}" \
            --output "$flair_singleframe" \
            --op first_frame \
            --verbose

        # Reorient T1 and FLAIR to RAS+ before passing to BaMoS.
        # BaMoS assumes RAS orientation (see BaMoS_WMH_080526_local.sh line 5).
        t1_ras="$scratch_path/t1_ras.nii.gz"
        flair_ras="$scratch_path/flair_ras.nii.gz"
        echo "Reorienting T1 to RAS+: $t1_singleframe -> $t1_ras"
        python /leukoquant/leukoquant/utils/image_utils.py \
            --input  "$t1_singleframe" \
            --output "$t1_ras" \
            --verbose
        echo "Reorienting FLAIR to RAS+: $flair_singleframe -> $flair_ras"
        python /leukoquant/leukoquant/utils/image_utils.py \
            --input  "$flair_singleframe" \
            --output "$flair_ras" \
            --verbose

        # BAMOS_ID is a fixed token used only for ephemeral scratch filenames inside
        # the BaMoS shell script - it is never surfaced in final output filenames.
        BAMOS_ID="subject"

        bamos_cmd="{params.bamos_script} \
            $BAMOS_ID \
            $flair_ras \
            $t1_ras \
            {params.gif_path_sing} \
            $bamos_results_dir \
            {params.jump_start} \
            {params.opt} \
            {params.space} \
            $scratch_path"

        echo "Running $bamos_cmd"
        $bamos_cmd

        # BaMoS writes its outputs into $scratch_path/BaMoS_$BAMOS_ID/
        scratch_bamos_dir="$scratch_path/BaMoS_$BAMOS_ID"
        echo "Checking scratch BaMoS dir: $scratch_bamos_dir"
        if [ -d "$scratch_bamos_dir" ]; then
            echo "Copying essential BaMoS outputs to {params.bamos_essential_final}"
            mkdir -p "{params.bamos_essential_final}"

            # Associative array: pattern -> output filename (easy to add more entries)
            declare -A file_mappings=(
                ["Correct_WS3WT3WC1Lesion_*_corr.nii.gz"]="Correct_WS3WT3WC1Lesion_corr.nii.gz"
                ["Connect_WS3WT3WC1Lesion_*_corr.nii.gz"]="Connect_WS3WT3WC1Lesion_corr.nii.gz"
                ["TxtLesion_WS3WT3WC1Lesion_*_corr.txt"]="TxtLesion_WS3WT3WC1Lesion_corr.txt"
                ["LesionMahal_T1FLAIR_BiASM*.nii.gz"]="LesionMahal_T1FLAIR_BiASM.nii.gz"
            )

            for pattern in "${{!file_mappings[@]}}"; do
                output_name="${{file_mappings[$pattern]}}"
                for f in "$scratch_bamos_dir"/$pattern; do
                    [ -f "$f" ] && cp "$f" "{params.bamos_essential_final}/$output_name"
                done
            done

            echo "Essential files copied successfully"
        
            if [ "{params.keep_intermediate}" = "True" ]; then
                echo "Copying intermediate BaMoS outputs to {params.bamos_intermediate_final}"
                mkdir -p "{params.bamos_intermediate_final}"

                # Array of file patterns to exclude from intermediate copy (easy to add more patterns)
                exclude_patterns=(
                    "Correct_WS3WT3WC1Lesion_*"
                    "Connect_WS3WT3WC1Lesion_*"
                    "TxtLesion_WS3WT3WC1Lesion_*"
                    "LesionMahal_T1FLAIR_BiASM*"
                )

                find "$scratch_bamos_dir" -maxdepth 1 -type f \
                    $(printf '! -name "%s" ' "${{exclude_patterns[@]}}") \
                    -exec cp {{}} "{params.bamos_intermediate_final}/" \\;

                echo "Intermediate files copied successfully"
            fi
        else
            echo "ERROR: BaMoS output directory not found: $scratch_bamos_dir" >&2
            echo "Contents of scratch_path:" >&2
            ls -la "$scratch_path" >&2 || true
            exit 1
        fi

        echo "BaMoS processing completed"
        echo "Date: $(date)"
        echo "Date: $(date)" >&2
        """

rule run_bamos_correction:
    input:
        bamos_lesion_file=f"{OUTPUT_DIR}/{{subject}}/{TOOL_NAME}/intermediate/essential/Correct_WS3WT3WC1Lesion_corr.nii.gz",
        bamos_connectivity_file=f"{OUTPUT_DIR}/{{subject}}/{TOOL_NAME}/intermediate/essential/Connect_WS3WT3WC1Lesion_corr.nii.gz",
        bamos_label_file=f"{OUTPUT_DIR}/{{subject}}/{TOOL_NAME}/intermediate/essential/TxtLesion_WS3WT3WC1Lesion_corr.txt",
        gif_parcellation_file=lambda wildcards: _gif_parcellation_path(wildcards.subject),
    output:
        output_bamos_correction=f"{OUTPUT_DIR}/{{subject}}/{TOOL_NAME}/outputs/CorrectLesion.nii.gz",
    params:
        bamos_corrections_script=BAMOS_CORRECTIONS_SCRIPT,
        gif_parc_filename=lambda wildcards: gif_parc_filename_map[wildcards.subject],
        bamos_essential_sing=lambda wildcards: f"{OUTPUT_DIR_SING}/{wildcards.subject}/{TOOL_NAME}/intermediate/essential",
        bamos_outputs_sing=lambda wildcards: f"{OUTPUT_DIR_SING}/{wildcards.subject}/{TOOL_NAME}/outputs",
        gif_path_sing=lambda wildcards: gif_sing_map[wildcards.subject] if gif_map.get(wildcards.subject) else f"{OUTPUT_DIR_SING}/{wildcards.subject}/gif/outputs",
        keep_intermediate=KEEP_INTERMEDIATE,
        log_file=lambda wildcards: f"{OUTPUT_DIR_SING}/{wildcards.subject}/{TOOL_NAME}/logs/bamos_correction.log",
        error_file=lambda wildcards: f"{OUTPUT_DIR_SING}/{wildcards.subject}/{TOOL_NAME}/logs/bamos_correction_error.log",
    resources:
        mem_mb=4*1024,
        time="05:00:00",
        name="BaMoS_Corr",
        workdir=lambda wildcards: f"{OUTPUT_DIR}/{wildcards.subject}/{TOOL_NAME}",
    shell:
        """
        set -euo pipefail

        source /leukoquant/leukoquant/utils/bash_utils.sh
        log_dir="$(dirname "{params.log_file}")"
        mkdir -p "$log_dir"
        if [ ! -w "$log_dir" ]; then
            echo "FATAL: bamos log directory not writable: $log_dir" >&2
            exit 1
        fi
        exec > "{params.log_file}"
        exec 2> "{params.error_file}"

        echo "Date: $(date)"
        echo "Date: $(date)" >&2

        enable_line_buffering

        chmod +x {params.bamos_corrections_script}

        mkdir -p "{params.bamos_outputs_sing}"
        mkdir -p "{params.bamos_essential_sing}"

        echo "Date: $(date)"
        echo "Date: $(date)" >&2
        echo "Running BaMoS correction for subject {wildcards.subject}"

        bamos_lesion_file="{params.bamos_essential_sing}/Correct_WS3WT3WC1Lesion_corr.nii.gz"
        bamos_connectivity_file="{params.bamos_essential_sing}/Connect_WS3WT3WC1Lesion_corr.nii.gz"
        bamos_label_file="{params.bamos_essential_sing}/TxtLesion_WS3WT3WC1Lesion_corr.txt"
        gif_parcellation_file="{params.gif_path_sing}/{params.gif_parc_filename}"

        bamos_corr_cmd="python {params.bamos_corrections_script} \
            -les $bamos_lesion_file \
            -connect $bamos_connectivity_file \
            -label $bamos_label_file \
            -parc $gif_parcellation_file \
            -corr choroid cortex sheet"

        echo "Running $bamos_corr_cmd"
        $bamos_corr_cmd

        corr_output="{params.bamos_essential_sing}/CorrectLesion.nii.gz"
        if [ ! -f "$corr_output" ]; then
            echo "ERROR: correction output not found: $corr_output" >&2
            exit 1
        fi
        echo "Moving correction output to {params.bamos_outputs_sing}"
        mv "$corr_output" "{params.bamos_outputs_sing}/CorrectLesion.nii.gz"

        # Clean up essential intermediate files once correction is done
        #if [ "{params.keep_intermediate}" = "False" ]; then
        #    rm -rf "{params.bamos_essential_sing}"
        #fi

        echo "BaMoS correction completed for subject {wildcards.subject}"
        echo "Date: $(date)"
        echo "Date: $(date)" >&2
        """
