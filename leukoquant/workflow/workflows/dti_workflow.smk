"""
Snakemake workflow to run DTI fitting.
Supports one or more subjects via the {subject} wildcard.

Output structure: {output_dir}/{subject}/dti/outputs/fa.nii.gz  (and other metric maps)
                  {output_dir}/{subject}/dti/intermediate/metadata.json
                  {output_dir}/{subject}/dti/logs/dti.log
"""
import os
import re
import sys

subjects = config.get("subjects", [])
if not subjects:
    subjects = [config.get("subject", os.path.basename(config.get("output_dir", "unknown")))]
if isinstance(subjects, str):
    subjects = [subjects]

dwi_paths      = config.get("dwi_paths",      [config.get("dwi_path", "")])
if isinstance(dwi_paths, str):
    dwi_paths = [dwi_paths]

dwi_paths_sing = config.get("dwi_paths_singularity", [config.get("dwi_path_singularity", "")])
if isinstance(dwi_paths_sing, str):
    dwi_paths_sing = [dwi_paths_sing]

mask_paths      = config.get("brain_mask_paths",      [config.get("brain_mask_path", "")])
if isinstance(mask_paths, str):
    mask_paths = [mask_paths]

mask_paths_sing = config.get("brain_mask_paths_singularity", [config.get("brain_mask_path_singularity", "")])
if isinstance(mask_paths_sing, str):
    mask_paths_sing = [mask_paths_sing]

bvecs_paths      = config.get("bvecs_paths",      [config.get("bvecs_path", "")])
if isinstance(bvecs_paths, str):
    bvecs_paths = [bvecs_paths]

bvecs_paths_sing = config.get("bvecs_paths_singularity", [config.get("bvecs_path_singularity", "")])
if isinstance(bvecs_paths_sing, str):
    bvecs_paths_sing = [bvecs_paths_sing]

bvals_paths      = config.get("bvals_paths",      [config.get("bvals_path", "")])
if isinstance(bvals_paths, str):
    bvals_paths = [bvals_paths]

bvals_paths_sing = config.get("bvals_paths_singularity", [config.get("bvals_path_singularity", "")])
if isinstance(bvals_paths_sing, str):
    bvals_paths_sing = [bvals_paths_sing]

n = len(subjects)
def _pad(lst, n): return lst + [""] * (n - len(lst))
mask_paths       = _pad(mask_paths,       n)
mask_paths_sing  = _pad(mask_paths_sing,  n)
bvecs_paths      = _pad(bvecs_paths,      n)
bvecs_paths_sing = _pad(bvecs_paths_sing, n)
bvals_paths      = _pad(bvals_paths,      n)
bvals_paths_sing = _pad(bvals_paths_sing, n)

dwi_map        = dict(zip(subjects, dwi_paths))
dwi_sing_map   = dict(zip(subjects, dwi_paths_sing))
mask_map       = dict(zip(subjects, mask_paths))
mask_sing_map  = dict(zip(subjects, mask_paths_sing))
bvecs_map      = dict(zip(subjects, bvecs_paths))
bvecs_sing_map = dict(zip(subjects, bvecs_paths_sing))
bvals_map      = dict(zip(subjects, bvals_paths))
bvals_sing_map = dict(zip(subjects, bvals_paths_sing))

OUTPUT_DIR      = config.get("output_dir")
OUTPUT_DIR_SING = config.get("output_dir_singularity", "/output")
SCRATCH_DIR_SING = config.get("scratch_dir_singularity", "/scratch0")
KEEP_INTERMEDIATE = config.get("keep_intermediate", False)
TOOL_NAME       = "dti"

LEUKOQUANT_DIR = config.get("leukoquant_dir")

sys.path.insert(0, os.path.dirname(LEUKOQUANT_DIR))
from leukoquant.utils.container_utils import ensure_container, FREESURFER_SIF_HF_PATH

SKULL_STRIP = config.get("skull_strip", False)

if SKULL_STRIP:
    CONTAINER_SIF = LEUKOQUANT_DIR + "/workflow/containers/freesurfer_unified_container.sif"
    ensure_container(CONTAINER_SIF, filename=FREESURFER_SIF_HF_PATH)
else:
    CONTAINER_SIF = LEUKOQUANT_DIR + "/workflow/containers/miniconda_unified_container.sif"
    ensure_container(CONTAINER_SIF)

container:
    CONTAINER_SIF

wildcard_constraints:
    subject="|".join(re.escape(s) for s in subjects),

rule all:
    input:
        [f"{OUTPUT_DIR}/{s}/{TOOL_NAME}/outputs/fa.nii.gz" for s in subjects]

rule dti:
    input:
        dwi=lambda wildcards: dwi_map[wildcards.subject],
        brain_mask=lambda wildcards: (
            [mask_map[wildcards.subject]] if mask_map.get(wildcards.subject) else []
        ),
    output:
        dti_fa=f"{OUTPUT_DIR}/{{subject}}/{TOOL_NAME}/outputs/fa.nii.gz",
        dti_md=f"{OUTPUT_DIR}/{{subject}}/{TOOL_NAME}/outputs/md.nii.gz",
        dti_rd=f"{OUTPUT_DIR}/{{subject}}/{TOOL_NAME}/outputs/rd.nii.gz",
        dti_ad=f"{OUTPUT_DIR}/{{subject}}/{TOOL_NAME}/outputs/ad.nii.gz",
    params:
        dwi_sing=lambda wildcards: dwi_sing_map[wildcards.subject],
        mask_sing=lambda wildcards: mask_sing_map.get(wildcards.subject, ""),
        bvecs_sing=lambda wildcards: bvecs_sing_map.get(wildcards.subject, ""),
        bvals_sing=lambda wildcards: bvals_sing_map.get(wildcards.subject, ""),
        # Container path to the per-subject tool directory
        tool_dir_sing=lambda wildcards: f"{OUTPUT_DIR_SING}/{wildcards.subject}/{TOOL_NAME}",
        keep_intermediate=KEEP_INTERMEDIATE,
        skull_strip=SKULL_STRIP,
        log_file=lambda wildcards: config.get("log_file", f"{OUTPUT_DIR_SING}/{wildcards.subject}/{TOOL_NAME}/logs/dti.log"),
        error_file=lambda wildcards: config.get("error_file", f"{OUTPUT_DIR_SING}/{wildcards.subject}/{TOOL_NAME}/logs/dti_error.log"),
        scratch_dir_singularity=SCRATCH_DIR_SING,
    resources:
        mem_mb=8*1024,
        time="24:00:00",
        scratch_size=0.25*1024,
        name="dtifit",
        workdir=lambda wildcards: f"{OUTPUT_DIR}/{wildcards.subject}/{TOOL_NAME}",
    shell:
        """
        source /leukoquant/leukoquant/utils/bash_utils.sh
        mkdir -p "$(dirname "{params.log_file}")"
        exec > "{params.log_file}"
        exec 2> "{params.error_file}"

        echo "Date: $(date)"
        echo "Date: $(date)" >&2

        enable_line_buffering

        MY_JOB_ID="$(get_job_id)"
        scratch_work_dir=""
        trap 'scratch_cleanup "$scratch_work_dir"' EXIT INT TERM

        tool_dir="{params.tool_dir_sing}"
        outputs_dir="$tool_dir/outputs"
        output_intermediate_dir="$tool_dir/intermediate"
        scratch_work_dir="{params.scratch_dir_singularity}/$USER/$MY_JOB_ID/dti"
        mkdir -p "$outputs_dir" "$scratch_work_dir"

        bvecs_arg=""
        if [ -n "{params.bvecs_sing}" ]; then
            bvecs_arg="--bvecs {params.bvecs_sing}"
        fi
        bvals_arg=""
        if [ -n "{params.bvals_sing}" ]; then
            bvals_arg="--bvals {params.bvals_sing}"
        fi

        # Prepare DWI (writes data.nii.gz, bvecs, bvals, metadata.json to scratch_work_dir)
        echo "Starting input preparation..."
        dwi_prepare_cmd="python /leukoquant/leukoquant/utils/dwi_utils.py --dwi {params.dwi_sing} $bvecs_arg $bvals_arg --outdir $scratch_work_dir"
        echo "Command: $dwi_prepare_cmd"
        $dwi_prepare_cmd

        dwi_output="$scratch_work_dir/data.nii.gz"
        brain_mask_output="$scratch_work_dir/brain_mask.nii.gz"
        bvecs_output="$scratch_work_dir/bvecs"
        bvals_output="$scratch_work_dir/bvals"

        mask_arg=""
        if [ -n "{params.mask_sing}" ]; then
            cp "{params.mask_sing}" "$brain_mask_output"
            mask_arg="--mask $brain_mask_output"
        elif [ "{params.skull_strip}" = "True" ]; then
            echo "Running skull stripping with mri_synthstrip..."
            mri_synthstrip -i "$dwi_output" -o "$scratch_work_dir/dwi_stripped.nii.gz" -m "$brain_mask_output"
            mask_arg="--mask $brain_mask_output"
            echo "Skull stripping complete"
        fi

        # Run DTI fitting - metric maps are written to outputs/.
        save_intermediate_flag=""
        if [ "{params.keep_intermediate}" = "True" ]; then
            save_intermediate_flag="--save-intermediate"
        fi
        echo "Running DTI fitting..."
        dtifit_cmd="python /leukoquant/leukoquant/utils/reconst_dti.py --dwi $dwi_output --bvecs $bvecs_output --bvals $bvals_output $mask_arg --outdir $outputs_dir $save_intermediate_flag"
        echo "Command: $dtifit_cmd"
        $dtifit_cmd

        # Move intermediate files to intermediate/ if requested; scratch is cleaned by trap on exit.
        # masked_dwi.nii.gz is only written by reconst_dti.py when --save-intermediate is passed.
        intermediates=("$scratch_work_dir/metadata.json" "$outputs_dir/masked_dwi.nii.gz")
        if [ "{params.keep_intermediate}" = "True" ]; then
            mkdir -p "$output_intermediate_dir"
            for f in "${{intermediates[@]}}"; do
                [ -f "$f" ] && mv "$f" "$output_intermediate_dir/" || true
            done
        fi
        """