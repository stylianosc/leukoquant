"""
Snakemake workflow to run FreeSurfer `recon-all` for one or more subjects.

Output structure: {output_dir}/{subject}/recon-all/outputs/mri/...
                  {output_dir}/{subject}/recon-all/logs/...
"""

import os
import re
import sys

subjects = config.get("subject_id")
if isinstance(subjects, str):
    subjects = [subjects]

t1_files = config.get("t1_file")
if isinstance(t1_files, str):
    t1_files = [t1_files]

t1_files_singularity = config.get("t1_file_singularity")
if isinstance(t1_files_singularity, str):
    t1_files_singularity = [t1_files_singularity]

output_dir             = config.get("output_dir")
output_dir_singularity = config.get("output_dir_singularity", "/output")
keep_intermediate      = config.get("keep_intermediate", False)

CONTAINER_NAME = config.get("container_name", "freesurfer_unified_container")
LEUKOQUANT_PARENT_DIR  = config.get("leukoquant_parent_dir", "")
CONTAINER_SIF = os.path.join(
    LEUKOQUANT_PARENT_DIR,
    "leukoquant/workflow/containers/",
    CONTAINER_NAME + ".sif",
)

sys.path.insert(0, LEUKOQUANT_PARENT_DIR)
from leukoquant.utils.container_utils import ensure_container
ensure_container(CONTAINER_SIF)

TOOL_NAME = "recon-all"

# Per-subject lookup dicts - same pattern as gif_workflow.smk
t1_map      = dict(zip(subjects, t1_files))
t1_sing_map = dict(zip(subjects, t1_files_singularity))

container:
    CONTAINER_SIF

wildcard_constraints:
    subject="|".join(re.escape(s) for s in subjects),

rule all:
    input:
        # ThalamicNuclei is intentionally NOT listed here: tracula's own rule
        # already declares it as a hard input (see thalamic_nuclei=... in
        # tracula_workflow.smk) and will correctly pull it in whenever tracula
        # itself still needs to run. Listing it here too would force it (and
        # therefore a full recon-all rerun to backfill aseg.mgz) for every
        # subject unconditionally, even ones whose tracula/metrics are already
        # complete and have no remaining need for it.
        expand(f"{output_dir}/{{subject}}/{TOOL_NAME}/outputs/mri/aparc+aseg.mgz", subject=subjects),

rule recon_all:
    input:
        t1=lambda wildcards: t1_map[wildcards.subject],
    output:
        aparc_aseg=f"{output_dir}/{{subject}}/{TOOL_NAME}/outputs/mri/aparc+aseg.mgz",
        # aseg, wmparc, norm, talairach.xfm: the files segmentThalamicNuclei.sh
        # explicitly checks for before running. Declaring them here (and as inputs
        # to thalamic_seg) lets Snakemake enforce the dependency precisely.
        aseg=f"{output_dir}/{{subject}}/{TOOL_NAME}/outputs/mri/aseg.mgz",
        wmparc=f"{output_dir}/{{subject}}/{TOOL_NAME}/outputs/mri/wmparc.mgz",
        norm=f"{output_dir}/{{subject}}/{TOOL_NAME}/outputs/mri/norm.mgz",
        talairach=f"{output_dir}/{{subject}}/{TOOL_NAME}/outputs/mri/transforms/talairach.xfm",
        # rh/lh.pial: Snakemake re-runs recon-all when the pial stage never
        # completed. rh.pial survives both pruned outputs (plain file after
        # thalamic_seg cp -Lf) and non-pruned outputs (symlink to rh.pial.T1).
        rh_pial=f"{output_dir}/{{subject}}/{TOOL_NAME}/outputs/surf/rh.pial",
        lh_pial=f"{output_dir}/{{subject}}/{TOOL_NAME}/outputs/surf/lh.pial",
    params:
        t1_file_singularity=lambda wildcards: t1_sing_map[wildcards.subject],
        recon_sd=lambda wildcards: f"{output_dir_singularity}/{wildcards.subject}/{TOOL_NAME}",
        recon_subject_name="outputs",
        log_file=lambda wildcards:   f"{output_dir_singularity}/{wildcards.subject}/{TOOL_NAME}/logs/recon_all.log",
        error_file=lambda wildcards: f"{output_dir_singularity}/{wildcards.subject}/{TOOL_NAME}/logs/recon_all_error.log",
    resources:
        mem_mb=12*1024,
        time="168:00:00",
        name="recon_all",
        scratch_size=1*1024,
        workdir=lambda wildcards: f"{output_dir}/{wildcards.subject}/{TOOL_NAME}",
    shell:
        """
        # Log files go to NFS (small sequential writes); all recon-all I/O
        # goes to local scratch to avoid hammering the SAN.
        source /leukoquant/leukoquant/utils/bash_utils.sh
        mkdir -p "$(dirname "{params.log_file}")"
        exec > "{params.log_file}"
        exec 2> "{params.error_file}"

        echo "Date: $(date)"
        echo "Date: $(date)" >&2

        enable_line_buffering

        export FS_LICENSE=/license/license.txt

        recon_sd="{params.recon_sd}"
        recon_name="{params.recon_subject_name}"
        recon_subject_dir="$recon_sd/$recon_name"

        MY_JOB_ID="$(get_job_id)"
        scratch_sd="/scratch0/$USER/$MY_JOB_ID/recon_all_{wildcards.subject}"
        scratch_subject_dir="$scratch_sd/$recon_name"

        # Always clean up scratch on exit, success or failure.
        trap 'scratch_cleanup "$scratch_sd"' EXIT INT TERM

        echo "Starting recon-all for subject {wildcards.subject}"
        echo "  T1:     {input.t1}"
        echo "  Job ID: $MY_JOB_ID"
        echo "  Scratch SUBJECTS_DIR: $scratch_sd"
        echo "  NFS output dir:       $recon_sd"

        mkdir -p "$scratch_sd"
        rm -rf "$scratch_subject_dir"
        export TMPDIR="$scratch_sd/tmp"
        mkdir -p "$TMPDIR"

        # Some raw T1 acquisitions are exported as 4D with a duplicate second
        # frame instead of a normal 3D volume -- recon-all rejects this
        # outright ("input(s) cannot have multiple frames!"). Collapse to the
        # first frame into scratch before running (no-op if already
        # single-volume -- see image_utils.py:select_first_frame); the
        # original NFS T1 is never touched.
        t1_singleframe="$scratch_sd/t1_singleframe.nii.gz"
        echo "Selecting first frame of T1: {params.t1_file_singularity} -> $t1_singleframe"
        python /leukoquant/leukoquant/utils/image_utils.py \
            --input  "{params.t1_file_singularity}" \
            --output "$t1_singleframe" \
            --op first_frame \
            --verbose

        recon-all -i "$t1_singleframe" -s "$recon_name" -sd "$scratch_sd" -all

        # Copy the whitelist files from scratch to NFS.
        # -L dereferences the full symlink chain at the kernel level (handles
        # double symlinks like rh.pial → rh.pial.T1 in FreeSurfer 7+).
        # -r handles any directory entries; -u skips files already up to date.
        echo "Copying whitelisted outputs from scratch to NFS..."
        mkdir -p "$recon_subject_dir"

        files_to_copy=(
            "mri/aparc+aseg.mgz"
            "mri/aseg.mgz"
            "mri/wmparc.mgz"
            "mri/orig.mgz"
            "mri/T1.mgz"
            "mri/norm.mgz"
            "mri/brainmask.mgz"
            "mri/brain.mgz"
            "mri/wm.mgz"
            "mri/ribbon.mgz"
            "mri/aseg.presurf.mgz"
            "mri/transforms/talairach.xfm"
            "surf/lh.white" "surf/rh.white"
            "surf/lh.pial"  "surf/rh.pial"
            "surf/lh.sphere.reg" "surf/rh.sphere.reg"
            "surf/lh.thickness"  "surf/rh.thickness"
            "label/lh.cortex.label" "label/rh.cortex.label"
            "stats/brainvol.stats"
        )

        copy_ok=true
        for file in "${{files_to_copy[@]}}"; do
            src="$scratch_subject_dir/$file"
            dst="$recon_subject_dir/$file"
            if [ -e "$src" ]; then
                mkdir -p "$(dirname "$dst")"
                cp -ruLf "$src" "$dst"
            else
                echo "WARNING: expected file missing on scratch: $file" >&2
                copy_ok=false
            fi
        done

        if [ "$copy_ok" = false ]; then
            echo "ERROR: one or more whitelist files were missing from scratch - aborting." >&2
            exit 1
        fi

        echo "Copy complete. Scratch will be cleaned up by trap on exit."
        """

rule thalamic_seg:
    input:
        aparc_aseg=f"{output_dir}/{{subject}}/{TOOL_NAME}/outputs/mri/aparc+aseg.mgz",
        aseg=f"{output_dir}/{{subject}}/{TOOL_NAME}/outputs/mri/aseg.mgz",
        wmparc=f"{output_dir}/{{subject}}/{TOOL_NAME}/outputs/mri/wmparc.mgz",
        norm=f"{output_dir}/{{subject}}/{TOOL_NAME}/outputs/mri/norm.mgz",
        talairach=f"{output_dir}/{{subject}}/{TOOL_NAME}/outputs/mri/transforms/talairach.xfm",
    output:
        thalamic_nuclei=f"{output_dir}/{{subject}}/{TOOL_NAME}/outputs/mri/ThalamicNuclei.v13.T1.FSvoxelSpace.mgz",
    params:
        recon_sd=lambda wildcards: f"{output_dir_singularity}/{wildcards.subject}/{TOOL_NAME}",
        recon_subject_name="outputs",
        log_file=lambda wildcards:   f"{output_dir_singularity}/{wildcards.subject}/{TOOL_NAME}/logs/thalamic_seg.log",
        error_file=lambda wildcards: f"{output_dir_singularity}/{wildcards.subject}/{TOOL_NAME}/logs/thalamic_seg_error.log",
    resources:
        mem_mb=8*1024,
        time="01:00:00",
        name="thalamic_seg",
        workdir=lambda wildcards: f"{output_dir}/{wildcards.subject}/{TOOL_NAME}",
    shell:
        """
        set -euo pipefail

        source /leukoquant/leukoquant/utils/bash_utils.sh
        log_dir="$(dirname "{params.log_file}")"
        mkdir -p "$log_dir"
        if [ ! -w "$log_dir" ]; then
            echo "FATAL: thalamic_seg log directory not writable: $log_dir" >&2
            exit 1
        fi
        exec > "{params.log_file}"
        exec 2> "{params.error_file}"

        echo "Date: $(date)"
        echo "Date: $(date)" >&2
        enable_line_buffering

        export FS_LICENSE=/license/license.txt
        export SUBJECTS_DIR="{params.recon_sd}"

        # recon-all's NFS copy step only whitelists mri/surf/label/stats files
        # (to save space/inodes) and never copies scripts/ -- but
        # segmentThalamicNuclei.sh unconditionally writes its own lock file to
        # scripts/IsRunningThalamicNuclei_mainFreeSurferT1, so without this the
        # tool fails immediately with "No such file or directory" every time,
        # for every subject (confirmed: 850/850 failures in one EPAD batch).
        mkdir -p "{params.recon_sd}/{params.recon_subject_name}/scripts"

        # Snakemake only schedules this rule when ThalamicNuclei output is
        # missing - any IsRunning* lock left by a previous attempt is therefore
        # stale (the prior run crashed before clearing it) and would otherwise
        # cause segmentThalamicNuclei.sh to refuse to start.
        rm -f "{params.recon_sd}/{params.recon_subject_name}/scripts/IsRunningThalamicNuclei_mainFreeSurferT1"

        echo "Starting thalamic segmentation for subject {wildcards.subject}"
        segmentThalamicNuclei.sh {params.recon_subject_name} {params.recon_sd}
        """
