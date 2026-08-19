"""
Snakemake workflow for GIF (Geodesic Information Flow) processing.
Supports one or more subjects via the {subject} wildcard.

Output structure: {output_dir}/{subject}/gif/outputs/{image_id}_NeuroMorph_Parcellation.nii.gz
                  {output_dir}/{subject}/gif/intermediate/...
                  {output_dir}/{subject}/gif/logs/gif.log
"""
import os
import re
import sys

subjects        = config.get("subjects", [config.get("subject_id", "")])
if isinstance(subjects, str):
    subjects = [subjects]

t1_files        = config.get("t1_files", [config.get("t1_file", "")])
if isinstance(t1_files, str):
    t1_files = [t1_files]

flair_files     = config.get("flair_files", [config.get("flair_file", "")])
if isinstance(flair_files, str):
    flair_files = [flair_files]

# image_ids: stem of the primary input image (T1 if provided, else FLAIR).
# Passed in config as "image_ids"; legacy callers may still use "t1_ids".
image_ids       = config.get("image_ids", config.get("t1_ids", [config.get("t1_id", "")]))
if isinstance(image_ids, str):
    image_ids = [image_ids]

t1_files_sing   = config.get("t1_files_singularity", [config.get("t1_file_singularity", "")])
if isinstance(t1_files_sing, str):
    t1_files_sing = [t1_files_sing]

flair_files_sing = config.get("flair_files_singularity", [config.get("flair_file_singularity", "")])
if isinstance(flair_files_sing, str):
    flair_files_sing = [flair_files_sing]

mask_files      = config.get("mask_files", [config.get("mask_file", "")])
if isinstance(mask_files, str):
    mask_files = [mask_files]

mask_files_sing = config.get("mask_files_singularity", [config.get("mask_file_singularity", "")])
if isinstance(mask_files_sing, str):
    mask_files_sing = [mask_files_sing]

# Normalise list lengths
n = len(subjects)
def _pad(lst, n): return lst + [""] * (n - len(lst))
t1_files       = _pad(t1_files,       n)
flair_files    = _pad(flair_files,    n)
t1_files_sing  = _pad(t1_files_sing,  n)
flair_files_sing = _pad(flair_files_sing, n)
mask_files     = _pad(mask_files,     n)
mask_files_sing = _pad(mask_files_sing, n)
image_ids      = _pad(image_ids,      n)

# Per-subject lookup dicts
t1_map         = dict(zip(subjects, t1_files))
flair_map      = dict(zip(subjects, flair_files))
image_id_map   = dict(zip(subjects, image_ids))
t1_sing_map    = dict(zip(subjects, t1_files_sing))
flair_sing_map = dict(zip(subjects, flair_files_sing))
mask_map       = dict(zip(subjects, mask_files))
mask_sing_map  = dict(zip(subjects, mask_files_sing))

OUTPUT_DIR      = config.get("output_dir")
OUTPUT_DIR_SING = config.get("output_dir_singularity", "/output")
KEEP_INTERMEDIATE = config.get("keep_intermediate", False)
GIF_DIR         = "/leukoquant/leukoquant/external/gif"
GIF_SCRIPT      = GIF_DIR + "/GIF_111125.sh"
TOOL_NAME       = "gif"

# GIF database paths inside the container.
# Both dbs are directories containing db.xml.  The GIF installation is bound
# at /GIF, so db_mideface/ and db_FLAIR_mideface/ are reachable without extra bind mounts.
# T1 processing uses db_mideface/db.xml; FLAIR-only uses db_FLAIR_mideface/db.xml.
T1_DB_SING    = config.get("t1_db_singularity",    "/GIF/db_mideface/db.xml")
FLAIR_DB_SING = config.get("flair_db_singularity", "/GIF/db_FLAIR_mideface/db.xml")

# Host-side counterpart of the above, only set by callers that resolve GIF's
# home directory without validating its contents in Python (currently
# bamos_processor.py, forwarded through bamos_workflow.smk's config_gif --
# BaMoS subjects may need GIF run for real, unlike ones with externally
# pre-computed --gif-results-dir results). When present, used to declare a
# real Snakemake input on rule run_gif so a missing database surfaces as a
# normal missing-input error instead of failing deep inside the container.
# Empty for callers (gif_processor.py, process_all_processor.py) that
# already validate the database in Python before Snakemake ever runs.
GIF_HOME_HOST = config.get("gif_home_host", "")
T1_DB_HOST    = f"{GIF_HOME_HOST}/db_mideface/db.xml"       if GIF_HOME_HOST else ""
FLAIR_DB_HOST = f"{GIF_HOME_HOST}/db_FLAIR_mideface/db.xml" if GIF_HOME_HOST else ""

LEUKOQUANT_PARENT_DIR = config.get("leukoquant_parent_dir")

# Per-subject output directories: callers (e.g. bamos_workflow) may override these
# to point at pre-existing GIF results directories so GIF does not re-run.
gif_output_dirs = config.get("gif_output_dirs", [])
if isinstance(gif_output_dirs, str):
    gif_output_dirs = [gif_output_dirs]
# Pad to match subjects list; subjects without an override use the standard path.
while len(gif_output_dirs) < len(subjects):
    gif_output_dirs.append("")

gif_output_dir_map = dict(zip(subjects, gif_output_dirs))

def _subject_gif_output_dir(subject):
    """Return the GIF outputs directory for a subject.

    When an external pre-existing GIF directory was provided (e.g. via
    gif_results_paths in process-bamos), use that directly so the GIF module
    targets the already-existing parcellation files and skips re-running GIF.
    """
    override = gif_output_dir_map.get(subject, "")
    if override:
        return override
    return f"{OUTPUT_DIR}/{subject}/{TOOL_NAME}/outputs"

sys.path.insert(0, LEUKOQUANT_PARENT_DIR)
from leukoquant.utils.container_utils import ensure_container

CONTAINER_SIF = os.path.join(
    LEUKOQUANT_PARENT_DIR,
    "leukoquant/workflow/containers/miniconda_unified_container.sif",
)
ensure_container(CONTAINER_SIF)

container:
    CONTAINER_SIF

# Reverse lookups (for resolving subject from image_id wildcard)
image_id_to_subject      = {iid: subj for subj, iid in image_id_map.items()}
image_id_to_t1_file      = {iid: t1_map[subj]      for subj, iid in image_id_map.items()}
image_id_to_t1_file_sing = {iid: t1_sing_map[subj]  for subj, iid in image_id_map.items()}
image_id_to_mask         = {iid: mask_map[subj]     for subj, iid in image_id_map.items()}
image_id_to_mask_sing    = {iid: mask_sing_map[subj] for subj, iid in image_id_map.items()}

wildcard_constraints:
    subject="|".join(re.escape(s) for s in subjects),
    image_id="|".join(re.escape(iid) for iid in image_ids),


target_files = [
    f"{_subject_gif_output_dir(s)}/{iid}_NeuroMorph_Parcellation.nii.gz"
    for s, iid in zip(subjects, image_ids)
] + [
    f"{_subject_gif_output_dir(s)}/{iid}_NeuroMorph_Brain.nii.gz"
    for s, iid in zip(subjects, image_ids)
]

rule all:
    input:
        target_files

def _resolve_subject(wildcards):
    """Return the subject for a given set of wildcards (resolves via image_id when needed)."""
    if hasattr(wildcards, "subject") and wildcards.subject:
        return wildcards.subject
    return image_id_to_subject[wildcards.image_id]

rule run_gif:
    input:
        # Collect only the files that are actually provided (non-empty paths).
        t1=lambda wildcards: (
            [t1_map[_resolve_subject(wildcards)]]
            if t1_map.get(_resolve_subject(wildcards))
            else []
        ),
        flair=lambda wildcards: (
            [flair_map[_resolve_subject(wildcards)]]
            if flair_map.get(_resolve_subject(wildcards))
            else []
        ),
        mask=lambda wildcards: (
            [mask_map[_resolve_subject(wildcards)]]
            if mask_map.get(_resolve_subject(wildcards))
            else []
        ),
        # Only a real dependency when a caller passed gif_home_host (see
        # GIF_HOME_HOST above) - empty otherwise, so callers that already
        # validate the database themselves in Python are unaffected. When
        # present, picks the same T1-vs-FLAIR database the shell command
        # below selects at runtime.
        gif_db=lambda wildcards: (
            [T1_DB_HOST if t1_map.get(_resolve_subject(wildcards)) else FLAIR_DB_HOST]
            if GIF_HOME_HOST
            else []
        ),
    output:
        parcellation=f"{OUTPUT_DIR}/{{subject}}/{TOOL_NAME}/outputs/{{image_id}}_NeuroMorph_Parcellation.nii.gz",
        brain=f"{OUTPUT_DIR}/{{subject}}/{TOOL_NAME}/outputs/{{image_id}}_NeuroMorph_Brain.nii.gz"
    params:
        subject=lambda wildcards: _resolve_subject(wildcards),
        t1_sing=lambda wildcards: t1_sing_map.get(_resolve_subject(wildcards), ""),
        flair_sing=lambda wildcards: flair_sing_map.get(_resolve_subject(wildcards), ""),
        mask_sing=lambda wildcards: mask_sing_map.get(_resolve_subject(wildcards), ""),
        gif_script=GIF_SCRIPT,
        # Database paths: T1-mode uses the T1 db; FLAIR-only mode uses the FLAIR db.
        t1_db_sing=T1_DB_SING,
        flair_db_sing=FLAIR_DB_SING,
        # Singularity path to the per-subject tool directory (outputs/ is a subdirectory)
        tool_dir_sing=lambda wildcards: f"{OUTPUT_DIR_SING}/{_resolve_subject(wildcards)}/{TOOL_NAME}",
        keep_intermediate=KEEP_INTERMEDIATE,
        log_file=lambda wildcards: f"{OUTPUT_DIR_SING}/{_resolve_subject(wildcards)}/{TOOL_NAME}/logs/gif.log",
        error_file=lambda wildcards: f"{OUTPUT_DIR_SING}/{_resolve_subject(wildcards)}/{TOOL_NAME}/logs/gif_error.log",
        scratch_dir_sing=config.get("scratch_dir_singularity", "/scratch0"),
    threads:
        1
    resources:
        mem_mb=12*1024,
        time="168:00:00",
        scratch_size=0.5*1024,
        name="gif_processing",
        workdir=lambda wildcards: f"{OUTPUT_DIR}/{_resolve_subject(wildcards)}/{TOOL_NAME}",
    shell:
        """
        source /leukoquant/leukoquant/utils/bash_utils.sh
        mkdir -p "$(dirname "{params.log_file}")"
        exec > "{params.log_file}"
        exec 2> "{params.error_file}"

        echo "Date: $(date)"
        echo "Date: $(date)" >&2

        enable_line_buffering

        # Prevent output buffering to avoid null-byte corruption in log files
        # caused by pre-allocated file buffers (Singularity/NFS interaction)
        export PYTHONUNBUFFERED=1
        export RUST_LOG=info

        MY_JOB_ID="$(get_job_id)"
        scratch_intermediate=""
        trap 'scratch_cleanup "$scratch_intermediate"' EXIT INT TERM

        export OMP_NUM_THREADS={threads}

        out_dir="{params.tool_dir_sing}/outputs"
        output_intermediate_dir="{params.tool_dir_sing}/intermediate"
        scratch_intermediate="{params.scratch_dir_sing}/$USER/$MY_JOB_ID/gif_{params.subject}"
        mkdir -p "$out_dir" "$scratch_intermediate"

        # Reorient the primary input (T1 if available, else FLAIR) to RAS+.
        # First collapse to a single frame (no-op if already single-volume --
        # see image_utils.py:select_first_frame) -- some raw acquisitions are
        # exported as 4D with a duplicate second frame instead of a normal 3D
        # volume, which GIF cannot handle.
        primary_singleframe="$scratch_intermediate/{wildcards.image_id}_singleframe.nii.gz"
        primary_ras="$scratch_intermediate/{wildcards.image_id}.nii.gz"
        if [ -n "{params.t1_sing}" ]; then
            echo "Selecting first frame of T1: {params.t1_sing} -> $primary_singleframe"
            python /leukoquant/leukoquant/utils/image_utils.py \
                --input  "{params.t1_sing}" \
                --output "$primary_singleframe" \
                --op first_frame \
                --verbose
        else
            echo "Selecting first frame of FLAIR: {params.flair_sing} -> $primary_singleframe"
            python /leukoquant/leukoquant/utils/image_utils.py \
                --input  "{params.flair_sing}" \
                --output "$primary_singleframe" \
                --op first_frame \
                --verbose
        fi
        echo "Reorienting to RAS+: $primary_singleframe -> $primary_ras"
        python /leukoquant/leukoquant/utils/image_utils.py \
            --input  "$primary_singleframe" \
            --output "$primary_ras" \
            --verbose

        # Determine primary input: T1 preferred; fall back to FLAIR.
        # Use the RAS-reoriented copy produced above.
        primary_img="$primary_ras"

        mask_arg=""
        if [ -n "{params.mask_sing}" ]; then
            mask_arg="--mask {params.mask_sing}"
        fi

        # Select the GIF database: FLAIR-only runs use the FLAIR db;
        # T1 or T1+FLAIR runs use the T1 db.
        if [ -z "{params.t1_sing}" ]; then
            db_arg="--db {params.flair_db_sing}"
        else
            db_arg="--db {params.t1_db_sing}"
        fi

        # NOTE: when both T1 and FLAIR are provided, only T1 is passed to GIF.
        # FLAIR is declared as a Snakemake input (and bound into the container)
        # so it is available at {params.flair_sing} for future use - e.g. to
        # concatenate T1+FLAIR into a 4D input once GIF supports multi-channel.
        bash {params.gif_script} \
            --filename "$primary_img" \
            --output-folder "$scratch_intermediate" \
            --threads {threads} \
            $db_arg \
            $mask_arg

        # Move files from scratch intermediate to outputs (files needed by downstream tools like BaMoS)
        for pattern in "*_NeuroMorph_Parcellation.nii.gz" "*_NeuroMorph_Segmentation.nii.gz" "*_NeuroMorph_prior.nii.gz" "*TIV*" "*_NeuroMorph_Brain.nii.gz"; do
            file=$(find "$scratch_intermediate" -name "$pattern" | head -n 1)
            if [ -n "$file" ]; then
                mv "$file" "$out_dir/"
            fi
        done

        # Copy intermediate to output if requested
        if [ "{params.keep_intermediate}" = "True" ]; then
            mkdir -p "$output_intermediate_dir"
            cp -r "$scratch_intermediate"/. "$output_intermediate_dir/" 2>/dev/null || true
        fi

        """