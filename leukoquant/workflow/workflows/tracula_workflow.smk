"""
Snakemake workflow to run TRACULA (trac-all) for one or more subjects.

If recon-all outputs (aparc+aseg.mgz) are absent for any subject, the
recon_all_workflow module is activated and recon-all is run automatically
before TRACULA - mirroring the pattern used by bamos_workflow for GIF.

The caller (tracula_processor.py) determines per-subject whether recon-all
is needed and passes the resolved paths via config.

Output structure: {output_dir}/{subject}/tracula/outputs/dpath/...
                  {output_dir}/{subject}/tracula/outputs/dmri/...
                  {output_dir}/{subject}/tracula/intermediate/dmrirc
                  {output_dir}/{subject}/tracula/logs/
"""


import os
import re
import sys
from pathlib import Path

# ── subject / file lists ──────────────────────────────────────────────────────
subjects = config.get("subjects", [config.get("subject", "")])
if isinstance(subjects, str):
    subjects = [subjects]

# DWI / bvec / bval are List[List[str]] - one inner list per subject, each
# inner list holding one or more scan paths for that subject's session.
dwi_imgs      = config.get("dwi_imgs",              [])
dwi_imgs_sing = config.get("dwi_paths_singularity", [])
bvecs_paths      = config.get("bvecs_paths",             [])
bvecs_paths_sing = config.get("bvecs_paths_singularity", [])
bvals_paths      = config.get("bvals_paths",             [])
bvals_paths_sing = config.get("bvals_paths_singularity", [])

dwi_map        = dict(zip(subjects, dwi_imgs))       # {subject: [path, ...]}
dwi_sing_map   = dict(zip(subjects, dwi_imgs_sing))
bvecs_map      = dict(zip(subjects, bvecs_paths))
bvecs_sing_map = dict(zip(subjects, bvecs_paths_sing))
bvals_map      = dict(zip(subjects, bvals_paths))
bvals_sing_map = dict(zip(subjects, bvals_paths_sing))

# ── recon-all dependency resolution ──────────────────────────────────────────
# recon_dirs: host path to the recon-all outputs/ dir per subject.
#   If the caller found existing outputs → the supplied (or default) dir.
#   If recon-all must be run → the standard output path where it will write.
# recon_dirs_sing: corresponding container-side paths (always /output/…).
# needs_recon: set of subjects for which recon-all must run.
recon_dirs      = config.get("recon_dirs",             {})   # {subject: host_path}
recon_dirs_sing = config.get("recon_dirs_singularity", {})   # {subject: sing_path}
needs_recon     = set(config.get("needs_recon", []))          # subjects needing recon-all

# T1 inputs for recon-all (only populated when needs_recon is non-empty)
t1_files      = config.get("t1_files",             [])
t1_files_sing = config.get("t1_files_singularity", [])

# ── shared config ─────────────────────────────────────────────────────────────
OUTPUT_DIR       = config.get("output_dir")
OUTPUT_DIR_SING  = config.get("output_dir_singularity", "/output")
SCRATCH_DIR_SING = config.get("scratch_dir_singularity", "/scratch0")
KEEP_INTERMEDIATE = config.get("keep_intermediate", False)

LEUKOQUANT_PARENT_DIR = config.get("leukoquant_parent_dir")
MAPPING_FILE     = config.get("mapping_file", "")
# Accept either a list (new) or a single string (legacy) for backwards compat.
PARCELLATIONS    = config.get("parcellations", [config.get("parcellation", "freesurfer")])
if isinstance(PARCELLATIONS, str):
    PARCELLATIONS = [p.strip() for p in PARCELLATIONS.split(",") if p.strip()]
BRAIN_PARCELLATION = config.get("brain_parcellation", "")

# Optional dependency hook when running in master workflow
FS_DONE_DEP = config.get("fs_done_file", [])

CONTAINER_NAME = config.get("container_name", "freesurfer_unified_container")
CONTAINER_SIF  = os.path.join(
    LEUKOQUANT_PARENT_DIR, "leukoquant/workflow/containers/", CONTAINER_NAME + ".sif")

sys.path.insert(0, LEUKOQUANT_PARENT_DIR)
from leukoquant.utils.container_utils import ensure_container
ensure_container(CONTAINER_SIF)

# ── recon-all module (activated only when at least one subject needs it) ──────
RUN_RECON = len(needs_recon) > 0

if RUN_RECON:
    # Build per-subject recon-all config exactly as freesurfer_processor would.
    recon_subjects     = [s for s in subjects if s in needs_recon]
    recon_t1_map       = dict(zip(subjects, t1_files))
    recon_t1_sing_map  = dict(zip(subjects, t1_files_sing))

    config_recon = {
        "subject_id":            recon_subjects,
        "t1_file":               [recon_t1_map[s]      for s in recon_subjects],
        "t1_file_singularity":   [recon_t1_sing_map[s] for s in recon_subjects],
        # recon-all writes into OUTPUT_DIR/{subject}/recon-all/outputs
        "output_dir":            OUTPUT_DIR,
        "output_dir_singularity": OUTPUT_DIR_SING,
        "keep_intermediate":     KEEP_INTERMEDIATE,
        "container_name":        CONTAINER_NAME,
        "leukoquant_parent_dir": LEUKOQUANT_PARENT_DIR,
    }
    module recon_all:
        snakefile: "recon_all_workflow.smk"
        config: config_recon
    use rule * from recon_all as recon_all_*

# ── per-subject parcellation conversion maps ──────────────────────────────────
# Host and container paths where each subject's converted parcellation lands.
CONVERTED_PARCELLATION_MAP = {
    s: f"{OUTPUT_DIR}/{s}/atlas_conversion/outputs/converted_atlas.mgz"
    for s in subjects
}
CONVERTED_PARCELLATION_MAP_SING = {
    s: f"{OUTPUT_DIR_SING}/{s}/atlas_conversion/outputs/converted_atlas.mgz"
    for s in subjects
}

# GIF brain mask maps (only for GIF parcellation)
BRAIN_MASK_MAP      = config.get("brain_mask_map",      {})  # {subject: host_path}
BRAIN_MASK_SING_MAP = config.get("brain_mask_sing_map", {})  # {subject: container_path}

# ── input parcellation maps (per-subject) ─────────────────────────────────────
# Two modes:
#   1. Multi-subject (called from process_all): brain_parcellation_map provides
#      explicit per-subject dicts; no single-subject restriction.
#   2. Single-subject standalone (process-tracula --brain-parcellation): uses the
#      legacy single-value brain_parcellation key; single-subject enforced below.
_bp_map_cfg  = config.get("brain_parcellation_map",  {})
_bp_sing_cfg = config.get("brain_parcellation_sing_map", {})

if _bp_map_cfg:
    BRAIN_PARCELLATION_MAP      = _bp_map_cfg
    BRAIN_PARCELLATION_SING_MAP = _bp_sing_cfg
else:
    # Standalone: bind the single supplied file at a fixed container path.
    bp_suffix = "".join(Path(BRAIN_PARCELLATION).suffixes) if BRAIN_PARCELLATION else ""
    BRAIN_PARCELLATION_MAP = (
        {s: BRAIN_PARCELLATION for s in subjects} if BRAIN_PARCELLATION else {}
    )
    BRAIN_PARCELLATION_SING_MAP = (
        {s: f"/input/input_parcellation{bp_suffix}" for s in subjects}
        if BRAIN_PARCELLATION else {}
    )

# Derive the container-side path for the mapping CSV.
# Files under leukoquant_parent_dir are reachable at /leukoquant/<rel_path>;
# files explicitly bound at /input/parcellation_mapping* are also accepted.
_mf_sing_cfg = config.get("mapping_file_singularity", "")
if not _mf_sing_cfg and MAPPING_FILE:
    _mf_path = Path(MAPPING_FILE)
    _lq_path = Path(LEUKOQUANT_PARENT_DIR)
    try:
        _rel = _mf_path.relative_to(_lq_path)
        _mf_sing_cfg = f"/leukoquant/{_rel}"
    except ValueError:
        _mf_sing_cfg = "/input/parcellation_mapping" + "".join(_mf_path.suffixes)
MAPPING_FILE_SING = _mf_sing_cfg

_non_fs_parcellations = [p for p in PARCELLATIONS if p != "freesurfer"]
if _non_fs_parcellations:
    # Standalone single-path brain_parcellation cannot serve multiple subjects.
    # When using process-all or the processor builds the map, _bp_map_cfg is non-empty.
    if not _bp_map_cfg and len(subjects) > 1:
        raise ValueError(
            "Custom parcellation via --brain-parcellation supports single-subject runs only. "
            "For multi-subject parcellation, use process-all or ensure brain_parcellation_map is provided."
        )

    # atlas_converter_workflow.smk's multi-subject wildcard mode mirrors the
    # pattern used by gif_workflow.smk and bamos_workflow.smk: pass subjects
    # and per-subject input/output maps; the workflow handles {subject} wildcards
    # internally.  output_subdir places converted_atlas.mgz exactly where
    # CONVERTED_PARCELLATION_MAP expects it.
    _atlas_config = {
        "subjects":                  subjects,
        "input_parcellation_map":    BRAIN_PARCELLATION_MAP,
        "input_singularity_map":     BRAIN_PARCELLATION_SING_MAP,
        "mapping_file":              MAPPING_FILE,
        "mapping_file_singularity":  MAPPING_FILE_SING,
        "output_dir":                OUTPUT_DIR,
        "output_dir_singularity":    OUTPUT_DIR_SING,
        # Output: {OUTPUT_DIR}/{subject}/atlas_conversion/outputs/converted_atlas.mgz
        # Consistent with standalone use; matches CONVERTED_PARCELLATION_MAP exactly.
        "output_subdir":             "atlas_conversion/outputs",
        "leukoquant_parent_dir":     LEUKOQUANT_PARENT_DIR,
        "container_name":            CONTAINER_NAME,
    }

    # Use only the primary non-freesurfer parcellation's mapping for atlas conversion.
    # Each atlas converter run produces converted_atlas.mgz which the non-fs TRACULA
    # run reads; multiple non-fs parcellations would each need a separate converter run.
    module atlas_converter:
        snakefile: "atlas_converter_workflow.smk"
        config: _atlas_config

    use rule * from atlas_converter as other_*

wildcard_constraints:
    subject="|".join(re.escape(s) for s in subjects),

# ── helpers ───────────────────────────────────────────────────────────────────
def _recon_aparc_aseg(subject):
    """Return the host path to aparc+aseg.mgz for this subject."""
    if subject in needs_recon:
        return f"{OUTPUT_DIR}/{subject}/recon-all/outputs/mri/aparc+aseg.mgz"
    return f"{recon_dirs[subject]}/mri/aparc+aseg.mgz"

def _recon_thalamic_nuclei(subject):
    """Return the host path to ThalamicNuclei.v13.T1.FSvoxelSpace.mgz for this subject."""
    if subject in needs_recon:
        return f"{OUTPUT_DIR}/{subject}/recon-all/outputs/mri/ThalamicNuclei.v13.T1.FSvoxelSpace.mgz"
    return f"{recon_dirs[subject]}/mri/ThalamicNuclei.v13.T1.FSvoxelSpace.mgz"

def _recon_rh_pial(subject):
    """Return the host path to rh.pial for this subject.

    rh.pial survives both pruned outputs (where it becomes a plain file via
    cp -Lf in thalamic_seg) and non-pruned outputs (where it is a symlink to
    rh.pial.T1).  Declaring it as a tracula input lets Snakemake re-run
    recon-all automatically when the pial stage never completed.
    """
    if subject in needs_recon:
        return f"{OUTPUT_DIR}/{subject}/recon-all/outputs/surf/rh.pial"
    return f"{recon_dirs[subject]}/surf/rh.pial"

def _recon_lh_pial(subject):
    """Return the host path to lh.pial for this subject (see _recon_rh_pial)."""
    if subject in needs_recon:
        return f"{OUTPUT_DIR}/{subject}/recon-all/outputs/surf/lh.pial"
    return f"{recon_dirs[subject]}/surf/lh.pial"

def _recon_dir_sing(subject):
    """Container-side path to the recon-all outputs/ dir for this subject."""
    return recon_dirs_sing.get(subject, f"{OUTPUT_DIR_SING}/{subject}/recon-all/outputs")

# ── rule all ──────────────────────────────────────────────────────────────────
rule all:
    input:
        [f"{OUTPUT_DIR}/{s}/tracula-{parc}/outputs/dpath/merged_avg16_syn_bbr.mgz"
         for s in subjects for parc in PARCELLATIONS],
        [f"{OUTPUT_DIR}/{s}/tracula-{parc}/outputs/dmri/nodif_brain_mask.nii.gz"
         for s in subjects for parc in PARCELLATIONS],
        [CONVERTED_PARCELLATION_MAP[s] for s in subjects] if _non_fs_parcellations else [],

# ── trac_all rules - one per parcellation ────────────────────────────────────
# Each rule has only {subject} as a wildcard so that SGE executor submits a
# separate N-subject array per parcellation.  Downstream per-subject rules
# (dti, noddi) can then depend on the primary-parcellation array with
# -hold_jid_ad without hitting SGE's requirement that both arrays have
# identical task counts (which would fail when N_parcellations > 1).
for _parc in PARCELLATIONS:
    rule:
        name: f"trac_all_{_parc}"
        input:
            dwi_input=lambda wildcards: dwi_map[wildcards.subject],
            container=CONTAINER_SIF,
            # Wire recon-all outputs as explicit inputs so Snakemake builds the
            # DAG correctly at --immediate-submit / DAG-build time.  Both
            # aparc+aseg and thalamic nuclei are consumed by tracula's dmrirc.
            aparc_aseg=lambda wildcards: _recon_aparc_aseg(wildcards.subject),
            thalamic_nuclei=lambda wildcards: _recon_thalamic_nuclei(wildcards.subject),
            rh_pial=lambda wildcards: _recon_rh_pial(wildcards.subject),
            lh_pial=lambda wildcards: _recon_lh_pial(wildcards.subject),
            # Atlas conversion output is only required for non-freesurfer parcellations.
            converted_parcellation=lambda wildcards, p=_parc: (
                [CONVERTED_PARCELLATION_MAP[wildcards.subject]] if p != "freesurfer" else []
            ),
            # GIF brain mask is only required for GIF parcellation.
            brain_mask=lambda wildcards, p=_parc: (
                [BRAIN_MASK_MAP[wildcards.subject]] if p != "freesurfer" and BRAIN_MASK_MAP else []
            ),
            fs_done=FS_DONE_DEP,
        output:
            tract_merged=f"{OUTPUT_DIR}/{{subject}}/tracula-{_parc}/outputs/dpath/merged_avg16_syn_bbr.mgz",
            brain_mask_dwi=f"{OUTPUT_DIR}/{{subject}}/tracula-{_parc}/outputs/dmri/nodif_brain_mask.nii.gz",
            dwi_preproc=f"{OUTPUT_DIR}/{{subject}}/tracula-{_parc}/outputs/dmri/dwi.nii.gz",
            bvecs_preproc=f"{OUTPUT_DIR}/{{subject}}/tracula-{_parc}/outputs/dmri/dwi.bvecs",
            bvals_preproc=f"{OUTPUT_DIR}/{{subject}}/tracula-{_parc}/outputs/dmri/dwi.bvals",
        params:
            subject=lambda wildcards: wildcards.subject,
            # parcellation is a fixed string per rule - not a wildcard.
            parcellation=_parc,
            bvecs_sing=lambda wildcards: bvecs_sing_map.get(wildcards.subject, ""),
            bvals_sing=lambda wildcards: bvals_sing_map.get(wildcards.subject, ""),
            dwi_sing=lambda wildcards: dwi_sing_map[wildcards.subject],
            converted_parcellation_sing=lambda wildcards: CONVERTED_PARCELLATION_MAP_SING[wildcards.subject],
            brain_mask_sing=lambda wildcards: BRAIN_MASK_SING_MAP.get(wildcards.subject, ""),
            freesurfer_dir_singularity=lambda wildcards: _recon_dir_sing(wildcards.subject),
            tool_dir_sing=lambda wildcards, p=_parc: f"{OUTPUT_DIR_SING}/{wildcards.subject}/tracula-{p}",
            keep_intermediate=KEEP_INTERMEDIATE,
            log_file=lambda wildcards, p=_parc: f"{OUTPUT_DIR_SING}/{wildcards.subject}/tracula-{p}/logs/tracula_all.log",
            error_file=lambda wildcards, p=_parc: f"{OUTPUT_DIR_SING}/{wildcards.subject}/tracula-{p}/logs/tracula_all_error.log",
        container:
            CONTAINER_SIF,
        resources:
            mem_mb=20*1024,
            time="168:00:00",
            scratch_size=3*1024,
            name=f"tracula_{_parc}",
            workdir=lambda wildcards, p=_parc: f"{OUTPUT_DIR}/{wildcards.subject}/tracula-{p}",
        shell:
            """
            source /leukoquant/leukoquant/utils/bash_utils.sh
            log_dir="$(dirname "{params.log_file}")"
            mkdir -p "$log_dir"
            if [ ! -w "$log_dir" ]; then
                echo "FATAL: tracula log directory not writable: $log_dir" >&2
                exit 1
            fi
            exec > "{params.log_file}"
            exec 2> "{params.error_file}"

            echo "Date: $(date)"
            echo "Date: $(date)" >&2

            enable_line_buffering
    
            MY_JOB_ID="$(get_job_id)"
            scratch_path=""

            trap 'scratch_cleanup "$scratch_path"' EXIT INT TERM
    
            export FS_LICENSE=/license/license.txt
            unset SGE_ROOT
    
            # Add FSL to PATH if not already available
            if ! command -v bet &>/dev/null; then
                for _fsl in /usr/local/fsl /usr/share/fsl/5.0 /usr/share/fsl; do
                    if [ -f "$_fsl/bin/bet" ]; then
                        export FSLDIR="$_fsl"
                        source "$FSLDIR/etc/fslconf/fsl.sh" 2>/dev/null || true
                        export PATH="$FSLDIR/bin:$PATH"
                        break
                    fi
                done
            fi
    
            set -euo pipefail
    
            MY_JOB_ID="$(get_job_id)"
    
            echo "Detected Job ID: $MY_JOB_ID"
    
            scratch_path="/scratch0/$USER/$MY_JOB_ID/tracula_{params.subject}_{params.parcellation}"
            mkdir -p $scratch_path
            echo "Using scratch directory: $scratch_path"
    
            # Point TMPDIR and FS_TMPDIR at a directory inside per-job scratch so
            # tracula's transient files don't bloat the NFS output area. Without
            # this, tracula's trac-preproc invokes fs_temp_file with SGE's default
            # /tmp/<jobid>.<queue> - not bound into the container - fs_temp_file
            # fails, the error message propagates into $slistfile in trac-preproc,
            # tcsh flags the spaces as "Ambiguous", and dmri_train is invoked with
            # --slist missing its file argument. Both TMPDIR and FS_TMPDIR are
            # set: fs_temp_file's help documents FS_TMPDIR as the override but
            # empirically TMPDIR is also honoured; setting both is belt-and-braces.
            export TMPDIR="$scratch_path/tmp"
            export FS_TMPDIR="$TMPDIR"
            mkdir -p "$TMPDIR"
            # Verify the directory really is writable from the same shell context
            # that trac-preproc will inherit - if a previous run failed here, we
            # want a precise diagnostic instead of fs_temp_file's generic message.
            if ! ( touch "$TMPDIR/.tracula_tmp_canary" && rm -f "$TMPDIR/.tracula_tmp_canary" ); then
                echo "FATAL: tracula TMPDIR not writable: $TMPDIR" >&2
                ls -ld "$TMPDIR" "$scratch_path" /scratch0 >&2 || true
                exit 1
            fi
    
            tool_dir="{params.tool_dir_sing}"
            outputs_dir="$tool_dir/outputs"
            intermediate_dir="$tool_dir/intermediate"
            mkdir -p "$outputs_dir" "$intermediate_dir"
    
            # Build bash arrays from the (possibly space-separated) Snakemake params.
            # When a subject has N DWI scans, {params.dwi_sing} expands to N
            # space-separated container paths; the arrays capture all of them.
            echo "DEBUG: params.dwi_sing = '{params.dwi_sing}'"
            dwi_list=( {params.dwi_sing} )
            bvecs_list=( {params.bvecs_sing} )
            bvals_list=( {params.bvals_sing} )
            echo "DEBUG: dwi_list has ${{#dwi_list[@]}} elements"
            if [ ${{#dwi_list[@]}} -eq 0 ]; then
                echo "ERROR: dwi_list is empty! The bind mounts may not have been created correctly." >&2
                exit 1
            fi

            # Prepare every DWI scan into its own subdirectory, then copy the
            # standardised outputs (data.nii.gz / bvecs / bvals) into a shared
            # flat directory with indexed names.  TRACULA requires all scans to
            # share a single dcmroot, so a flat layout with indexed filenames is
            # the correct approach.
            dwi_output_directory="$scratch_path/dwi_extracted"
            mkdir -p "$dwi_output_directory"

            dwi_prepared=()
            bvecs_prepared=()
            bvals_prepared=()

            for i in "${{!dwi_list[@]}}"; do
                scan_outdir="$dwi_output_directory/scan_${{i}}"
                mkdir -p "$scan_outdir"

                dwi_file="${{dwi_list[$i]}}"
                echo "Preparing DWI scan ${{i}}: $dwi_file"

                # Check if file exists before processing
                if [ ! -f "$dwi_file" ]; then
                    echo "ERROR: DWI file not found at: $dwi_file" >&2
                    echo "" >&2
                    echo "Bind mount diagnostics:" >&2
                    echo "  Looking for /input_dwi_* directories..." >&2
                    if ls -d /input_dwi_* >/dev/null 2>&1; then
                        for dir in /input_dwi_*; do
                            echo "  Found: $dir" >&2
                            ls -lh "$dir" >&2 2>/dev/null || true
                        done
                    else
                        echo "  ERROR: No /input_dwi_* directories found!" >&2
                        echo "  This suggests the bind mounts were not created by Apptainer." >&2
                        echo "  Check that --apptainer-args are being passed correctly to Snakemake." >&2
                    fi
                    echo "" >&2
                    echo "Mount points in container:" >&2
                    mount | grep -E 'input_|output|/scratch' || echo "  (no relevant mounts found)"  >&2
                    exit 1
                fi

                bvecs_arg=""
                if [ -n "${{bvecs_list[$i]:-}}" ]; then
                    bvecs_arg="--bvecs ${{bvecs_list[$i]}}"
                fi
                bvals_arg=""
                if [ -n "${{bvals_list[$i]:-}}" ]; then
                    bvals_arg="--bvals ${{bvals_list[$i]}}"
                fi

                python /leukoquant/leukoquant/utils/dwi_utils.py \
                    --dwi "$dwi_file" $bvecs_arg $bvals_arg \
                    --outdir "$scan_outdir"

                cp "$scan_outdir/data.nii.gz" "$dwi_output_directory/data_${{i}}.nii.gz"
                cp "$scan_outdir/bvecs"       "$dwi_output_directory/bvecs_${{i}}"
                cp "$scan_outdir/bvals"       "$dwi_output_directory/bvals_${{i}}"
                # Copy JSON metadata if present; tracula_utils reads it per-scan.
                if [ -f "$scan_outdir/metadata.json" ]; then
                    cp "$scan_outdir/metadata.json" "$dwi_output_directory/metadata_${{i}}.json"
                fi

                dwi_prepared+=( "$dwi_output_directory/data_${{i}}.nii.gz" )
                bvecs_prepared+=( "$dwi_output_directory/bvecs_${{i}}" )
                bvals_prepared+=( "$dwi_output_directory/bvals_${{i}}" )
            done

            echo "Prepared ${{#dwi_prepared[@]}} DWI scan(s) for subject {params.subject}."

            dwi_output="${{dwi_prepared[*]}}"
            bvecs_output="${{bvecs_prepared[*]}}"
            bvals_output="${{bvals_prepared[*]}}"
    
            recon_dir_singularity="{params.freesurfer_dir_singularity}"
    
            # Check if Thalamic Nuclei segmentation exists (recon_dir points directly at outputs/)
            thalamic_seg="$recon_dir_singularity/mri/ThalamicNuclei.v13.T1.FSvoxelSpace.mgz"
            usethalnuc_value=0
            if [ -f "$thalamic_seg" ]; then
                usethalnuc_value=1
                echo "Thalamic Nuclei segmentation found for subject {params.subject} at $thalamic_seg."
            else
                echo "Thalamic Nuclei segmentation NOT found for subject {params.subject} at $thalamic_seg."
            fi
    
            # Pre-flight: verify pial surfaces are present and reachable before
            # doing any heavy DWI work. In FreeSurfer 7+, rh.pial is a symlink
            # → rh.pial.T1; in pruned recon-all outputs (after thalamic_seg the
            # pruner dereferences it with cp -Lf), rh.pial is a plain file and
            # rh.pial.T1 no longer exists. Using -f (follows symlinks) catches
            # both a genuinely missing file and a broken symlink, which indicates
            # the pial stage of recon-all never completed. trac-all -path calls
            # mris_anatomical_stats at the very end, so without this guard a
            # broken pial surface wastes 8+ hours of tractography before failing.
            for _pial in "surf/rh.pial" "surf/lh.pial"; do
                if [ ! -f "$recon_dir_singularity/$_pial" ]; then
                    echo "ERROR: pial surface missing or broken symlink: $recon_dir_singularity/$_pial" >&2
                    echo "       recon-all likely crashed before completing the pial stage." >&2
                    echo "       Re-run recon-all for this subject before running tracula." >&2
                    ls -la "$recon_dir_singularity/surf/" 2>/dev/null >&2 || true
                    exit 1
                fi
            done
            echo "Pre-flight: rh.pial and lh.pial verified for subject {params.subject}."

            recon_scratch_dir="$scratch_path/recon-all"
            echo "PARCELLATION variable is set to: {params.parcellation}"
            if [ "{params.parcellation}" = "freesurfer" ]; then
                echo "Using recon-all parcellation for subject {params.subject}."
                mkdir -p "$recon_scratch_dir/{params.subject}"
                # recon_dir_singularity points at outputs/ which contains mri/, surf/, etc.
                cp -r "$recon_dir_singularity/." "$recon_scratch_dir/{params.subject}/"
            else
                echo "Using custom parcellation from: {params.converted_parcellation_sing} for subject {params.subject}."
                mkdir -p "$recon_scratch_dir/{params.subject}"
    
                files_to_copy=(
                    "surf/lh.white" "surf/rh.white" "surf/lh.thickness" "surf/rh.thickness"
                    "surf/rh.sphere.reg" "surf/lh.sphere.reg" "surf/rh.pial" "surf/lh.pial"
                    "label/lh.cortex.label" "label/rh.cortex.label"
                    "mri/orig.mgz" "mri/T1.mgz"
                    "mri/norm.mgz" "mri/wm.mgz" "mri/transforms/talairach.xfm"
                    "mri/aseg.presurf.mgz" "mri/ribbon.mgz" "mri/ThalamicNuclei.v13.T1.FSvoxelSpace.mgz"
                    "stats/brainvol.stats"
                )
    
                (
                    cd "$recon_dir_singularity" || exit
                    for pattern in "${{files_to_copy[@]}}"; do
                        for file in $pattern; do
                            if [ -e "$file" ]; then
                                # cp -L dereferences any symlink chain at copy time, writing a
                                # plain file to scratch regardless of whether the source is a
                                # symlink into the container's own overlayFS layer.  The old
                                # readlink -f approach failed when the resolved absolute path
                                # lay outside the container's bind mounts, silently skipping the
                                # file and causing late errors in mris_anatomical_stats.
                                # --parents preserves the relative path (e.g. surf/rh.pial) so
                                # no separate mkdir -p is needed.
                                cp -ruLf --parents "$file" "$recon_scratch_dir/{params.subject}/"
                                chmod u+rw "$recon_scratch_dir/{params.subject}/$file"
                            else
                                echo "Warning: Source file $file not found for {params.subject}"
                            fi
                        done
                    done
                )

                # Verify that the surface files tracula needs for tract-endpoint statistics
                # are present in scratch.  They are used late (mris_anatomical_stats), so
                # a missing file would otherwise waste 20+ hours of tractography before failing.
                for _critical in "surf/rh.pial" "surf/lh.pial" "surf/rh.white" "surf/lh.white"; do
                    if [ ! -f "$recon_scratch_dir/{params.subject}/$_critical" ]; then
                        echo "ERROR: required FreeSurfer file missing from scratch: $_critical" >&2
                        echo "       expected source: $recon_dir_singularity/$_critical" >&2
                        ls -la "$recon_dir_singularity/surf/" 2>/dev/null >&2 || true
                        exit 1
                    fi
                done
                echo "Critical FreeSurfer surface files verified in scratch."
    
                target_parcellation_path="$recon_scratch_dir/{params.subject}/mri/aparc+aseg.mgz"
                echo "Copying custom parcellation from {params.converted_parcellation_sing} to $target_parcellation_path"
                if [ -f "{params.converted_parcellation_sing}" ]; then
                    cp -rf "{params.converted_parcellation_sing}" "$target_parcellation_path"

                    # Resample GIF atlas to match T1 space dimensions first
                    # GIF atlas may have different dimensions (e.g., 176x256x256) than T1 (256x256x256)
                    # This causes tracula to fail with 0 streamlines. Resample to T1 space.
                    t1_file="$recon_scratch_dir/{params.subject}/mri/T1.mgz"
                    echo "Resampling GIF atlas to match T1 space..."
                    mri_convert "$target_parcellation_path" "$target_parcellation_path" --like "$t1_file"
                    echo "Resampled aparc+aseg.mgz to T1 space"

                    # Apply GIF brain mask to T1.mgz to produce brain.mgz and brainmask.mgz
                    # Now done after resampling aparc+aseg, so brain.mgz is created in correct space
                    echo "Applying GIF brain mask to T1.mgz..."
                    brain_file="$recon_scratch_dir/{params.subject}/mri/brain.mgz"
                    brainmask_file="$recon_scratch_dir/{params.subject}/mri/brainmask.mgz"

                    mri_mask "$t1_file" "{params.brain_mask_sing}" "$brain_file"
                    cp "$brain_file" "$brainmask_file"
                    echo "Created brain.mgz and brainmask.mgz from GIF brain mask applied to T1.mgz"
                else
                    echo "Error: Converted parcellation file {params.converted_parcellation_sing} not found. Exiting."
                    exit 1
                fi
            fi
    
            # dmrirc lives in scratch; only copied to intermediate/ if keep_intermediate=True
            dmrirc_path="$scratch_path/dmrirc"
            is_shelled_flag_file="$scratch_path/is_shelled_flag"

            echo "Writing dmrirc at $dmrirc_path (${{#dwi_prepared[@]}} scan(s))"
            python /leukoquant/leukoquant/utils/tracula_utils.py \
            --subject {params.subject} \
            --subjects-dir "$recon_scratch_dir" \
            --outdir "$scratch_path" \
            --dwi $dwi_output \
            --bvecs $bvecs_output \
            --bvals $bvals_output \
            --usethalnuc $usethalnuc_value \
            --doeddy 2 \
            --dob0 0 \
            --dmrirc-out "$dmrirc_path" \
            --is-shelled-flag-file "$is_shelled_flag_file"
    
            trac_preproc_file="/usr/local/freesurfer/bin/trac-preproc"
            bedpostx_mgh_file="/usr/local/freesurfer/bin/bedpostx_mgh"
            fsl_sub_mgh_file="/usr/local/freesurfer/bin/fsl_sub_mgh"
    
            local_bin="$scratch_path/bin"
            mkdir -p "$local_bin"
    
            if [ -f "$trac_preproc_file" ]; then
                cp "$trac_preproc_file" "$local_bin/trac-preproc"
                chmod +x "$local_bin/trac-preproc"
                sed -i 's/set cmd = eddy_openmp/set cmd = ( eddy diffusion )/g' "$local_bin/trac-preproc"
                grep "set cmd = ( eddy diffusion )" "$local_bin/trac-preproc" || echo "Modification not found in trac-preproc!"
                sed -i 's/\\$esp \\* .001/\\$esp/g' "$local_bin/trac-preproc"
                grep "\\$esp \\* .001" "$local_bin/trac-preproc" && echo "Modification found in trac-preproc!"
    
                # BIDS-style subject IDs contain a slash (e.g.
                # "sub-XXX/sub-XXX_ses-YY"). trac-preproc builds the dmri_train
                # slist file name via:
                #     fs_temp_file --suffix .$ntrain.$subj.txt
                # which makes the suffix contain that slash; mkstemp then fails
                # with "Could not create temporary file in <TMPDIR>" because
                # the intermediate directory tmp.XXXXXX.<ntrain>.sub-XXX/
                # does not exist. The $subj part is decorative - fs_temp_file
                # already injects a random 6-char unique token - so drop it.
                sed -i 's/--suffix \\.\\$ntrain\\.\\$subj\\.txt/--suffix .$ntrain.txt/g' "$local_bin/trac-preproc"
                grep -F -- '--suffix .$ntrain.txt' "$local_bin/trac-preproc" >/dev/null || echo "WARNING: trac-preproc slash-suffix patch did not apply"

                cp "$bedpostx_mgh_file" "$local_bin/bedpostx_mgh"
                chmod +x "$local_bin/bedpostx_mgh"
                sed -i 's|#!/bin/sh|#!/bin/bash|g' "$local_bin/bedpostx_mgh"
                grep "#!/bin/bash" "$local_bin/bedpostx_mgh" || echo "Modification not found in bedpostx_mgh!"
    
                cp "$fsl_sub_mgh_file" "$local_bin/fsl_sub_mgh"
                chmod +x "$local_bin/fsl_sub_mgh"

                # Patch orientLAS to fix PSL/PSR DWI handling (see utils/patch_orientlas.py).
                orient_las_src="/usr/local/freesurfer/bin/orientLAS"
                if [ -f "$orient_las_src" ]; then
                    cp "$orient_las_src" "$local_bin/orientLAS"
                    chmod +x "$local_bin/orientLAS"
                    python3 /leukoquant/leukoquant/utils/patch_orientlas.py "$local_bin/orientLAS"
                else
                    echo "Warning: orientLAS not found at $orient_las_src; PSL/PSR DWI will fail."
                fi

                export PATH="$local_bin:$PATH"
            else
                echo "Warning: trac-preproc file not found at $trac_preproc_file."
            fi

            # is_shelled_flag written by tracula_utils.py; 1=shelled, 0=not shelled.
            # When not shelled, tracula_utils already lowered doeddy to 1 (eddy_correct).
            is_shelled_flag="$(cat "$is_shelled_flag_file" 2>/dev/null || echo 1)"
            echo "Info: is_shelled_flag=$is_shelled_flag (doeddy adjusted in tracula_utils.py if 0)"

            debug_dir="$tool_dir/debug"

            echo "Running trac-all -prep"
            trac-all -c $dmrirc_path -no-isrunning -prep

            echo "Running trac-all -bedp"
            trac-all -c $dmrirc_path -no-isrunning -bedp

            echo "Running trac-all -path"
            trac-all -c $dmrirc_path -no-isrunning -path

            echo "Copying results to outputs directory..."
            # Copy only the files we actually use to save disk space and inode count
            folder_paths=(
                "dpath/merged_avg16_syn_bbr.mgz"
                "dpath/*avg16*/path.pd.trk"
                "dpath/*avg16*/path.pd.nii.gz"
                "dmri/nodif_brain_mask.nii.gz"
                "dmri/dwi.nii.gz"
                "dmri/dwi.bvals"
                "dmri/dwi.bvecs"
            )

            # Tracula outputs to $scratch_path/{params.subject}/
            tracula_subject_dir="$scratch_path/{params.subject}/"

            (
                cd "$tracula_subject_dir" || exit
                for pattern in "${{folder_paths[@]}}"; do
                    for file in $pattern; do
                        if [ -e "$file" ]; then
                            cp -ruL --parents "$file" "$outputs_dir/"
                            echo "Copied $file to outputs directory for subject {params.subject}."
                        fi
                    done
                done
            )

            echo "Tracula processing completed for subject {params.subject}."

            # Copy dmrirc to intermediate/ only if requested
            if [ "{params.keep_intermediate}" = "True" ]; then
                mkdir -p "$intermediate_dir"
                cp "$dmrirc_path" "$intermediate_dir/dmrirc"
                echo "Copied dmrirc to intermediate directory for subject {params.subject}."
            fi
            """
