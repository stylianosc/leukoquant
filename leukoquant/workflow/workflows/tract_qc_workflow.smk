"""
Snakemake workflow for tractography QC.

This workflow:
1. Reads subject IDs (single or batch)
2. Locates tractography outputs for each subject
3. Runs QC utility script for each subject
4. Generates QC report outputs

Output structure: {output_dir}/{subject}/tract_qc/outputs/qc_report.csv
                  {output_dir}/{subject}/tract_qc/logs/

Uses Singularity/Conda containers for isolated execution.
"""

import os
import re
import sys
import glob as _pyglob
from pathlib import Path
from leukoquant.utils.z_score_utils import translate_path


# ============================================================================
# Configuration & Setup
# ============================================================================
LEUKOQUANT_PARENT_DIR = config.get("leukoquant_parent_dir")
if not LEUKOQUANT_PARENT_DIR:
    raise ValueError("leukoquant_parent_dir not provided in config")

sys.path.insert(0, LEUKOQUANT_PARENT_DIR)
from leukoquant.utils.container_utils import ensure_container

# Standalone `process-tract-qc` invokes Snakemake with --configfile=<...>/tract_qc_config.yaml,
# which Snakemake's own CLI already merges into `config` before this file runs -- so `config`
# already carries every key that used to be re-read from disk here. Re-parsing a same-named
# file from the working directory on top of that was redundant for that path, and would be
# dangerous for process_all's module-based invocation the same way the equivalent pattern in
# metrics_workflow.smk was: process_all's workdir is the dataset's shared, persistent
# output_dir, so a stray tract_qc_config.yaml left behind by an unrelated standalone run would
# silently hijack every later process_all run sharing that directory.
cfg = config

SUBJECTS              = cfg.get("subjects", [])
TRACTOGRAPHY_PATH     = cfg.get("tractography_path", "")
OUTPUT_DIR            = cfg.get("output_dir", "./tract_qc_output")
VALIDATE_TRACTOGRAPHY = cfg.get("validate_tractography_glob", True)
# Accept list (new) or single string (legacy) for backwards compat.
PARCELLATIONS         = cfg.get("parcellations", [cfg.get("parcellation", "freesurfer")])
if isinstance(PARCELLATIONS, str):
    PARCELLATIONS = [p.strip() for p in PARCELLATIONS.split(",") if p.strip()]

# When True (set by process-all), outputs go to tract_qc-{parcellation}/.
# When False (standalone), outputs go to tract_qc/ with no suffix.
USE_PARCELLATION_SUFFIX = cfg.get("use_parcellation_suffix", False)

if not SUBJECTS:
    raise ValueError("No subjects provided in config")
if not TRACTOGRAPHY_PATH:
    raise ValueError("tractography_path not provided in config")

CONTAINER_SIF = os.path.join(
    LEUKOQUANT_PARENT_DIR,
    "leukoquant/workflow/containers/",
    cfg.get("container_name", "miniconda_unified_container") + ".sif"
)
ensure_container(CONTAINER_SIF)

BIND_MAP  = cfg.get("singularity_binds", {})
QC_SCRIPT = "/leukoquant/leukoquant/utils/tract_qc_util.py"

def _resolve_tractography_path(subject, parcellation=None):
    """Substitute {subject} and {parcellation} placeholders in TRACTOGRAPHY_PATH."""
    path = TRACTOGRAPHY_PATH
    if "{subject}" in path:
        path = path.replace("{subject}", subject)
    if parcellation and "{parcellation}" in path:
        path = path.replace("{parcellation}", parcellation)
    return path

def _resolve_tract_files(subject, parcellation=None):
    """Resolve tractography files on the host for a given subject at parse time."""
    _TRACTOGRAPHY_PATH = _resolve_tractography_path(subject, parcellation)
    _tract_parts = _TRACTOGRAPHY_PATH.split(":")
    _tract_base  = _tract_parts[0].rstrip("/")
    _tract_glob  = ":".join(_tract_parts[1:]).lstrip("/") if len(_tract_parts) > 1 else ""
    _subject_in_base = subject in str(Path(_tract_base).resolve())
    if _subject_in_base:
        _tract_pattern = os.path.join(_tract_base, _tract_glob) if _tract_glob else _tract_base
    else:
        _tract_pattern = os.path.join(_tract_base, subject, _tract_glob) if _tract_glob else os.path.join(_tract_base, subject)
    _tract_files = sorted(_pyglob.glob(_tract_pattern))
    if not _tract_files:
        print(
            f"[tract_qc] No tractography files found for subject {subject} "
            f"with pattern: {_tract_pattern}",
            file=sys.stderr,
        )
    return _tract_files

def _sif_tractography_path(subject, parcellation=None):
    """Return the container-mounted tractography path for a given subject."""
    _TRACTOGRAPHY_PATH = _resolve_tractography_path(subject, parcellation)
    _tract_parts = _TRACTOGRAPHY_PATH.split(":")
    _tract_base  = _tract_parts[0].rstrip("/")
    _tract_glob  = ":".join(_tract_parts[1:]).lstrip("/") if len(_tract_parts) > 1 else ""
    _sif_tract_base = translate_path(str(Path(_tract_base).resolve()), BIND_MAP)
    return f"{_sif_tract_base}:{_tract_glob}" if _tract_glob else _sif_tract_base

def _sif_output_dir(subject, parcellation=None):
    """Return the container-mounted output directory for a given subject."""
    if USE_PARCELLATION_SUFFIX:
        _parc = parcellation or PARCELLATIONS[0]
        folder = f"{OUTPUT_DIR}/{subject}/tract_qc-{_parc}/outputs"
    else:
        folder = f"{OUTPUT_DIR}/{subject}/tract_qc/outputs"
    return translate_path(str(Path(folder).resolve()), BIND_MAP)

# ============================================================================
# Rules - output folder naming depends on USE_PARCELLATION_SUFFIX:
#   True  (process-all): tract_qc-{parcellation}/  fans out per subject × parcellation
#   False (standalone):  tract_qc/                 no parcellation suffix
# ============================================================================

if USE_PARCELLATION_SUFFIX:
    # {subject: {parcellation: host_path}} - populated by process_all so the
    # for-loop rules can declare a Snakemake DAG edge to the tracula outputs.
    # Empty when running standalone (no DAG edge needed in that context).
    _tracula_map = cfg.get("tracula_tract_merged_map", {})

    wildcard_constraints:
        subject="|".join(re.escape(s) for s in SUBJECTS),

    rule all:
        input:
            [f"{OUTPUT_DIR}/{subject}/tract_qc-{parc}/outputs/qc_report.csv"
             for subject in SUBJECTS for parc in PARCELLATIONS]

    # One rule per parcellation so each SGE array has only {subject} as a wildcard.
    # This gives separate N-subject arrays per parcellation, matching the pattern
    # used by tracula and metrics and satisfying SGE's -hold_jid_ad requirement.
    for _parc in PARCELLATIONS:
        rule:
            name: f"tract_qc_{_parc}"
            input:
                tractography=lambda wildcards, p=_parc: (
                    _resolve_tract_files(wildcards.subject, p)
                    if VALIDATE_TRACTOGRAPHY else []
                ),
                # DAG edge to the tracula output for this parcellation.
                # Empty list in standalone context (no process_all dependency map).
                tract_merged=lambda wildcards, p=_parc, m=_tracula_map: (
                    [m[wildcards.subject][p]]
                    if wildcards.subject in m and p in m[wildcards.subject] else []
                ),
            output:
                qc_report=f"{OUTPUT_DIR}/{{subject}}/tract_qc-{_parc}/outputs/qc_report.csv",
            container:
                CONTAINER_SIF,
            params:
                sif_tractography_path=lambda wildcards, p=_parc: _sif_tractography_path(wildcards.subject, p),
                sif_output_dir=lambda wildcards, p=_parc: _sif_output_dir(wildcards.subject, p),
                qc_script=QC_SCRIPT,
                log_file=lambda wildcards, p=_parc: translate_path(
                    str(Path(f"{OUTPUT_DIR}/{wildcards.subject}/tract_qc-{p}/logs/tract_qc_log.txt").resolve()),
                    BIND_MAP,
                ),
                error_file=lambda wildcards, p=_parc: translate_path(
                    str(Path(f"{OUTPUT_DIR}/{wildcards.subject}/tract_qc-{p}/logs/tract_qc_error.txt").resolve()),
                    BIND_MAP,
                ),
            resources:
                mem_mb=4*1024,
                time="1:00:00",
                name=f"tractqc_{_parc}",
                workdir=lambda wildcards, p=_parc: f"{OUTPUT_DIR}/{wildcards.subject}/tract_qc-{p}",
            shell:
                """
                source /leukoquant/leukoquant/utils/bash_utils.sh
                mkdir -p "$(dirname "{params.log_file}")"
                exec > "{params.log_file}"
                exec 2> "{params.error_file}"

                echo "Date: $(date)"
                echo "Date: $(date)" >&2

                mkdir -p "{params.sif_output_dir}"
                python {params.qc_script} \
                    --subject {wildcards.subject} \
                    --tractography-path {params.sif_tractography_path} \
                    --output {params.sif_output_dir} \
                    --skip-subject-dir
                """

else:
    wildcard_constraints:
        subject="|".join(re.escape(s) for s in SUBJECTS),

    rule all:
        input:
            [f"{OUTPUT_DIR}/{subject}/tract_qc/outputs/qc_report.csv"
             for subject in SUBJECTS]

    rule tract_qc:
        input:
            tractography=lambda wildcards: (
                _resolve_tract_files(wildcards.subject)
                if VALIDATE_TRACTOGRAPHY else []
            ),
        output:
            qc_report=f"{OUTPUT_DIR}/{{subject}}/tract_qc/outputs/qc_report.csv",
        container:
            CONTAINER_SIF
        params:
            sif_tractography_path=lambda wildcards: _sif_tractography_path(wildcards.subject),
            sif_output_dir=lambda wildcards: _sif_output_dir(wildcards.subject),
            qc_script=QC_SCRIPT,
            log_file=lambda wildcards: translate_path(str(Path(f"{OUTPUT_DIR}/{wildcards.subject}/tract_qc/logs/tract_qc_log.txt").resolve()), BIND_MAP),
            error_file=lambda wildcards: translate_path(str(Path(f"{OUTPUT_DIR}/{wildcards.subject}/tract_qc/logs/tract_qc_error.txt").resolve()), BIND_MAP),
        resources:
            mem_mb=4*1024,
            time="1:00:00",
            name="tractqc",
            workdir=lambda wildcards: f"{OUTPUT_DIR}/{wildcards.subject}/tract_qc",
        shell:
            """
            source /leukoquant/leukoquant/utils/bash_utils.sh
            mkdir -p "$(dirname "{params.log_file}")"
            exec > "{params.log_file}"
            exec 2> "{params.error_file}"

            echo "Date: $(date)"
            echo "Date: $(date)" >&2
            
            mkdir -p "{params.sif_output_dir}"
            python {params.qc_script} \
                --subject {wildcards.subject} \
                --tractography-path {params.sif_tractography_path} \
                --output {params.sif_output_dir} \
                --skip-subject-dir
            """
