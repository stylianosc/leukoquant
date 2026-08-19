"""
Snakemake workflow for computing metrics along tracts, lesions, and WMH regions.

This workflow:
1. Reads subject IDs (single or batch)
2. Locates tractography outputs for each subject
3. Processes metrics (DTI, NODDI, etc.) along tracts and in ROIs
4. Generates CSV outputs with metric statistics

Output structure: {output_dir}/{subject}/metrics-{parcellation}/outputs/metrics/whole_brain_metrics.csv
                  {output_dir}/{subject}/metrics-{parcellation}/logs/

The parcellation suffix is always applied (defaults to 'freesurfer' when not specified).
Uses Singularity containers for isolated execution with bind-mounted paths.
Python helpers (path translation) live in leukoquant/utils/z_score_utils.py.
"""

import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional

# ============================================================================
# Configuration & Setup
# ============================================================================

LEUKOQUANT_PARENT_DIR = config.get("leukoquant_parent_dir")

if not LEUKOQUANT_PARENT_DIR:
    raise ValueError("leukoquant_parent_dir not provided in config")

sys.path.insert(0, LEUKOQUANT_PARENT_DIR)

from leukoquant.utils.z_score_utils import translate_path
from leukoquant.utils.container_utils import ensure_container

# Standalone `process-metrics` invokes Snakemake with --configfile=<...>/metrics_config.yaml,
# which Snakemake's own CLI already merges into `config` before this file runs -- so `config`
# already carries every key that used to be re-read from disk here. Re-parsing a same-named
# file from the working directory on top of that was redundant for that path, and actively
# dangerous for process_all's module-based invocation: process_all's workdir is the dataset's
# shared, persistent output_dir, so any metrics_config.yaml ever left behind by an unrelated
# standalone `process-metrics` run permanently and silently hijacked every later process_all
# run sharing that directory, collapsing SUBJECTS to whatever that old run had targeted.
cfg = config

# Extract configuration
SUBJECTS          = cfg.get("subjects", [])
TRACTOGRAPHY_PATH = cfg.get("tractography_path", "")
LESION_PATH       = cfg.get("lesion_path", None)
T1_PATH           = cfg.get("t1_path", None)
DWI_PATH          = cfg.get("dwi_path", None)
METRICS           = cfg.get("metrics", {})
OUTPUT_DIR        = cfg.get("output_dir", "./metrics_output")
SINGULARITY_BINDS = cfg.get("singularity_binds", {})
METRIC_MOUNTS     = cfg.get("metric_mounts", {})
VERBOSE           = cfg.get("verbose", False)
# Accept list (new) or single string (legacy) for backwards compat.
PARCELLATIONS     = cfg.get("parcellations", [cfg.get("parcellation", "freesurfer")])
if isinstance(PARCELLATIONS, str):
    PARCELLATIONS = [p.strip() for p in PARCELLATIONS.split(",") if p.strip()]

# When True, metric base paths already point to the exact directory containing
# the files (no per-subject subfolder), so the subject ID is not appended.
# This is used when process_all passes fully-qualified per-subject metric paths.
SKIP_SUBJECT_DIR_METRICS = cfg.get("skip_subject_dir_metrics", False)

if not SUBJECTS:
    raise ValueError("No subjects provided in config")

if not TRACTOGRAPHY_PATH:
    raise ValueError("tractography_path not provided in config")

# Container setup
CONTAINER_SIF = os.path.join(
    LEUKOQUANT_PARENT_DIR,
    "leukoquant/workflow/containers/",
    cfg.get("container_name", "freesurfer_unified_container") + ".sif"
)
ensure_container(CONTAINER_SIF)

BIND_MAP = cfg.get("singularity_binds", {})
METRICS_SCRIPT_SIF = "/leukoquant/leukoquant/utils/metrics_calc.sh"


def _get_lesion_names():
    """Parse LESION_PATH to extract lesion set names (the key before '=' in each CSV token).

    Lesion names are subject-independent - they come from the config key, not the path.
    E.g. "wmh=/path/{subject}/wmh.nii.gz,pvs=/path/{subject}/pvs.nii.gz" → ["wmh", "pvs"]
    """
    if not LESION_PATH:
        return []
    names = []
    for token in [t.strip() for t in LESION_PATH.split(",") if t.strip()]:
        name = token.split("=", 1)[0].strip() if "=" in token else "lesion"
        names.append(name)
    return names

LESION_NAMES = _get_lesion_names()


# ============================================================================
# Per-subject metrics extraction
# ============================================================================

T1_MAP  = cfg.get("t1_map",  {})
DWI_MAP = cfg.get("dwi_map", {})

def _subst(path_str, subject):
    if path_str and "{subject}" in str(path_str):
        return str(path_str).replace("{subject}", subject)
    return path_str

def _sif_tractography_path(subject, parcellation=None):
    _TRACTOGRAPHY_PATH = _subst(TRACTOGRAPHY_PATH, subject)
    if parcellation and "{parcellation}" in _TRACTOGRAPHY_PATH:
        _TRACTOGRAPHY_PATH = _TRACTOGRAPHY_PATH.replace("{parcellation}", parcellation)
    _tract_base   = _TRACTOGRAPHY_PATH.split(":", 1)[0]
    _tract_suffix = _TRACTOGRAPHY_PATH.split(":", 1)[1] if ":" in _TRACTOGRAPHY_PATH else ""
    _sif_tract_base = translate_path(str(Path(_tract_base).resolve()), BIND_MAP)
    return f"{_sif_tract_base}:{_tract_suffix}" if _tract_suffix else _sif_tract_base

def _sif_output_dir(subject, parcellation=None):
    _parc = parcellation or PARCELLATIONS[0]
    folder = f"{OUTPUT_DIR}/{subject}/metrics-{_parc}/outputs"
    return translate_path(str(Path(folder).resolve()), BIND_MAP)

def _sif_logs_dir(subject, parcellation=None):
    _parc = parcellation or PARCELLATIONS[0]
    folder = f"{OUTPUT_DIR}/{subject}/metrics-{_parc}/logs"
    return translate_path(str(Path(folder).resolve()), BIND_MAP)

def _sif_qc_report(subject, parcellation=None):
    """Return the container path for the tract_qc report for this subject."""
    _parc = parcellation or PARCELLATIONS[0]
    qc_folder = f"tract_qc-{_parc}"
    qc_report_host = str(Path(f"{OUTPUT_DIR}/{subject}/{qc_folder}/outputs/qc_report.csv").resolve())
    return translate_path(qc_report_host, BIND_MAP)

def _sif_t1_path(subject):
    _subject_t1_path = T1_MAP.get(subject, T1_PATH) if T1_MAP else T1_PATH
    if not _subject_t1_path:
        return None
    t1_base = _subject_t1_path.split(":")[0]
    t1_glob = ":".join(_subject_t1_path.split(":")[1:]) if ":" in _subject_t1_path else ""
    t1_subject_path = f"{t1_base}/{subject}:{t1_glob}" if t1_glob.strip(":") else t1_base
    sif = translate_path(str(Path(t1_subject_path.split(":")[0]).resolve()), BIND_MAP)
    return f"{sif}:{t1_subject_path.split(':', 1)[1]}" if ":" in t1_subject_path else sif

def _sif_dwi_path(subject):
    _subject_dwi_path = DWI_MAP.get(subject, DWI_PATH) if DWI_MAP else DWI_PATH
    if not _subject_dwi_path:
        return None
    dwi_base = _subject_dwi_path.split(":")[0]
    dwi_glob = ":".join(_subject_dwi_path.split(":")[1:]) if ":" in _subject_dwi_path else ""
    dwi_subject_path = f"{dwi_base}/{subject}:{dwi_glob}" if dwi_glob.strip(":") else dwi_base
    sif = translate_path(str(Path(dwi_subject_path.split(":")[0]).resolve()), BIND_MAP)
    return f"{sif}:{dwi_subject_path.split(':', 1)[1]}" if ":" in dwi_subject_path else sif

def _sif_lesion_path(subject):
    _LESION_PATH = _subst(LESION_PATH, subject) if LESION_PATH else LESION_PATH
    if not _LESION_PATH:
        return None
    _lesion_entries = []
    for lesion_token in [t.strip() for t in _LESION_PATH.split(",") if t.strip()]:
        if "=" in lesion_token:
            lesion_name, lesion_spec = lesion_token.split("=", 1)
            lesion_name = lesion_name.strip()
        else:
            lesion_name = "lesion"
            lesion_spec = lesion_token
        lesion_parts = lesion_spec.split(":")
        lesion_base  = lesion_parts[0]
        lesion_glob  = lesion_parts[1] if len(lesion_parts) > 1 else ""
        lesion_space = lesion_parts[2] if len(lesion_parts) > 2 else "t1"
        _subject_in_base = subject in str(Path(lesion_base).resolve())
        if (lesion_glob and subject in lesion_glob) or _subject_in_base:
            lesion_subject_base = lesion_base
        else:
            lesion_subject_base = f"{lesion_base}/{subject}"
        sif = translate_path(str(Path(lesion_subject_base).resolve()), BIND_MAP)
        spec = sif
        if lesion_glob:
            spec = f"{spec}:{lesion_glob}"
        if lesion_space:
            spec = f"{spec}:{lesion_space}"
        if lesion_name:
            spec = f"{lesion_name}={spec}"
        _lesion_entries.append(spec)
    return ",".join(_lesion_entries) if _lesion_entries else None

def _metrics_spec_str(subject):
    _METRICS = {k: _subst(v, subject) for k, v in METRICS.items()} if METRICS else {}
    sif_metrics = {}
    for metric_name, metric_spec in _METRICS.items():
        if ":" in metric_spec:
            parts       = metric_spec.split(":")
            metric_base = parts[0]
            metric_glob = parts[1] if len(parts) > 1 else ""
            _subject_in_base = subject in str(Path(metric_base).resolve())
            if SKIP_SUBJECT_DIR_METRICS or (metric_glob and subject in metric_glob) or _subject_in_base:
                metric_subject_path = metric_base
            else:
                metric_subject_path = f"{metric_base}/{subject}"
            metric_base_abs = str(Path(metric_subject_path).resolve())
            sif_base = METRIC_MOUNTS.get(metric_name) or translate_path(metric_base_abs, BIND_MAP)
            sif_metrics[metric_name] = f"{sif_base}:{':'.join(parts[1:])}" if len(parts) > 1 else sif_base
        else:
            _subject_in_base = subject in str(Path(metric_spec).resolve())
            metric_subject_path = metric_spec if _subject_in_base else f"{metric_spec}/{subject}"
            sif_metrics[metric_name] = translate_path(str(Path(metric_subject_path).resolve()), BIND_MAP)
    # z-score maps (populated by process_all via cfg["z_score_map"], see
    # z_score_outputs() in process_all_workflow.smk) -- added as extra named
    # metrics ("dti_fa_z_score" etc, derived from each file's own basename)
    # so the metrics script summarizes them into skeleton_metrics.csv /
    # tract_level_metrics.csv alongside the raw dti/noddi maps, instead of
    # them existing only as an (until now unused) DAG dependency edge.
    # metrics_calc.sh's --metrics parser requires "base:glob:space" per entry
    # (see parse_spec in metrics_calc.sh) -- a bare file path has no glob, so
    # it was silently dropped with a "no glob pattern" warning. z-score maps
    # are registered to T1 space (confirmed via freeview: same shape/affine
    # as native T1), so space="t1" here, matching how T1/lesion entries skip
    # the DWI->T1 registration step.
    for z_path in (_z_score_map.get(subject) or []):
        map_name = Path(z_path).name
        for ext in (".nii.gz", ".nii"):
            if map_name.endswith(ext):
                map_name = map_name[: -len(ext)]
                break
        z_resolved = Path(z_path).resolve()
        sif_z_dir = translate_path(str(z_resolved.parent), BIND_MAP)
        sif_metrics[map_name] = f"{sif_z_dir}:{z_resolved.name}:t1"
    return ",".join([f"{k}={v}" for k, v in sif_metrics.items()])

# ============================================================================
# Rules - always use metrics-{parcellation}/ for consistent output paths
# across both standalone and process-all invocations.
# One rule per parcellation so each SGE array has only {subject} as a wildcard.
# This ensures N-task arrays (one per subject) across all dependent rules,
# which is required for SGE's -hold_jid_ad task-aligned dependency mode.
# ============================================================================

# Cross-module dependency maps (populated by process_all; empty dicts in standalone).
# These wire tracula/bamos/dti/noddi → metrics DAG edges without requiring
# use-rule overrides (which cannot be used inside a Python for loop).
_bamos_map     = cfg.get("bamos_correction_map", {})
_dti_fa_map    = cfg.get("dti_fa_map", {})
_dti_md_map    = cfg.get("dti_md_map", {})
_dti_ad_map    = cfg.get("dti_ad_map", {})
_dti_rd_map    = cfg.get("dti_rd_map", {})
_noddi_odi_map = cfg.get("noddi_odi_map", {})
_thalamic_map  = cfg.get("thalamic_nuclei_map", {})
_tract_qc_map  = cfg.get("tract_qc_report_map", {})  # {subject: {parcellation: path}}
_z_score_map   = cfg.get("z_score_map", {})          # {subject: [paths]}

wildcard_constraints:
    subject="|".join(re.escape(s) for s in SUBJECTS),

rule all:
    input:
        # Always-present per-subject outputs
        [f"{OUTPUT_DIR}/{subject}/metrics-{parc}/outputs/metrics/whole_brain_metrics.csv"
         for subject in SUBJECTS for parc in PARCELLATIONS] +
        [f"{OUTPUT_DIR}/{subject}/metrics-{parc}/outputs/metrics/tract_level_metrics.csv"
         for subject in SUBJECTS for parc in PARCELLATIONS] +
        [f"{OUTPUT_DIR}/{subject}/metrics-{parc}/outputs/metrics/skeleton_metrics.csv"
         for subject in SUBJECTS for parc in PARCELLATIONS] +
        # Per-lesion outputs (empty list when no lesion path is configured)
        [f"{OUTPUT_DIR}/{subject}/metrics-{parc}/outputs/metrics/{ln}_lesion_metrics.csv"
         for subject in SUBJECTS for parc in PARCELLATIONS for ln in LESION_NAMES] +
        [f"{OUTPUT_DIR}/{subject}/metrics-{parc}/outputs/metrics/{ln}_tract_aggregated_lesion_metrics.csv"
         for subject in SUBJECTS for parc in PARCELLATIONS for ln in LESION_NAMES]

for _parc in PARCELLATIONS:
    rule:
        name: f"extract_metrics_{_parc}"
        input:
            # Cross-module dependencies - empty lists in standalone context.
            # process_all populates these maps so Snakemake builds the full DAG.
            bamos_correction=lambda wc, m=_bamos_map: ([m[wc.subject]] if wc.subject in m else []),
            dti_fa=lambda wc, m=_dti_fa_map: ([m[wc.subject]] if wc.subject in m else []),
            dti_md=lambda wc, m=_dti_md_map: ([m[wc.subject]] if wc.subject in m else []),
            dti_ad=lambda wc, m=_dti_ad_map: ([m[wc.subject]] if wc.subject in m else []),
            dti_rd=lambda wc, m=_dti_rd_map: ([m[wc.subject]] if wc.subject in m else []),
            noddi_odi=lambda wc, m=_noddi_odi_map: ([m[wc.subject]] if wc.subject in m else []),
            thalamic_nuclei=lambda wc, m=_thalamic_map: ([m[wc.subject]] if wc.subject in m else []),
            tract_qc_report=lambda wc, p=_parc, m=_tract_qc_map: (
                [m[wc.subject][p]] if wc.subject in m and p in m[wc.subject] else []
            ),
            z_scores=lambda wc, m=_z_score_map: (m.get(wc.subject) or []),
        output:
            csv=f"{OUTPUT_DIR}/{{subject}}/metrics-{_parc}/outputs/metrics/whole_brain_metrics.csv",
            tract_level_metrics=f"{OUTPUT_DIR}/{{subject}}/metrics-{_parc}/outputs/metrics/tract_level_metrics.csv",
            skeleton_metrics=f"{OUTPUT_DIR}/{{subject}}/metrics-{_parc}/outputs/metrics/skeleton_metrics.csv",
            lesion_metrics=[
                f"{OUTPUT_DIR}/{{subject}}/metrics-{_parc}/outputs/metrics/{ln}_lesion_metrics.csv"
                for ln in LESION_NAMES
            ],
            lesion_tract_aggregated=[
                f"{OUTPUT_DIR}/{{subject}}/metrics-{_parc}/outputs/metrics/{ln}_tract_aggregated_lesion_metrics.csv"
                for ln in LESION_NAMES
            ],
        container:
            CONTAINER_SIF,
        params:
            tractography_path_sif=lambda wildcards, p=_parc: _sif_tractography_path(wildcards.subject, p),
            lesion_path_sif=lambda wildcards: _sif_lesion_path(wildcards.subject),
            t1_path_sif=lambda wildcards: _sif_t1_path(wildcards.subject),
            dwi_path_sif=lambda wildcards: _sif_dwi_path(wildcards.subject),
            metrics_sif_spec=lambda wildcards: _metrics_spec_str(wildcards.subject),
            output_dir_sif=lambda wildcards, p=_parc: _sif_output_dir(wildcards.subject, p),
            verbose=VERBOSE,
            metrics_script_sif=METRICS_SCRIPT_SIF,
            qc_report_sif=lambda wildcards, p=_parc: _sif_qc_report(wildcards.subject, p),
            log_file=lambda wildcards, p=_parc: f"{_sif_logs_dir(wildcards.subject, p)}/metrics_log.txt",
            error_file=lambda wildcards, p=_parc: f"{_sif_logs_dir(wildcards.subject, p)}/metrics_error.txt",
        resources:
            # Bumped from 16GB (2026-08-16): mri_synthstrip (CPU-only CNN
            # skull-strip, register_dwi_to_t1.sh) spikes to ~4.7GB in one
            # allocation and intermittently OOM'd 24/816 ADNI3 subjects at
            # 16GB. 18GB is a modest first try, not a large blanket bump.
            mem_mb=18 * 1024,
            disk_mb=2 * 1024,
            time="24:00:00",
            name=f"metrics_{_parc}",
            workdir=lambda wildcards, p=_parc: f"{OUTPUT_DIR}/{wildcards.subject}/metrics-{p}",
        shell:
            """
            source /leukoquant/leukoquant/utils/bash_utils.sh
            mkdir -p "$(dirname "{params.log_file}")"
            exec > "{params.log_file}"
            exec 2> "{params.error_file}"

            echo "Date: $(date)"
            echo "Date: $(date)" >&2

            export PYTHONUNBUFFERED=1
            export FS_LICENSE=/license/license.txt
            [ -f /usr/local/fsl/etc/fslconf/fsl.sh ] && . /usr/local/fsl/etc/fslconf/fsl.sh

            mkdir -p {params.output_dir_sif}

            VERBOSE_FLAG=""
            if [ "{params.verbose}" = "True" ]; then
                VERBOSE_FLAG="--verbose"
            fi

            bash {params.metrics_script_sif} \
                --subject {wildcards.subject} \
                --tractography-path {params.tractography_path_sif} \
                --output-dir {params.output_dir_sif} \
                --lesion-path "{params.lesion_path_sif}" \
                --t1-path "{params.t1_path_sif}" \
                --dwi "{params.dwi_path_sif}" \
                --metrics "{params.metrics_sif_spec}" \
                --qc-report "{params.qc_report_sif}" \
                $VERBOSE_FLAG \
                --output-csv {params.output_dir_sif}
            """

