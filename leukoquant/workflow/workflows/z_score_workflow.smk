"""
Snakemake workflow to compute Z-scores for target subjects against a healthy cohort.

For each target subject a single rule runs z_score_calc.sh, which:
  - registers every healthy T1 to target space once (reused for all metrics),
  - then processes all metrics in one pass (resample, merge, GLM, Z-score).

Python helpers (path resolution) live in leukoquant/utils/z_score_utils.py.
Requires config key "leukoquant_parent_dir" pointing to the leukoquant parent directory.
"""
import os
import re
import sys
from pathlib import Path

LEUKOQUANT_PARENT_DIR = config.get("leukoquant_parent_dir")

sys.path.insert(0, LEUKOQUANT_PARENT_DIR)

from leukoquant.utils.z_score_utils import get_t1_path, get_metric_path, get_dwi_path, get_bval_path, translate_path
from leukoquant.utils.container_utils import ensure_container


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
if isinstance(config.get("healthy_subjects_list"), list):
    HEALTHY_SUBJECTS = config["healthy_subjects_list"]
else:
    with open(config["healthy_subjects_list"], "r") as _f:
        HEALTHY_SUBJECTS = [line.strip() for line in _f if line.strip()]

num_healthy = len(HEALTHY_SUBJECTS)
if num_healthy == 0:
    raise ValueError("Healthy subjects list is empty. Please check the path and contents of 'healthy_subjects_list' in the config.")

if isinstance(config.get("target_subjects_list"), list):
    TARGET_SUBJECTS = config["target_subjects_list"]
else:
    with open(config["target_subjects_list"], "r") as _f:
        TARGET_SUBJECTS = [line.strip() for line in _f if line.strip()]

# Drop any requested metric that isn't available for every HEALTHY subject --
# the healthy cohort's data for a given metric is loaded once and shared
# across every target subject's z-score computation (a single reference
# distribution), so it genuinely must be complete or the metric can't be
# computed for anyone -- this really is a hard, permanent requirement (not a
# "still running" one), since a NODDI/DTI run for a healthy subject that
# never produced a given map isn't something a rerun of THIS DAG can fix.
_requested_metrics = list(config.get("metrics", {}).keys())
METRICS = []
_healthy_gaps = {}
for _m in _requested_metrics:
    _missing = [_s for _s in HEALTHY_SUBJECTS if not os.path.isfile(get_metric_path(_s, _m, config))]
    if _missing:
        _healthy_gaps[_m] = _missing
    else:
        METRICS.append(_m)
if _healthy_gaps:
    _lines = "\n".join(
        f"  - '{_m}': missing for {len(_ss)}/{len(HEALTHY_SUBJECTS)} healthy subject(s), e.g. {_ss[:5]}"
        for _m, _ss in _healthy_gaps.items()
    )
    print(
        f"WARNING: dropping {len(_healthy_gaps)} z-score metric(s) -- incomplete for the healthy "
        f"reference cohort (this blocks the metric for every subject, since the healthy "
        f"cohort is a single shared reference distribution):\n{_lines}",
        file=sys.stderr,
    )
if not METRICS:
    raise ValueError(
        "No requested z-score metric is available for every healthy subject in this "
        "run -- the healthy reference cohort itself is incomplete, nothing left to compute."
    )

METRIC_SPACE = config.get("metric_space", "t1")
OUTPUT_SPACE = config.get("output_space", "t1")
OUTPUT_DIR   = config.get("output_dir", "outputs/z_scores")
TOOL_NAME    = "z-score"

# Build filename suffix from covariates and poly terms
COVARIATES = config.get("covariates", [])
POLY_TERMS = config.get("polynomial_terms", [])
FILENAME_SUFFIX = ""
if COVARIATES and len(COVARIATES) > 0:
    cov_safe = "_".join(COVARIATES)
    FILENAME_SUFFIX += f"_cov_{cov_safe}"
if POLY_TERMS and len(POLY_TERMS) > 0:
    poly_safe = "_".join(POLY_TERMS)
    FILENAME_SUFFIX += f"_poly_{poly_safe}"

# Validate DWI requirements early (before job submission)
if METRIC_SPACE == "dwi" or OUTPUT_SPACE == "dwi":
    if "dwi_pattern" not in config and "dwi_map" not in config:
        raise ValueError(
            "ERROR: metric_space/output_space requires DWI but neither dwi_pattern "
            "nor dwi_map is present in config."
        )

CONTAINER_SIF = os.path.join(
    LEUKOQUANT_PARENT_DIR,
    "leukoquant/workflow/containers/",
    config.get("container_name", "freesurfer_unified_container") + ".sif"
)
ensure_container(CONTAINER_SIF)

BIND_MAP        = config.get("singularity_binds", {})
CALC_SCRIPT_SIF = "/leukoquant/leukoquant/utils/z_score_calc.sh"

# ---------------------------------------------------------------------------
# Pre-compute constants shared across all target subjects
# (healthy subjects' T1s and metrics are the same for every target)
# ---------------------------------------------------------------------------

# Healthy T1 host paths and their container equivalents
_h_ht1s     = [get_t1_path(hs, config) for hs in HEALTHY_SUBJECTS]
_sif_ht1s   = " ".join(translate_path(p, BIND_MAP) for p in _h_ht1s)

# Healthy metric host paths (grouped metric-first, then per healthy subject)
_h_hmetrics   = [get_metric_path(hs, m, config) for m in METRICS for hs in HEALTHY_SUBJECTS]
_sif_hmetrics = " ".join(translate_path(p, BIND_MAP) for p in _h_hmetrics)

# Demographics and healthy-subjects list (identical for every target)
_demographics_path = config.get("demographics_csv", "")
_sif_demographics  = translate_path(str(Path(_demographics_path).resolve()), BIND_MAP) if _demographics_path else ""

_healthy_list_path = config["healthy_subjects_list"] if not isinstance(config["healthy_subjects_list"], list) else ""
_sif_healthy_list  = translate_path(str(Path(_healthy_list_path).resolve()), BIND_MAP) if _healthy_list_path else ""

# DWI paths for healthy subjects (only when metric/output space is DWI)
_h_healthy_dwis = []
_h_healthy_bvals = []
if METRIC_SPACE == "dwi":
    _h_healthy_dwis  = [get_dwi_path(hs, config) for hs in HEALTHY_SUBJECTS]
    if "bval_pattern" in config or "bval_map" in config:
        _h_healthy_bvals = [get_bval_path(hs, config) for hs in HEALTHY_SUBJECTS]

# Per-target lookup dicts built at module load time so lambda functions stay fast
_target_t1_map = {s: get_t1_path(s, config) for s in TARGET_SUBJECTS}

# Target metrics declared unconditionally (not existence-checked here) -- for
# process-all callers (deterministic paths from upstream dti/noddi rules in
# this SAME DAG), this lets Snakemake's own dependency resolution (hold_jid
# chaining via the sge executor plugin) wait for dti/noddi to finish before
# z-score runs, rather than a Python pre-check treating "doesn't exist yet"
# as "never will" and skipping the subject before Snakemake gets a chance to
# schedule the wait. Only exclude a target subject if get_metric_path itself
# raises -- which only happens for the standalone glob-based CLI caller
# (resolve_single_match finds zero/multiple matches, a structural
# impossibility no rerun of this DAG can fix); process-all's deterministic
# colon-path format never raises here regardless of current file existence.
_target_metrics_map = {}
_targets_unresolvable = []
for _s in TARGET_SUBJECTS:
    try:
        _target_metrics_map[_s] = [get_metric_path(_s, _m, config) for _m in METRICS]
    except Exception as _e:
        _targets_unresolvable.append((_s, str(_e)))
if _targets_unresolvable:
    print(
        f"WARNING: excluding {len(_targets_unresolvable)} target subject(s) -- metric path "
        f"could not be resolved at all (e.g. {_targets_unresolvable[0][0]}: {_targets_unresolvable[0][1]})",
        file=sys.stderr,
    )
    _unresolvable_set = {s for s, _ in _targets_unresolvable}
    TARGET_SUBJECTS = [s for s in TARGET_SUBJECTS if s not in _unresolvable_set]

_target_dwi_map  = {}
_target_bval_map = {}
if METRIC_SPACE == "dwi" or OUTPUT_SPACE == "dwi":
    _target_dwi_map  = {s: get_dwi_path(s, config) for s in TARGET_SUBJECTS}
    if "bval_pattern" in config or "bval_map" in config:
        _target_bval_map = {s: get_bval_path(s, config) for s in TARGET_SUBJECTS}

# Per-job /scratch0 sizing (see the z_score rule's `resources:` below).
#
# Two additive components live in the SAME per-job $TMP_DIR (see
# z_score_calc.sh's trap/cleanup comment: "target's own T1/DWI prep, ALL
# healthy->target registration output, design matrix, xfm files" -- the
# earlier version of this formula only counted the first component,
# confirmed as the real cause of repeated "No space left on device"
# failures on ADNI3 (2026-08-18) despite the reservation mechanism itself
# working correctly (SGE's tscratch consumable, mapped from this value,
# does get decremented -- the requested amount was just too small):
#
#   1. Target's own raw-DWI prep: prepare_dwi_input() copies one subject's
#      raw DWI into scratch at a time (cleaned up immediately after, see
#      dwi_utils.py), so this scales with the single LARGEST raw DWI file
#      this run could ever touch, not with healthy-cohort size.
#   2. Healthy-cohort registration: z_score_calc.sh registers EVERY healthy
#      subject's metric map into target space, one .nii.gz per
#      (metric, healthy subject) pair (see the printf "%04d" indexed
#      $METRIC_TMP/$IDX/ loop), plus a merged 4D concatenation per metric
#      for the GLM step -- both scale with num_healthy * n_metrics, and
#      neither is covered by (1) since they never touch DWI files at all.
#      Sized off real registered/merged z-score output files on disk
#      (~4MB compressed at typical DTI resolution, confirmed 2026-08-18),
#      doubled for margin since T1-space registration can differ in
#      resolution from the DWI-space files this was measured against.
def _file_size_mb(path):
    try:
        return os.path.getsize(path) / (1024 * 1024)
    except OSError:
        return 0

_max_dwi_mb = max(
    [_file_size_mb(p) for p in list(_target_dwi_map.values()) + _h_healthy_dwis] or [1024]
)
_target_dwi_scratch_mb = 2 * 1024 + 3 * _max_dwi_mb
_per_registered_metric_mb = 8  # ~4MB real file size, doubled for margin
_n_metrics = max(len(METRICS), 1)
# x2: per-healthy registered files (IDX/*.nii.gz) plus the merged 4D
# concatenation per metric, both scaling with num_healthy * n_metrics.
_healthy_registration_scratch_mb = 2 * num_healthy * _n_metrics * _per_registered_metric_mb
# Capped at 15GB (not the old fixed 10GB) -- confirmed 2026-08-18 that real
# node-local /scratch0 capacity varies widely across the cluster, from
# ~15GB to 160G+ depending on host. Kept at the smallest real node size
# rather than raised further, so a job can never request more than the
# smallest node has in total and get stuck permanently unschedulable there
# -- a single pathological job claiming a whole small node's scratch is an
# acceptable tradeoff against that.
ZSCORE_SCRATCH_MB = min(
    int(_target_dwi_scratch_mb + _healthy_registration_scratch_mb), 15 * 1024
)

# Per-job memory (see the z_score rule's `resources:` below). fsl_glm
# regresses across the full merged healthy-cohort 4D volume when covariates
# are configured, so its real memory need does grow with num_healthy -- but
# scaling the request up per-dataset (tried 2026-08-13, reverted) makes
# larger-cohort datasets (EPAD, 72 healthy) request 32G+ jobs that sit
# queued far longer on a busy shared cluster than a 16G request would.
# Fixed at 16G for now (confirmed real OOM was at 12G) while the actual fix
# -- bounding fsl_glm's memory independent of cohort size, e.g. chunked
# GLM instead of one whole-volume call -- is investigated separately.
ZSCORE_MEM_MB = 16 * 1024


def _sif_target_t1(subject):
    return translate_path(_target_t1_map[subject], BIND_MAP)

def _sif_target_metrics(subject):
    return " ".join(translate_path(p, BIND_MAP) for p in _target_metrics_map[subject])

def _sif_output_dir(subject):
    return translate_path(str(Path(f"{OUTPUT_DIR}/{subject}/{TOOL_NAME}/outputs").resolve()), BIND_MAP)

# Dataset-level (not per-subject) -- shared healthy-cohort DWI/T1 prep
# cache, reused across every target subject's run. See z_score_calc.sh's
# prepare_dwi_inputs/prepare_t1_inputs for why this is safe to share: those
# steps never reference the target at all.
_sif_shared_cache_dir = translate_path(
    str(Path(f"{OUTPUT_DIR}/_zscore_healthy_cache").resolve()), BIND_MAP
)

def _sif_target_dwis(subject):
    """Space-separated DWI paths in container space: target first, then healthy."""
    if not _target_dwi_map:
        return ""
    paths = [_target_dwi_map[subject]] + _h_healthy_dwis
    return " ".join(translate_path(p, BIND_MAP) for p in paths)

def _sif_target_bvals(subject):
    if not _target_bval_map:
        return ""
    paths = [_target_bval_map[subject]] + _h_healthy_bvals
    return " ".join(translate_path(p, BIND_MAP) for p in paths)


# ---------------------------------------------------------------------------
# Rules
# ---------------------------------------------------------------------------

wildcard_constraints:
    subject="|".join(re.escape(s) for s in TARGET_SUBJECTS),

rule all:
    input:
        [
            f"{OUTPUT_DIR}/{subject}/{TOOL_NAME}/outputs/{metric}_z_score{FILENAME_SUFFIX}.nii.gz"
            for subject in TARGET_SUBJECTS
            for metric in METRICS
        ]

rule z_score:
    input:
        target_t1      = lambda wildcards: _target_t1_map[wildcards.subject],
        demographics   = [_demographics_path] if _demographics_path else [],
        healthy_t1s    = _h_ht1s,
        target_metrics = lambda wildcards: _target_metrics_map[wildcards.subject],
        healthy_metrics = _h_hmetrics,
        dwis  = lambda wildcards: (
            [_target_dwi_map[wildcards.subject]] + _h_healthy_dwis
            if METRIC_SPACE == "dwi" else
            ([_target_dwi_map[wildcards.subject]] if OUTPUT_SPACE == "dwi" and _target_dwi_map else [])
        ),
        bvals = lambda wildcards: (
            [_target_bval_map[wildcards.subject]] + _h_healthy_bvals
            if METRIC_SPACE == "dwi" and _target_bval_map else []
        ),
    output:
        # Static full METRICS list (Snakemake doesn't allow output: to be a
        # per-wildcard function, unlike input:/params:) -- every target
        # subject is expected to eventually produce all of these, with
        # Snakemake's own DAG/hold_jid chaining waiting on dti/noddi as
        # needed (see _target_metrics_map above).
        z_scores = [
            f"{OUTPUT_DIR}/{{subject}}/{{TOOL_NAME}}/outputs/{m}_z_score{FILENAME_SUFFIX}.nii.gz"
            for m in METRICS
        ],
    params:
        target_id           = lambda wildcards: wildcards.subject,
        covariates          = ",".join(COVARIATES) if COVARIATES else "",
        poly_terms          = ",".join(POLY_TERMS) if POLY_TERMS else "",
        metric_names        = ",".join(METRICS),
        metric_space        = METRIC_SPACE,
        output_space        = OUTPUT_SPACE,
        skip_skullstrip_t1  = config.get("skip_skullstrip_t1", False),
        skip_skullstrip_dwi = config.get("skip_skullstrip_dwi", False),
        calc_script_sif     = CALC_SCRIPT_SIF,
        target_t1_sif       = lambda wildcards: _sif_target_t1(wildcards.subject),
        demographics_sif    = _sif_demographics,
        healthy_list_sif    = _sif_healthy_list,
        output_dir_sif      = lambda wildcards: _sif_output_dir(wildcards.subject),
        shared_cache_sif    = _sif_shared_cache_dir,
        target_metrics_sif  = lambda wildcards: _sif_target_metrics(wildcards.subject),
        healthy_t1s_sif     = _sif_ht1s,
        healthy_metrics_sif = _sif_hmetrics,
        dwis_sif            = lambda wildcards: _sif_target_dwis(wildcards.subject),
        bvals_sif           = lambda wildcards: _sif_target_bvals(wildcards.subject),
        threads             = config.get("threads", 4),
        log_file            = lambda wildcards: f"{OUTPUT_DIR}/{wildcards.subject}/{TOOL_NAME}/logs/z_score.log",
        error_file          = lambda wildcards: f"{OUTPUT_DIR}/{wildcards.subject}/{TOOL_NAME}/logs/z_score_error.log",
    container:
        CONTAINER_SIF,
    resources:
        # Scales with healthy-cohort size -- see ZSCORE_MEM_MB above.
        mem_mb = ZSCORE_MEM_MB,
        # Not previously requested, so each task fell back to whatever
        # small default node-local /tmp offers. prepare_dwi_input() copies
        # one subject's raw DWI into a scratch workspace at a time
        # (cleaned up immediately after -- see dwi_utils.py), so this does
        # NOT need to scale with healthy-cohort size, only with the size of
        # the single largest DWI file this run could ever touch --
        # computed once above as ZSCORE_SCRATCH_MB (real file sizes, 2GB
        # base + 3x headroom, capped at 10GB so one job can't claim most of
        # a node's scratch partition, confirmed only ~15GB total and shared
        # with every other job on that node). A fixed 3GB here previously
        # exhausted node-local /tmp on a 46-subject healthy cohort because
        # nothing was cleaning up between subjects at all (separate leak,
        # also fixed) -- this is deliberately lenient headroom on top of
        # that real fix, not a substitute for it.
        scratch_size=ZSCORE_SCRATCH_MB,
        time   = f"{num_healthy}:00:00",
        workdir=lambda wildcards: f"{OUTPUT_DIR}/{wildcards.subject}/{TOOL_NAME}",
        name   = lambda wildcards: f"z_score_{wildcards.subject.replace('/', '_')}",
        # Each task's container does its first mkdir/write to the shared NFS
        # SAN within the same second under --immediate-submit; submitting a
        # large array fully unthrottled can overwhelm the NFS server and cause
        # a burst of transient "Read-only file system" errors on nearly every
        # task's very first write (observed: 352/353 tasks died this way on a
        # single ADNI3 run). No limit by default (matches prior behavior);
        # set config["task_concurrency"] / --task-concurrency to cap
        # concurrently-running array tasks via SGE's -tc if this recurs.
        sge_task_concurrency = config.get("task_concurrency"),
    shell:
        """
        source /leukoquant/leukoquant/utils/bash_utils.sh
        mkdir -p "$(dirname '{params.log_file}')"
        exec > '{params.log_file}'
        exec 2> '{params.error_file}'

        # Point TMPDIR at per-job scratch (/scratch0), never the node's
        # local /tmp -- same pattern as tracula_workflow.smk. Without this,
        # dwi_utils.py's prepare_dwi_input() (via tempfile.mkdtemp(), which
        # honours TMPDIR) falls back to node-local /tmp for every healthy
        # subject's DWI scratch copy; on a busy shared node that's a small,
        # contended partition, and 46 healthy subjects processed
        # sequentially in one job exhausted it ("No space left on device")
        # even after fixing the separate cleanup leak (see dwi_utils.py).
        MY_JOB_ID="$(get_job_id)"
        scratch_path=""
        trap 'scratch_cleanup "$scratch_path"' EXIT INT TERM
        # Isolated single-quoted assignment (never embedded inside a longer
        # double-quoted string) -- matches the only substitution form
        # confirmed to survive Snakemake's own templating intact.
        ZSCORE_TARGET_ID='{params.target_id}'
        scratch_path="/scratch0/$USER/$MY_JOB_ID/zscore_$ZSCORE_TARGET_ID"
        mkdir -p "$scratch_path"
        export TMPDIR="$scratch_path/tmp"
        mkdir -p "$TMPDIR"

        # Every param substitution below is single-quoted -- confirmed by
        # capturing the real, fully-rendered dollar-CMD from a live run
        # (dumped to disk immediately before eval) that Snakemake's own
        # param substitution silently STRIPS surrounding double-quote
        # characters (a double-quoted reference renders as a bare,
        # unquoted value -- true for every double-quoted entry, not just
        # empty ones), while single-quoted references survive substitution
        # intact, including as a real empty '' token when the value is
        # empty. A double-quoted empty param (e.g. demographics_sif when
        # --demographics-csv isn't passed) therefore rendered as literally
        # nothing, collapsing its whole line and shifting every subsequent
        # flag/value pair left by one position -- confirmed in practice: an
        # empty demographics_sif caused --healthy-ids to be silently
        # consumed as --demographics's value, leaving the real healthy-ids
        # path as an unrecognized stray argument ("Unknown argument:
        # <path>") with no output ever produced. An earlier fix attempted to
        # solve this with double-quotes, which looked correct in the source
        # but never actually survived Snakemake's substitution -- single
        # quotes are the only form that does.
        CMD="bash {params.calc_script_sif} \\
            --target-t1       '{params.target_t1_sif}' \\
            --target-id       '{params.target_id}' \\
            --demographics    '{params.demographics_sif}' \\
            --healthy-ids     '{params.healthy_list_sif}' \\
            --output-dir      '{params.output_dir_sif}' \\
            --shared-cache-dir '{params.shared_cache_sif}' \\
            --metrics         '{params.metric_names}' \\
            --covariates      '{params.covariates}' \\
            --poly-terms      '{params.poly_terms}' \\
            --threads         '{params.threads}' \\
            --metric-space    '{params.metric_space}' \\
            --output-space    '{params.output_space}' \\
            --healthy-t1s     {params.healthy_t1s_sif} \\
            --target-metrics  {params.target_metrics_sif} \\
            --healthy-metrics {params.healthy_metrics_sif}"

        if [ -n '{params.dwis_sif}' ]; then
            CMD="$CMD --dwi-paths {params.dwis_sif}"
            if [ -n '{params.bvals_sif}' ]; then
                CMD="$CMD --bval-paths {params.bvals_sif}"
            fi
        fi

        if [ '{params.skip_skullstrip_t1}' = "True" ]; then
            CMD="$CMD --skip-skullstrip-t1"
        fi

        if [ '{params.skip_skullstrip_dwi}' = "True" ]; then
            CMD="$CMD --skip-skullstrip-dwi"
        fi

        eval "$CMD"
        """
