#!/usr/bin/env python
"""
Z-score processing for leukoquant.

Runs the `z_score_workflow.smk` Snakemake workflow to compute Z-scores
for target subjects against a healthy cohort.
"""

import glob
import logging
import re
import subprocess
import sys
from pathlib import Path
import os
from typing import Dict, Iterable, List, Optional, Tuple

import yaml

from leukoquant.utils.container_utils import (
    ensure_container,
    FREESURFER_SIF_FILENAME,
    FREESURFER_SIF_HF_PATH,
)
from leukoquant.utils.external_utils import check_sge_plugin
from leukoquant.utils.snakemake_utils import add_forcerun_args, load_yaml_config, first_truthy

logger = logging.getLogger(__name__)


def _normalize_list_value(value: Optional[object]) -> Optional[list[str]]:
    if value is None:
        return None
    if isinstance(value, str):
        parts = [item.strip() for item in value.split(",") if item.strip()]
        return parts
    if isinstance(value, (list, tuple)):
        return [str(item).strip() for item in value if str(item).strip()]
    raise ValueError("Expected a string or list for list-like config values")


def _parse_polynomial_terms(value: Optional[object]) -> Optional[list[str]]:
    """Parse polynomial terms with optional power notation.
    
    Supports formats:
    - "age" -> ["age_2"] (power defaults to 2)
    - "age:2" -> ["age_2"] (squared using colon)
    - "age^2" -> ["age_2"] (squared using caret)
    - "age:2,age:3" -> ["age_2", "age_3"] (multiple powers)
    - ["age", "age:2", "age^3"] -> ["age_2", "age_2", "age_3"] (mixed formats)
    
    Args:
        value: Polynomial term specification (string, list, or None)
    
    Returns:
        List of formatted polynomial term names (underscore notation) or None if input is None
    """
    if value is None:
        return None
    
    # Normalize input to a list
    if isinstance(value, str):
        terms = [item.strip() for item in value.split(",") if item.strip()]
    elif isinstance(value, (list, tuple)):
        terms = [str(item).strip() for item in value if str(item).strip()]
    else:
        raise ValueError("Expected a string or list for polynomial terms")
    
    # Parse each term
    parsed_terms: list[str] = []
    for term in terms:
        # Support both ":" and "^" as power separators
        if ":" in term:
            name, power_str = term.split(":", 1)
            name = name.strip()
        elif "^" in term:
            name, power_str = term.split("^", 1)
            name = name.strip()
        else:
            # No power specified, default to 2
            name = term.strip()
            power_str = "2"
        
        try:
            power = int(power_str.strip())
            if power < 1:
                raise ValueError(f"Polynomial power must be >= 1, got {power} for term '{term}'")
            parsed_terms.append(f"{name}_{power}")
        except ValueError as e:
            if "invalid literal" in str(e):
                raise ValueError(f"Invalid polynomial power in term '{term}': must be an integer")
            raise
    
    return parsed_terms if parsed_terms else None


def _parse_metrics_spec(metrics: Optional[object]) -> Dict[str, str]:
    if metrics is None:
        return {}
    if isinstance(metrics, dict):
        return {str(k): str(v) for k, v in metrics.items()}
    if isinstance(metrics, (list, tuple)):
        parts: list[str] = []
        for item in metrics:
            if isinstance(item, str):
                parts.extend([seg for seg in item.split(",") if seg.strip()])
            else:
                raise ValueError("Metric list must contain only strings")
        metrics = parts
    if isinstance(metrics, str):
        metrics = [seg for seg in metrics.split(",") if seg.strip()]

    if isinstance(metrics, list):
        parsed: Dict[str, str] = {}
        for item in metrics:
            if "=" not in item:
                raise ValueError("Metric mapping must use name=pattern format")
            name, pattern = item.split("=", 1)
            parsed[name.strip()] = pattern.strip()
        return parsed

    raise ValueError("Unsupported metrics format")


def _normalize_subject_pattern(pattern_spec: str, label: str) -> str:
    #logger.info(f"Normalizing pattern for {label}: {pattern_spec}")
    if ":" in pattern_spec:
        base, glob_part = pattern_spec.split(":", 1)
        #logger.info(f"Parsed {label} pattern - base: '{base}', glob_part: '{glob_part}', pattern_spec: '{pattern_spec}'")
        base = base.strip().rstrip("/\\")
        glob_part = glob_part.strip().lstrip("/\\")
        #logger.info(f"Parsed {label} pattern - base: '{base}', glob_part: '{glob_part}', pattern_spec: '{pattern_spec}'")
        if not base or not glob_part:
            raise ValueError(f"Invalid {label} pattern. Expected base:glob format.")
        full_path = Path(base) / "{subject}" / Path(glob_part)
        full_path = full_path.resolve()
        #logger.info(f"Path base {Path(base)}, Path glob_part: '{Path(glob_part)}'")

        #logger.info(f"Label: {label}, path: {str(full_path)}")
        return str(full_path)
    elif "{subject}" in pattern_spec:
        # Already-normalized pattern (e.g. re-fed from a z_score_config.yaml
        # this same function previously produced, which stores the expanded
        # {subject}-placeholder path rather than the raw base:glob input
        # format). Idempotent passthrough so --config-yaml re-entry works.
        return pattern_spec
    else:
        #logger.info(f"No base specified for {label} pattern.")
        raise ValueError(f"Invalid {label} pattern. Expected base:glob format.")


def _read_subject_list(list_path: str) -> list[str]:
    with open(list_path, "r") as f:
        return [line.strip() for line in f if line.strip()]


def _validate_glob_matches(pattern: str, subjects: Iterable[str], label: str) -> None:
    for subject in subjects:
        candidate = pattern.format(subject=subject)
        matches = sorted(glob.glob(candidate))
        if not matches:
            raise ValueError(f"No {label} file found for subject '{subject}' with pattern: {candidate}. Full path: {Path(candidate).resolve()}")
        if len(matches) > 1:
            raise ValueError(
                f"Multiple {label} files found for subject '{subject}' with pattern: {candidate}. "
                f"Matches: {matches}"
            )


def _subject_has_match(pattern: str, subject: str) -> bool:
    """Like _validate_glob_matches for one subject, but returns False instead
    of raising when the file is missing (multiple matches still count as an
    error condition -- ambiguous data should never be silently dropped)."""
    candidate = pattern.format(subject=subject)
    matches = glob.glob(candidate)
    if len(matches) > 1:
        raise ValueError(
            f"Multiple files found for subject '{subject}' with pattern: {candidate}. Matches: {matches}"
        )
    return len(matches) == 1


def _drop_metrics_missing_for_any_subject(
    metrics: Dict[str, str], subjects: Iterable[str]
) -> Dict[str, str]:
    """Drop any metric from `metrics` that isn't available for every subject,
    rather than crashing the whole run. Subjects are typically a small,
    intentionally-curated cohort (healthy averages + targets), so a metric
    that's missing for even one of them (e.g. an older NODDI fit that never
    produced ndi/fwf, only odi) is dropped for the whole run with a warning,
    keeping every metric that *is* available for all subjects."""
    subjects = list(subjects)
    kept: Dict[str, str] = {}
    for name, pattern in metrics.items():
        missing = [s for s in subjects if not _subject_has_match(pattern, s)]
        if missing:
            logger.warning(
                f"Dropping metric '{name}' from this z-score run: missing for "
                f"{len(missing)}/{len(subjects)} subject(s) (e.g. {missing[:5]}). "
                f"Pattern: {pattern}"
            )
        else:
            kept[name] = pattern
    if not kept:
        raise ValueError(
            "No requested metric is available for every subject in this run -- "
            "nothing left to compute Z-scores for."
        )
    return kept


def _discover_subjects_from_pattern(pattern: str) -> List[str]:
    """Discover every subject ID that has a match for a {subject}-templated
    pattern, by globbing with {subject} replaced by a wildcard and recovering
    the {subject} value from each match via regex. Handles both single-level
    subject IDs (e.g. "subj-002-s-6007") and multi-level ones (e.g. DPUK's
    "NEW002/NEW002_PETMR_V2").

    Tries bounded depths first (1 path segment, then 2) since those are a
    single non-recursive glob() call each and cover every known leukoquant
    dataset's subject-ID convention; only falls back to an unbounded
    recursive "**" glob (which walks the entire tree and can be very slow on
    a large NFS-backed output directory) if neither bounded attempt matches
    anything.
    """
    if "{subject}" not in pattern:
        raise ValueError(f"Pattern must contain '{{subject}}' for auto-discovery: {pattern}")

    before, after = pattern.split("{subject}", 1)

    def _glob_to_regex(segment: str) -> str:
        segment = re.escape(segment)
        return segment.replace(r"\*", "[^/]*").replace(r"\?", "[^/]")

    def _try(wildcard: str, recursive: bool) -> List[str]:
        # For a bounded wildcard (e.g. "*/*"), the {subject} capture must
        # span exactly that many path segments -- match "[^/]*" per segment
        # joined by literal "/", not a greedy cross-segment ".+?" (reserved
        # for the unbounded "**" fallback).
        if wildcard == "**":
            subject_re = "(.+?)"
        else:
            subject_re = "(" + "/".join("[^/]*" for _ in wildcard.split("/")) + ")"
        matcher = re.compile("^" + _glob_to_regex(before) + subject_re + _glob_to_regex(after) + "$")
        glob_pattern = pattern.replace("{subject}", wildcard)
        found: set = set()
        for match_path in glob.glob(glob_pattern, recursive=recursive):
            m = matcher.match(match_path)
            if m:
                found.add(m.group(1))
        return sorted(found)

    for wildcard, recursive in (("*", False), ("*/*", False), ("**", True)):
        result = _try(wildcard, recursive)
        if result:
            return result
    return []


# ---------------------------------------------------------------------------
# Singularity bind-map helpers
# ---------------------------------------------------------------------------

def _base_dir_from_pattern(pattern: str) -> str:
    """Return the absolute host base directory – the path component before {subject}."""
    idx = pattern.find("{subject}")
    if idx == -1:
        return str(Path(pattern).parent.resolve())
    prefix = pattern[:idx].rstrip("/")
    return str(Path(prefix).resolve()) if prefix else "/"


def _build_bind_map(named_dirs: Dict[str, str]) -> Dict[str, str]:
    """Build an ordered mapping  host_abs_dir → /mnt/<role>  from {role: host_dir}.

    Identical host directories that appear under multiple roles share a single
    container mount (the first role name encountered wins).

    Returns: {"/cluster/data/t1": "/mnt/t1-images", ...}
    """
    seen_hosts: Dict[str, str] = {}   # host_abs_path -> container_path
    used_containers: set = set()
    for role, host_dir in named_dirs.items():
        abs_dir = str(Path(host_dir).resolve())
        if abs_dir in seen_hosts:
            continue
        candidate = f"/mnt/{role}"
        # Avoid container-path collisions (unlikely but safe)
        counter = 2
        while candidate in used_containers:
            candidate = f"/mnt/{role}-{counter}"
            counter += 1
        used_containers.add(candidate)
        seen_hosts[abs_dir] = candidate
    return seen_hosts


def _translate_path(path: str, bind_map: Dict[str, str]) -> str:
    """Translate an absolute host path to its container-mounted equivalent.

    Uses longest-prefix match so nested directories resolve correctly.
    Returns the original path unchanged if no bind covers it.
    """
    for host_dir, container_dir in sorted(bind_map.items(), key=lambda x: -len(x[0])):
        if path == host_dir:
            return container_dir
        if path.startswith(host_dir + "/"):
            return container_dir + path[len(host_dir):]
    return path


class ZScoreProcessor:
    """Processor to run the z-score Snakemake workflow."""

    def __init__(self):
        self.current_dir = Path(__file__).parent
        self.leukoquant_dir = self.current_dir.parent
        self.leukoquant_parent_dir = self.leukoquant_dir.parent

        self.workflow_dir = self.leukoquant_dir / "workflow" / "workflows"

    def run_zscore(self,
                   healthy_subjects_list: str,
                   metrics: Dict[str, str],
                   output_dir: str,
                   t1_pattern: str,
                   demographics_csv: str,
                   target_subjects_list: Optional[str] = None,
                   covariates: Optional[List[str]] = None,
                   polynomial_terms: Optional[List[str]] = None,
                   metric_space: str = 't1',
                   output_space: str = 't1',
                   dwi_pattern: Optional[str] = None,
                   bval_pattern: Optional[str] = None,
                   skip_skullstrip_t1: bool = False,
                   skip_skullstrip_dwi: bool = False,
                   scheduler: str = "local",
                   cores: int = 1,
                   task_concurrency: Optional[int] = None,
                   force_rules: Optional[List[str]] = None,
                   verbose: bool = False) -> Tuple[str, str]:
        """Run the z-score Snakemake workflow.

        Args:
            healthy_subjects_list: path to file listing healthy cohort subjects
            target_subjects_list:  path to file listing target subjects. If not
                                    provided, every subject with a T1 matching
                                    `t1_pattern` is used as a target, minus the
                                    healthy cohort.
            metrics:               mapping of metric name -> {subject} pattern
            output_dir:            workflow output directory
            t1_pattern:            {subject} pattern to locate subject T1s
            demographics_csv:      CSV file with demographics for GLM
            covariates:            list of covariate column names (optional)
            polynomial_terms:      list of polynomial terms with optional power notation (default power is 2, e.g., ["age", "age:2", "bmi^3"], optional)
            metric_space:          't1' or 'dwi' (default: 't1')
            output_space:          't1' or 'dwi' (default: 't1')
            dwi_pattern:           {subject} pattern to locate DWI images (required if metric_space='dwi' or output_space='dwi')
            bval_pattern:          {subject} pattern to locate bval files (optional)
            skip_skullstrip_t1:    skip skull stripping for target/healthy T1s
            skip_skullstrip_dwi:   skip skull stripping for DWI b0 inputs
            scheduler:             "local" or "sge"
            cores:                 number of cores for Snakemake
            task_concurrency:      max concurrently-running SGE array tasks (-tc).
                                    Each task's container writes to the shared NFS
                                    output tree on its very first mkdir; submitting
                                    hundreds of tasks unthrottled under
                                    --immediate-submit can overwhelm the NFS server
                                    and cause a burst of transient write failures.
                                    Defaults to 20 if not set (see z_score_workflow.smk).

        Returns:
            (job_id, results_dir)
        """
        wf_file = self.workflow_dir / "z_score_workflow.smk"


        out_path = Path(output_dir)
        out_path.mkdir(parents=True, exist_ok=True)

        # ------------------------------------------------------------------
        # Resolve subject lists (auto-discover targets if not given), then
        # validate every subject has exactly one matching file on host
        # ------------------------------------------------------------------
        healthy_subjects = _read_subject_list(healthy_subjects_list)

        if target_subjects_list:
            target_subjects = _read_subject_list(target_subjects_list)
        else:
            logger.info(
                "No target subject list provided -- auto-discovering every "
                "subject with a T1 matching the given pattern, excluding the "
                "healthy cohort."
            )
            discovered = _discover_subjects_from_pattern(t1_pattern)
            target_subjects = sorted(set(discovered) - set(healthy_subjects))
            if not target_subjects:
                raise ValueError(
                    f"Auto-discovery found no target subjects (pattern: {t1_pattern}, "
                    f"{len(discovered)} total subjects found, {len(healthy_subjects)} "
                    "excluded as healthy)."
                )
            target_subjects_list = str(out_path / "auto_discovered_targets.txt")
            with open(target_subjects_list, "w") as f:
                f.write("\n".join(target_subjects) + "\n")
            logger.info(
                f"Auto-discovered {len(target_subjects)} target subjects "
                f"(from {len(discovered)} total, {len(healthy_subjects)} excluded as "
                f"healthy) -> {target_subjects_list}"
            )

        all_subjects = sorted(set(healthy_subjects + target_subjects))

        # Drop any requested metric that isn't available for every subject in
        # this run, rather than crashing (e.g. an older NODDI fit that only
        # produced odi.nii.gz, never ndi/fwf).
        metrics = _drop_metrics_missing_for_any_subject(metrics, all_subjects)

        _validate_glob_matches(t1_pattern, all_subjects, "T1")
        for metric_name, metric_pattern in metrics.items():
            _validate_glob_matches(metric_pattern, all_subjects, f"metric '{metric_name}'")
        
        # Validate DWI patterns if DWI inputs are required
        if metric_space == 'dwi' or output_space == 'dwi':
            if not dwi_pattern:
                raise ValueError("dwi_pattern is required when metric_space='dwi' or output_space='dwi'")
            dwi_subjects = all_subjects if metric_space == 'dwi' else target_subjects
            _validate_glob_matches(dwi_pattern, dwi_subjects, "DWI")
            if bval_pattern:
                _validate_glob_matches(bval_pattern, dwi_subjects, "bval")

        # ------------------------------------------------------------------
        # Build Singularity bind map
        # named_dirs: role -> absolute host directory
        # Roles sharing the same host directory are collapsed into one mount.
        # ------------------------------------------------------------------
        named_dirs: Dict[str, str] = {
            "healthy-list": str(Path(healthy_subjects_list).parent.resolve()),
            "target-list":  str(Path(target_subjects_list).parent.resolve()),
            "t1-images":    _base_dir_from_pattern(t1_pattern),
            "demographics": str(Path(demographics_csv).parent.resolve()),
            "output":       str(out_path.resolve()),
        }
        for name, pattern in metrics.items():
            named_dirs[f"metric-{name.lower()}"] = _base_dir_from_pattern(pattern)
        
        # Add DWI and bval directories if needed
        if (metric_space == 'dwi' or output_space == 'dwi') and dwi_pattern:
            named_dirs["dwi-images"] = _base_dir_from_pattern(dwi_pattern)
        if bval_pattern:
            named_dirs["bval-files"] = _base_dir_from_pattern(bval_pattern)

        logger.info(f"Named directories for Singularity binds: {named_dirs}")

        bind_map = _build_bind_map(named_dirs)
        # Always mount the leukoquant repo at /leukoquant inside the container
        bind_map[str(self.leukoquant_parent_dir.resolve())] = "/leukoquant"

        # ------------------------------------------------------------------
        # Ensure container is present (downloaded from HuggingFace if missing)
        # ------------------------------------------------------------------
        containers_folder = self.leukoquant_dir / "workflow" / "containers"
        sif_path = str(containers_folder / FREESURFER_SIF_FILENAME)
        ensure_container(sif_path, filename=FREESURFER_SIF_HF_PATH)

        # ------------------------------------------------------------------
        # Write Snakemake config YAML
        # Host paths are used by Snakemake for dependency tracking.
        # singularity_binds lets the workflow translate any host path to its
        # container-mounted equivalent at runtime (inside params: lambdas).
        # ------------------------------------------------------------------
        config_file = str(out_path.resolve() / "z_score_config.yaml")
        config_dict: Dict[str, object] = {
            "healthy_subjects_list": str(Path(healthy_subjects_list).resolve()),
            "target_subjects_list":  str(Path(target_subjects_list).resolve()),
            "output_dir":            str(out_path.resolve()),
            "demographics_csv":      str(Path(demographics_csv).resolve()),
            "t1_pattern":            t1_pattern,
            "metrics":               metrics,
            "metric_space":          metric_space,
            "output_space":          output_space,
            "skip_skullstrip_t1":    bool(skip_skullstrip_t1),
            "skip_skullstrip_dwi":   bool(skip_skullstrip_dwi),
            "singularity_binds":     bind_map,
            "leukoquant_parent_dir":      str(self.leukoquant_parent_dir.resolve()),
            "container_name":        "freesurfer_unified_container",
            "conda_env_name":        "z_scores_env",
        }
        if dwi_pattern:
            config_dict["dwi_pattern"] = dwi_pattern
        if bval_pattern:
            config_dict["bval_pattern"] = bval_pattern
        if covariates:
            config_dict["covariates"] = covariates
        if polynomial_terms:
            config_dict["polynomial_terms"] = polynomial_terms
        if task_concurrency is not None:
            config_dict["task_concurrency"] = int(task_concurrency)

        with open(config_file, "w") as f:
            yaml.safe_dump(config_dict, f, sort_keys=False)

        # ------------------------------------------------------------------
        # Build Snakemake command with auto-generated singularity bind args
        # ------------------------------------------------------------------
        bind_str = ",".join(f"{host}:{container}" for host, container in bind_map.items())
        singularity_bind = f"--bind {bind_str}"

        snakemake_cmd = [
            "snakemake",
            "-s", str(wf_file),
            "--directory", str(out_path),
            "--cores", str(cores),
            "--jobs", "unlimited",
            "--software-deployment-method", "apptainer",
            "--apptainer-args", singularity_bind,
            f"--configfile={config_file}",
        ]

        if scheduler == "sge":
            logger.info("Using SGE scheduler for z-score")
            snakemake_cmd.extend([
                "--immediate-submit",
                "--notemp",
                # Cold NFS mounts on some compute nodes exceed the default 5s.
                "--latency-wait", "60",
                "--max-jobs-per-timespan", "75000/1s",
                "--executor", "sge",
            ])

        add_forcerun_args(snakemake_cmd, force_rules)
        snakemake_cmd.append("all")

        #logger.info(f"Singularity binds: {bind_str}")

        env = os.environ.copy()
        # Isolate the Snakemake source cache per output directory to avoid NFS
        # race conditions when multiple instances run in parallel on a shared filesystem.
        env["SNAKEMAKE_SOURCECACHE_PATH"] = str(out_path / ".snakemake" / "source-cache")
        try:
            result = subprocess.run(
                    snakemake_cmd,
                    check=True,
                    cwd=str(self.workflow_dir),
                    env=env,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
            )
            if verbose and result.stdout:
                    print(result.stdout)
            return ("0", str(out_path.resolve()))
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"Z-score workflow failed: {e}\nError Details:\n{e.stderr}"
            )


def apply_zscore(
    healthy_subjects_list: Optional[str] = None,
    target_subjects_list: Optional[str] = None,
    metrics: Optional[Dict[str, str]] = None,
    output_dir: Optional[str] = None,
    t1_pattern: Optional[str] = None,
    demographics_csv: Optional[str] = None,
    covariates: Optional[str] = None,
    polynomial_terms: Optional[str] = None,
    metric_space: str = 't1',
    output_space: str = 't1',
    dwi_pattern: Optional[str] = None,
    bval_pattern: Optional[str] = None,
    skip_skullstrip_t1: bool = False,
    skip_skullstrip_dwi: bool = False,
    scheduler: Optional[str] = None,
    cores: Optional[int] = None,
    task_concurrency: Optional[int] = None,
    force_rules: Optional[List[str]] = None,
    verbose: bool = False,
    config_yaml: Optional[str] = None,
) -> dict:
    """Run the Z-score workflow and return a summary dict.

    CLI / caller arguments take priority over values in ``config_yaml``
    when both are provided.
    """
    if verbose:
        logging.basicConfig(level=logging.INFO)

    # --- Load YAML config (provides defaults; caller args override) --------
    yaml_cfg = load_yaml_config(config_yaml)

    # Each line: try the caller value first, then one or more YAML keys.
    healthy_subjects_list = first_truthy(healthy_subjects_list,
                                   yaml_cfg.get("healthy_subjects_list"),
                                   yaml_cfg.get("healthy_list"))
    target_subjects_list  = first_truthy(target_subjects_list,
                                   yaml_cfg.get("target_subjects_list"),
                                   yaml_cfg.get("target_list"))
    output_dir            = first_truthy(output_dir,       yaml_cfg.get("output_dir"))
    demographics_csv      = first_truthy(demographics_csv, yaml_cfg.get("demographics_csv"))
    t1_pattern_raw        = first_truthy(t1_pattern,       yaml_cfg.get("t1_pattern"),
                                   yaml_cfg.get("t1_path"))
    dwi_pattern_raw       = first_truthy(dwi_pattern,      yaml_cfg.get("dwi_pattern"),
                                   yaml_cfg.get("dwi_path"))
    bval_pattern_raw      = first_truthy(bval_pattern,     yaml_cfg.get("bval_pattern"),
                                   yaml_cfg.get("bval_path"))
    metric_space          = first_truthy(metric_space,     yaml_cfg.get("metric_space")) or 't1'
    output_space          = first_truthy(output_space,     yaml_cfg.get("output_space")) or 't1'

    # CLI flags are boolean (default False). YAML can also set defaults; CLI True wins.
    skip_skullstrip_t1_val = bool(yaml_cfg.get("skip_skullstrip_t1", False)) or bool(skip_skullstrip_t1)
    skip_skullstrip_dwi_val = bool(yaml_cfg.get("skip_skullstrip_dwi", False)) or bool(skip_skullstrip_dwi)

    # metrics: caller dict wins; fall back to YAML
    if not metrics:
        metrics = _parse_metrics_spec(yaml_cfg.get("metrics")) or None

    # scheduler / cores: caller wins; fall back to YAML; then hard defaults
    final_scheduler = first_truthy(scheduler, yaml_cfg.get("scheduler")) or "local"
    cores_raw       = first_truthy(cores,     yaml_cfg.get("cores"))
    final_cores     = int(cores_raw) if cores_raw is not None else 1
    task_conc_raw   = first_truthy(task_concurrency, yaml_cfg.get("task_concurrency"))
    final_task_concurrency = int(task_conc_raw) if task_conc_raw is not None else None

    # --- Validate required parameters -------------------------------------
    # target_subjects_list is deliberately NOT required here -- when omitted,
    # ZScoreProcessor.run_zscore() auto-discovers targets as every subject
    # matching t1_pattern minus the healthy cohort.
    missing = [
        name for name, val in [
            ("healthy_subjects_list", healthy_subjects_list),
            ("output_dir",            output_dir),
            ("t1_pattern",            t1_pattern_raw),
            ("demographics_csv",      demographics_csv),
            ("metrics",               metrics),
        ] if not val
    ]
    if missing:
        raise ValueError(f"Missing required parameters: {', '.join(missing)}")

    # Validate DWI pattern if DWI inputs are required
    if (metric_space == 'dwi' or output_space == 'dwi') and not dwi_pattern_raw:
        raise ValueError("--dwi-pattern is required when --metric-space=dwi or --output-space=dwi")

    # --- Normalise patterns and lists ------------------------------------
    t1_final      = _normalize_subject_pattern(str(t1_pattern_raw), "T1")
    metrics_final = {
        str(name): _normalize_subject_pattern(str(pattern), f"metric '{name}'")
        for name, pattern in (metrics or {}).items()
    }
    dwi_final = None
    if dwi_pattern_raw:
        dwi_final = _normalize_subject_pattern(str(dwi_pattern_raw), "DWI")
    bval_final = None
    if bval_pattern_raw:
        bval_final = _normalize_subject_pattern(str(bval_pattern_raw), "bval")

    # covariates / poly_terms: caller value first, then YAML
    covariates_raw  = covariates       if covariates       is not None else yaml_cfg.get("covariates")
    poly_terms_raw  = polynomial_terms if polynomial_terms is not None else yaml_cfg.get("polynomial_terms")
    covariates_list = _normalize_list_value(covariates_raw)
    poly_terms_list = _parse_polynomial_terms(poly_terms_raw)

    proc = ZScoreProcessor()
    try:
        job_id, results_dir = proc.run_zscore(
            healthy_subjects_list=str(healthy_subjects_list),
            target_subjects_list=str(target_subjects_list) if target_subjects_list else None,
            metrics=metrics_final,
            output_dir=str(output_dir),
            t1_pattern=t1_final,
            demographics_csv=str(demographics_csv),
            covariates=covariates_list,
            polynomial_terms=poly_terms_list,
            metric_space=metric_space,
            output_space=output_space,
            dwi_pattern=dwi_final,
            bval_pattern=bval_final,
            skip_skullstrip_t1=bool(skip_skullstrip_t1_val),
            skip_skullstrip_dwi=bool(skip_skullstrip_dwi_val),
            scheduler=final_scheduler,
            cores=final_cores,
            task_concurrency=final_task_concurrency,
            force_rules=force_rules,
            verbose=verbose,
        )
        print("✅ Z-score job submitted successfully")
        return {
            "success":     True,
            "job_id":      job_id,
            "results_dir": results_dir,
            "output_dir":  str(Path(str(output_dir)).resolve()),
            "scheduler":   final_scheduler,
        }
    except Exception as e:
        print(f"❌ Z-score job failed:\n{e}", file=sys.stderr)
        return {
            "success":    False,
            "error":      str(e),
            "output_dir": str(Path(str(output_dir)).resolve()),
        }
