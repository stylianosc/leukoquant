#!/usr/bin/env python
"""
Metrics processor for leukoquant.

Runs the `metrics_workflow.smk` Snakemake workflow to compute metrics
along tracts, lesions, and other ROIs for a subject or batch of subjects.
"""

import glob
import logging
import subprocess
import sys
import yaml
from pathlib import Path
import os
from typing import Dict, List, Optional, Tuple

from leukoquant.utils.subject_utils import read_subjects
from leukoquant.utils.snakemake_utils import add_forcerun_args, load_yaml_config, first_truthy
from leukoquant.utils.container_utils import (
    ensure_container,
    FREESURFER_SIF_FILENAME,
    FREESURFER_SIF_HF_PATH,
)
from leukoquant.utils.external_utils import _resolve_fs_license, _resolve_gif_home

logger = logging.getLogger(__name__)


def _base_dir_from_pattern(pattern: str) -> str:
    """Return the absolute host base directory – the path component before the glob."""
    # Pattern format: /path/to/base:/glob/pattern
    if ":" not in pattern:
        return str(Path(pattern).parent.resolve())
    
    base_part = pattern.split(":")[0]
    return str(Path(base_part).resolve())


def _resolve_pattern_spec(pattern_spec: str) -> str:
    """Resolve only the base path in a base:glob[:space] spec.

    Keeps the suffix untouched while making the base absolute.
    """
    if ":" not in pattern_spec:
        return str(Path(pattern_spec).resolve())

    base, rest = pattern_spec.split(":", 1)
    return f"{Path(base).resolve()}:{rest}"


def _parse_lesion_entries(spec: str) -> List[Tuple[str, str, str, str]]:
    """Parse potentially named, comma-separated lesion spec.

    Format per entry: ``[name=]base:glob[:space]``
    ``name`` is separated by ``=`` (optional, default ``lesion``).
    ``space`` is the third ``:``-separated component (optional, default ``t1``).

    Returns a list of ``(name, abs_base, glob_pat, space)`` tuples.
    """
    results: List[Tuple[str, str, str, str]] = []
    for token in [t.strip() for t in spec.split(",") if t.strip()]:
        if "=" in token:
            name, path_spec = token.split("=", 1)
            name = name.strip()
            path_spec = path_spec.strip()
        else:
            name = "lesion"
            path_spec = token
        parts = path_spec.split(":")
        base = parts[0].strip()
        glob_pat = parts[1].strip() if len(parts) > 1 else ""
        space = parts[2].strip() if len(parts) > 2 else "t1"
        if not base:
            raise ValueError(f"Lesion spec '{token}': base cannot be empty.")
        results.append((name, str(Path(base).resolve()), glob_pat, space))
    return results


def _resolve_lesion_spec(spec: str) -> str:
    """Resolve base paths in a named, multi-entry lesion spec.

    Preserves name prefixes and appends path suffixes unchanged; only
    the base directory component is made absolute.
    """
    resolved_tokens: List[str] = []
    for token in [t.strip() for t in spec.split(",") if t.strip()]:
        if "=" in token:
            name, path_spec = token.split("=", 1)
            base, _, rest = path_spec.partition(":")
            resolved_path = f"{Path(base.strip()).resolve()}:{rest}" if rest else str(Path(base.strip()).resolve())
            resolved_tokens.append(f"{name.strip()}={resolved_path}")
        else:
            resolved_tokens.append(_resolve_pattern_spec(token))
    return ",".join(resolved_tokens)


def _parse_metric_spec(metric_name: str, spec: str, default_space: str = "dwi") -> Tuple[str, str, str]:
    """Parse a metric spec of the form 'base:glob:space' and return (base, glob, space).
    
    Args:
        metric_name: Name of the metric (for error messages)
        spec: Specification string in format 'base:glob:space' or 'base:glob'
    
    Returns:
        Tuple of (base_dir, glob_pattern, space)
        - base_dir: absolute path to the base directory
        - glob_pattern: glob pattern relative to base
        - space: space designation ('dwi', 't1', 'atlas', or default)
    """
    if ":" not in spec:
        raise ValueError(f"Metric '{metric_name}' spec must use format 'base:glob' or 'base:glob:space'")
    
    parts = spec.split(":")
    if len(parts) < 2:
        raise ValueError(f"Metric '{metric_name}' spec must have at least base and glob parts")
    
    base = parts[0].strip()
    glob_pat = parts[1].strip()
    space = parts[2].strip() if len(parts) > 2 else default_space
    
    if not base or not glob_pat:
        raise ValueError(f"Metric '{metric_name}': base and glob parts cannot be empty")
    
    # Validate space
    valid_spaces = {"dwi", "t1", "atlas"}
    if space not in valid_spaces:
        raise ValueError(
            f"Metric '{metric_name}': invalid space '{space}'. Must be one of: {valid_spaces}"
        )
    
    abs_base = str(Path(base).resolve())
    return abs_base, glob_pat, space


def _build_bind_plan(named_dirs: Dict[str, str]) -> Tuple[Dict[str, str], Dict[str, str], List[Tuple[str, str]]]:
    """Build bind structures from {role: host_dir}.

    Returns:
    - default_bind_map: host_abs_dir -> first container mount (for generic translate_path)
    - role_mounts: role -> container mount (always unique per role)
    - bind_pairs: list of (host_abs_dir, container_mount), includes duplicate hosts
      so different roles can mount the same host at different container paths.
    """
    default_bind_map: Dict[str, str] = {}
    role_mounts: Dict[str, str] = {}
    bind_pairs: List[Tuple[str, str]] = []
    used_containers: set = set()

    for role, host_dir in named_dirs.items():
        abs_dir = str(Path(host_dir).resolve())

        candidate = f"/mnt/{role}"
        counter = 2
        while candidate in used_containers:
            candidate = f"/mnt/{role}_{counter}"
            counter += 1

        used_containers.add(candidate)
        role_mounts[role] = candidate
        bind_pairs.append((abs_dir, candidate))

        # Keep first mapping for generic path translation logic.
        if abs_dir not in default_bind_map:
            default_bind_map[abs_dir] = candidate

    return default_bind_map, role_mounts, bind_pairs


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


class MetricsProcessor:
    """Processor to run the metrics Snakemake workflow."""

    def __init__(self):
        self.current_dir = Path(__file__).parent
        self.leukoquant_dir = self.current_dir.parent
        self.leukoquant_parent_dir = self.leukoquant_dir.parent
        self.workflow_dir = self.leukoquant_dir / "workflow" / "workflows"

    def run_metrics(self,
                    subject_input: str,
                    tractography_path: str,
                    t1_path: Optional[str],
                    dwi_path: Optional[str],
                    lesion_path: Optional[str],
                    metrics: Optional[Dict[str, str]],
                    tract_mode: str,
                    output_dir: str,
                    scheduler: str = "local",
                    cores: int = 1,
                    keep_intermediate: bool = False,
                    parcellation: str = "freesurfer",
                    force_rules: Optional[List[str]] = None,
                    verbose: bool = False) -> Tuple[str, str]:
        """Run the metrics Snakemake workflow.

        Args:
            subject_input: Single subject ID or path to file listing subjects
            tractography_path: Path to tractography root directory
            dwi_path: DWI pattern in base:glob:space format (optional)
            lesion_path: Lesion pattern in base:glob format (optional)
            metrics: Dict of metric_name -> base:glob:space patterns
            output_dir: Workflow output directory
            scheduler: "local" or "sge"
            cores: Number of cores for Snakemake
            force_rules: Rule name(s) to force-rerun (e.g. ["extract_metrics_gif"])
                even if their outputs already exist and are up to date. Anything
                upstream that's already complete is left untouched.
            verbose: Enable verbose output

        Returns:
            (job_id, results_dir)
        """
        wf_file = self.workflow_dir / "metrics_workflow.smk"

        out_path = Path(output_dir)
        out_path.mkdir(parents=True, exist_ok=True)

        # ------------------------------------------------------------------
        # Determine subjects (single or batch)
        # ------------------------------------------------------------------
        subjects = read_subjects(subject_input)
        subjects_file = str(Path(subject_input).resolve()) if Path(subject_input).is_file() else None

        if not subjects:
            raise ValueError("No subjects provided")

        # ------------------------------------------------------------------
        # Build Singularity bind map
        # named_dirs: role -> absolute host directory
        # ------------------------------------------------------------------
        # Extract base directories from patterns (handle base:glob:space format)
        tractography_base = tractography_path.split(":")[0] if ":" in tractography_path else tractography_path
        
        named_dirs: Dict[str, str] = {
            "tractography": str(Path(tractography_base).resolve()),
            "output": str(out_path.resolve()),
        }

        # Add T1 path if provided
        if t1_path:
            t1_base, _, _ = _parse_metric_spec("t1", t1_path, default_space="t1")
            named_dirs["t1"] = t1_base

        # Add DWI path if provided
        if dwi_path:
            dwi_base, _, _ = _parse_metric_spec("dwi", dwi_path)
            named_dirs["dwi"] = dwi_base

        # Add lesion path(s) if provided - supports named multi-entry specs
        if lesion_path:
            for lesion_name, abs_base, _, _ in _parse_lesion_entries(lesion_path):
                role = f"lesion-{lesion_name}"
                named_dirs[role] = abs_base

        # Add metric directories
        if metrics:
            for metric_name, metric_spec in metrics.items():
                metric_base, _, _ = _parse_metric_spec(metric_name, metric_spec)
                named_dirs[f"metric-{metric_name}"] = metric_base

        #logger.info(f"Named directories for Singularity binds: {named_dirs}")

        bind_map, role_mounts, bind_pairs = _build_bind_plan(named_dirs)

        # Always mount the leukoquant repo at /leukoquant inside the container.
        repo_host = str(self.leukoquant_parent_dir.resolve())
        bind_map[repo_host] = "/leukoquant"
        bind_pairs.append((repo_host, "/leukoquant"))

        fs_license_dir = _resolve_fs_license(self.leukoquant_dir, verbose=verbose)
        fs_license_singularity = "/license/license.txt"
        host_fs_license_dir = str(fs_license_dir.resolve())
        bind_map[host_fs_license_dir] = fs_license_singularity
        bind_pairs.append((host_fs_license_dir, fs_license_singularity))

        # Bind scratch space
        scratch_host = Path("/scratch0") if Path("/scratch0").exists() else out_path / ".scratch"
        scratch_host.mkdir(parents=True, exist_ok=True)
        scratch_host_abs = str(scratch_host.resolve())
        scratch_dir_singularity = "/scratch0"
        
        bind_map[scratch_host_abs] = scratch_dir_singularity
        bind_pairs.append((scratch_host_abs, scratch_dir_singularity))

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
        # container-mounted equivalent at runtime.
        # ------------------------------------------------------------------
        config_file = str(out_path.resolve() / "metrics_config.yaml")
        config_dict: Dict[str, object] = {
            "subjects": subjects,
            "tractography_path": _resolve_pattern_spec(tractography_path),
            "tract_mode": tract_mode,
            "output_dir": str(out_path.resolve()),
            "singularity_binds": bind_map,
            "metric_mounts": {k.replace("metric-", "", 1): v for k, v in role_mounts.items() if k.startswith("metric-")},
            "keep_intermediate": keep_intermediate,
            "leukoquant_parent_dir": str(self.leukoquant_parent_dir.resolve()),
            "container_name": "freesurfer_unified_container",
            "conda_env_name": "metrics_env",
            "verbose": verbose,
            "parcellations": [p.strip() for p in parcellation.split(",") if p.strip()],
        }

        if subjects_file:
            config_dict["subjects_file"] = subjects_file
        
        if t1_path:
            config_dict["t1_path"] = _resolve_pattern_spec(t1_path)

        if dwi_path:
            config_dict["dwi_path"] = _resolve_pattern_spec(dwi_path)
        
        if lesion_path:
            config_dict["lesion_path"] = _resolve_lesion_spec(lesion_path)

        if metrics:
            # Resolve metric paths to absolute paths while preserving glob patterns and space
            resolved_metrics = {}
            for metric_name, metric_spec in metrics.items():
                resolved_metrics[metric_name] = _resolve_pattern_spec(metric_spec)
            config_dict["metrics"] = resolved_metrics

        with open(config_file, "w") as f:
            yaml.dump(config_dict, f, default_flow_style=False)

        #logger.info(f"Config written to {config_file}")

        # ------------------------------------------------------------------
        # Build Snakemake command with auto-generated singularity bind args
        # ------------------------------------------------------------------
        bind_str = ",".join(f"{host}:{container}" for host, container in bind_pairs)
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
            logger.info("Using SGE scheduler for metrics")
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
                f"Metrics workflow failed: {e}\nError Details:\n{e.stderr}"
            )


def apply_metrics(
    subject: Optional[str] = None,
    tractography_path: Optional[str] = None,
    t1_path: Optional[str] = None,
    dwi_path: Optional[str] = None,
    lesion_path: Optional[str] = None,
    metrics: Optional[Dict[str, str]] = None,
    tract_mode: Optional[str] = None,
    output_dir: Optional[str] = None,
    scheduler: Optional[str] = None,
    cores: Optional[int] = None,
    keep_intermediate: bool = False,
    parcellation: str = "freesurfer",
    force_rules: Optional[List[str]] = None,
    verbose: bool = False,
    config_yaml: Optional[str] = None,
) -> dict:
    """Run the metrics workflow and return a summary dict.

    CLI / caller arguments take priority over values in ``config_yaml`` when
    both are provided.

    Args:
        subject: Subject ID or path to file with subject IDs
        tractography_path: Path to tractography root directory
        dwi_path: DWI pattern (optional)
        lesion_path: Lesion pattern (optional)
        metrics: Dict of metric_name -> base:glob:space patterns
        tract_mode: Tractography mode passed through to the workflow config
        output_dir: Output directory
        scheduler: "local" or "sge"
        cores: Number of cores
        parcellation: Comma-separated parcellation(s) to compute metrics for
            (e.g. "gif", "freesurfer", or "gif,freesurfer"). Determines which
            extract_metrics_{parcellation} rule(s) the workflow builds.
        force_rules: Rule name(s) to force-rerun even if outputs exist
        verbose: Verbose output
        config_yaml: Optional path to a YAML file providing defaults for any
            of the above (CLI/caller arguments always win over this).

    Returns:
        Dictionary with success status and results
    """
    # --- Load YAML config (provides defaults; caller args override) --------
    yaml_cfg = load_yaml_config(config_yaml)
    subject           = first_truthy(subject,           yaml_cfg.get("subject"))
    tractography_path = first_truthy(tractography_path, yaml_cfg.get("tractography_path"))
    t1_path           = first_truthy(t1_path,           yaml_cfg.get("t1_path"))
    dwi_path          = first_truthy(dwi_path,          yaml_cfg.get("dwi_path"))
    lesion_path       = first_truthy(lesion_path,       yaml_cfg.get("lesion_path"))
    tract_mode        = first_truthy(tract_mode,        yaml_cfg.get("tract_mode")) or "tractography-atlas"
    output_dir        = first_truthy(output_dir,        yaml_cfg.get("output_dir"))
    parcellation      = first_truthy(parcellation,      yaml_cfg.get("parcellation")) or "freesurfer"
    keep_intermediate = bool(yaml_cfg.get("keep_intermediate", False)) or bool(keep_intermediate)
    if not metrics:
        metrics = yaml_cfg.get("metrics") or None

    if verbose:
        logging.basicConfig(level=logging.INFO)

    # Validate required parameters
    missing = [
        name for name, val in [
            ("subject", subject),
            ("tractography_path", tractography_path),
            ("output_dir", output_dir),
        ] if not val
    ]
    if missing:
        raise ValueError(f"Missing required parameters: {', '.join(missing)}")

    if not metrics:
        logger.warning("No metrics provided")

    # Apply defaults
    final_scheduler = first_truthy(scheduler, yaml_cfg.get("scheduler")) or "local"
    cores_raw = first_truthy(cores, yaml_cfg.get("cores"))
    final_cores = int(cores_raw) if cores_raw is not None else 1

    proc = MetricsProcessor()
    try:
        job_id, results_dir = proc.run_metrics(
            subject_input=subject,
            tractography_path=tractography_path,
            t1_path=t1_path,
            dwi_path=dwi_path,
            lesion_path=lesion_path,
            metrics=metrics,
            tract_mode=tract_mode,
            output_dir=output_dir,
            scheduler=final_scheduler,
            cores=final_cores,
            keep_intermediate=keep_intermediate,
            parcellation=parcellation,
            force_rules=force_rules,
            verbose=verbose,
        )
        print("✅ Metrics calculation job submitted successfully")
        return {
            "success": True,
            "job_id": job_id,
            "results_dir": results_dir,
            "output_dir": str(Path(output_dir).resolve()),
            "scheduler": final_scheduler,
        }
    except Exception as e:
        print(f"❌ Metrics calculation job failed:\n{e}", file=sys.stderr)
        return {
            "success": False,
            "error": str(e),
            "output_dir": str(Path(output_dir).resolve()) if output_dir else None,
        }
