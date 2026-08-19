#!/usr/bin/env python3
"""
Tractography QC processor for leukoquant.

Runs the `tract_qc_workflow.smk` Snakemake workflow to perform QC on tractography outputs for a subject or batch of subjects.
"""
import os
from pathlib import Path
from typing import Optional, Dict, List, Tuple
import logging
import subprocess
import sys
import yaml

from leukoquant.utils.subject_utils import read_subjects
from leukoquant.utils.external_utils import check_sge_plugin
from leukoquant.utils.container_utils import (
    ensure_container,
    MINICONDA_SIF_FILENAME,
    MINICONDA_SIF_HF_PATH,
)
from leukoquant.utils.snakemake_utils import add_forcerun_args, load_yaml_config, first_truthy

logger = logging.getLogger(__name__)

class TractQCProcessor:
    """Processor to run the tract QC Snakemake workflow."""
    def __init__(self):
        self.current_dir = Path(__file__).parent
        self.leukoquant_dir = self.current_dir.parent
        self.leukoquant_parent_dir = self.leukoquant_dir.parent
        self.workflow_dir = self.leukoquant_dir / "workflow" / "workflows"

    def run_tract_qc(self,
                     subject_input: str,
                     tractography_path: str,
                     output_dir: str,
                     scheduler: str = "local",
                     cores: int = 1,
                     force_rules: Optional[List[str]] = None,
                     verbose: bool = False) -> Tuple[str, str]:
        """Run the tract QC Snakemake workflow.

        Args:
            subject_input: Single subject ID or path to file listing subjects
            tractography_path: Path to tractography root directory
            output_dir: Workflow output directory
            scheduler: "local" or "sge"
            cores: Number of cores for Snakemake
            verbose: Enable verbose output

        Returns:
            (job_id, results_dir)
        """

        wf_file = self.workflow_dir / "tract_qc_workflow.smk"
        if os.path.isdir(output_dir):
            out_path = Path(output_dir)
        else:
            out_path = Path(output_dir).parent
        out_path.mkdir(parents=True, exist_ok=True)

        # Determine subjects (single or batch)
        subjects = read_subjects(subject_input)
        subjects_file = str(Path(subject_input).resolve()) if Path(subject_input).is_file() else None
        if not subjects:
            raise ValueError("No subjects provided")

        # Handle colon-separated tractography_path (base:glob:space)
        def _resolve_pattern_spec(pattern_spec: str) -> str:
            if ":" not in pattern_spec:
                return str(Path(pattern_spec).resolve())
            base, glob = pattern_spec.split(":", 1)
            return f"{Path(base).resolve()}:{glob}"

        tractography_base = tractography_path.split(":")[0] if ":" in tractography_path else tractography_path
        named_dirs = {
            "tractography": str(Path(tractography_base).resolve()),
            "output": str(out_path.resolve()),
        }
        # Always mount the leukoquant repo at /leukoquant inside the container.
        repo_host = str(self.leukoquant_parent_dir.resolve())
        bind_pairs = [
            (named_dirs["tractography"], "/mnt/tractography"),
            (named_dirs["output"], "/mnt/output"),
            (repo_host, "/leukoquant"),
        ]
        bind_map = {
            named_dirs["tractography"]: "/mnt/tractography",
            named_dirs["output"]: "/mnt/output",
            repo_host: "/leukoquant",
        }

        # Ensure container is present (downloaded from HuggingFace if missing)
        containers_folder = self.leukoquant_dir / "workflow" / "containers"
        sif_path = str(containers_folder / MINICONDA_SIF_FILENAME)
        ensure_container(sif_path, filename=MINICONDA_SIF_HF_PATH)

        # Build config dict
        config_file = str(out_path.resolve() / "tract_qc_config.yaml")
        config_dict = {
            "subjects": subjects,
            "tractography_path": _resolve_pattern_spec(tractography_path),
            "output_dir": str(out_path.resolve()),
            "singularity_binds": bind_map,
            "leukoquant_parent_dir": repo_host,
            "container_name": "miniconda_unified_container",
            "conda_env_name": "metrics_env",
        }
        if subjects_file:
            config_dict["subjects_file"] = subjects_file

        with open(config_file, "w") as f:
            yaml.safe_dump(config_dict, f)

        # Build Snakemake command with singularity binds
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
            logger.info("Using SGE scheduler for tract_qc")
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
                f"TractQC workflow failed: {e}\nError Details:\n{e.stderr}"
            )

def apply_tract_qc(
    subject: Optional[str] = None,
    tractography_path: Optional[str] = None,
    output_dir: Optional[str] = None,
    scheduler: str = "local",
    cores: int = 1,
    force_rules: Optional[List[str]] = None,
    verbose: bool = False,
    config_yaml: Optional[str] = None,
) -> dict:
    """Run the tract QC workflow and return a summary dict.

    CLI / caller arguments take priority over values in ``config_yaml`` when
    both are provided.
    """
    yaml_cfg = load_yaml_config(config_yaml)
    subject           = first_truthy(subject,           yaml_cfg.get("subject"))
    tractography_path = first_truthy(tractography_path, yaml_cfg.get("tractography_path"))
    output_dir        = first_truthy(output_dir,        yaml_cfg.get("output_dir"))
    scheduler         = first_truthy(scheduler,         yaml_cfg.get("scheduler")) or "local"
    cores_raw         = first_truthy(cores,             yaml_cfg.get("cores"))
    cores             = int(cores_raw) if cores_raw is not None else 1

    if not subject:
        raise ValueError("subject is required (via argument or config_yaml)")
    if not tractography_path:
        raise ValueError("tractography_path is required (via argument or config_yaml)")
    if not output_dir:
        raise ValueError("output_dir is required (via argument or config_yaml)")

    try:
        proc = TractQCProcessor()
        job_id, results_dir = proc.run_tract_qc(
            subject_input=subject,
            tractography_path=tractography_path,
            output_dir=output_dir,
            scheduler=scheduler,
            cores=cores,
            force_rules=force_rules,
            verbose=verbose,
        )
        print("✅ Tract Quality Control job submitted successfully")
        return {"success": True, "job_id": job_id, "results_dir": results_dir}
    except Exception as e:
        print(f"❌ Tract Quality Control job failed:\n{e}", file=sys.stderr)
        return {"success": False, "error": str(e), "subject": subject, "output_dir": output_dir}
