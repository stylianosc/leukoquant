#!/usr/bin/env python
"""
FreeSurfer recon-all processing for leukoquant.

Provides a processor class to run FreeSurfer's `recon-all` via a Snakemake workflow.
Supports single or multiple subjects via subject ID or text file.
"""

import logging
import os
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Tuple

from leukoquant.utils.subject_utils import read_subjects, resolve_subject_pattern
from leukoquant.utils.bind_utils import dir_level_bind_files, consolidate_bind_entries
from leukoquant.utils.container_utils import (
    ensure_container,
    FREESURFER_SIF_FILENAME,
    FREESURFER_SIF_HF_PATH,
)
from leukoquant.utils.external_utils import _resolve_fs_license
from leukoquant.utils.snakemake_utils import add_forcerun_args, load_yaml_config, first_truthy

logger = logging.getLogger(__name__)

class FreeSurferProcessor:
    """Processor to run FreeSurfer `recon-all` through Snakemake."""

    def __init__(self, external_dir: Optional[str] = None):
        self.current_dir = Path(__file__).parent
        self.leukoquant_dir = self.current_dir.parent
        self.leukoquant_parent_dir = self.leukoquant_dir.parent

        if external_dir is None:
            external_dir = self.leukoquant_dir / "external"

        self.external_dir = Path(external_dir)
        self.workflow_dir = self.leukoquant_dir / "workflow" / "workflows"

    def run_recon_all(self,
                      subject_input: str,
                      t1_pattern: str,
                      output_dir: str = ".",
                      scheduler: str = "local",
                      cores: int = 1,
                      keep_intermediate: bool = False,
                      force_rules: Optional[List[str]] = None,
                      verbose: bool = False) -> Tuple[str, str]:
        """Run FreeSurfer `recon-all` for one or more subjects using Snakemake.

        Args:
            subject_input: bare subject ID or path to a text file (one ID per line)
            t1_pattern: explicit path or {subject} glob pattern for T1 files
            output_dir: directory to place workflow outputs and config
            scheduler: "local" or "sge" to use existing SGE profile
            cores: number of cores to request for Snakemake

        Returns:
            Tuple of (job_id, results_dir). job_id is "0" for local runs.
        """
        subjects = read_subjects(subject_input)
        t1_files = [resolve_subject_pattern(t1_pattern, s) for s in subjects]

        for f in t1_files:
            if not Path(f).exists():
                raise FileNotFoundError(f"T1 file not found: {f}")

        fs_license_dir = _resolve_fs_license(self.leukoquant_dir, verbose=verbose)

        wf_file = self.workflow_dir / "recon_all_workflow.smk"
        if not wf_file.exists():
            raise FileNotFoundError(f"Reconstruction workflow not found: {wf_file}")

        out_path = Path(output_dir)
        out_path.mkdir(parents=True, exist_ok=True)

        fs_license_singularity = Path("/license/license.txt")

        containers_folder = self.leukoquant_dir / "workflow" / "containers"

        # Ensure container is present (downloaded from HuggingFace if missing)
        sif_path = str(containers_folder / FREESURFER_SIF_FILENAME)
        ensure_container(sif_path, filename=FREESURFER_SIF_HF_PATH)

        # Directory-level binds + consolidation: keeps --bind short for large cohorts.
        t1_files_singularity, t1_binds = dir_level_bind_files(t1_files, "input")
        t1_binds, remap = consolidate_bind_entries(t1_binds)
        t1_files_singularity = [remap(p) for p in t1_files_singularity]

        # Setup scratch directory for intermediate files
        scratch_host = Path("/scratch0") if Path("/scratch0").exists() else out_path / ".scratch"
        scratch_host.mkdir(parents=True, exist_ok=True)
        scratch_dir_singularity = "/scratch0"

        config_file = str(out_path.absolute()) + "/recon_all_config.yaml"
        config_dict = {
            "subject_id": subjects,
            "t1_file": t1_files,
            "t1_file_singularity": t1_files_singularity,
            "output_dir": str(out_path.absolute()),
            "output_dir_singularity": "/output",
            "scratch_dir_singularity": scratch_dir_singularity,
            "keep_intermediate": keep_intermediate,
            "container_name": "freesurfer_unified_container",
            "leukoquant_parent_dir": str(self.leukoquant_parent_dir.absolute()),
        }
        import yaml as _yaml
        with open(config_file, "w") as f:
            _yaml.dump(config_dict, f, default_flow_style=False)

        bind_parts = (
            t1_binds
            + [f"{str(self.leukoquant_parent_dir.absolute())}:/leukoquant"]
            + [f"{str(out_path.absolute())}:/output"]
            + [f"{str(scratch_host.absolute())}:{scratch_dir_singularity}"]
            + [f"{str(fs_license_dir.absolute())}:{str(fs_license_singularity.absolute())}"]
        )
        singularity_bind = "--bind " + ",".join(bind_parts)

        snakemake_cmd = [
            "snakemake",
            "-s", str(wf_file),
            "--directory", str(out_path),
            "--cores", str(cores),
            "--jobs", "unlimited",
            "--immediate-submit", "--notemp",
            # Cold NFS mounts on some compute nodes exceed the default 5s.
            "--latency-wait", "60",
            "--software-deployment-method", "apptainer",
            "--apptainer-args", singularity_bind,
            f"--configfile={config_file}",
        ]

        if scheduler == "sge":
            snakemake_cmd.extend([
                "--max-jobs-per-timespan", "75000/1s",
                "--executor", "sge",
            ])
        else:
            logger.info("Using local scheduler")

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
            return "0", str(out_path)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"recon-all failed with exit code {e.returncode}\n{e.stderr}"
            )


def apply_recon_all(subject_input: Optional[str] = None,
                    t1_pattern: Optional[str] = None,
                    output_dir: Optional[str] = None,
                    scheduler: str = "local",
                    cores: int = 1,
                    keep_intermediate: bool = False,
                    force_rules: Optional[List[str]] = None,
                    verbose: bool = False,
                    config_yaml: Optional[str] = None) -> dict:
    """Run recon-all for one or more subjects and return a summary dict.

    CLI / caller arguments take priority over values in ``config_yaml`` when
    both are provided.
    """
    yaml_cfg = load_yaml_config(config_yaml)
    subject_input = first_truthy(subject_input, yaml_cfg.get("subject_input"), yaml_cfg.get("subject"))
    t1_pattern    = first_truthy(t1_pattern,    yaml_cfg.get("t1_pattern"), yaml_cfg.get("t1"))
    output_dir    = first_truthy(output_dir,    yaml_cfg.get("output_dir"))
    scheduler     = first_truthy(scheduler,     yaml_cfg.get("scheduler")) or "local"
    cores_raw     = first_truthy(cores,         yaml_cfg.get("cores"))
    cores         = int(cores_raw) if cores_raw is not None else 1
    keep_intermediate = bool(yaml_cfg.get("keep_intermediate", False)) or bool(keep_intermediate)

    if not subject_input:
        raise ValueError("subject_input is required (via argument or config_yaml)")
    if not t1_pattern:
        raise ValueError("t1_pattern is required (via argument or config_yaml)")
    if not output_dir:
        raise ValueError("output_dir is required (via argument or config_yaml)")

    if verbose:
        logging.basicConfig(level=logging.INFO)
    else:
        logging.basicConfig(level=logging.WARNING)

    proc = FreeSurferProcessor()
    try:
        job_id, results_dir = proc.run_recon_all(
            subject_input=subject_input,
            t1_pattern=t1_pattern,
            output_dir=output_dir,
            scheduler=scheduler,
            cores=cores,
            keep_intermediate=keep_intermediate,
            force_rules=force_rules,
            verbose=verbose,
        )
        print("✅ FreeSurfer recon-all job submitted successfully")
        return {
            "success": True,
            "job_id": job_id,
            "results_dir": results_dir,
            "subject_input": subject_input,
            "t1_pattern": t1_pattern,
            "output_dir": output_dir,
            "scheduler": scheduler,
        }
    except Exception as e:
        print(f"❌ FreeSurfer recon-all job failed:\n{e}", file=sys.stderr)
        return {
            "success": False,
            "error": str(e),
            "subject_input": subject_input,
            "t1_pattern": t1_pattern,
            "output_dir": output_dir,
        }
