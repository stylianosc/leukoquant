#!/usr/bin/env python
"""
NODDI processing for leukoquant.
Supports single or multiple subjects via subject ID or text file.
"""

import logging
import subprocess
import sys
from pathlib import Path
import os
from typing import List, Optional, Tuple

import yaml

from leukoquant.utils.subject_utils import read_subjects, resolve_subject_pattern
from leukoquant.utils.external_utils import check_sge_plugin, _resolve_fs_license
from leukoquant.utils.bind_utils import dir_level_bind_files, consolidate_bind_entries
from leukoquant.utils.snakemake_utils import add_forcerun_args, load_yaml_config, first_truthy
from leukoquant.utils.container_utils import (
    ensure_container,
    MINICONDA_SIF_FILENAME,
    MINICONDA_SIF_HF_PATH,
    FREESURFER_SIF_FILENAME,
    FREESURFER_SIF_HF_PATH,
)

logger = logging.getLogger(__name__)



class NODDIProcessor:

    def __init__(self):
        self.current_dir = Path(__file__).parent
        self.leukoquant_dir = self.current_dir.parent
        self.leukoquant_parent_dir = self.leukoquant_dir.parent
        self.workflow_dir = self.leukoquant_dir / "workflow" / "workflows"

    def run_noddi(self,
                  subject_input: str,
                  dwi_pattern: str,
                  mask_pattern: str,
                  bvecs_pattern: str,
                  bvals_pattern: str,
                  output_dir: str,
                  scheduler: str = "local",
                  cores: int = 1,
                  keep_intermediate: bool = False,
                  skull_strip: bool = False,
                  force_rules: Optional[List[str]] = None,
                  verbose: bool = False) -> Tuple[str, str]:
        """Run NODDI fitting for one or more subjects."""
        wf_file = self.workflow_dir / "noddi_workflow.smk"

        subjects = read_subjects(subject_input)
        dwi_files = [resolve_subject_pattern(dwi_pattern, s) for s in subjects]

        mask_files: List[str] = []
        for s in subjects:
            try:
                mask_files.append(resolve_subject_pattern(mask_pattern, s) if mask_pattern else "")
            except FileNotFoundError:
                mask_files.append("")

        bvecs_files: List[str] = []
        for s in subjects:
            try:
                bvecs_files.append(resolve_subject_pattern(bvecs_pattern, s) if bvecs_pattern else "")
            except FileNotFoundError:
                bvecs_files.append("")

        bvals_files: List[str] = []
        for s in subjects:
            try:
                bvals_files.append(resolve_subject_pattern(bvals_pattern, s) if bvals_pattern else "")
            except FileNotFoundError:
                bvals_files.append("")

        out_path = Path(output_dir)
        out_path.mkdir(parents=True, exist_ok=True)

        dwi_sing,   dwi_binds   = dir_level_bind_files(dwi_files,   "input_dwi")
        mask_sing,  mask_binds  = dir_level_bind_files(mask_files,  "input_mask")
        bvecs_sing, bvecs_binds = dir_level_bind_files(bvecs_files, "input_bvecs")
        bvals_sing, bvals_binds = dir_level_bind_files(bvals_files, "input_bvals")

        # Setup scratch directory for intermediate files
        scratch_host = Path("/scratch0") if Path("/scratch0").exists() else out_path / ".scratch"
        scratch_host.mkdir(parents=True, exist_ok=True)
        scratch_dir_singularity = "/scratch0"

        # Consolidate all per-subject input binds to a single common-ancestor entry
        # per filesystem root, keeping --bind well within ARG_MAX for large cohorts.
        input_binds = dwi_binds + mask_binds + bvecs_binds + bvals_binds
        input_binds, remap = consolidate_bind_entries(input_binds)
        dwi_sing   = [remap(p) for p in dwi_sing]
        mask_sing  = [remap(p) for p in mask_sing]
        bvecs_sing = [remap(p) for p in bvecs_sing]
        bvals_sing = [remap(p) for p in bvals_sing]

        bind_parts = (
            input_binds
            + [f"{str(self.leukoquant_parent_dir.absolute())}:/leukoquant"]
            + [f"{str(out_path.absolute())}:/output"]
            + [f"{str(scratch_host.absolute())}:{scratch_dir_singularity}"]
        )

        if skull_strip:
            fs_license_path = _resolve_fs_license(self.leukoquant_dir)
            bind_parts.append(f"{str(fs_license_path.absolute())}:/license/license.txt")

        singularity_bind = "--bind " + ",".join(bind_parts)

        config_dict = {
            "subjects": subjects,
            "dwi_paths": dwi_files,
            "dwi_paths_singularity": dwi_sing,
            "brain_mask_paths": mask_files,
            "brain_mask_paths_singularity": mask_sing,
            "bvecs_paths": bvecs_files,
            "bvecs_paths_singularity": bvecs_sing,
            "bvals_paths": bvals_files,
            "bvals_paths_singularity": bvals_sing,
            "output_dir": str(out_path.absolute()),
            "output_dir_singularity": "/output",
            "scratch_dir_singularity": scratch_dir_singularity,
            "keep_intermediate": keep_intermediate,
            "skull_strip": skull_strip,
            "leukoquant_dir": str(self.leukoquant_dir.absolute()),
        }

        config_file = str(out_path.absolute()) + "/noddi_config.yaml"
        with open(config_file, "w") as f:
            yaml.dump(config_dict, f, default_flow_style=False)

        # Ensure container is present (downloaded from HuggingFace if missing).
        # Skull-strip mode requires the FreeSurfer container (mri_synthstrip);
        # standard fitting uses the miniconda container.
        containers_folder = self.leukoquant_dir / "workflow" / "containers"
        if skull_strip:
            sif_path = str(containers_folder / FREESURFER_SIF_FILENAME)
            ensure_container(sif_path, filename=FREESURFER_SIF_HF_PATH)
        else:
            sif_path = str(containers_folder / MINICONDA_SIF_FILENAME)
            ensure_container(sif_path, filename=MINICONDA_SIF_HF_PATH)

        snakemake_cmd = [
            "snakemake", "-s", str(wf_file),
            "--directory", str(out_path),
            "--cores", str(cores),
            "--jobs", "unlimited",
            "--software-deployment-method", "apptainer",
            "--apptainer-args", singularity_bind,
            f"--configfile={config_file}",
        ]

        if scheduler == "sge":
            snakemake_cmd.extend([
                "--immediate-submit", "--notemp",
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
            return "0", str(out_path.absolute())
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"NODDI failed: {e}\nError Details:\n{e.stderr}"
            )


def apply_noddi(subject_input: Optional[str] = None,
                dwi_pattern: Optional[str] = None,
                mask_pattern: Optional[str] = None,
                bvecs_pattern: Optional[str] = None,
                bvals_pattern: Optional[str] = None,
                output_dir: Optional[str] = None,
                scheduler: str = "local",
                cores: int = 1,
                keep_intermediate: bool = False,
                skull_strip: bool = False,
                force_rules: Optional[List[str]] = None,
                verbose: bool = False,
                config_yaml: Optional[str] = None) -> dict:
    """Run NODDI fitting and return a summary dict.

    CLI / caller arguments take priority over values in ``config_yaml`` when
    both are provided.
    """
    yaml_cfg = load_yaml_config(config_yaml)
    subject_input = first_truthy(subject_input, yaml_cfg.get("subject_input"), yaml_cfg.get("subject"))
    dwi_pattern   = first_truthy(dwi_pattern,   yaml_cfg.get("dwi_pattern"), yaml_cfg.get("dwi"))
    mask_pattern  = first_truthy(mask_pattern,  yaml_cfg.get("mask_pattern"))
    bvecs_pattern = first_truthy(bvecs_pattern, yaml_cfg.get("bvecs_pattern"), yaml_cfg.get("bvecs")) or ""
    bvals_pattern = first_truthy(bvals_pattern, yaml_cfg.get("bvals_pattern"), yaml_cfg.get("bvals")) or ""
    output_dir    = first_truthy(output_dir,    yaml_cfg.get("output_dir"))
    scheduler     = first_truthy(scheduler,     yaml_cfg.get("scheduler")) or "local"
    cores_raw     = first_truthy(cores,         yaml_cfg.get("cores"))
    cores         = int(cores_raw) if cores_raw is not None else 1
    keep_intermediate = bool(yaml_cfg.get("keep_intermediate", False)) or bool(keep_intermediate)
    skull_strip   = bool(yaml_cfg.get("skull_strip", False)) or bool(skull_strip)

    if not subject_input:
        raise ValueError("subject_input is required (via argument or config_yaml)")
    if not dwi_pattern:
        raise ValueError("dwi_pattern is required (via argument or config_yaml)")
    if not output_dir:
        raise ValueError("output_dir is required (via argument or config_yaml)")

    if verbose:
        logging.basicConfig(level=logging.INFO)
    proc = NODDIProcessor()
    try:
        job_id, results_dir = proc.run_noddi(
            subject_input=subject_input,
            dwi_pattern=dwi_pattern,
            mask_pattern=mask_pattern,
            bvecs_pattern=bvecs_pattern,
            bvals_pattern=bvals_pattern,
            output_dir=output_dir,
            scheduler=scheduler,
            cores=cores,
            keep_intermediate=keep_intermediate,
            skull_strip=skull_strip,
            force_rules=force_rules,
            verbose=verbose,
        )
        print("✅ NODDI job submitted successfully")
        return {"success": True, "job_id": job_id, "results_dir": results_dir,
                "subject_input": subject_input, "output_dir": str(Path(output_dir).absolute())}
    except Exception as e:
        print(f"❌ NODDI job failed:\n{e}", file=sys.stderr)
        return {"success": False, "error": str(e),
                "subject_input": subject_input, "output_dir": str(Path(output_dir).absolute())}
