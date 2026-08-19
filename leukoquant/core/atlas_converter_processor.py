#!/usr/bin/env python
"""
Atlas conversion processing for leukoquant.

Provides a processor class to run atlas conversion via a Snakemake workflow.
The workflow converts a brain parcellation file (e.g. GIF output) from one
label space to another (e.g. FreeSurfer aparc+aseg) using a CSV mapping file.

Output path: {output_dir}/{subject}/outputs/converted_atlas.mgz
"""

import logging
import subprocess
import sys
from pathlib import Path
import os
from typing import List, Optional

import yaml

from leukoquant.utils.container_utils import (
    ensure_container,
    MINICONDA_SIF_FILENAME,
    MINICONDA_SIF_HF_PATH,
)
from leukoquant.utils.external_utils import check_sge_plugin
from leukoquant.utils.snakemake_utils import add_forcerun_args, load_yaml_config, first_truthy

logger = logging.getLogger(__name__)


class AtlasConverterProcessor:
    """Processor to run atlas conversion through Snakemake."""

    def __init__(self, external_dir: Optional[str] = None):
        self.current_dir = Path(__file__).parent
        self.leukoquant_dir = self.current_dir.parent
        self.leukoquant_parent_dir = self.leukoquant_dir.parent

        if external_dir is None:
            external_dir = self.leukoquant_dir / "external"

        self.external_dir = Path(external_dir)
        self.workflow_dir = self.leukoquant_dir / "workflow" / "workflows"

    def run_atlas_conversion(self,
                             subject: str,
                             input_parcellation: str,
                             mapping_file: str,
                             output_dir: Optional[str] = None,
                             validate: bool = True,
                             verbose: bool = False,
                             scheduler: str = "local",
                             cores: int = 1,
                             force_rules: Optional[List[str]] = None) -> tuple[str, str]:
        """Run atlas conversion using Snakemake.

        Args:
            subject:            Subject ID. Used as the output subdirectory name.
            input_parcellation: Path to the input parcellation NIfTI file.
            mapping_file:       Path to the CSV label mapping file.
            output_dir:         Root output directory. The converted atlas is
                                written to {output_dir}/{subject}/outputs/converted_atlas.mgz.
                                Defaults to the parent directory of input_parcellation.
            validate:           Run post-conversion label validation.
            verbose:            Print Snakemake output to stdout.
            scheduler:          'local' or 'sge'.
            cores:              Number of Snakemake cores.

        Returns:
            Tuple of (job_id_str, output_subject_dir).
        """
        input_path = Path(input_parcellation).absolute()
        mapping_path = Path(mapping_file).absolute()

        if not input_path.exists():
            raise FileNotFoundError(f"Input parcellation not found: {input_parcellation}")
        if not mapping_path.exists():
            raise FileNotFoundError(f"Mapping file not found: {mapping_file}")

        wf_file = self.workflow_dir / "atlas_converter_workflow.smk"
        if not wf_file.exists():
            raise FileNotFoundError(f"Atlas converter workflow not found: {wf_file}")

        # Default output root: parent directory of the input file
        if output_dir is None:
            output_dir = str(input_path.parent)
        out_root = Path(output_dir).absolute()
        out_root.mkdir(parents=True, exist_ok=True)

        # Container-side paths
        input_sing = f"/input/{input_path.name}"

        # If the mapping CSV is inside the leukoquant package tree it is already
        # accessible inside the container via the /leukoquant bind; otherwise bind
        # it at a dedicated mount point.
        try:
            rel = mapping_path.relative_to(self.leukoquant_parent_dir.absolute())
            mapping_sing = f"/leukoquant/{rel}"
            extra_mapping_bind = ""
        except ValueError:
            mapping_sing = f"/input/parcellation_mapping{mapping_path.suffix}"
            extra_mapping_bind = f"{str(mapping_path)}:{mapping_sing}"

        singularity_bind = (
            f"--bind {str(self.leukoquant_parent_dir.absolute())}:/leukoquant,"
            f"{str(input_path)}:{input_sing},"
            f"{str(out_root)}:/output"
        )
        if extra_mapping_bind:
            singularity_bind += f",{extra_mapping_bind}"

        config_dict = {
            "subjects":                 [subject],
            "input_parcellation_map":   {subject: str(input_path)},
            "input_singularity_map":    {subject: input_sing},
            "mapping_file":             str(mapping_path),
            "mapping_file_singularity": mapping_sing,
            "output_dir":               str(out_root),
            "output_dir_singularity":   "/output",
            # Matches the TOOL_NAME default in atlas_converter_workflow.smk so
            # logs land at {subject}/atlas_conversion/logs/ and output at
            # {subject}/atlas_conversion/outputs/converted_atlas.mgz
            "output_subdir":            "atlas_conversion/outputs",
            "leukoquant_parent_dir":    str(self.leukoquant_parent_dir.absolute()),
            "validate":                 validate,
            "verbose":                  verbose,
        }

        config_file = str(out_root / f"atlas_config_{subject}.yaml")
        with open(config_file, "w") as f:
            yaml.dump(config_dict, f, default_flow_style=False)

        # Ensure container is present (downloaded from HuggingFace if missing)
        containers_folder = self.workflow_dir / "containers"
        sif_path = str(containers_folder / MINICONDA_SIF_FILENAME)
        ensure_container(sif_path, filename=MINICONDA_SIF_HF_PATH)

        snakemake_cmd = [
            "snakemake",
            "-s", str(wf_file),
            "--directory", str(out_root),
            "--cores", str(cores),
            "--jobs", "unlimited",
            "--software-deployment-method", "apptainer",
            "--apptainer-args", singularity_bind,
            f"--configfile={config_file}",
        ]

        if scheduler == "sge":
            logger.info("Using SGE scheduler")
            snakemake_cmd.extend([
                "--immediate-submit", "--notemp",
                "--latency-wait", "60",
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
        env["SNAKEMAKE_SOURCECACHE_PATH"] = str(out_root / ".snakemake" / "source-cache")
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
            subject_out_dir = str(out_root / subject / "atlas_conversion" / "outputs")
            return "0", subject_out_dir
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"Atlas conversion failed: {e}\nError Details:\n{e.stderr}"
            )


def apply_atlas_conversion(subject: Optional[str] = None,
                           input_parcellation: Optional[str] = None,
                           mapping_file: Optional[str] = None,
                           output_dir: Optional[str] = None,
                           validate: bool = True,
                           verbose: bool = False,
                           scheduler: str = "local",
                           cores: int = 1,
                           force_rules: Optional[List[str]] = None,
                           config_yaml: Optional[str] = None) -> dict:
    """Convenience function to run atlas conversion and return a summary dict.

    CLI / caller arguments take priority over values in ``config_yaml`` when
    both are provided.
    """
    yaml_cfg = load_yaml_config(config_yaml)
    subject            = first_truthy(subject,            yaml_cfg.get("subject"))
    input_parcellation = first_truthy(input_parcellation, yaml_cfg.get("input_parcellation"))
    mapping_file        = first_truthy(mapping_file,       yaml_cfg.get("mapping_file"))
    output_dir          = first_truthy(output_dir,         yaml_cfg.get("output_dir"))
    scheduler           = first_truthy(scheduler,          yaml_cfg.get("scheduler")) or "local"
    cores_raw           = first_truthy(cores,              yaml_cfg.get("cores"))
    cores               = int(cores_raw) if cores_raw is not None else 1

    if not subject:
        raise ValueError("subject is required (via argument or config_yaml)")
    if not input_parcellation:
        raise ValueError("input_parcellation is required (via argument or config_yaml)")
    if not mapping_file:
        raise ValueError("mapping_file is required (via argument or config_yaml)")

    if verbose:
        logging.basicConfig(level=logging.INFO)

    proc = AtlasConverterProcessor()
    try:
        job_id, results_dir = proc.run_atlas_conversion(
            subject=subject,
            input_parcellation=input_parcellation,
            mapping_file=mapping_file,
            output_dir=output_dir,
            validate=validate,
            verbose=verbose,
            scheduler=scheduler,
            cores=cores,
            force_rules=force_rules,
        )
        print("✅ Atlas conversion job submitted successfully")
        return {
            "success": True,
            "job_id": job_id,
            "results_dir": results_dir,
            "subject": subject,
            "input_parcellation": input_parcellation,
            "mapping_file": mapping_file,
            "output_dir": output_dir,
            "validate": validate,
            "scheduler": scheduler,
        }
    except Exception as e:
        print(f"❌ Atlas conversion job failed:\n{e}", file=sys.stderr)
        return {
            "success": False,
            "error": str(e),
            "subject": subject,
            "input_parcellation": input_parcellation,
            "mapping_file": mapping_file,
            "output_dir": output_dir,
        }
