#!/usr/bin/env python
"""
GIF (Geodesic Information Flow) processing for leukoquant.
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
from leukoquant.utils.bind_utils import dir_level_bind_files, consolidate_bind_entries
from leukoquant.utils.snakemake_utils import add_forcerun_args, load_yaml_config, first_truthy
from leukoquant.utils.container_utils import (
    ensure_container,
    MINICONDA_SIF_FILENAME,
    MINICONDA_SIF_HF_PATH,
)
from leukoquant.utils.external_utils import _resolve_gif_home

logger = logging.getLogger(__name__)
# Set root logger to WARNING by default (only show errors/warnings unless --verbose)
logging.getLogger().setLevel(logging.WARNING)



class GIFProcessor:
    """Processor for GIF neuroimaging tool."""

    def __init__(self, external_dir: Optional[str] = None):
        self.current_dir = Path(__file__).parent
        self.leukoquant_dir = self.current_dir.parent
        self.leukoquant_parent_dir = self.leukoquant_dir.parent

        if external_dir is None:
            external_dir = self.leukoquant_dir / "external"
        self.external_dir = Path(external_dir)
        self.gif_dir = self.external_dir / "gif"

        self._validate_directories()
        self._set_executable_permissions()

    def _validate_directories(self):
        if not self.gif_dir.exists():
            raise FileNotFoundError(f"Required directory not found: {self.gif_dir}")
        for script in ["GIF_111125.sh"]:
            sp = self.gif_dir / script
            if not sp.exists():
                raise FileNotFoundError(f"Required script not found: {sp}")

    def _set_executable_permissions(self):
        for script_path in self.gif_dir.iterdir():
            if script_path.exists():
                script_path.chmod(0o755)

    def run_gif(self,
                subject_input: str,
                output_dir: str,
                t1_pattern: Optional[str] = None,
                flair_pattern: Optional[str] = None,
                mask_pattern: Optional[str] = None,
                flair_db: Optional[str] = None,
                scheduler: str = "local",
                cores: int = 1,
                keep_intermediate: bool = False,
                force_rules: Optional[List[str]] = None,
                verbose: bool = False) -> Tuple[str, str]:
        """Run GIF segmentation for one or more subjects.

        At least one of t1_pattern or flair_pattern must be provided.
        When both are given, T1 is used as the primary input and FLAIR is ignored
        (the GIF script does not support multi-channel input).
        When only one is given, that image is used as the primary input and the
        output ID is derived from its filename.

        flair_db: optional path to a directory containing a FLAIR-specific GIF
            database (db.xml). Used when the primary input is FLAIR (i.e. no T1
            is provided). Defaults to <GIF_HOME>/db_FLAIR_mideface/, which must contain
            db.xml. Ignored when --t1 is provided (T1 runs always use
            <GIF_HOME>/db_mideface/db.xml).
        """
        if not t1_pattern and not flair_pattern:
            raise ValueError("At least one of --t1 or --flair must be provided for GIF.")

        subjects = read_subjects(subject_input)

        def _resolve_optional_list(pattern: Optional[str]) -> List[str]:
            if not pattern:
                return [""] * len(subjects)
            result = []
            for s in subjects:
                try:
                    result.append(resolve_subject_pattern(pattern, s))
                except FileNotFoundError:
                    result.append("")
            return result

        t1_files    = _resolve_optional_list(t1_pattern)
        flair_files = _resolve_optional_list(flair_pattern)

        # The primary image drives the output filename (t1 preferred over flair).
        primary_files = [
            t1 if t1 else flair
            for t1, flair in zip(t1_files, flair_files)
        ]
        image_ids = [Path(f).name.split(".", 1)[0] for f in primary_files]

        mask_files: List[str] = []
        if mask_pattern:
            for s in subjects:
                try:
                    mask_files.append(resolve_subject_pattern(mask_pattern, s))
                except FileNotFoundError:
                    mask_files.append("")
        else:
            mask_files = [""] * len(subjects)

        out_path = Path(output_dir)
        out_path.mkdir(parents=True, exist_ok=True)

        for s in subjects:
            (out_path / s).mkdir(parents=True, exist_ok=True)

        # Per gif_workflow.smk's per-subject db selection (a subject uses the
        # T1 db if it has a T1 file, and only falls back to the FLAIR db when
        # it doesn't - see db_arg in gif_workflow.smk), a single batch can mix
        # both, need only one, or need neither exclusively. Gate which
        # database(s) we require/auto-download on what this batch actually
        # uses, so a T1-only run never pulls the ~830MB FLAIR tarball (and
        # vice versa).
        needs_t1_db    = any(t1_files)
        needs_flair_db = any(not t1 for t1 in t1_files)

        gif_software_path = _resolve_gif_home(self.gif_dir, ensure_t1_db=needs_t1_db)

        t1_sing,    t1_binds    = dir_level_bind_files(t1_files,    "input_t1")
        flair_sing, flair_binds = dir_level_bind_files(flair_files, "input_flair")
        mask_sing,  mask_binds  = dir_level_bind_files(mask_files,  "input_mask")

        # Resolve the FLAIR db container path.
        # flair_db is a directory expected to contain db.xml.
        # Default: db_FLAIR_mideface/ inside the GIF installation tree, which is already
        # covered by the /GIF bind - no extra mount needed.
        # If a directory outside the GIF tree is given, bind it separately.
        flair_db_binds: List[str] = []
        flair_db_container_dir: str = ""
        flair_db_dir = Path(flair_db).absolute() if flair_db else gif_software_path / "db_FLAIR_mideface"
        if needs_flair_db:
            if not (flair_db_dir / "db.xml").exists() and not flair_db:
                # Only auto-download the default location - an explicitly-provided
                # --gif-flair-db path is the caller's own directory and should
                # never be silently overwritten/populated by us.
                from leukoquant.utils.container_utils import (
                    ensure_gif_db,
                    GIF_DB_FLAIR_FILENAME,
                )
                ensure_gif_db(flair_db_dir, GIF_DB_FLAIR_FILENAME)
            if not flair_db_dir.exists():
                raise FileNotFoundError(
                    f"FLAIR GIF database directory not found: {flair_db_dir}\n"
                    "Auto-download from Hugging Face was attempted and did not produce it "
                    "(check HF_TOKEN is set if the GIF database repo is still private).\n"
                    "Provide a directory containing db.xml via --gif-flair-db, "
                    "or place the FLAIR database at <GIF_HOME>/db_FLAIR_mideface/."
                )
        try:
            # Inside the GIF installation tree → reachable via the /GIF bind
            rel = flair_db_dir.relative_to(gif_software_path.absolute())
            flair_db_sing = f"/GIF/{rel}/db.xml"
        except ValueError:
            # Outside the GIF tree → needs its own bind mount; included in consolidation
            flair_db_container_dir = "/input_flair_db"
            flair_db_sing = f"{flair_db_container_dir}/db.xml"
            flair_db_binds = [f"{str(flair_db_dir)}:{flair_db_container_dir}"]

        # Setup scratch directory for intermediate files
        scratch_host = Path("/scratch0") if Path("/scratch0").exists() else out_path / ".scratch"
        scratch_host.mkdir(parents=True, exist_ok=True)
        scratch_dir_singularity = "/scratch0"

        # Consolidate all per-subject input binds (plus optional flair_db if external).
        input_binds = t1_binds + flair_binds + mask_binds + flair_db_binds
        input_binds, remap = consolidate_bind_entries(input_binds)
        t1_sing    = [remap(p) for p in t1_sing]
        flair_sing = [remap(p) for p in flair_sing]
        mask_sing  = [remap(p) for p in mask_sing]
        if flair_db_container_dir:
            flair_db_sing = remap(flair_db_container_dir) + "/db.xml"

        bind_parts = (
            input_binds
            + [f"{str(self.leukoquant_parent_dir.absolute())}:/leukoquant"]
            + [f"{str(out_path.absolute())}:/output"]
            + [f"{str(gif_software_path.absolute())}:/GIF"]
            + [f"{str(scratch_host.absolute())}:{scratch_dir_singularity}"]
        )
        singularity_bind = "--bind " + ",".join(bind_parts)

        workflow_dir = self.leukoquant_dir / "workflow" / "workflows"
        snakefile = workflow_dir / "gif_workflow.smk"
        if not snakefile.exists():
            raise FileNotFoundError(f"Snakefile not found: {snakefile}")

        config_dict = {
            "subjects": subjects,
            "t1_files": t1_files,
            "flair_files": flair_files,
            "image_ids": image_ids,
            "t1_files_singularity": t1_sing,
            "flair_files_singularity": flair_sing,
            "mask_files": mask_files,
            "mask_files_singularity": mask_sing,
            "output_dir": str(out_path.absolute()),
            "output_dir_singularity": "/output",
            "scratch_dir_singularity": scratch_dir_singularity,
            "keep_intermediate": keep_intermediate,
            "leukoquant_parent_dir": str(self.leukoquant_parent_dir.absolute()),
            # Database paths inside the container (t1_db always uses the default)
            "flair_db_singularity": flair_db_sing,
        }

        config_file = str(out_path.absolute()) + "/gif_config.yaml"
        with open(config_file, "w") as f:
            yaml.dump(config_dict, f, default_flow_style=False)

        # Ensure container is present (downloaded from HuggingFace if missing)
        containers_folder = self.leukoquant_dir / "workflow" / "containers"
        sif_path = str(containers_folder / MINICONDA_SIF_FILENAME)
        ensure_container(sif_path, filename=MINICONDA_SIF_HF_PATH)

        snakemake_cmd = [
            "snakemake", "-s", str(snakefile),
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
                "--executor", "sge"
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
                    env=env,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
            )
            if verbose and result.stdout:
                    print(result.stdout)
            return "0", str(out_path)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"GIF processing failed: {e}\nError Details:\n{e.stderr}")


def apply_gif(subject_input: Optional[str] = None,
              output_dir: Optional[str] = None,
              t1_pattern: Optional[str] = None,
              flair_pattern: Optional[str] = None,
              mask_pattern: Optional[str] = None,
              flair_db: Optional[str] = None,
              scheduler: str = "local",
              cores: int = 1,
              keep_intermediate: bool = False,
              force_rules: Optional[List[str]] = None,
              verbose: bool = False,
              config_yaml: Optional[str] = None) -> dict:
    """Run GIF segmentation and return a summary dict.

    CLI / caller arguments take priority over values in ``config_yaml`` when
    both are provided.
    """
    yaml_cfg = load_yaml_config(config_yaml)
    subject_input = first_truthy(subject_input, yaml_cfg.get("subject_input"), yaml_cfg.get("subject"))
    output_dir    = first_truthy(output_dir,    yaml_cfg.get("output_dir"))
    t1_pattern    = first_truthy(t1_pattern,    yaml_cfg.get("t1_pattern"), yaml_cfg.get("t1"))
    flair_pattern = first_truthy(flair_pattern, yaml_cfg.get("flair_pattern"), yaml_cfg.get("flair"))
    mask_pattern  = first_truthy(mask_pattern,  yaml_cfg.get("mask_pattern"))
    flair_db      = first_truthy(flair_db,      yaml_cfg.get("flair_db"))
    scheduler     = first_truthy(scheduler,     yaml_cfg.get("scheduler")) or "local"
    cores_raw     = first_truthy(cores,         yaml_cfg.get("cores"))
    cores         = int(cores_raw) if cores_raw is not None else 1
    keep_intermediate = bool(yaml_cfg.get("keep_intermediate", False)) or bool(keep_intermediate)

    if not subject_input:
        raise ValueError("subject_input is required (via argument or config_yaml)")
    if not output_dir:
        raise ValueError("output_dir is required (via argument or config_yaml)")
    if not t1_pattern and not flair_pattern:
        raise ValueError("at least one of t1_pattern or flair_pattern is required (via argument or config_yaml)")

    if verbose:
        logging.basicConfig(level=logging.INFO)
    try:
        processor = GIFProcessor()
        job_id, results_dir = processor.run_gif(
            subject_input=subject_input,
            t1_pattern=t1_pattern,
            flair_pattern=flair_pattern,
            output_dir=output_dir,
            mask_pattern=mask_pattern,
            flair_db=flair_db,
            scheduler=scheduler,
            cores=cores,
            keep_intermediate=keep_intermediate,
            force_rules=force_rules,
            verbose=verbose,
        )
        print("✅ GIF job submitted successfully")
        return {"success": True, "job_id": job_id, "results_dir": results_dir,
                "subject_input": subject_input, "output_dir": output_dir}
    except Exception as e:
        print(f"❌ GIF job failed:\n{e}", file=sys.stderr)
        return {"success": False, "error": str(e),
                "subject_input": subject_input, "output_dir": output_dir}
