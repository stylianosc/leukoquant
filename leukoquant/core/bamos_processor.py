#!/usr/bin/env python
"""
BaMoS (Bayesian Model Selection) processing for leukoquant.
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
from leukoquant.utils.container_utils import (
    ensure_container,
    MINICONDA_SIF_FILENAME,
    MINICONDA_SIF_HF_PATH,
)
from leukoquant.utils.external_utils import _find_gif_home_dir
from leukoquant.utils.bind_utils import dir_level_bind_files, consolidate_bind_entries
from leukoquant.utils.snakemake_utils import add_forcerun_args, load_yaml_config, first_truthy

logger = logging.getLogger(__name__)



class BaMoSProcessor:

    def __init__(self):
        self.current_dir = Path(__file__).parent
        self.leukoquant_dir = self.current_dir.parent
        self.leukoquant_parent_dir = self.leukoquant_dir.parent

    def run_bamos(self,
                  subject_input: str,
                  flair_pattern: str,
                  t1_pattern: str,
                  output_dir: str,
                  gif_results_pattern: Optional[str] = None,
                  space: int = 1,
                  opt: str = "TA",
                  jump_start: int = 0,
                  scheduler: str = "local",
                  cores: int = 1,
                  keep_intermediate: bool = False,
                  force_rules: Optional[List[str]] = None,
                  verbose: bool = False) -> Tuple[str, str]:
        """Run BaMoS lesion detection for one or more subjects.

        Args:
            subject_input: bare subject ID or path to a text file (one ID per line)
            flair_pattern: explicit path or {subject} glob pattern for FLAIR files
            t1_pattern: explicit path or {subject} glob pattern for T1 files
            gif_results_pattern: explicit path or {subject} glob pattern for GIF result directories
            output_dir: base output directory (per-subject subdirs created automatically)
            space: 1=T1, 2=FLAIR
            opt: optimisation parameter (default "TA")
            jump_start: 0=full run, 1=skip initial step
            scheduler: "local" or "sge"
            cores: Snakemake cores

        Returns:
            (job_id, results_dir)
        """
        subjects = read_subjects(subject_input)
        flair_files = [resolve_subject_pattern(flair_pattern, s) for s in subjects]
        t1_files = [resolve_subject_pattern(t1_pattern, s) for s in subjects]
        # A subject whose gif_results_pattern doesn't match anything falls
        # back to "" (no external results) rather than raising -- a single
        # subject's missing/mis-named directory shouldn't abort the whole
        # batch. Same underlying resolve_subject_pattern() glob call either
        # way, just with all=True so a zero-match result is returned instead
        # of raised.
        gif_results_paths = [
            (resolve_subject_pattern(gif_results_pattern, s, all=True) or [""])[0]
            if gif_results_pattern else ""
            for s in subjects
        ]

        # Best-effort only: BaMoS may run entirely against externally-supplied
        # --gif-results-dir results, in which case a full GIF install isn't
        # needed at all. Resolved once here (not per subject) both for the
        # /GIF bind below and to decide, just below, which of the subjects
        # lacking external results can actually have internal GIF run for
        # them.
        gif_software_path = _find_gif_home_dir(self.leukoquant_dir / "external" / "gif")

        # Subjects with neither external GIF results nor a usable GIF
        # database can't be processed at all -- exclude just those subjects
        # (with a clear warning) rather than failing the whole batch or
        # silently attempting internal GIF with no database to run it
        # against. This is a single database-existence check, not a
        # per-subject one: there is only one GIF install for the whole run,
        # so its usability is the same answer for every subject that needs
        # it. BaMoS always requires T1 (see the required flair_pattern/
        # t1_pattern args below), so any subject falling back to internal
        # GIF always uses the T1 database, never the FLAIR-only one.
        needs_internal_gif = [s for s, p in zip(subjects, gif_results_paths) if not p]
        if needs_internal_gif and gif_software_path:
            # Auto-download, gated strictly on actually needing internal GIF:
            # if every subject supplies external --gif-results-dir output,
            # BaMoS never touches the GIF database at all and shouldn't force
            # a ~900MB download just because a bundled/external GIF install
            # happens to be configured. T1-only, per the comment above (BaMoS
            # always requires T1) -- the FLAIR db is never needed here.
            db_mideface_dir = gif_software_path / "db_mideface"
            if not (db_mideface_dir / "db.xml").exists():
                from leukoquant.utils.container_utils import (
                    ensure_gif_db,
                    GIF_DB_T1_FILENAME,
                )
                ensure_gif_db(db_mideface_dir, GIF_DB_T1_FILENAME)
        if needs_internal_gif:
            t1_db_usable = bool(
                gif_software_path and (gif_software_path / "db_mideface" / "db.xml").exists()
            )
            if not t1_db_usable:
                for s in needs_internal_gif:
                    print(
                        f"WARNING: skipping subject '{s}': no matching --gif-results-dir "
                        "and no usable GIF database found for internal GIF.",
                        file=sys.stderr,
                    )
                keep = [s not in needs_internal_gif for s in subjects]
                subjects          = [s for s, k in zip(subjects, keep) if k]
                flair_files       = [v for v, k in zip(flair_files, keep) if k]
                t1_files          = [v for v, k in zip(t1_files, keep) if k]
                gif_results_paths = [v for v, k in zip(gif_results_paths, keep) if k]
                if not subjects:
                    raise ValueError(
                        "No subjects remain after excluding those with neither "
                        "external GIF results nor a usable GIF database."
                    )

        out_path = Path(output_dir)
        out_path.mkdir(parents=True, exist_ok=True)

        for s in subjects:
            (out_path / s / "bamos").mkdir(parents=True, exist_ok=True)

        flair_sing, flair_binds = dir_level_bind_files(flair_files, "input_flair")
        t1_sing, t1_binds     = dir_level_bind_files(t1_files,    "input_t1")

        # GIF results directories: already directory-level, one per subject.
        gif_container_paths = [f"/gif_input_{i}" for i in range(len(gif_results_paths))]
        gif_binds = [
            f"{p}:{gif_container_paths[i]}"
            for i, p in enumerate(gif_results_paths) if p
        ]

        scratch_dir = Path("/scratch0") if Path("/scratch0").exists() else out_path / ".scratch"
        scratch_dir.mkdir(parents=True, exist_ok=True)
        scratch_dir_sing = "/scratch0"

        # Consolidate all per-subject input binds to a single common-ancestor entry
        # per filesystem root, keeping --bind well within ARG_MAX for large cohorts.
        input_binds = flair_binds + t1_binds + gif_binds
        input_binds, remap = consolidate_bind_entries(input_binds)
        flair_sing     = [remap(p) for p in flair_sing]
        t1_sing        = [remap(p) for p in t1_sing]
        gif_sing_paths = [remap(p) for p in gif_container_paths]

        bind_parts = (
            input_binds
            + [f"{str(self.leukoquant_parent_dir.absolute())}:/leukoquant"]
            + [f"{str(out_path.absolute())}:/bamos_output"]
            + ([f"{str(gif_software_path.absolute())}:/GIF"] if gif_software_path else [])
            + [f"{str(scratch_dir.absolute())}:{scratch_dir_sing}"]
        )
        singularity_bind = "--bind " + ",".join(bind_parts)

        workflow_dir = self.leukoquant_dir / "workflow" / "workflows"
        snakefile = workflow_dir / "bamos_workflow.smk"
        if not snakefile.exists():
            raise FileNotFoundError(f"Snakefile not found: {snakefile}")

        config_dict = {
            "subjects": subjects,
            "flair_files": flair_files,
            "flair_files_singularity": flair_sing,
            "t1_files": t1_files,
            "t1_files_singularity": t1_sing,
            "gif_results_paths": gif_results_paths,
            "gif_results_paths_singularity": gif_sing_paths,
            "output_dir": str(out_path.absolute()),
            "output_dir_singularity": "/bamos_output",
            "scratch_dir_singularity": scratch_dir_sing,
            "keep_intermediate": keep_intermediate,
            "space": space,
            "opt": opt,
            "jump_start": jump_start,
            "leukoquant_parent_dir": str(self.leukoquant_parent_dir.absolute()),
            **({"gif_home_host": str(gif_software_path.absolute())} if gif_software_path else {}),
        }

        config_file = str(out_path.absolute()) + "/bamos_config.yaml"
        with open(config_file, "w") as f:
            yaml.dump(config_dict, f, default_flow_style=False)

        # Ensure container is present (downloaded from HuggingFace if missing)
        containers_folder = self.leukoquant_dir / "workflow" / "containers"
        sif_path = str(containers_folder / MINICONDA_SIF_FILENAME)
        ensure_container(sif_path, filename=MINICONDA_SIF_HF_PATH)

        snakemake_cmd = [
            "snakemake", "-s", str(snakefile),
            "--directory", str(out_path),
            "--cores", "all",
            "--jobs", "unlimited",
            "--max-jobs-per-timespan", "75000/1s",
            "--drop-metadata",
            "--immediate-submit", "--notemp",
            # Cold NFS mounts on some compute nodes exceed the default 5s.
            "--latency-wait", "60",
            "--software-deployment-method", "apptainer",
            "--apptainer-args", singularity_bind,
            f"--configfile={config_file}",
        ]

        if scheduler == "sge":
            snakemake_cmd.extend([
                "--executor", "sge"
            ])

        add_forcerun_args(snakemake_cmd, force_rules)
        # Explicitly target BaMoS's own `rule all` so that imported module rules
        # (e.g. gif_all from the GIF module) don't become the default target.
        snakemake_cmd.append("all")

        env = os.environ.copy()
        # Isolate the Snakemake source cache per output directory to avoid NFS
        # race conditions when multiple instances run in parallel on a shared filesystem.
        env["SNAKEMAKE_SOURCECACHE_PATH"] = str(out_path / ".snakemake" / "source-cache")
        try:
            result = subprocess.run(
                    snakemake_cmd,
                    check=True,
                    cwd=str(workflow_dir),
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
                f"BaMoS processing failed: {e}\nError Details:\n{e.stderr}"
            )


def apply_bamos(subject_input: Optional[str] = None,
                flair_pattern: Optional[str] = None,
                t1_pattern: Optional[str] = None,
                output_dir: Optional[str] = None,
                gif_results_pattern: Optional[str] = None,
                scheduler: str = "local",
                cores: int = 1,
                keep_intermediate: bool = False,
                force_rules: Optional[List[str]] = None,
                verbose: bool = False,
                config_yaml: Optional[str] = None) -> dict:
    """Run BaMoS and return a summary dict.

    CLI / caller arguments take priority over values in ``config_yaml`` when
    both are provided.
    """
    yaml_cfg = load_yaml_config(config_yaml)
    subject_input = first_truthy(subject_input, yaml_cfg.get("subject_input"), yaml_cfg.get("subject"))
    flair_pattern = first_truthy(flair_pattern, yaml_cfg.get("flair_pattern"), yaml_cfg.get("flair"))
    t1_pattern    = first_truthy(t1_pattern,    yaml_cfg.get("t1_pattern"), yaml_cfg.get("t1"))
    output_dir    = first_truthy(output_dir,    yaml_cfg.get("output_dir"))
    gif_results_pattern = first_truthy(gif_results_pattern, yaml_cfg.get("gif_results_pattern"))
    scheduler     = first_truthy(scheduler,     yaml_cfg.get("scheduler")) or "local"
    cores_raw     = first_truthy(cores,         yaml_cfg.get("cores"))
    cores         = int(cores_raw) if cores_raw is not None else 1
    keep_intermediate = bool(yaml_cfg.get("keep_intermediate", False)) or bool(keep_intermediate)

    if not subject_input:
        raise ValueError("subject_input is required (via argument or config_yaml)")
    if not flair_pattern:
        raise ValueError("flair_pattern is required (via argument or config_yaml)")
    if not t1_pattern:
        raise ValueError("t1_pattern is required (via argument or config_yaml)")
    if not output_dir:
        raise ValueError("output_dir is required (via argument or config_yaml)")

    if verbose:
        logging.basicConfig(level=logging.INFO)
    try:
        processor = BaMoSProcessor()
        job_id, results_dir = processor.run_bamos(
            subject_input=subject_input,
            flair_pattern=flair_pattern,
            t1_pattern=t1_pattern,
            output_dir=output_dir,
            gif_results_pattern=gif_results_pattern,
            scheduler=scheduler,
            cores=cores,
            keep_intermediate=keep_intermediate,
            force_rules=force_rules,
            verbose=verbose,
        )
        print("✅ BaMoS job submitted successfully")
        return {"success": True, "job_id": job_id, "results_dir": results_dir,
                "subject_input": subject_input, "output_dir": output_dir}
    except Exception as e:
        print(f"❌ BaMoS job failed:\n{e}", file=sys.stderr)
        return {"success": False, "error": str(e),
                "subject_input": subject_input, "output_dir": output_dir}
