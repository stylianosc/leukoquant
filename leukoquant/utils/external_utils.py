#!/usr/bin/env python
"""
Utilities for resolving bundled external tool paths (GIF, FreeSurfer license)
and validating runtime dependencies such as scheduler executor plugins.

Each helper checks for a bundled installation inside the leukoquant package
first, then falls back to the corresponding environment variable.
"""

import importlib.util
import os
import subprocess
import yaml
from pathlib import Path
from typing import List, Optional


def check_sge_plugin(scheduler: str) -> None:
    """Raise a clear error if SGE scheduler is requested but the plugin is not installed."""
    if scheduler.lower() != "sge":
        return
    if importlib.util.find_spec("snakemake_executor_plugin_sge") is None:
        raise RuntimeError(
            "SGE scheduler selected but 'snakemake-executor-plugin-sge' is not installed.\n"
            "Install it with:\n"
            "  pip install snakemake-executor-plugin-sge\n"
            "Or if it is hosted on GitHub:\n"
            "  pip install git+https://github.com/stylianosc/snakemake-executor-plugin-sge.git"
        )


def _resolve_gif_home(gif_dir: Path, ensure_t1_db: bool = True) -> Path:
    """Resolve the GIF software path.

    Checks for a bundled GIF installation at external/gif/ first. The directory
    must contain bin/ (binaries) since external/gif/ is bound to /GIF inside the
    Singularity container.
    Falls back to the GIF_HOME environment variable if the bundled copy is absent.

    ensure_t1_db: whether to also require/auto-download db_mideface/db.xml (the
        T1 GIF database) at the returned path. Defaults to True, matching every
        current caller except gif_processor.py's run_gif(), which passes False
        when no subject in its batch actually has a T1 image (per
        gif_workflow.smk's per-subject db selection, a pure-FLAIR batch never
        touches db_mideface at all, so there is no reason to require or
        download it there).
    """
    if gif_dir.exists():
        gif_home = gif_dir / "GIF"
        if ensure_t1_db:
            db_mideface_dir = gif_home / "db_mideface"
            db_xml = db_mideface_dir / "db.xml"
            if not db_xml.exists():
                from leukoquant.utils.container_utils import (
                    ensure_gif_db,
                    GIF_DB_T1_FILENAME,
                )
                ensure_gif_db(db_mideface_dir, GIF_DB_T1_FILENAME)
            if not db_xml.exists():
                raise FileNotFoundError(
                    f"Bundled GIF directory found at {gif_dir} but db_mideface/db.xml is missing.\n"
                    f"Expected: {db_xml}\n"
                    "Auto-download from Hugging Face was attempted and did not produce it "
                    "(check HF_TOKEN is set if the GIF database repo is still private).\n"
                    "Please add the GIF database files or set GIF_HOME to an external installation."
                )
        # Bind external/gif/GIF → /GIF so db is at /GIF/db_mideface/db.xml inside the container
        return gif_home

    # Fall back to GIF_HOME environment variable
    gif_home = os.environ.get("GIF_HOME", "").strip()
    if not gif_home:
        raise EnvironmentError(
            "GIF installation not found.\n"
            f"Either place the GIF files (including db_mideface/db.xml) in: {gif_dir / 'GIF'}\n"
            "Or set the GIF_HOME environment variable:\n"
            "  export GIF_HOME=/path/to/gif/installation"
        )
    gif_software_path = Path(gif_home).expanduser()
    if not gif_software_path.exists():
        raise FileNotFoundError(
            f"GIF_HOME path does not exist: {gif_software_path}\n"
            "Please set GIF_HOME to a valid directory."
        )
    return gif_software_path


def _find_gif_home_dir(gif_dir: Path) -> Optional[Path]:
    """Best-effort GIF home directory lookup, without validating its contents.

    Unlike _resolve_gif_home(), this never raises: it returns the bundled
    external/gif/GIF directory if present, else GIF_HOME if set and it
    exists, else None. Whether db_mideface/db.xml actually exists inside the
    returned directory is left for Snakemake to check as a rule input, so
    callers that don't always need GIF (e.g. BaMoS with externally-supplied
    results) aren't forced to have a complete GIF install just to bind it.
    """
    gif_home = gif_dir / "GIF"
    if gif_home.exists():
        return gif_home

    env_gif_home = os.environ.get("GIF_HOME", "").strip()
    if env_gif_home:
        gif_software_path = Path(env_gif_home).expanduser()
        if gif_software_path.exists():
            return gif_software_path

    return None


def _resolve_fs_license(leukoquant_dir: Path, verbose: bool = False) -> Path:
    """Resolve the FreeSurfer license file path.

    Checks leukoquant/licenses/license.txt first. Falls back to the FS_LICENSE
    environment variable (which may point to the file or its parent directory)
    if the bundled file is absent.
    """
    bundled = leukoquant_dir / "licenses" / "license.txt"
    if bundled.exists():
        if verbose:
            print(f"[leukoquant] Using FreeSurfer license: {bundled.absolute()}", flush=True)
        return bundled

    fs_license = os.environ.get("FS_LICENSE", "").strip()
    if not fs_license:
        raise EnvironmentError(
            "FreeSurfer license not found.\n"
            f"Either place license.txt in: {bundled.parent}\n"
            "Or set the FS_LICENSE environment variable:\n"
            "  export FS_LICENSE=/path/to/license.txt"
        )
    fs_license_path = Path(fs_license).expanduser()
    if not fs_license_path.exists():
        raise FileNotFoundError(
            f"FS_LICENSE path does not exist: {fs_license_path}\n"
            "Please set FS_LICENSE to a valid license.txt file or directory."
        )
    if verbose:
        print(f"[leukoquant] Using FreeSurfer license: {fs_license_path.absolute()}", flush=True)
    return fs_license_path


def run_snakemake_per_subject(
    cmd_base: List[str],
    subjects: List[str],
    output_dir: str,
    scheduler: str = "local",
    verbose: bool = False
) -> None:
    """Run Snakemake once per subject with subject-specific working directories.

    When using SGE scheduler, this ensures each subject has isolated logs and
    working directories at {subject}/.snakemake/ and {subject}/.sge_logs/.
    This makes debugging much easier when processing multiple subjects.

    For local scheduler, runs all subjects in a single Snakemake call.

    Parameters
    ----------
    cmd_base : List[str]
        Base Snakemake command (without subject-specific modifications)
    subjects : List[str]
        List of subject IDs to process
    output_dir : str
        Output directory where subject folders are located
    scheduler : str
        Scheduler type ('local' or 'sge')
    verbose : bool
        Enable verbose output
    """
    if scheduler.lower() != "sge":
        # For local scheduler, run all subjects in one call
        result = subprocess.run(cmd_base, check=True)
        return

    # For SGE: run one Snakemake job per subject with separate working directories
    output_path = Path(output_dir).absolute()
    failed_subjects = []

    for subject in subjects:
        subject_dir = output_path / subject
        subject_dir.mkdir(parents=True, exist_ok=True)

        # Create subject-specific working directory
        # SGE logs will be in .snakemake/sge_logs/ by default
        subject_workdir = subject_dir / ".snakemake"

        # Filter config to only this subject
        cmd = cmd_base.copy()

        # Add subject-specific working directory
        cmd.extend(["--directory", str(subject_dir)])

        # Modify config file if present to only include this subject
        # Find and update --configfile entries
        for i, arg in enumerate(cmd):
            if arg == "--configfile" and i + 1 < len(cmd):
                config_file = Path(cmd[i + 1])
                if config_file.exists():
                    # Create subject-specific config
                    import yaml
                    with open(config_file, "r") as f:
                        config = yaml.safe_load(f)
                    config["subjects"] = [subject]
                    subject_config = subject_dir / f"{subject}_config.yaml"
                    with open(subject_config, "w") as f:
                        yaml.dump(config, f, default_flow_style=False)
                    cmd[i + 1] = str(subject_config)
                break

        if verbose:
            print(f"\n{'='*60}")
            print(f"Processing subject: {subject}")
            print(f"Working directory: {subject_dir}")
            print(f"SGE logs: {subject_dir}/.snakemake/sge_logs/")
            print(f"{'='*60}\n")

        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            failed_subjects.append(subject)
            print(f"\n❌ Processing failed for subject: {subject}", flush=True)
            if verbose:
                print(f"Error details: {e}", flush=True)

    # Report final status
    if failed_subjects:
        print(f"\n⚠️  Failed subjects: {', '.join(failed_subjects)}", flush=True)
        raise RuntimeError(
            f"SGE processing failed for {len(failed_subjects)} subject(s): {', '.join(failed_subjects)}"
        )
    else:
        print(f"\n✅ All {len(subjects)} subject(s) processed successfully", flush=True)
