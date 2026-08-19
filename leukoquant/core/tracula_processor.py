#!/usr/bin/env python
"""
TRACULA processing for leukoquant.
Supports single or multiple subjects via subject ID or text file.

If recon-all outputs are absent for any subject the recon_all_workflow module
is activated automatically - mirroring how bamos_workflow handles missing GIF
results.  The caller may supply --freesurfer-recon-dir to point at an existing
SUBJECTS_DIR; if omitted (or if aparc+aseg.mgz is missing there), recon-all
will be run into the standard output path and --t1 is required.
"""

import logging
import subprocess
import sys
from pathlib import Path
import os
from typing import Dict, List, Optional, Set, Tuple

import yaml

from leukoquant.utils.subject_utils import read_subjects, resolve_subject_pattern
from leukoquant.utils.container_utils import (
    ensure_container,
    FREESURFER_SIF_FILENAME,
    FREESURFER_SIF_HF_PATH,
)
from leukoquant.utils.external_utils import _resolve_fs_license
from leukoquant.utils.bind_utils import (
    dir_level_bind_files,
    dir_level_bind_file_lists,
    consolidate_bind_entries,
)
from leukoquant.utils.snakemake_utils import add_forcerun_args, load_yaml_config, first_truthy

logger = logging.getLogger(__name__)

_RECON_SENTINEL = "aparc+aseg.mgz"


def _is_nifti(path: Path) -> bool:
    """Return True only for files that already carry bvec/bval sidecars.

    .zip archives and DICOM files are converted by dcm2niix inside the
    container; their sidecars don't exist on the host until after conversion.
    """
    name = path.name.lower()
    return name.endswith('.nii.gz') or name.endswith('.nii')



def _discover_dwi_with_sidecars(
    dwi_pattern: str,
    subjects: List[str],
    bvecs_pattern: Optional[str],
    bvals_pattern: Optional[str],
) -> Tuple[List[List[str]], List[List[str]], List[List[str]]]:
    """Return (dwi_lists, bvec_lists, bval_lists) - one inner list per subject.

    For each subject all DWI files matching dwi_pattern are collected.  Each
    DWI must have a same-stem .bvec and .bval sidecar in the same directory;
    DWIs lacking either sidecar are skipped with a warning.

    Special case: when only one DWI matches AND explicit bvecs/bvals patterns
    are provided, those patterns are honoured instead of auto-discovery (backward
    compatibility with callers that supply separate --bvecs / --bvals arguments).

    Raises:
        FileNotFoundError: if no DWI file at all matches the pattern for a subject.
        ValueError: if every matched DWI for a subject is missing bvec or bval.
    """
    all_dwi: List[List[str]] = []
    all_bvec: List[List[str]] = []
    all_bval: List[List[str]] = []

    for s in subjects:
        # Collect every DWI file that matches - all=True returns [] on no match.
        candidates: List[str] = resolve_subject_pattern(dwi_pattern, s, all=True)  # type: ignore[assignment]
        if not candidates:
            raise FileNotFoundError(
                f"No DWI files found for subject '{s}' with pattern: {dwi_pattern}"
            )

        # Single DWI + explicit bvecs/bvals patterns → honour those patterns
        # exactly as the old code did (first match or empty string on missing).
        if len(candidates) == 1 and (bvecs_pattern or bvals_pattern):
            dwi = candidates[0]
            try:
                bvec: str = resolve_subject_pattern(bvecs_pattern, s) if bvecs_pattern else ""  # type: ignore[assignment]
            except FileNotFoundError:
                bvec = ""
            try:
                bval: str = resolve_subject_pattern(bvals_pattern, s) if bvals_pattern else ""  # type: ignore[assignment]
            except FileNotFoundError:
                bval = ""
            all_dwi.append([dwi])
            all_bvec.append([bvec])
            all_bval.append([bval])
            continue

        # Multi-DWI path (or single DWI without explicit sidecar patterns):
        # auto-discover bvec/bval alongside each DWI by matching the file stem.
        valid_dwi, valid_bvec, valid_bval = [], [], []
        for dwi_path in candidates:
            dwi_p = Path(dwi_path)
            # Strip compound extensions (.nii.gz → stem without .nii or .gz).
            stem = dwi_p.name.split(".")[0]

            if _is_nifti(dwi_p):
                bvec_p = dwi_p.with_name(stem + ".bvec")
                bval_p = dwi_p.with_name(stem + ".bval")
                if not bvec_p.exists() or not bval_p.exists():
                    logger.warning(
                        "Skipping DWI '%s' for subject '%s': "
                        "missing bvec (%s exists=%s) or bval (%s exists=%s).",
                        dwi_path, s,
                        bvec_p, bvec_p.exists(),
                        bval_p, bval_p.exists(),
                    )
                    continue
                valid_dwi.append(str(dwi_p.absolute()))
                valid_bvec.append(str(bvec_p.absolute()))
                valid_bval.append(str(bval_p.absolute()))
            else:
                # .zip / DICOM: no host-side sidecars until dcm2niix runs in-container.
                valid_dwi.append(str(dwi_p.absolute()))
                valid_bvec.append("")
                valid_bval.append("")

        if not valid_dwi:
            raise ValueError(
                f"No valid DWI files remain for subject '{s}' after filtering: "
                f"all {len(candidates)} matched DWI(s) are missing a bvec or bval sidecar. "
                f"Pattern: {dwi_pattern}"
            )

        all_dwi.append(valid_dwi)
        all_bvec.append(valid_bvec)
        all_bval.append(valid_bval)

    return all_dwi, all_bvec, all_bval


def _resolve_recon_outputs_dir(recon_dir: str, subject: str) -> str:
    """Return the absolute host path to recon-all outputs/ for this subject.

    Accepts two forms:
      - Pattern form:  contains {subject}, e.g. "/data/{subject}/recon-all/outputs"
        → substituted directly, no extra path segments appended.
      - Root form:     no {subject} placeholder, e.g. "/data/recon_multi"
        → appends /{subject}/recon-all/outputs (legacy behaviour).

    Always returns an absolute path so Snakemake does not warn about relative paths.
    """
    if "{subject}" in recon_dir:
        resolved = recon_dir.format(subject=subject)
    else:
        resolved = str(Path(recon_dir) / subject / "recon-all" / "outputs")
    return str(Path(resolved).resolve())


def _recon_exists(recon_dir: Optional[str], subject: str) -> bool:
    """Return True if aparc+aseg.mgz is present for this subject under recon_dir."""
    if not recon_dir:
        return False
    p = Path(_resolve_recon_outputs_dir(recon_dir, subject)) / "mri" / _RECON_SENTINEL
    return p.exists()


class TraculaProcessor:

    def __init__(self):
        self.current_dir = Path(__file__).parent
        self.leukoquant_dir = self.current_dir.parent
        self.leukoquant_parent_dir = self.leukoquant_dir.parent
        self.workflow_dir = self.leukoquant_dir / "workflow" / "workflows"

    def run_tracula(self,
                    subject_input: str,
                    output_dir: str,
                    dwi_pattern: str,
                    freesurfer_recon_dir: Optional[str] = None,
                    t1_pattern: Optional[str] = None,
                    bvecs_pattern: Optional[str] = None,
                    bvals_pattern: Optional[str] = None,
                    scheduler: str = "local",
                    cores: int = 1,
                    parcellation: str = "freesurfer",
                    brain_parcellation: Optional[str] = None,
                    mapping_file: Optional[str] = None,
                    keep_intermediate: bool = False,
                    force_rules: Optional[List[str]] = None,
                    verbose: bool = False) -> Tuple[str, str]:
        """Run TRACULA for one or more subjects.

        Args:
            subject_input:        Bare subject ID or path to a text file (one ID per line).
            output_dir:           Directory for TRACULA outputs and config.
            dwi_pattern:          Explicit path or {subject} glob pattern for DWI files.
            freesurfer_recon_dir: Path to existing recon-all outputs. Accepts two forms:
                                  - Pattern: contains {subject}, e.g. "/data/{subject}/recon-all/outputs"
                                  - Root dir: no {subject}, e.g. "/data/recon_multi"
                                    (appends /{subject}/recon-all/outputs automatically).
                                  Optional - if omitted or outputs are missing, recon-all is run
                                  automatically and t1_pattern is required.
            t1_pattern:           Explicit path or {subject} glob pattern for T1 files.
                                  Required when any subject is missing recon-all outputs.
            bvecs_pattern:        Optional bvecs file or {subject} pattern.
            bvals_pattern:        Optional bvals file or {subject} pattern.
        """
        wf_file = self.workflow_dir / "tracula_workflow.smk"

        subjects = read_subjects(subject_input)

        # Parse parcellations from comma-separated string or list
        if isinstance(parcellation, str):
            parcellations = [p.strip() for p in parcellation.split(",") if p.strip()]
        else:
            parcellations = list(parcellation)
        non_fs_parcellations = [p for p in parcellations if p != "freesurfer"]

        # Validate non-gif non-freesurfer parcellations require explicit mapping
        if non_fs_parcellations and any(p != "gif" for p in non_fs_parcellations) and not mapping_file:
            non_gif_other = [p for p in non_fs_parcellations if p != "gif"]
            raise ValueError(
                f"--mapping-file is required when --parcellation is '{non_gif_other[0]}' "
                f"(it is only auto-resolved for 'gif')."
            )

        out_path = Path(output_dir)
        out_path.mkdir(parents=True, exist_ok=True)

        # ── determine which subjects need recon-all ──────────────────────────
        # Per-subject: decide the host path to the recon-all outputs/ dir.
        # If the user supplied freesurfer_recon_dir and aparc+aseg.mgz exists
        # there, use that dir.  Otherwise fall back to the standard output path
        # where the recon_all module will produce results.
        recon_dirs: Dict[str, str] = {}        # {subject: host outputs/ path}
        needs_recon: Set[str] = set()

        for s in subjects:
            if _recon_exists(freesurfer_recon_dir, s):
                recon_dirs[s] = _resolve_recon_outputs_dir(freesurfer_recon_dir, s)  # type: ignore[arg-type]  # None case excluded by _recon_exists
            else:
                # Will be produced by the recon_all module
                recon_dirs[s] = str(out_path / s / "recon-all" / "outputs")
                needs_recon.add(s)

        if needs_recon and not t1_pattern:
            missing = sorted(needs_recon)
            raise ValueError(
                f"recon-all outputs not found for subject(s): {missing}\n"
                "--t1 is required to run recon-all automatically."
            )

        # ── resolve file paths ───────────────────────────────────────────────
        # DWI: collect ALL matching files per subject; filter by bvec/bval presence.
        dwi_files_per_subject, bvec_files_per_subject, bval_files_per_subject = (
            _discover_dwi_with_sidecars(dwi_pattern, subjects, bvecs_pattern, bvals_pattern)
        )

        t1_files: List[str] = []
        if needs_recon:
            for s in subjects:
                t1_files.append(
                    resolve_subject_pattern(t1_pattern, s) if s in needs_recon else ""
                )
            for s in needs_recon:
                f = t1_files[subjects.index(s)]
                if not Path(f).exists():
                    raise FileNotFoundError(f"T1 file not found for subject '{s}': {f}")

        # ── per-subject brain parcellation maps ──────────────────────────────
        # Build per-subject maps so that tracula_workflow.smk can resolve each
        # subject's parcellation input independently, supporting multi-subject runs.
        brain_parcellation_map: Dict[str, str] = {}       # {subject: host_path}
        brain_parcellation_sing_map: Dict[str, str] = {}  # {subject: container_path}

        if non_fs_parcellations:
            for s in subjects:
                if "gif" in non_fs_parcellations:
                    if brain_parcellation:
                        # Explicit file supplied: bind at /input/ inside the container.
                        bp_abs = str(Path(brain_parcellation).absolute())
                        bp_suffix = "".join(Path(bp_abs).suffixes)
                        brain_parcellation_map[s] = bp_abs
                        brain_parcellation_sing_map[s] = f"/input/input_parcellation{bp_suffix}"
                    else:
                        # Auto-resolve from GIF output under out_path.
                        # out_path is bound as /output, so no extra bind entry is needed.
                        gif_out_dir = out_path / s / "gif" / "outputs"
                        gif_matches = sorted(gif_out_dir.glob("*_NeuroMorph_Parcellation.nii.gz"))
                        if not gif_matches:
                            raise FileNotFoundError(
                                f"GIF parcellation for subject '{s}' not found at "
                                f"{gif_out_dir}/*_NeuroMorph_Parcellation.nii.gz. "
                                f"Run 'process-gif' first, or supply --brain-parcellation explicitly."
                            )
                        if len(gif_matches) > 1:
                            logger.warning(
                                f"Multiple GIF parcellation files found for subject '{s}'; "
                                f"using: {gif_matches[0]}"
                            )
                        gif_parc = gif_matches[0]
                        brain_parcellation_map[s] = str(gif_parc)
                        brain_parcellation_sing_map[s] = (
                            f"/output/{s}/gif/outputs/{gif_parc.name}"
                        )
                else:
                    # Other non-freesurfer parcellation: brain_parcellation is required.
                    # Validated above; bind at a fixed container path.
                    bp_abs = str(Path(brain_parcellation).absolute())
                    bp_suffix = "".join(Path(bp_abs).suffixes)
                    brain_parcellation_map[s] = bp_abs
                    brain_parcellation_sing_map[s] = f"/input/input_parcellation{bp_suffix}"

        # ── per-subject GIF brain mask maps ──────────────────────────────────
        # GIF brain mask is needed for custom (non-freesurfer) parcellations.
        brain_mask_map: Dict[str, str] = {}
        brain_mask_sing_map: Dict[str, str] = {}

        if "gif" in non_fs_parcellations:
            for s in subjects:
                # Glob for the NeuroMorph_Brain file from GIF outputs.
                # out_path is bound as /output, so no extra bind entry is needed.
                gif_out_dir = out_path / s / "gif" / "outputs"
                mask_matches = sorted(gif_out_dir.glob("*_NeuroMorph_Brain.nii.gz"))
                if not mask_matches:
                    raise FileNotFoundError(
                        f"GIF brain mask for subject '{s}' not found at "
                        f"{gif_out_dir}/*_NeuroMorph_Brain.nii.gz. "
                        f"Run 'process-gif' first."
                    )
                if len(mask_matches) > 1:
                    logger.warning(
                        f"Multiple GIF brain mask files found for subject '{s}'; "
                        f"using: {mask_matches[0]}"
                    )
                gif_brain = mask_matches[0]
                brain_mask_map[s] = str(gif_brain)
                brain_mask_sing_map[s] = f"/output/{s}/gif/outputs/{gif_brain.name}"

        # ── container / license ──────────────────────────────────────────────
        scratch_dir = Path("/scratch0") if Path("/scratch0").exists() else out_path / ".scratch"
        scratch_dir.mkdir(parents=True, exist_ok=True)
        scratch_dir_singularity = "/scratch0"

        fs_license_dir = _resolve_fs_license(self.leukoquant_dir, verbose=verbose)
        fs_license_singularity = "/license/license.txt"

        containers_folder = self.leukoquant_dir / "workflow" / "containers"
        sif_path = str(containers_folder / FREESURFER_SIF_FILENAME)
        ensure_container(sif_path, filename=FREESURFER_SIF_HF_PATH)

        # ── singularity bind mounts ──────────────────────────────────────────
        # DWI, bvec, bval are nested lists (one inner list per subject).
        # Directory-level binds automatically cover JSON sidecars in the same
        # directory - no explicit JSON bind loop needed.
        dwi_sing,   dwi_binds   = dir_level_bind_file_lists(dwi_files_per_subject,  "input_dwi")
        bvecs_sing, bvecs_binds = dir_level_bind_file_lists(bvec_files_per_subject, "input_bvecs")
        bvals_sing, bvals_binds = dir_level_bind_file_lists(bval_files_per_subject, "input_bvals")

        # T1 binds for subjects that need recon-all (flat list, empty for others).
        t1_sing: List[str]
        t1_binds: List[str]
        t1_sing, t1_binds = dir_level_bind_files(t1_files, "input_t1")

        # Bind each pre-existing recon-all outputs directory.
        recon_dirs_sing: Dict[str, str] = {}
        recon_binds: List[str] = []
        for i, s in enumerate(subjects):
            if s not in needs_recon:
                host_path = recon_dirs[s]
                sing_path = f"/recon_input_{i}"
                recon_dirs_sing[s] = sing_path
                recon_binds.append(f"{host_path}:{sing_path}")
            else:
                # Will be produced inside /output - no extra bind needed.
                recon_dirs_sing[s] = f"/output/{s}/recon-all/outputs"

        # Consolidate all per-subject input binds to a single common-ancestor entry
        # per filesystem root, keeping --bind well within ARG_MAX for large cohorts.
        input_binds = dwi_binds + bvecs_binds + bvals_binds + t1_binds + recon_binds
        input_binds, remap = consolidate_bind_entries(input_binds)
        dwi_sing   = [[remap(p) for p in subj] for subj in dwi_sing]
        bvecs_sing = [[remap(p) for p in subj] for subj in bvecs_sing]
        bvals_sing = [[remap(p) for p in subj] for subj in bvals_sing]
        t1_sing    = [remap(p) for p in t1_sing]
        recon_dirs_sing = {s: remap(path) for s, path in recon_dirs_sing.items()}

        bind_parts = (
            input_binds
            + [f"{str(self.leukoquant_parent_dir.absolute())}:/leukoquant"]
            + [f"{str(out_path.absolute())}:/output"]
            + [f"{str(scratch_dir.absolute())}:{scratch_dir_singularity}"]
            + [f"{str(fs_license_dir.absolute())}:{fs_license_singularity}"]
        )

        # Bind explicit brain_parcellation for non-gif parcellations (gif auto-resolved
        # paths live under out_path which is already bound as /output).
        if (brain_parcellation and non_fs_parcellations
                and any(p != "gif" for p in non_fs_parcellations)):
            bp_path = Path(brain_parcellation)
            bp_suffix = "".join(bp_path.suffixes)
            bind_parts.append(f"{str(bp_path.absolute())}:/input/input_parcellation{bp_suffix}")

        if mapping_file:
            mf_path = Path(mapping_file)
            mf_suffix = "".join(mf_path.suffixes)
            bind_parts.append(f"{str(mf_path.absolute())}:/input/parcellation_mapping{mf_suffix}")
        # gif bundled mapping is accessible via /leukoquant mount - no extra bind needed.

        singularity_bind = "--bind " + ",".join(bind_parts)

        # ── snakemake config ─────────────────────────────────────────────────
        # DWI / bvec / bval are stored as List[List[str]] - one inner list per
        # subject, potentially containing multiple scans per session.
        config_dict = {
            "subjects":                subjects,
            "dwi_imgs":                dwi_files_per_subject,
            "dwi_paths_singularity":   dwi_sing,
            "bvecs_paths":             bvec_files_per_subject,
            "bvecs_paths_singularity": bvecs_sing,
            "bvals_paths":             bval_files_per_subject,
            "bvals_paths_singularity": bvals_sing,
            # Per-subject recon-all info consumed by the workflow
            "recon_dirs":              recon_dirs,
            "recon_dirs_singularity":  recon_dirs_sing,
            "needs_recon":             list(needs_recon),
            # T1 inputs for subjects that need recon-all (empty list when not needed)
            "t1_files":                t1_files,
            "t1_files_singularity":    t1_sing,
            "output_dir":              str(out_path.absolute()),
            "output_dir_singularity":  "/output",
            "scratch_dir_singularity": scratch_dir_singularity,
            "keep_intermediate":       keep_intermediate,
            "leukoquant_parent_dir":   str(self.leukoquant_parent_dir.absolute()),
            "container_name":          "freesurfer_unified_container",
        }

        # Always pass the full list; the workflow fans out via {parcellation} wildcard.
        config_dict["parcellations"] = parcellations

        # Pass per-subject brain parcellation maps for non-freesurfer parcellations.
        if brain_parcellation_map:
            config_dict["brain_parcellation_map"]      = brain_parcellation_map
            config_dict["brain_parcellation_sing_map"] = brain_parcellation_sing_map

        # Pass per-subject GIF brain mask maps.
        if brain_mask_map:
            config_dict["brain_mask_map"]      = brain_mask_map
            config_dict["brain_mask_sing_map"] = brain_mask_sing_map

        if mapping_file:
            mf_path = Path(mapping_file)
            mf_suffix = "".join(mf_path.suffixes)
            config_dict["mapping_file"] = str(mf_path.absolute())
            config_dict["mapping_file_singularity"] = f"/input/parcellation_mapping{mf_suffix}"
        elif "gif" in non_fs_parcellations:
            # Bundled gif→freesurfer mapping; accessible via /leukoquant mount.
            mapping_default = self.leukoquant_parent_dir / "mappings/gif_to_freesurfer.csv"
            config_dict["mapping_file"] = str(mapping_default.absolute())
            # mapping_file_singularity is auto-computed by the workflow from the relative path.

        config_file = str(out_path.absolute() / "tracula_config.yaml")
        with open(config_file, "w") as f:
            yaml.dump(config_dict, f, default_flow_style=False)

        # Validate DWI paths are in config for all subjects
        for i, s in enumerate(subjects):
            if i < len(dwi_sing):
                if not dwi_sing[i]:
                    logger.warning(f"Subject {s}: DWI singularity paths are empty!")

        # ── snakemake command ────────────────────────────────────────────────
        snakemake_cmd = [
            "snakemake", "-s", str(wf_file),
            "--directory", str(out_path),
            "--cores", str(cores),
            "--jobs", "unlimited",
            "--drop-metadata",
            "--immediate-submit", "--notemp",
            # Some compute nodes mount NFS lazily; the default 5 s wait is too
            # short for cold mounts and yields spurious "missing input" errors.
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

        add_forcerun_args(snakemake_cmd, force_rules)
        # Explicitly target tracula's own `rule all` so that imported module
        # rules (recon_all_all from the recon-all module) don't become the
        # default target when the module is activated.
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
                f"TRACULA failed: {e}\nError Details:\n{e.stderr}"
            )


def apply_tracula(subject_input: Optional[str] = None,
                  output_dir: Optional[str] = None,
                  dwi_pattern: Optional[str] = None,
                  freesurfer_recon_dir: Optional[str] = None,
                  t1_pattern: Optional[str] = None,
                  bvecs_pattern: Optional[str] = None,
                  bvals_pattern: Optional[str] = None,
                  scheduler: str = "local",
                  cores: int = 1,
                  keep_intermediate: bool = False,
                  force_rules: Optional[List[str]] = None,
                  verbose: bool = False,
                  parcellation: str = "freesurfer",
                  brain_parcellation: Optional[str] = None,
                  mapping_file: Optional[str] = None,
                  config_yaml: Optional[str] = None) -> dict:
    """Run TRACULA and return a summary dict.

    CLI / caller arguments take priority over values in ``config_yaml`` when
    both are provided.
    """
    yaml_cfg = load_yaml_config(config_yaml)
    subject_input = first_truthy(subject_input, yaml_cfg.get("subject_input"), yaml_cfg.get("subject"))
    output_dir    = first_truthy(output_dir,    yaml_cfg.get("output_dir"))
    dwi_pattern   = first_truthy(dwi_pattern,   yaml_cfg.get("dwi_pattern"), yaml_cfg.get("dwi"))
    freesurfer_recon_dir = first_truthy(freesurfer_recon_dir, yaml_cfg.get("freesurfer_recon_dir"))
    t1_pattern    = first_truthy(t1_pattern,    yaml_cfg.get("t1_pattern"), yaml_cfg.get("t1"))
    bvecs_pattern = first_truthy(bvecs_pattern, yaml_cfg.get("bvecs_pattern"), yaml_cfg.get("bvecs"))
    bvals_pattern = first_truthy(bvals_pattern, yaml_cfg.get("bvals_pattern"), yaml_cfg.get("bvals"))
    scheduler     = first_truthy(scheduler,     yaml_cfg.get("scheduler")) or "local"
    cores_raw     = first_truthy(cores,         yaml_cfg.get("cores"))
    cores         = int(cores_raw) if cores_raw is not None else 1
    keep_intermediate = bool(yaml_cfg.get("keep_intermediate", False)) or bool(keep_intermediate)
    parcellation  = first_truthy(parcellation,  yaml_cfg.get("parcellation")) or "freesurfer"
    brain_parcellation = first_truthy(brain_parcellation, yaml_cfg.get("brain_parcellation"))
    mapping_file  = first_truthy(mapping_file,  yaml_cfg.get("mapping_file"))

    if not subject_input:
        raise ValueError("subject_input is required (via argument or config_yaml)")
    if not output_dir:
        raise ValueError("output_dir is required (via argument or config_yaml)")
    if not dwi_pattern:
        raise ValueError("dwi_pattern is required (via argument or config_yaml)")

    if verbose:
        logging.basicConfig(level=logging.INFO)
    proc = TraculaProcessor()
    try:
        job_id, results_dir = proc.run_tracula(
            subject_input=subject_input,
            output_dir=output_dir,
            dwi_pattern=dwi_pattern,
            freesurfer_recon_dir=freesurfer_recon_dir,
            t1_pattern=t1_pattern,
            bvecs_pattern=bvecs_pattern,
            bvals_pattern=bvals_pattern,
            scheduler=scheduler,
            cores=cores,
            keep_intermediate=keep_intermediate,
            parcellation=parcellation,
            brain_parcellation=brain_parcellation,
            mapping_file=mapping_file,
            force_rules=force_rules,
            verbose=verbose,
        )
        print("✅ TRACULA job submitted successfully")
        return {"success": True, "job_id": job_id, "results_dir": results_dir,
                "subject_input": subject_input, "output_dir": str(Path(output_dir).absolute())}
    except Exception as e:
        print(f"❌ TRACULA job failed:\n{e}", file=sys.stderr)
        return {"success": False, "error": str(e),
                "subject_input": subject_input, "output_dir": str(Path(output_dir).absolute())}
