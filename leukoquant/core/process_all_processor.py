#!/usr/bin/env python3
"""
process_all processor - runs the full pipeline for one or more subjects.

All subjects are submitted in a single Snakemake invocation so that DAG
construction, Singularity image discovery, and conda environment setup are
paid only once. The smk workflow fans out per-subject jobs via wildcards
(the same pattern used by every other individual tool processor).
"""
import logging
import os
import sys
from pathlib import Path
from typing import List, Optional, Tuple

import yaml
from snakemake.api import SnakemakeApi

# snakemake.settings.types is the canonical location in Snakemake 8+/9+.
# A small number of intermediate 8.x builds re-exported from snakemake.settings;
# try the canonical path first so we always get the most specific module.
try:
    from snakemake.settings.types import (
        ConfigSettings,
        DAGSettings,
        DeploymentMethod,
        DeploymentSettings,
        ExecutionSettings,
        GroupSettings,
        OutputSettings,
        RemoteExecutionSettings,
        ResourceSettings,
        SchedulingSettings,
        SharedFSUsage,
        StorageSettings,
    )
except ImportError:
    from snakemake.settings import (  # type: ignore[no-redef]
        ConfigSettings,
        DAGSettings,
        DeploymentMethod,
        DeploymentSettings,
        ExecutionSettings,
        GroupSettings,
        OutputSettings,
        RemoteExecutionSettings,
        ResourceSettings,
        SchedulingSettings,
        SharedFSUsage,
        StorageSettings,
    )

from snakemake_executor_plugin_sge import ExecutorSettings as SgeExecutorSettings

from leukoquant.utils.subject_utils import read_subjects, resolve_subject_pattern
from leukoquant.utils.external_utils import _resolve_gif_home, _resolve_fs_license
from leukoquant.utils.bind_utils import consolidate_bind_entries
from leukoquant.utils.snakemake_utils import load_yaml_config, first_truthy

logger = logging.getLogger(__name__)


def _parse_covariates(value: Optional[str]) -> List[str]:
    """Split a comma-separated covariate string into a list, or return [] if None."""
    if not value:
        return []
    return [c.strip() for c in value.split(",") if c.strip()]


def _parse_poly_terms(value: Optional[str]) -> List[str]:
    """Convert CLI polynomial terms (e.g. 'age:2,bmi:3') to underscore notation ('age_2,bmi_3').

    Accepts colon (``age:2``), caret (``age^2``), or bare name (``age`` → ``age_2``).
    Returns a list of formatted strings suitable for z_score_utils.generate_design_matrix.
    """
    if not value:
        return []
    result = []
    for term in value.split(","):
        term = term.strip()
        if not term:
            continue
        if ":" in term:
            name, power = term.split(":", 1)
        elif "^" in term:
            name, power = term.split("^", 1)
        else:
            name, power = term, "2"
        result.append(f"{name.strip()}_{power.strip()}")
    return result


def _resolve_optional(pattern: Optional[str], subject: str) -> str:
    if not pattern:
        return ""
    try:
        return resolve_subject_pattern(pattern, subject)
    except FileNotFoundError:
        return ""


def _is_nifti(path: Path) -> bool:
    """Return True only for files that already carry bvec/bval sidecars.

    .zip archives and DICOM files are converted by dcm2niix inside the
    container; their sidecars don't exist on the host until after conversion.
    """
    name = path.name.lower()
    return name.endswith('.nii.gz') or name.endswith('.nii')


class ProcessAllProcessor:
    """Processor to run the entire pipeline end-to-end for one or more subjects."""

    def __init__(self, external_dir: Optional[str] = None):
        self.current_dir = Path(__file__).parent
        self.leukoquant_dir = self.current_dir.parent
        self.leukoquant_parent_dir = self.leukoquant_dir.parent

        if external_dir is None:
            external_dir = self.leukoquant_dir / "external"
        self.external_dir = Path(external_dir)
        self.workflow_dir = self.leukoquant_dir / "workflow" / "workflows"

    @staticmethod
    def _register_file_bind(
        host_file: Path,
        subject_idx: int,
        dir_map: dict,
        bind_entries: list,
    ) -> str:
        """Register a directory-level Singularity bind for host_file (if not yet registered).

        Rather than adding one bind entry per file, this method binds the parent directory
        once and derives the container file path from it. The total bind count stays
        proportional to unique host directories rather than unique files - critical for
        large cohorts where per-file binds exhaust Apptainer's file-descriptor budget and
        cause a "bad file descriptor" error.

        Returns the in-container path for host_file.
        """
        host_dir_str = str(host_file.parent)
        if host_dir_str not in dir_map:
            container_dir = f"/input_{subject_idx}_{len(dir_map)}"
            dir_map[host_dir_str] = container_dir
            bind_entries.append(f"{host_dir_str}:{container_dir}")
        return f"{dir_map[host_dir_str]}/{host_file.name}"

    def run_process_all(self,
                        subject_input: str,
                        t1_pattern: str,
                        flair_pattern: str,
                        dwi_pattern: str,
                        bvecs_pattern: str,
                        bvals_pattern: str,
                        mask_pattern: Optional[str],
                        output_dir: str,
                        scheduler: str = "local",
                        cores: int = 1,
                        skip_zscore: bool = False,
                        healthy_list: Optional[str] = None,
                        demographics: Optional[str] = None,
                        covariates: Optional[str] = None,
                        poly_terms: Optional[str] = None,
                        parcellation: str = "freesurfer",
                        force_rules: Optional[List[str]] = None,
                        verbose: bool = False) -> Tuple[str, str]:
        """Run process-all for one or more subjects in a single Snakemake call.

        All subjects share one DAG-build pass and one container-discovery step,
        which removes the per-subject Snakemake overhead that a for-loop approach
        would incur.

        force_rules: Rule name(s) to force-rerun (e.g. ["extract_metrics_gif"])
            even if their outputs already exist and are up to date. Snakemake
            scopes this to the named rule(s) plus anything downstream, so
            subjects whose earlier stages (recon-all, DTI, GIF, tractography,
            ...) are already complete have those stages correctly skipped --
            only the forced rule(s) recompute. This is how to re-run just the
            metrics stage for a subset of already-processed subjects without
            standing up a separate process-metrics invocation.
        """
        subjects = read_subjects(subject_input)

        # If healthy list is provided, ensure those subjects are in the processing queue
        # so Snakemake can build their metrics before z-score calculation uses them.
        healthy_list_resolved = None
        healthy_subs: List[str] = []
        if healthy_list:
            healthy_list_resolved = str(Path(healthy_list).resolve())
            healthy_subs = read_subjects(healthy_list_resolved)
            for hs in healthy_subs:
                if hs not in subjects:
                    subjects.append(hs)

        # target_subjects = all subjects minus the healthy cohort
        healthy_set = set(healthy_subs)
        target_subjects = [s for s in subjects if s not in healthy_set]

        # ── Validate required environment variables ──────────────────────────
        fs_license_dir = _resolve_fs_license(self.leukoquant_dir, verbose=verbose)
        gif_software_path = _resolve_gif_home(self.leukoquant_dir / "external" / "gif")

        base_out = Path(output_dir).absolute()
        base_out.mkdir(parents=True, exist_ok=True)

        # ── Resolve per-subject inputs ───────────────────────────────────────
        t1_files:    List[str] = []
        flair_files: List[str] = []
        # DWI, bvec, bval are nested: one inner list per subject, each inner
        # list containing all valid scan paths for that session.
        dwi_files:   List[List[str]] = []
        bvecs_list:  List[List[str]] = []
        bvals_list:  List[List[str]] = []
        mask_list:   List[str] = []

        # Singularity container paths - T1/FLAIR/mask are one per subject;
        # DWI/bvec/bval mirror the nested structure (one inner list per subject).
        t1_sing_paths:    List[str] = []
        flair_sing_paths: List[str] = []
        dwi_sing_paths:   List[List[str]] = []
        bvecs_sing_paths: List[List[str]] = []
        bvals_sing_paths: List[List[str]] = []
        mask_sing_paths:  List[str] = []

        # Bind entries: "host_path:container_path"
        bind_entries: List[str] = []

        # z_score_workflow.smk (shared with the standalone process-zscore CLI)
        # translates demographics_csv / healthy_subjects_list to their in-container
        # path via config["singularity_binds"] -- register their parent directories
        # now so they're covered by the same consolidation pass as the T1/DWI binds,
        # or z_score_calc.sh fails with FileNotFoundError trying to read the raw
        # host path from inside the container.
        zscore_dir_map: dict[str, str] = {}
        if demographics:
            self._register_file_bind(Path(demographics).resolve(), "zscore", zscore_dir_map, bind_entries)
        if healthy_list:
            self._register_file_bind(Path(healthy_list).resolve(), "zscore", zscore_dir_map, bind_entries)

        for i, subject in enumerate(subjects):
            # T1 and FLAIR: take the first match (single file per session).
            t1    = resolve_subject_pattern(t1_pattern,    subject)
            flair = resolve_subject_pattern(flair_pattern, subject)
            mask  = _resolve_optional(mask_pattern, subject)

            t1_p    = Path(t1).absolute()
            flair_p = Path(flair).absolute()

            # Directory-level bind map for this subject: str(host_dir) → container_dir.
            # Each unique host directory is bound once; all files within it share the same
            # container directory path. This keeps the total Singularity bind count
            # proportional to unique directories rather than unique files, preventing the
            # Apptainer "bad file descriptor" error that appears when many per-file bind
            # entries exhaust the process's open file-descriptor limit on large cohorts.
            subj_dir_map: dict[str, str] = {}

            t1_sing    = self._register_file_bind(t1_p,    i, subj_dir_map, bind_entries)
            flair_sing = self._register_file_bind(flair_p, i, subj_dir_map, bind_entries)

            mask_sing = ""
            if mask:
                mask_p = Path(mask).absolute()
                mask_sing = self._register_file_bind(mask_p, i, subj_dir_map, bind_entries)

            # DWI: collect ALL matching files; filter by bvec/bval presence.
            dwi_candidates: List[str] = resolve_subject_pattern(dwi_pattern, subject, all=True)  # type: ignore[assignment]
            if not dwi_candidates:
                raise FileNotFoundError(
                    f"No DWI files found for subject '{subject}' with pattern: {dwi_pattern}"
                )

            subj_dwi:       List[str] = []
            subj_bvec:      List[str] = []
            subj_bval:      List[str] = []
            subj_dwi_sing:  List[str] = []
            subj_bvec_sing: List[str] = []
            subj_bval_sing: List[str] = []

            # Single DWI + explicit bvec/bval patterns → honour those patterns
            # (backward compatibility with callers that pass separate patterns).
            if len(dwi_candidates) == 1 and (bvecs_pattern or bvals_pattern):
                dwi_p = Path(dwi_candidates[0]).absolute()
                dwi_sing = self._register_file_bind(dwi_p, i, subj_dir_map, bind_entries)
                # JSON sidecar lives in the same BIDS directory as the DWI; registering it
                # reuses the existing directory bind entry rather than adding a new one.
                json_p = dwi_p.with_name(dwi_p.name.split(".")[0] + ".json")
                if json_p.exists():
                    self._register_file_bind(json_p, i, subj_dir_map, bind_entries)

                bvecs = _resolve_optional(bvecs_pattern, subject)
                bvals = _resolve_optional(bvals_pattern,  subject)

                bvec_sing = ""
                if bvecs:
                    bvec_p = Path(bvecs).absolute()
                    bvec_sing = self._register_file_bind(bvec_p, i, subj_dir_map, bind_entries)

                bval_sing = ""
                if bvals:
                    bval_p = Path(bvals).absolute()
                    bval_sing = self._register_file_bind(bval_p, i, subj_dir_map, bind_entries)

                subj_dwi.append(str(dwi_p))
                subj_bvec.append(bvecs)
                subj_bval.append(bvals)
                subj_dwi_sing.append(dwi_sing)
                subj_bvec_sing.append(bvec_sing)
                subj_bval_sing.append(bval_sing)

            else:
                # Multi-DWI (or single DWI without explicit sidecar patterns):
                # auto-discover bvec/bval by matching the file stem.
                for dwi_candidate in dwi_candidates:
                    dwi_p = Path(dwi_candidate).absolute()
                    stem = dwi_p.name.split(".")[0]

                    if _is_nifti(dwi_p):
                        bvec_p = dwi_p.with_name(stem + ".bvec")
                        bval_p = dwi_p.with_name(stem + ".bval")
                        if not bvec_p.exists() or not bval_p.exists():
                            logger.warning(
                                "Skipping DWI '%s' for subject '%s': "
                                "missing bvec (%s exists=%s) or bval (%s exists=%s).",
                                dwi_candidate, subject,
                                bvec_p, bvec_p.exists(),
                                bval_p, bval_p.exists(),
                            )
                            continue
                        # JSON sidecar: in BIDS, it lives in the same directory as the DWI,
                        # so registering it typically reuses the existing bind entry.
                        json_p = dwi_p.with_name(stem + ".json")
                        if json_p.exists():
                            self._register_file_bind(json_p, i, subj_dir_map, bind_entries)
                        bvec_sing = self._register_file_bind(bvec_p, i, subj_dir_map, bind_entries)
                        bval_sing = self._register_file_bind(bval_p, i, subj_dir_map, bind_entries)
                        bvec_host = str(bvec_p)
                        bval_host = str(bval_p)
                    else:
                        # .zip / DICOM: dcm2niix runs inside the container and produces
                        # bvec/bval there. No host-side sidecars exist yet.
                        bvec_sing = ""
                        bval_sing = ""
                        bvec_host = ""
                        bval_host = ""

                    dwi_sing = self._register_file_bind(dwi_p, i, subj_dir_map, bind_entries)

                    subj_dwi.append(str(dwi_p))
                    subj_bvec.append(bvec_host)
                    subj_bval.append(bval_host)
                    subj_dwi_sing.append(dwi_sing)
                    subj_bvec_sing.append(bvec_sing)
                    subj_bval_sing.append(bval_sing)

                if not subj_dwi:
                    raise ValueError(
                        f"No valid DWI files remain for subject '{subject}' after filtering: "
                        f"all {len(dwi_candidates)} matched DWI(s) are missing a bvec or bval "
                        f"sidecar. Pattern: {dwi_pattern}"
                    )

            t1_files.append(t1)
            flair_files.append(flair)
            dwi_files.append(subj_dwi)
            bvecs_list.append(subj_bvec)
            bvals_list.append(subj_bval)
            mask_list.append(mask)

            t1_sing_paths.append(t1_sing)
            flair_sing_paths.append(flair_sing)
            dwi_sing_paths.append(subj_dwi_sing)
            bvecs_sing_paths.append(subj_bvec_sing)
            bvals_sing_paths.append(subj_bval_sing)
            mask_sing_paths.append(mask_sing)

        # ── Consolidate per-subject input binds to a single common-ancestor bind ──
        # Reduces O(subjects × modalities) directory entries to O(1) by finding
        # the longest common ancestor of all subject input directories.
        # consolidate_bind_entries handles multiple filesystem roots (e.g. Windows
        # drive letters) by producing one /input_N mount per root group.
        bind_entries, _remap = consolidate_bind_entries(bind_entries)
        t1_sing_paths    = [_remap(p) for p in t1_sing_paths]
        flair_sing_paths = [_remap(p) for p in flair_sing_paths]
        dwi_sing_paths   = [[_remap(p) for p in subj] for subj in dwi_sing_paths]
        bvecs_sing_paths = [[_remap(p) for p in subj] for subj in bvecs_sing_paths]
        bvals_sing_paths = [[_remap(p) for p in subj] for subj in bvals_sing_paths]
        mask_sing_paths  = [_remap(p) for p in mask_sing_paths]
        # host_dir -> post-consolidation container_dir, for z_score_workflow.smk's
        # translate_path() calls (demographics_csv / healthy_subjects_list).
        zscore_singularity_binds = {
            host_dir: _remap(container_dir)
            for host_dir, container_dir in zscore_dir_map.items()
        }

        # ── Shared bind entries (constant across subjects) ───────────────────
        bind_entries += [
            f"{str(self.leukoquant_parent_dir.absolute())}:/leukoquant",
            f"{str(base_out.absolute())}:/output",
            f"{str(fs_license_dir.absolute())}:/license/license.txt",
            f"{str(gif_software_path.absolute())}:/GIF",
        ]

        scratch_dir = Path("/scratch0") if Path("/scratch0").exists() else base_out / ".scratch"
        scratch_dir.mkdir(parents=True, exist_ok=True)
        bind_entries.append(f"{str(scratch_dir.absolute())}:/scratch0")

        # Common-ancestor consolidation above reduces bind_entries to O(1) entries
        # (one per filesystem root plus shared dirs), so the --bind string is short
        # enough to pass directly via --apptainer-args without hitting ARG_MAX.
        singularity_bind = "--bind " + ",".join(bind_entries)

        # ── Write config ─────────────────────────────────────────────────────
        # z-score config keys (optional)
        demographics_resolved = str(Path(demographics).resolve()) if demographics else ""

        config_dict = {
            "subjects": subjects,
            "target_subjects": target_subjects,
            "t1_files": t1_files,
            "t1_files_singularity": t1_sing_paths,
            "flair_files": flair_files,
            "flair_files_singularity": flair_sing_paths,
            # DWI / bvec / bval are List[List[str]] - one inner list per subject.
            "dwi_files": dwi_files,
            "dwi_paths_singularity": dwi_sing_paths,
            "bvecs_list": bvecs_list,
            "bvecs_paths_singularity": bvecs_sing_paths,
            "bvals_list": bvals_list,
            "bvals_paths_singularity": bvals_sing_paths,
            "mask_files": mask_list,
            "mask_files_singularity": mask_sing_paths,
            "output_dir": str(base_out),
            "leukoquant_parent_dir": str(self.leukoquant_parent_dir.absolute()),
            "containers": {"freesurfer_unified_container": "docker://stylianosc/compsvd:freesurfer"},
            "singularity_binds": zscore_singularity_binds,
            "skip_zscore": skip_zscore,
            "healthy_subjects_list": healthy_list_resolved,
            "demographics_csv": demographics_resolved,
            "covariates": _parse_covariates(covariates),
            "polynomial_terms": _parse_poly_terms(poly_terms),
            "verbose": verbose,
            # Write as a list so process_all_workflow.smk can fan out via the
            # {parcellation} wildcard without string-splitting at workflow load time.
            "parcellations": [p.strip() for p in parcellation.split(",") if p.strip()]
                              if isinstance(parcellation, str) else list(parcellation),
        }

        config_file = str(base_out / "process_all_config.yaml")
        with open(config_file, "w") as f:
            yaml.dump(config_dict, f, default_flow_style=False)

        # ── Run via Snakemake Python API ─────────────────────────────────────
        # Using the API runs Snakemake in-process rather than via subprocess.
        # This bypasses the kernel's execve() argument-length limit (~128 KB
        # per argument on Linux), which is hit when many subjects produce a
        # very long --apptainer-args --bind string.  The bind string is passed
        # as a Python object in memory; it never touches the command line.
        wf_file = self.workflow_dir / "process_all_workflow.smk"

        # Isolate the source cache per output directory to prevent NFS races
        # when multiple pipeline instances run concurrently on a shared filesystem.
        os.environ["SNAKEMAKE_SOURCECACHE_PATH"] = str(
            base_out / ".snakemake" / "source-cache"
        )

        resource_settings = ResourceSettings(
            cores=cores,
            # sys.maxsize mirrors what the CLI does for --jobs unlimited
            # (see snakemake/cli.py parse_jobs: "unlimited" → sys.maxsize).
            nodes=sys.maxsize,
        )

        deployment_settings = DeploymentSettings(
            deployment_method=frozenset({DeploymentMethod.APPTAINER}),
            apptainer_args=singularity_bind,
        )

        if scheduler == "sge":
            execution_settings = ExecutionSettings(
                # Cold NFS mounts on some compute nodes exceed the default 5 s.
                latency_wait=60,
                keep_metadata=False,
                standalone=True,
                # A batch run covers hundreds/thousands of subjects; one
                # subject's genuine failure (bad raw data, an edge case)
                # must not cancel everyone else's legitimate in-flight work.
                # Confirmed necessary in practice: a single false-positive
                # SGE status report once cancelled ~680 remaining steps of an
                # EPAD run at 33% (fixed separately in
                # snakemake-executor-plugin-sge), and a genuine per-subject
                # failure would hit the same default cancel-everything
                # behavior without this.
                keep_going=True,
            )
            # immediate_submit lives in RemoteExecutionSettings (Snakemake 9).
            remote_execution_settings = RemoteExecutionSettings(
                immediate_submit=True,
            )
            # max_jobs_per_second replaces the old --max-jobs-per-timespan 75000/1s.
            scheduling_settings = SchedulingSettings(
                max_jobs_per_second=75_000,
            )
            # notemp is required alongside immediate_submit - it lives in
            # StorageSettings in Snakemake 9 (confirmed from settings.types source).
            #
            # SharedFSUsage.SOURCE_CACHE is deliberately EXCLUDED here. When it is
            # present (the default), Snakemake forwards a single
            # --runtime-source-cache-path pointing at *this* process's own
            # ephemeral SourceCache.runtime_cache directory to every spawned
            # per-job invocation (snakemake/spawn_jobs.py: general_args()), so
            # that jobs sharing the same filesystem can reuse already-cached
            # Snakefiles instead of re-caching them individually. That
            # optimisation assumes the submitting process stays alive for the
            # whole run. With immediate_submit=True it does not: this process
            # submits the entire SGE array job graph and exits immediately,
            # which garbage-collects its runtime_cache TemporaryDirectory and
            # deletes it. Any array task that hasn't started yet (routine on a
            # busy queue) then tries to _do_cache() a Snakefile into that
            # already-deleted shared directory and crashes with a
            # FileNotFoundError from sourcecache.py's os.replace(), well before
            # it even reaches its assigned rule. Excluding SOURCE_CACHE makes
            # each spawned job fall back to its own private, process-scoped
            # runtime cache (workflow.py: source_cache_path property routes to
            # self.snakemake_tmp_dir instead of the shared one), which
            # eliminates the shared directory - and the race - entirely.
            storage_settings = StorageSettings(
                notemp=True,
                shared_fs_usage=frozenset(
                    u for u in SharedFSUsage.all() if u != SharedFSUsage.SOURCE_CACHE
                ),
            )
            executor = "sge"
        else:
            execution_settings = ExecutionSettings(
                keep_metadata=False,
                standalone=True,
            )
            remote_execution_settings = RemoteExecutionSettings()
            scheduling_settings = SchedulingSettings()
            storage_settings = StorageSettings()
            executor = "local"

        try:
            # Redirect stdout and stderr at the OS file-descriptor level so
            # nothing from Snakemake or its plugins leaks to the terminal.
            # Python-level redirects (sys.stdout/stderr) are insufficient because
            # logging.StreamHandler captures the stream object at init time, not
            # the sys.stderr attribute, so it bypasses a simple attribute swap.
            # os.dup2() redirects the actual fd 1/2, catching all writes
            # regardless of whether they come from Python logging, C extensions,
            # or inherited fds in subprocesses. The Snakemake file log writes to
            # a separate named fd (not fd 1 or 2) so it is unaffected.
            output_settings = OutputSettings()
            _devnull_fd = os.open(os.devnull, os.O_WRONLY)
            _saved_stdout_fd = os.dup(1)
            _saved_stderr_fd = os.dup(2)
            os.dup2(_devnull_fd, 1)
            os.dup2(_devnull_fd, 2)
            os.close(_devnull_fd)
            try:
                with SnakemakeApi(output_settings) as snk:
                    workflow_api = snk.workflow(
                        snakefile=wf_file,
                        workdir=base_out,
                        config_settings=ConfigSettings(
                            configfiles=[Path(config_file)],
                        ),
                        resource_settings=resource_settings,
                        storage_settings=storage_settings,
                        deployment_settings=deployment_settings,
                    )
                    dag_api = workflow_api.dag(
                        dag_settings=DAGSettings(
                            targets=["all"],
                            forcerun=frozenset(force_rules or []),
                        ),
                    )
                    success = dag_api.execute_workflow(
                        executor=executor,
                        execution_settings=execution_settings,
                        scheduling_settings=scheduling_settings,
                        remote_execution_settings=remote_execution_settings,
                        group_settings=GroupSettings(),
                        executor_settings=SgeExecutorSettings() if executor == "sge" else None,
                    )
            finally:
                os.dup2(_saved_stdout_fd, 1)
                os.dup2(_saved_stderr_fd, 2)
                os.close(_saved_stdout_fd)
                os.close(_saved_stderr_fd)

            # With --immediate-submit (SGE), the aggregation rule `all` is a
            # local rule with no shell command and is always skipped, causing
            # execute_workflow() to return False even though all compute jobs
            # were submitted and completed successfully.  The CLI exits 0 in
            # the same situation, so we treat this as success for SGE mode.
            if not success and executor != "sge":
                raise RuntimeError("Snakemake workflow returned failure")
            return "0", str(base_out)
        except Exception as e:
            raise RuntimeError(f"process-all failed: {e}")


def apply_process_all(subject_input: Optional[str] = None,
                      t1_pattern: Optional[str] = None,
                      flair_pattern: Optional[str] = None,
                      dwi_pattern: Optional[str] = None,
                      bvecs_pattern: Optional[str] = None,
                      bvals_pattern: Optional[str] = None,
                      mask_pattern: Optional[str] = None,
                      output_dir: Optional[str] = None,
                      scheduler: str = "local",
                      cores: int = 1,
                      skip_zscore: bool = False,
                      healthy_list: Optional[str] = None,
                      demographics: Optional[str] = None,
                      covariates: Optional[str] = None,
                      poly_terms: Optional[str] = None,
                      parcellation: str = "freesurfer",
                      force_rules: Optional[List[str]] = None,
                      verbose: bool = False,
                      config_yaml: Optional[str] = None) -> dict:
    """Run the full process-all pipeline and return a summary dict.

    CLI / caller arguments take priority over values in ``config_yaml`` when
    both are provided.
    """
    yaml_cfg = load_yaml_config(config_yaml)
    subject_input = first_truthy(subject_input, yaml_cfg.get("subject_input"), yaml_cfg.get("subject"))
    t1_pattern    = first_truthy(t1_pattern,    yaml_cfg.get("t1_pattern"), yaml_cfg.get("t1"))
    flair_pattern = first_truthy(flair_pattern, yaml_cfg.get("flair_pattern"), yaml_cfg.get("flair"))
    dwi_pattern   = first_truthy(dwi_pattern,   yaml_cfg.get("dwi_pattern"), yaml_cfg.get("dwi"))
    bvecs_pattern = first_truthy(bvecs_pattern, yaml_cfg.get("bvecs_pattern"), yaml_cfg.get("bvecs"))
    bvals_pattern = first_truthy(bvals_pattern, yaml_cfg.get("bvals_pattern"), yaml_cfg.get("bvals"))
    mask_pattern  = first_truthy(mask_pattern,  yaml_cfg.get("mask_pattern"), yaml_cfg.get("mask"))
    output_dir    = first_truthy(output_dir,    yaml_cfg.get("output_dir"))
    scheduler     = first_truthy(scheduler,     yaml_cfg.get("scheduler")) or "local"
    cores_raw     = first_truthy(cores,         yaml_cfg.get("cores"))
    cores         = int(cores_raw) if cores_raw is not None else 1
    skip_zscore   = bool(yaml_cfg.get("skip_zscore", False)) or bool(skip_zscore)
    healthy_list  = first_truthy(healthy_list,  yaml_cfg.get("healthy_list"), yaml_cfg.get("healthy_subjects_list"))
    demographics  = first_truthy(demographics,  yaml_cfg.get("demographics"), yaml_cfg.get("demographics_csv"))
    covariates    = first_truthy(covariates,    yaml_cfg.get("covariates"))
    poly_terms    = first_truthy(poly_terms,    yaml_cfg.get("poly_terms"), yaml_cfg.get("polynomial_terms"))
    parcellation  = first_truthy(parcellation,  yaml_cfg.get("parcellation")) or "freesurfer"

    if not subject_input:
        raise ValueError("subject_input is required (via argument or config_yaml)")
    if not t1_pattern:
        raise ValueError("t1_pattern is required (via argument or config_yaml)")
    if not flair_pattern:
        raise ValueError("flair_pattern is required (via argument or config_yaml)")
    if not dwi_pattern:
        raise ValueError("dwi_pattern is required (via argument or config_yaml)")
    if not output_dir:
        raise ValueError("output_dir is required (via argument or config_yaml)")

    if verbose:
        logging.basicConfig(level=logging.INFO)
    else:
        logging.basicConfig(level=logging.WARNING, format="%(levelname)s: %(message)s")

    proc = ProcessAllProcessor()
    try:
        _, results_dir = proc.run_process_all(
            subject_input=subject_input,
            t1_pattern=t1_pattern,
            flair_pattern=flair_pattern,
            dwi_pattern=dwi_pattern,
            bvecs_pattern=bvecs_pattern,
            bvals_pattern=bvals_pattern,
            mask_pattern=mask_pattern,
            output_dir=output_dir,
            scheduler=scheduler,
            cores=cores,
            skip_zscore=skip_zscore,
            healthy_list=healthy_list,
            demographics=demographics,
            covariates=covariates,
            poly_terms=poly_terms,
            parcellation=parcellation,
            force_rules=force_rules,
            verbose=verbose,
        )
        print("✅ Full processing job submitted successfully", flush=True)
        return {"success": True, "results_dir": results_dir}
    except (EnvironmentError, FileNotFoundError) as e:
        print(f"❌ Full processing job failed:\n{e}", file=sys.stderr, flush=True)
        return {"success": False, "error": str(e)}
    except Exception as e:
        print(f"❌ Full processing job failed:\nUnexpected error: {e}", file=sys.stderr, flush=True)
        return {"success": False, "error": f"Unexpected error: {e}"}
