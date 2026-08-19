"""
Unified workflow to run FreeSurfer, GIF, BaMoS, Tractography, and Metrics for one or more subjects.

All modules share a single base output_dir. Each module writes to:
  {output_dir}/{subject_id}/{tool_name}/outputs/
  {output_dir}/{subject_id}/{tool_name}/intermediate/   (where applicable)
  {output_dir}/{subject_id}/{tool_name}/logs/

Supports multiple subjects in a single Snakemake invocation to avoid repeated DAG-build and
container-locate overhead. The processor passes all subjects as lists; this workflow fans them
out to the appropriate sub-workflow modules, which already handle per-subject wildcards.
"""

import os
import re
import yaml
from pathlib import Path

# ── Subject and per-subject input lists ─────────────────────────────────────
subjects = config.get("subjects", [config.get("subject_id")])
if isinstance(subjects, str):
    subjects = [subjects]
subjects = [s for s in subjects if s]  # drop any empty entries

# target_subjects: all subjects minus the healthy cohort.
# The processor derives this and passes it explicitly; fall back to all subjects
# when running standalone (no healthy list provided).
_healthy_set = set()
_healthy_list_cfg = config.get("healthy_subjects_list", "")
if _healthy_list_cfg:
    import pathlib as _pathlib
    _hfile = _pathlib.Path(_healthy_list_cfg)
    if _hfile.exists():
        _healthy_set = {line.strip() for line in _hfile.read_text().splitlines() if line.strip()}
target_subjects = [s for s in subjects if s not in _healthy_set]

t1_files = config.get("t1_files", [config.get("t1_file")])
if isinstance(t1_files, str):
    t1_files = [t1_files]

t1_files_singularity = config.get("t1_files_singularity", [config.get("t1_file_singularity")])
if isinstance(t1_files_singularity, str):
    t1_files_singularity = [t1_files_singularity]

flair_files = config.get("flair_files", [config.get("flair_file")])
if isinstance(flair_files, str):
    flair_files = [flair_files]

flair_files_singularity = config.get("flair_files_singularity", [config.get("flair_file_singularity", "")])
if isinstance(flair_files_singularity, str):
    flair_files_singularity = [flair_files_singularity]

dwi_files = config.get("dwi_files", [config.get("dwi_file")])
if isinstance(dwi_files, str):
    dwi_files = [dwi_files]

dwi_paths_singularity = config.get("dwi_paths_singularity", [config.get("dwi_path_singularity")])
if isinstance(dwi_paths_singularity, str):
    dwi_paths_singularity = [dwi_paths_singularity]

bvecs = config.get("bvecs_list", [config.get("bvecs", "")])
if isinstance(bvecs, str):
    bvecs = [bvecs]

bvecs_singularity = config.get("bvecs_paths_singularity", [config.get("bvecs_path_singularity", "")])
if isinstance(bvecs_singularity, str):
    bvecs_singularity = [bvecs_singularity]

bvals = config.get("bvals_list", [config.get("bvals", "")])
if isinstance(bvals, str):
    bvals = [bvals]

bvals_singularity = config.get("bvals_paths_singularity", [config.get("bvals_path_singularity", "")])
if isinstance(bvals_singularity, str):
    bvals_singularity = [bvals_singularity]

mask_files = config.get("mask_files", [config.get("mask_file", "")])
if isinstance(mask_files, str):
    mask_files = [mask_files]

mask_files_singularity = config.get("mask_files_singularity", [config.get("mask_file_singularity", "")])
if isinstance(mask_files_singularity, str):
    mask_files_singularity = [mask_files_singularity]

# Pad shorter lists to match number of subjects
_n = len(subjects)
def _pad(lst, n, default=""):
    return (lst + [default] * n)[:n]

t1_files              = _pad(t1_files, _n)
t1_files_singularity  = _pad(t1_files_singularity, _n)
flair_files           = _pad(flair_files, _n)
flair_files_singularity = _pad(flair_files_singularity, _n)
dwi_files             = _pad(dwi_files, _n)
dwi_paths_singularity = _pad(dwi_paths_singularity, _n)
bvecs                 = _pad(bvecs, _n)
bvecs_singularity     = _pad(bvecs_singularity, _n)
bvals                 = _pad(bvals, _n)
bvals_singularity     = _pad(bvals_singularity, _n)
mask_files            = _pad(mask_files, _n)
mask_files_singularity = _pad(mask_files_singularity, _n)

# ── Per-subject lookup dicts ─────────────────────────────────────────────────
t1_map          = dict(zip(subjects, t1_files))
t1_sing_map     = dict(zip(subjects, t1_files_singularity))
t1_id_map       = {s: t1_map[s].split("/")[-1].split(".", 1)[0] for s in subjects}
flair_map       = dict(zip(subjects, flair_files))
flair_sing_map  = dict(zip(subjects, flair_files_singularity))
dwi_map         = dict(zip(subjects, dwi_files))
dwi_sing_map    = dict(zip(subjects, dwi_paths_singularity))
bvecs_map       = dict(zip(subjects, bvecs))
bvecs_sing_map  = dict(zip(subjects, bvecs_singularity))
bvals_map       = dict(zip(subjects, bvals))
bvals_sing_map  = dict(zip(subjects, bvals_singularity))
mask_map        = dict(zip(subjects, mask_files))
mask_sing_map   = dict(zip(subjects, mask_files_singularity))

# ── Global settings ──────────────────────────────────────────────────────────
output_dir            = config.get("output_dir")
shared_conda_env_name = config.get("conda_env_name", "process_all_env")
keep_intermediate     = config.get("keep_intermediate", False)
leukoquant_parent_dir = config.get("leukoquant_parent_dir")
verbose               = config.get("verbose", False)
# PARCELLATIONS: list of parcellations to fan out across.
# Accepts a list (new) or a single string (legacy) for backwards compat.
PARCELLATIONS = config.get("parcellations", [config.get("parcellation", "freesurfer")])
if isinstance(PARCELLATIONS, str):
    PARCELLATIONS = [p.strip() for p in PARCELLATIONS.split(",") if p.strip()]

# DTI and NODDI are parcellation-independent: they read preprocessed DWI from
# TRACULA's dmri/ folder. Prefer freesurfer if available (most common); otherwise
# use the first parcellation in the list.
_primary_parc = "freesurfer" if "freesurfer" in PARCELLATIONS else PARCELLATIONS[0]

# ── Concrete cross-module dependency paths (per-subject helpers) ─────────────
def gif_parcellation(subject):
    return f"{output_dir}/{subject}/gif/outputs/{t1_id_map[subject]}_NeuroMorph_Parcellation.nii.gz"

def gif_brain_mask(subject):
    return f"{output_dir}/{subject}/gif/outputs/{t1_id_map[subject]}_NeuroMorph_Brain.nii.gz"

def recon_aparc_aseg(subject):
    return f"{output_dir}/{subject}/recon-all/outputs/mri/aparc+aseg.mgz"

def recon_thalamic(subject):
    return f"{output_dir}/{subject}/recon-all/outputs/mri/ThalamicNuclei.v13.T1.FSvoxelSpace.mgz"

def tracula_brain_mask(subject, parcellation=None):
    parc = parcellation if parcellation is not None else _primary_parc
    return f"{output_dir}/{subject}/tracula-{parc}/outputs/dmri/nodif_brain_mask.nii.gz"

def tracula_tract_merged(subject, parcellation=None):
    parc = parcellation if parcellation is not None else _primary_parc
    return f"{output_dir}/{subject}/tracula-{parc}/outputs/dpath/merged_avg16_syn_bbr.mgz"

def tracula_dwi(subject):
    # DTI/NODDI always read from the primary parcellation's TRACULA run.
    return f"{output_dir}/{subject}/tracula-{_primary_parc}/outputs/dmri/dwi.nii.gz"

def tracula_bvecs(subject):
    return f"{output_dir}/{subject}/tracula-{_primary_parc}/outputs/dmri/dwi.bvecs"

def tracula_bvals(subject):
    return f"{output_dir}/{subject}/tracula-{_primary_parc}/outputs/dmri/dwi.bvals"

def dti_fa(subject):
    return f"{output_dir}/{subject}/dti/outputs/fa.nii.gz"

def dti_md(subject):
    return f"{output_dir}/{subject}/dti/outputs/md.nii.gz"

def dti_ad(subject):
    return f"{output_dir}/{subject}/dti/outputs/ad.nii.gz"

def dti_rd(subject):
    return f"{output_dir}/{subject}/dti/outputs/rd.nii.gz"

def noddi_odi(subject):
    return f"{output_dir}/{subject}/noddi/outputs/odi.nii.gz"

def bamos_correction(subject):
    return f"{output_dir}/{subject}/bamos/outputs/CorrectLesion.nii.gz"

def z_score_outputs(subject):
    """Return z-score output paths for a target subject.

    Returns an empty list when z-scores are disabled for this run (no healthy
    cohort supplied or skip_zscore=True), and for healthy subjects themselves.
    The `_include_z_score` check must mirror the gate further below that loads
    the z_score module - otherwise we'd demand files no rule produces.

    Declared unconditionally (not existence-checked) -- dti/noddi are part of
    this same DAG for target subjects, so Snakemake's own dependency
    resolution (hold_jid chaining via the sge executor plugin) waits for them
    to finish before z-score runs, and z-score before metrics-gif/fs (which
    takes this function's return value as its own required input via
    z_score_map below) -- exactly the intended dti/noddi -> z-score -> metrics
    ordering, no Python pre-check needed.
    """
    if not _include_z_score:
        return []
    if subject not in target_subjects:
        return []
    if _zscore_demo_complete_subjects is not None and subject not in _zscore_demo_complete_subjects:
        return []
    # Build suffix matching what z_score_workflow.smk produces
    covs  = config.get("covariates", [])
    polys = config.get("polynomial_terms", [])
    suffix = ""
    if covs:
        suffix += "_cov_" + "_".join(covs)
    if polys:
        suffix += "_poly_" + "_".join(polys)
    z_metrics = config.get("z_score_metrics", {
        "dti_fa": None, "dti_md": None, "dti_ad": None, "dti_rd": None, "noddi_odi": None,
    })
    # "z-score" (hyphenated, singular) matches TOOL_NAME in
    # z_score_workflow.smk -- this previously said "z_scores" (underscore,
    # plural), which never matched any real output path, so this function's
    # returned paths never lined up with a file any rule actually produced.
    return [
        f"{output_dir}/{subject}/z-score/outputs/{m}_z_score{suffix}.nii.gz"
        for m in z_metrics
    ]

# ── Ensure key directories exist before SGE submission ──────────────────────
for _s in subjects:
    os.makedirs(f"{output_dir}/{_s}/recon-all/logs",      exist_ok=True)
    os.makedirs(f"{output_dir}/{_s}/gif/logs",            exist_ok=True)
    os.makedirs(f"{output_dir}/{_s}/bamos/outputs",       exist_ok=True)
    os.makedirs(f"{output_dir}/{_s}/bamos/intermediate",  exist_ok=True)
    os.makedirs(f"{output_dir}/{_s}/bamos/logs",          exist_ok=True)
    os.makedirs(f"{output_dir}/{_s}/dti/logs",            exist_ok=True)
    os.makedirs(f"{output_dir}/{_s}/noddi/logs",          exist_ok=True)
    for _parc in PARCELLATIONS:
        os.makedirs(f"{output_dir}/{_s}/tracula-{_parc}/logs",        exist_ok=True)
        os.makedirs(f"{output_dir}/{_s}/tract_qc-{_parc}/logs",       exist_ok=True)
        os.makedirs(f"{output_dir}/{_s}/metrics-{_parc}/logs",        exist_ok=True)
        os.makedirs(f"{output_dir}/{_s}/metrics-{_parc}/outputs",     exist_ok=True)

# ── Sub-module configs (lists of all subjects) ───────────────────────────────

config_recon = {
    "subject_id": subjects,
    "t1_file": t1_files,
    "t1_file_singularity": t1_files_singularity,
    "output_dir": output_dir,
    "output_dir_singularity": "/output",
    "container_name": "freesurfer_unified_container",
    "leukoquant_parent_dir": leukoquant_parent_dir,
    "conda_env_name": shared_conda_env_name,
    "containers": config.get("containers", {}),
    "keep_intermediate": keep_intermediate,
}

config_gif = {
    "subjects": subjects,
    "t1_files": t1_files,
    "image_ids": [t1_id_map[s] for s in subjects],
    "t1_files_singularity": t1_files_singularity,
    "mask_files": mask_files,
    "mask_files_singularity": mask_files_singularity,
    "output_dir": output_dir,
    "output_dir_singularity": "/output",
    "leukoquant_parent_dir": leukoquant_parent_dir,
    "conda_env_name": shared_conda_env_name,
    "containers": config.get("containers", {}),
    "keep_intermediate": keep_intermediate,
}

config_bamos = {
    "subjects": subjects,
    "flair_files": flair_files,
    # Use the FLAIR singularity paths produced by the processor (each FLAIR is
    # bound at /input_flair_{i}/<basename>). The previous synthesis of
    # /input/<basename> referenced a mount point that was never bound, causing
    # niftyreg inside BaMoS to fail with "failed to find header file".
    "flair_files_singularity": [flair_sing_map[s] for s in subjects],
    "t1_files": t1_files,
    "t1_files_singularity": t1_files_singularity,
    # GIF outputs/ dir for each subject - bamos reads parcellation from here
    "gif_results_paths": [f"{output_dir}/{s}/gif/outputs" for s in subjects],
    "gif_results_paths_singularity": [f"/output/{s}/gif/outputs" for s in subjects],
    "output_dir": output_dir,
    "output_dir_singularity": "/output",
    "space": 1,
    "jump_start": 0,
    "opt": "TA",
    "leukoquant_parent_dir": leukoquant_parent_dir,
    "conda_env_name": shared_conda_env_name,
    "containers": config.get("containers", {}),
    "keep_intermediate": keep_intermediate,
}

# The default GIF → FreeSurfer label mapping ships with the repo.
_gif_to_fs_mapping = f"{leukoquant_parent_dir}/mappings/gif_to_freesurfer.csv"

config_tracula = {
    "subjects": subjects,
    "dwi_imgs": dwi_files,
    "dwi_paths_singularity": dwi_paths_singularity,
    # recon-all outputs/ dir - tracula copies from here into scratch
    "freesurfer_recon_root": f"{output_dir}",
    "freesurfer_recon_root_singularity": "/output",
    "output_dir": output_dir,
    "output_dir_singularity": "/output",
    "bvecs_paths": bvecs,
    "bvecs_paths_singularity": bvecs_singularity,
    "bvals_paths": bvals,
    "bvals_paths_singularity": bvals_singularity,
    # Full list so tracula_workflow.smk creates one rule per parcellation.
    "parcellations": PARCELLATIONS,
    # Per-subject recon-all output directories; the per-parcellation trac_all rules
    # read aparc+aseg.mgz and ThalamicNuclei.mgz via _recon_aparc_aseg / _recon_thalamic_nuclei.
    # process_all always runs recon for all subjects, so needs_recon is empty
    # (disables tracula_workflow's embedded recon sub-workflow).
    "recon_dirs":             {s: f"{output_dir}/{s}/recon-all/outputs" for s in subjects},
    "recon_dirs_singularity": {s: f"/output/{s}/recon-all/outputs" for s in subjects},
    "needs_recon":            [],
    # Per-subject parcellation maps for gif - each subject's GIF output is already
    # in output_dir (bound as /output), so no extra bind entry is needed.
    "brain_parcellation_map": (
        {s: gif_parcellation(s) for s in subjects} if "gif" in PARCELLATIONS else {}
    ),
    "brain_parcellation_sing_map": (
        {
            s: f"/output/{s}/gif/outputs/{t1_id_map[s]}_NeuroMorph_Parcellation.nii.gz"
            for s in subjects
        }
        if "gif" in PARCELLATIONS else {}
    ),
    # GIF brain mask maps - passed to tracula workflow for custom parcellations
    "brain_mask_map": (
        {s: gif_brain_mask(s) for s in subjects} if "gif" in PARCELLATIONS else {}
    ),
    "brain_mask_sing_map": (
        {
            s: f"/output/{s}/gif/outputs/{t1_id_map[s]}_NeuroMorph_Brain.nii.gz"
            for s in subjects
        }
        if "gif" in PARCELLATIONS else {}
    ),
    "mapping_file": _gif_to_fs_mapping if "gif" in PARCELLATIONS else "",
    # leukoquant_parent_dir is bound as /leukoquant, so the mapping CSV is reachable
    # at a fixed container path without an extra bind entry.
    "mapping_file_singularity": (
        "/leukoquant/mappings/gif_to_freesurfer.csv" if "gif" in PARCELLATIONS else ""
    ),
    "container_name": "freesurfer_unified_container",
    "leukoquant_parent_dir": leukoquant_parent_dir,
    "conda_env_name": shared_conda_env_name,
    "containers": config.get("containers", {}),
    "keep_intermediate": keep_intermediate,
}

config_dti = {
    "subjects": subjects,
    "dwi_paths": [tracula_dwi(s) for s in subjects],
    "dwi_paths_singularity": [f"/output/{s}/tracula-{_primary_parc}/outputs/dmri/dwi.nii.gz" for s in subjects],
    "brain_mask_paths": [tracula_brain_mask(s) for s in subjects],
    "brain_mask_paths_singularity": [f"/output/{s}/tracula-{_primary_parc}/outputs/dmri/nodif_brain_mask.nii.gz" for s in subjects],
    "bvecs_paths": [tracula_bvecs(s) for s in subjects],
    "bvecs_paths_singularity": [f"/output/{s}/tracula-{_primary_parc}/outputs/dmri/dwi.bvecs" for s in subjects],
    "bvals_paths": [tracula_bvals(s) for s in subjects],
    "bvals_paths_singularity": [f"/output/{s}/tracula-{_primary_parc}/outputs/dmri/dwi.bvals" for s in subjects],
    "output_dir": output_dir,
    "output_dir_singularity": "/output",
    "leukoquant_dir": leukoquant_parent_dir + "/leukoquant",
    "conda_env_name": shared_conda_env_name,
    "containers": config.get("containers", {}),
    "keep_intermediate": keep_intermediate,
}

config_noddi = {
    "subjects": subjects,
    "dwi_paths": [tracula_dwi(s) for s in subjects],
    "dwi_paths_singularity": [f"/output/{s}/tracula-{_primary_parc}/outputs/dmri/dwi.nii.gz" for s in subjects],
    "brain_mask_paths": [tracula_brain_mask(s) for s in subjects],
    "brain_mask_paths_singularity": [f"/output/{s}/tracula-{_primary_parc}/outputs/dmri/nodif_brain_mask.nii.gz" for s in subjects],
    "bvecs_paths": [tracula_bvecs(s) for s in subjects],
    "bvecs_paths_singularity": [f"/output/{s}/tracula-{_primary_parc}/outputs/dmri/dwi.bvecs" for s in subjects],
    "bvals_paths": [tracula_bvals(s) for s in subjects],
    "bvals_paths_singularity": [f"/output/{s}/tracula-{_primary_parc}/outputs/dmri/dwi.bvals" for s in subjects],
    "output_dir": output_dir,
    "output_dir_singularity": "/output",
    "leukoquant_dir": leukoquant_parent_dir + "/leukoquant",
    "conda_env_name": shared_conda_env_name,
    "containers": config.get("containers", {}),
    "keep_intermediate": keep_intermediate,
}

config_tract_qc = {
    "leukoquant_parent_dir": leukoquant_parent_dir,
    "subjects": subjects,
    "parcellations": PARCELLATIONS,
    # {subject} is resolved per-wildcard; {parcellation} is substituted at params-expansion
    # time by _resolve_tractography_path(subject, p) inside each per-parcellation rule.
    "tractography_path": f"{output_dir}/{{subject}}/tracula-{{parcellation}}/outputs:dpath/*/path.pd.trk:dwi",
    "output_dir": output_dir,
    "validate_tractography_glob": False,
    # Enable parcellation-suffixed output folders (tract_qc-{parcellation}/).
    # Standalone process-tract-qc leaves this unset (defaults False → tract_qc/).
    "use_parcellation_suffix": True,
    "conda_env_name": shared_conda_env_name,
    "containers": config.get("containers", {}),
    "singularity_binds": {
        f"{output_dir}": "/output",
        f"{leukoquant_parent_dir}": "/leukoquant",
    },
    # Per-subject per-parcellation tracula tract_merged paths; the for-loop tract_qc
    # rules declare a DAG edge to tracula via this map (no use-rule override needed).
    "tracula_tract_merged_map": {
        s: {p: tracula_tract_merged(s, p) for p in PARCELLATIONS}
        for s in subjects
    },
}

# Build per-subject T1/DWI directory bind map for metrics container.
# metrics_calc.sh needs to read T1 and DWI files from the data directories,
# which are outside output_dir and leukoquant_parent_dir.  We mount each
# unique host directory to the same container path the processor already
# allocated (e.g. /input_t1_0/), deriving it from t1_files_singularity.
_metrics_binds: dict = {
    f"{output_dir}": "/output",
    f"{leukoquant_parent_dir}/leukoquant": "/leukoquant",
    f"{leukoquant_parent_dir}": "/leukoquant",
}
for _s, _t1_host, _t1_sif in zip(subjects, t1_files, t1_files_singularity):
    if _t1_host and _t1_sif:
        _metrics_binds[str(Path(_t1_host).parent)] = str(Path(_t1_sif).parent)
for _s, _dwi_host_list, _dwi_sif_list in zip(subjects, dwi_files, dwi_paths_singularity):
    # dwi_files is now a nested list (one inner list per subject).
    for _dwi_host, _dwi_sif in zip(_dwi_host_list or [], _dwi_sif_list or []):
        if _dwi_host and _dwi_sif:
            _metrics_binds[str(Path(_dwi_host).parent)] = str(Path(_dwi_sif).parent)

# Must be defined before config_metrics because z_score_outputs() references it,
# and config_metrics eagerly evaluates z_score_map via a dict comprehension.
_include_z_score = not config.get("skip_zscore", False) and config.get("healthy_subjects_list")

# ── One-time demographics-completeness check for z-score covariates ─────────
# z_score_utils.py's generate_design_matrix() already degrades gracefully
# per-subject (falls back to a simple, non-covariate-adjusted z-score) when a
# subject's row is missing or has a null value for a configured covariate --
# see the isna() checks around z_score_utils.py:285-291. That fallback is
# silent, though: a subject would get a different z-score methodology than
# the rest of the cohort with nothing in the DAG/job stats to show it. This
# makes that exclusion explicit and visible instead, and only applies when
# covariates/poly-terms are actually configured -- with none configured,
# z_score_calc.sh never builds a design matrix at all (see z_score_calc.sh's
# own "Demographics CSV not provided" skip), so there's nothing to check.
_zscore_covs  = config.get("covariates", [])
_zscore_polys = config.get("polynomial_terms", [])
_zscore_demo_csv = config.get("demographics_csv", "")
_zscore_required_cols = {c.strip().lower() for c in _zscore_covs}
for _p in _zscore_polys:
    # Poly terms arrive pre-formatted as "covariate_power" (e.g. "age_2") --
    # same rsplit("_", 1) parsing z_score_utils.py itself uses to recover the
    # base covariate name.
    _parts = _p.rsplit("_", 1)
    if len(_parts) == 2:
        _zscore_required_cols.add(_parts[0].strip().lower())

# None = check not applicable (no covariates/poly-terms configured for this
# run) -- z_score_outputs() below treats None as "don't filter by this".
_zscore_demo_complete_subjects = None
if _include_z_score and _zscore_required_cols and _zscore_demo_csv:
    import pandas as _pd

    _demo_df = _pd.read_csv(_zscore_demo_csv)
    _demo_df.columns = _demo_df.columns.str.strip().str.lower()
    _demo_id_col = "subject" if "subject" in _demo_df.columns else _demo_df.columns[0]
    _demo_df[_demo_id_col] = _demo_df[_demo_id_col].astype(str)
    _demo_required_cols_present = [c for c in _zscore_required_cols if c in _demo_df.columns]

    def _subject_has_complete_covariates(sid):
        _row = _demo_df[_demo_df[_demo_id_col] == sid]
        if _row.empty:
            return False
        if len(_demo_required_cols_present) < len(_zscore_required_cols):
            return False  # a required covariate column doesn't exist at all
        return not _row[_demo_required_cols_present].isna().any(axis=None)

    # The healthy cohort feeds every target subject's design matrix -- one
    # healthy subject with a missing covariate value would silently degrade
    # ALL target subjects to simple z-score (not just that one healthy
    # subject), so this fails loudly up front instead of letting that happen
    # invisibly across the whole run.
    _zscore_bad_healthy = [s for s in _healthy_set if not _subject_has_complete_covariates(s)]
    if _zscore_bad_healthy:
        raise ValueError(
            f"Healthy reference subjects missing required covariate values "
            f"{sorted(_zscore_required_cols)}: {sorted(_zscore_bad_healthy)}. "
            f"Fix the demographics CSV or remove them from --healthy-subjects "
            f"before running with --covariates/--poly-terms."
        )

    _zscore_demo_complete_subjects = {
        s for s in target_subjects if _subject_has_complete_covariates(s)
    }

config_metrics = {
    "leukoquant_parent_dir": leukoquant_parent_dir,
    "subjects": subjects,
    "parcellations": PARCELLATIONS,
    # Enable parcellation-suffixed output folders (metrics-{parcellation}/).
    # Standalone process-metrics leaves this unset (defaults False → metrics/).
    "use_parcellation_suffix": True,
    # {subject} is a Snakemake wildcard; {parcellation} is substituted at params-expansion
    # time by _sif_tractography_path(subject, p) inside each per-parcellation rule.
    "tractography_path": f"{output_dir}/{{subject}}/tracula-{{parcellation}}/outputs:dpath/*/path.pd.nii.gz:dwi",
    "lesion_path": f"wmh={output_dir}/{{subject}}/bamos/outputs:CorrectLesion.nii.gz:t1",
    # t1/dwi use singularity-translated paths because translate_path() needs the
    # bind map to convert host paths to container paths, and t1_map holds host paths.
    # Passing the already-translated sif paths avoids a second translation step and
    # ensures metrics_calc.sh receives paths valid inside the container.
    # Use TRACULA-preprocessed DWI (from _primary_parc) for consistency with DTI/NODDI fitting.
    "t1_path": None,
    "dwi_path": None,
    "t1_map":  {s: t1_sing_map[s]  for s in subjects},
    "dwi_map": {s: f"/output/{s}/tracula-{_primary_parc}/outputs/dmri/dwi.nii.gz" for s in subjects},
    "metrics": {
        "dti_fa":    f"{output_dir}/{{subject}}/dti/outputs:fa.nii.gz:dwi",
        "dti_md":    f"{output_dir}/{{subject}}/dti/outputs:md.nii.gz:dwi",
        "dti_ad":    f"{output_dir}/{{subject}}/dti/outputs:ad.nii.gz:dwi",
        "dti_rd":    f"{output_dir}/{{subject}}/dti/outputs:rd.nii.gz:dwi",
        "noddi_odi": f"{output_dir}/{{subject}}/noddi/outputs:odi.nii.gz:dwi",
        "noddi_ndi": f"{output_dir}/{{subject}}/noddi/outputs:ndi.nii.gz:dwi",
        "noddi_fwf": f"{output_dir}/{{subject}}/noddi/outputs:fwf.nii.gz:dwi",
    },
    "skip_subject_dir_metrics": False,
    "output_dir": output_dir,
    "verbose": verbose,
    "conda_env_name": shared_conda_env_name,
    "containers": config.get("containers", {}),
    "keep_intermediate": keep_intermediate,
    "singularity_binds": _metrics_binds,
    # Cross-module dependency maps: per-subject paths to upstream outputs.
    # metrics_workflow.smk's for-loop rules read these to declare Snakemake DAG
    # edges without requiring use-rule overrides (which cannot be used in a loop).
    "bamos_correction_map":  {s: bamos_correction(s) for s in subjects},
    "dti_fa_map":            {s: dti_fa(s) for s in subjects},
    "dti_md_map":            {s: dti_md(s) for s in subjects},
    "dti_ad_map":            {s: dti_ad(s) for s in subjects},
    "dti_rd_map":            {s: dti_rd(s) for s in subjects},
    "noddi_odi_map":         {s: noddi_odi(s) for s in subjects},
    "thalamic_nuclei_map":   {s: recon_thalamic(s) for s in subjects},
    "tract_qc_report_map":   {s: {p: f"{output_dir}/{s}/tract_qc-{p}/outputs/qc_report.csv" for p in PARCELLATIONS} for s in subjects},
    "z_score_map":           {s: z_score_outputs(s) for s in subjects},
}

# ── Z-score config (cross-subject, optional) ─────────────────────────────────

_t1_sing_map = dict(zip(subjects, t1_files_singularity))

# Build bind map for z-score: maps each T1's directory (host) → container mount point.
# The processor mounts each T1 as /input_t1_{i}/{filename}; we need the directory-level
# mapping so translate_path() can rewrite host paths to their container equivalents.
_zscore_binds = {
    output_dir: "/output",
    leukoquant_parent_dir: "/leukoquant",
    # demographics_csv / healthy_subjects_list container-mount mapping, computed
    # by ProcessAllProcessor.run_process_all() from the consolidated bind_entries
    # actually passed to --apptainer-args. Without this, translate_path() below
    # can't find a matching prefix and passes the raw host path straight into the
    # container, where it doesn't exist (FileNotFoundError from z_score_calc.sh).
    **config.get("singularity_binds", {}),
}
for _s, _t1_host, _t1_sif in zip(subjects, t1_files, t1_files_singularity):
    if _t1_host and _t1_sif:
        _host_dir = str(Path(_t1_host).parent)
        _sif_dir  = str(Path(_t1_sif).parent)
        _zscore_binds[_host_dir] = _sif_dir

config_z_score = {
    "leukoquant_parent_dir": leukoquant_parent_dir,
    # z_score_workflow.smk builds its OWN rule-all target list straight from
    # this list -- it has no knowledge of the demographics-completeness gate
    # above (that logic lives entirely in this file). Filtering here, at the
    # source, is what actually keeps demographics-incomplete subjects out of
    # the DAG; z_score_outputs()'s own check further down is a second,
    # independent safety net for the metrics-dependency wiring, not the
    # primary gate.
    "target_subjects_list": (
        [s for s in target_subjects if s in _zscore_demo_complete_subjects]
        if _zscore_demo_complete_subjects is not None
        else target_subjects
    ),
    "healthy_subjects_list": config.get("healthy_subjects_list", ""),
    "demographics_csv": config.get("demographics_csv", ""),
    "output_dir": output_dir,
    # t1_map covers all subjects (targets + healthy merged by the processor).
    # z_score_utils.get_t1_path() uses this when t1_pattern is absent.
    "t1_map": t1_map,
    # dwi_map: host-path equivalent of the TRACULA-preprocessed DWI used to fit
    # DTI/NODDI (see tracula_dwi() above). Required whenever metric_space="dwi"
    # (the default here) so healthy subjects' DWI b0s can be registered to their
    # own T1 before being chained into target-T1 space.
    "dwi_map": {s: tracula_dwi(s) for s in subjects},
    # noddi_ndi and noddi_fwf are only produced for multi-shell acquisitions.
    # Including them as Snakemake inputs would break single-shell runs at DAG time
    # because no rule produces them. z_score_calc.sh handles missing maps gracefully.
    "metrics": config.get("z_score_metrics", {
        "dti_fa":    f"{output_dir}/{{subject}}/dti/outputs:fa.nii.gz:dwi",
        "dti_md":    f"{output_dir}/{{subject}}/dti/outputs:md.nii.gz:dwi",
        "dti_ad":    f"{output_dir}/{{subject}}/dti/outputs:ad.nii.gz:dwi",
        "dti_rd":    f"{output_dir}/{{subject}}/dti/outputs:rd.nii.gz:dwi",
        "noddi_odi": f"{output_dir}/{{subject}}/noddi/outputs:odi.nii.gz:dwi",
    }),
    # DTI/NODDI are always computed in native diffusion space in this pipeline
    # (see the metrics dict above -- every entry is tagged ":dwi") -- "t1" is
    # NOT a safe default here, unlike the standalone process-zscore CLI where
    # metric_space is a free user choice. Defaulting to "t1" silently skipped
    # the DWI->T1 resample and fed native-space (e.g. 116x116x80) metrics
    # straight into a GLM built from T1-space (e.g. 176x240x256) registrations,
    # crashing at the final z-score step with a dimension mismatch.
    "metric_space": config.get("metric_space", "dwi"),
    "output_space": config.get("output_space", "t1"),
    "covariates": config.get("covariates", []),
    "polynomial_terms": config.get("polynomial_terms", []),
    "threads": config.get("z_score_threads", 4),
    "container_name": "freesurfer_unified_container",
    "conda_env_name": shared_conda_env_name,
    "containers": config.get("containers", {}),
    "singularity_binds": _zscore_binds,
}

# ── Module imports ────────────────────────────────────────────────────────────

# bamos_workflow imports gif_workflow internally, which causes gif_run_gif to
# appear under two namespaces (gif_run_gif and bamos_gif_run_gif) with identical
# outputs. Declare ruleorder so Snakemake resolves the ambiguity without error.
ruleorder: gif_run_gif > bamos_gif_run_gif

# FreeSurfer recon-all
module recon:
    snakefile: "recon_all_workflow.smk"
    config: config_recon

use rule * from recon as recon_*

# GIF
module gif:
    snakefile: "gif_workflow.smk"
    config: config_gif

use rule * from gif as gif_*

# BaMoS - override run_bamos to add per-subject GIF dependency
module bamos:
    snakefile: "bamos_workflow.smk"
    config: config_bamos

use rule * from bamos as bamos_*

use rule run_bamos from bamos as bamos_run_bamos with:
    input:
        flair=lambda wildcards: flair_map[wildcards.subject],
        t1=lambda wildcards: t1_map[wildcards.subject],
        gif_done=lambda wildcards: gif_parcellation(wildcards.subject),

use rule run_bamos_correction from bamos as bamos_run_bamos_correction with:
    input:
        bamos_lesion_file=lambda wildcards: f"{output_dir}/{wildcards.subject}/bamos/intermediate/essential/Correct_WS3WT3WC1Lesion_corr.nii.gz",
        bamos_connectivity_file=lambda wildcards: f"{output_dir}/{wildcards.subject}/bamos/intermediate/essential/Connect_WS3WT3WC1Lesion_corr.nii.gz",
        bamos_label_file=lambda wildcards: f"{output_dir}/{wildcards.subject}/bamos/intermediate/essential/TxtLesion_WS3WT3WC1Lesion_corr.txt",
        gif_parcellation_file=lambda wildcards: gif_parcellation(wildcards.subject),

# TRACULA - per-parcellation rules created by for-loop inside tracula_workflow.smk
module tracula:
    snakefile: "tracula_workflow.smk"
    config: config_tracula

use rule * from tracula as tracula_*

# DTI - override to enforce TRACULA brain mask dependency
module dti:
    snakefile: "dti_workflow.smk"
    config: config_dti

use rule * from dti as dti_*

use rule dti from dti as dti_dti with:
    input:
        dwi=lambda wildcards: tracula_dwi(wildcards.subject),
        brain_mask=lambda wildcards: tracula_brain_mask(wildcards.subject),

# NODDI - override to enforce TRACULA brain mask dependency
module noddi:
    snakefile: "noddi_workflow.smk"
    config: config_noddi

use rule * from noddi as noddi_*

use rule noddi from noddi as noddi_noddi with:
    input:
        dwi=lambda wildcards: tracula_dwi(wildcards.subject),
        brain_mask=lambda wildcards: tracula_brain_mask(wildcards.subject),

# TRACT_QC - per-parcellation rules created by for-loop inside tract_qc_workflow.smk
module tract_qc:
    snakefile: "tract_qc_workflow.smk"
    config: config_tract_qc

use rule * from tract_qc as tract_qc_*

# METRICS - per-parcellation rules created by for-loop inside metrics_workflow.smk;
# upstream DAG edges are declared via config dependency maps.
module metrics:
    snakefile: "metrics_workflow.smk"
    config: config_metrics

use rule * from metrics as metrics_*

# ── Z-score (optional cross-subject module) ──────────────────────────────────
if _include_z_score:
    for _s in subjects:
        os.makedirs(f"{output_dir}/{_s}/z_scores/logs",    exist_ok=True)
        os.makedirs(f"{output_dir}/{_s}/z_scores/outputs", exist_ok=True)

    module z_score:
        snakefile: "z_score_workflow.smk"
        config: config_z_score

    use rule * from z_score as z_score_*

# ── Terminal target ───────────────────────────────────────────────────────────
# gif is intentionally NOT listed here: bamos's own rule already declares GIF's
# parcellation output as a hard input (see gif_parcellation_file=... in
# bamos_workflow.smk's run_bamos rule, needed for lesion detection regardless
# of which parcellation is requested for tracula) and will correctly pull it
# in whenever bamos itself still needs to run. Listing it here too would force
# it unconditionally for every subject, including ones whose bamos/metrics are
# already complete and have no remaining need for it -- same redundant-target
# pattern already fixed for recon_all/ThalamicNuclei above.
rule all:
    input:
        rules.recon_all.input,
        rules.bamos_all.input,
        rules.dti_all.input,
        rules.noddi_all.input,
        rules.tract_qc_all.input,
        rules.tracula_all.input,
        rules.metrics_all.input,
        *([rules.z_score_all.input] if _include_z_score else []),
