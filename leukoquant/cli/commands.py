#!/usr/bin/env python3
"""
CLI command implementations for leukoquant.

This module contains the implementation logic for CLI commands.
"""

import sys
from pathlib import Path
from typing import List, Optional
import os


# Each standalone process-* command has one (or, for parcellation-templated
# workflows, one-per-parcellation) Snakemake rule that actually does the work.
# --force on these commands is a plain boolean -- callers shouldn't need to
# know Snakemake rule names -- so each command maps it to its own known
# rule(s) internally. process-all's --force stays rule-name-based since it
# spans multiple stages and there's no single obvious rule to force.
_FORCE_RULE_NAMES = {
    "gif": ["run_gif"],
    "bamos": ["run_bamos"],
    "recon_all": ["recon_all"],
    "noddi": ["noddi"],
    "atlas_conversion": ["convert_atlas"],
    "dti": ["dti"],
    "zscore": ["z_score"],
    "tract_qc": ["tract_qc"],
}


def _force_rules_for(processor: str, force: bool) -> Optional[List[str]]:
    """Return the rule name(s) to force for a fixed-rule-name processor."""
    if not force:
        return None
    return _FORCE_RULE_NAMES[processor]


def _force_rules_for_parcellation(rule_prefix: str, parcellation: str, force: bool) -> Optional[List[str]]:
    """Return per-parcellation rule name(s) to force (metrics, tracula)."""
    if not force:
        return None
    return [f"{rule_prefix}_{p.strip()}" for p in parcellation.split(",") if p.strip()]


from leukoquant.core.gif_processor import apply_gif
from leukoquant.core.bamos_processor import apply_bamos
from leukoquant.core.freesurfer_processor import apply_recon_all
from leukoquant.core.noddi_processor import apply_noddi
from leukoquant.core.atlas_converter_processor import apply_atlas_conversion
from leukoquant.core.dti_processor import apply_dti
from leukoquant.core.tracula_processor import apply_tracula
from leukoquant.core.zscore_processor import apply_zscore
from leukoquant.core.metrics_processor import apply_metrics
from leukoquant.core.tract_qc_processor import apply_tract_qc
from leukoquant.core.process_all_processor import apply_process_all

def process_all(subject: Optional[str] = None,
                t1: Optional[str] = None,
                flair: Optional[str] = None,
                dwi: Optional[str] = None,
                bvecs: Optional[str] = None,
                bvals: Optional[str] = None,
                output_dir: Optional[str] = None,
                mask: Optional[str] = None,
                verbose: bool = False,
                scheduler: str = "local",
                cores: int = 1,
                skip_zscore: bool = False,
                healthy_subset: Optional[str] = None,
                demographics_csv: Optional[str] = None,
                covariates: Optional[str] = None,
                poly_terms: Optional[str] = None,
                parcellation: str = "freesurfer",
                force: bool = False,
                config_yaml: Optional[str] = None) -> dict:
    """Run full pipeline process-all.

    force: Force-rerun metrics extraction (for the given --parcellation)
        even if outputs already exist. Snakemake scopes this to the metrics
        rule(s) plus anything downstream, so subjects whose earlier stages
        (recon-all/DTI/GIF/tractography) are already complete have those
        stages correctly skipped.
    """
    try:
        if verbose:
            print("Running full pipeline process-all")

        results = apply_process_all(
            subject_input=subject, t1_pattern=t1, flair_pattern=flair,
            dwi_pattern=dwi, bvecs_pattern=bvecs, bvals_pattern=bvals, mask_pattern=mask,
            output_dir=output_dir, scheduler=scheduler, cores=cores,
            skip_zscore=skip_zscore, healthy_list=healthy_subset, demographics=demographics_csv,
            covariates=covariates, poly_terms=poly_terms,
            parcellation=parcellation,
            force_rules=_force_rules_for_parcellation("metrics_extract_metrics", parcellation, force),
            verbose=verbose,
            config_yaml=config_yaml,
        )
        return results
    except Exception as e:
        return {"success": False, "error": str(e)}


def process_gif(subject: Optional[str] = None, output_dir: Optional[str] = None,
               t1: Optional[str] = None, flair: Optional[str] = None,
               verbose: bool = False, mask_file: Optional[str] = None,
               scheduler: str = "local", cores: int = 1,
               force: bool = False, config_yaml: Optional[str] = None) -> dict:
    """Process image(s) with GIF segmentation using Snakemake."""
    try:
        if verbose:
            print("Processing GIF")
            print(f"Subject: {subject}")
            if t1:
                print(f"T1 pattern: {t1}")
            if flair:
                print(f"FLAIR pattern: {flair}")
            print(f"Output directory: {output_dir}")

        results = apply_gif(
            subject_input=subject,
            t1_pattern=t1,
            flair_pattern=flair,
            output_dir=output_dir,
            mask_pattern=mask_file,
            scheduler=scheduler,
            cores=cores,
            force_rules=_force_rules_for("gif", force),
            verbose=verbose,
            config_yaml=config_yaml,
        )

        if not results.get("success"):
            sys.exit(1)

        return results

    except Exception as e:
        return {"success": False, "error": str(e), "output_dir": output_dir}

def process_bamos(subject: Optional[str] = None, flair: Optional[str] = None, t1: Optional[str] = None,
                 output_dir: Optional[str] = None, gif_results_dir: Optional[str] = None,
                 verbose: bool = False, scheduler: str = "local", cores: int = 1,
                 force: bool = False, config_yaml: Optional[str] = None) -> dict:
    """Process FLAIR and T1 images with BaMoS lesion detection and corrections."""
    try:
        if verbose:
            print("Processing BaMoS")
            print(f"Subject: {subject}")
            print(f"FLAIR pattern: {flair}")
            print(f"T1 pattern: {t1}")
            print(f"GIF results pattern: {gif_results_dir}")
            print(f"Output directory: {output_dir}")

        results = apply_bamos(
            subject_input=subject,
            flair_pattern=flair,
            t1_pattern=t1,
            output_dir=output_dir,
            gif_results_pattern=gif_results_dir,
            scheduler=scheduler,
            cores=cores,
            force_rules=_force_rules_for("bamos", force),
            verbose=verbose,
            config_yaml=config_yaml,
        )

        return results

    except Exception as e:
        return {"success": False, "error": str(e), "output_dir": output_dir}


def process_recon_all(subject: Optional[str] = None,
                      t1: Optional[str] = None,
                      output_dir: Optional[str] = None,
                      verbose: bool = False,
                      scheduler: str = "local",
                      cores: int = 1,
                      engine: str = "snakemake",
                      force: bool = False,
                      config_yaml: Optional[str] = None) -> dict:
    """Run FreeSurfer recon-all for one or more subjects via Snakemake."""
    try:
        if verbose:
            print("Running FreeSurfer recon-all")
            print(f"Subject: {subject}")
            print(f"T1 pattern: {t1}")
            print(f"Output dir: {output_dir}")
            print(f"Scheduler: {scheduler}")
            print(f"Cores: {cores}")
            print(f"Engine: {engine}")

        results = apply_recon_all(
            subject_input=subject,
            t1_pattern=t1,
            output_dir=output_dir,
            scheduler=scheduler,
            cores=cores,
            force_rules=_force_rules_for("recon_all", force),
            verbose=verbose,
            config_yaml=config_yaml,
        )

        return results

    except Exception as e:
        return {"success": False, "error": str(e), "t1_pattern": t1, "output_dir": output_dir}


def process_noddi(subject: Optional[str] = None,
                  dwi: Optional[str] = None,
                  output_dir: Optional[str] = None,
                  bvecs: Optional[str] = None,
                  bvals: Optional[str] = None,
                  mask_file: Optional[str] = None,
                  skull_strip: bool = False,
                  verbose: bool = False,
                  scheduler: str = "local",
                  cores: int = 1,
                  force: bool = False,
                  config_yaml: Optional[str] = None) -> dict:
    """Run NODDI fitting for one or more subjects."""
    try:
        if verbose:
            print("Running NODDI fitting")
            print(f"Subject: {subject}")
            print(f"DWI pattern: {dwi}")
            if bvecs: print(f"Bvecs pattern: {bvecs}")
            if bvals: print(f"Bvals pattern: {bvals}")
            if mask_file: print("Mask pattern: " + str(mask_file))
            if skull_strip: print("Skull stripping: enabled")
            print(f"Output dir: {output_dir}")

        results = apply_noddi(
            subject_input=subject,
            dwi_pattern=dwi,
            mask_pattern=mask_file,
            bvecs_pattern=bvecs or "",
            bvals_pattern=bvals or "",
            output_dir=output_dir,
            scheduler=scheduler,
            cores=cores,
            skull_strip=skull_strip,
            force_rules=_force_rules_for("noddi", force),
            verbose=verbose,
            config_yaml=config_yaml,
        )

        return results

    except Exception as e:
        return {"success": False, "error": str(e), "dwi": dwi, "output_dir": output_dir}


def process_atlas_conversion(subject: Optional[str] = None,
                           input_parcellation: Optional[str] = None,
                           mapping_file: Optional[str] = None,
                           output_dir: Optional[str] = None,
                           validate: bool = True,
                           verbose: bool = False,
                           scheduler: str = "local",
                           cores: int = 1,
                           force: bool = False,
                           config_yaml: Optional[str] = None) -> dict:
    """Run atlas conversion via Snakemake.

    Args:
        subject:            Subject ID used as the output subdirectory name.
        input_parcellation: Path to the input parcellation NIfTI file.
        mapping_file:       Path to the CSV label mapping file.
        output_dir:         Root output directory. Converted atlas is written to
                            {output_dir}/{subject}/outputs/converted_atlas.mgz.
        validate:           Run post-conversion label validation.
        verbose:            Enable verbose output.
        scheduler:          Execution scheduler ("local" or "sge").
        cores:              Number of cores for Snakemake.

    Returns:
        Dictionary with processing results.
    """
    try:
        if verbose:
            print("Running atlas conversion")
            print(f"Subject: {subject}")
            print(f"Input parcellation: {input_parcellation}")
            print(f"Mapping file: {mapping_file}")
            print(f"Output directory: {output_dir}")

        results = apply_atlas_conversion(
            subject=subject,
            input_parcellation=input_parcellation,
            mapping_file=mapping_file,
            output_dir=output_dir,
            validate=validate,
            verbose=verbose,
            scheduler=scheduler,
            cores=cores,
            force_rules=_force_rules_for("atlas_conversion", force),
            config_yaml=config_yaml,
        )

        return results

    except Exception as e:
        return {
            "success": False,
            "error": str(e),
            "subject": subject,
            "input_parcellation": input_parcellation,
            "mapping_file": mapping_file,
            "output_dir": output_dir,
        }


def process_dti(subject: Optional[str] = None,
                  dwi: Optional[str] = None,
                  output_dir: Optional[str] = None,
                  bvecs: Optional[str] = None,
                  bvals: Optional[str] = None,
                  mask_file: Optional[str] = None,
                  skull_strip: bool = False,
                  verbose: bool = False,
                  scheduler: str = "local",
                  cores: int = 1,
                  force: bool = False,
                  config_yaml: Optional[str] = None) -> dict:
    """Run DTI fitting for one or more subjects."""
    try:
        if verbose:
            print("Running DTI fitting")
            print(f"Subject: {subject}")
            print(f"DWI pattern: {dwi}")
            if bvecs: print(f"Bvecs pattern: {bvecs}")
            if bvals: print(f"Bvals pattern: {bvals}")
            if mask_file: print("Mask pattern: " + str(mask_file))
            if skull_strip: print("Skull stripping: enabled")
            print(f"Output dir: {output_dir}")

        results = apply_dti(
            subject_input=subject,
            dwi_pattern=dwi,
            mask_pattern=mask_file,
            bvecs_pattern=bvecs or "",
            bvals_pattern=bvals or "",
            output_dir=output_dir,
            scheduler=scheduler,
            cores=cores,
            skull_strip=skull_strip,
            force_rules=_force_rules_for("dti", force),
            verbose=verbose,
            config_yaml=config_yaml,
        )

        return results

    except Exception as e:
        return {"success": False, "error": str(e), "dwi": dwi, "output_dir": output_dir}


def process_zscore(healthy_subjects_list: Optional[str],
                   target_subjects_list: Optional[str],
                   metrics: Optional[dict],
                   output_dir: Optional[str],
                   t1_pattern: Optional[str],
                   demographics_csv: Optional[str],
                   covariates: Optional[str] = None,
                   polynomial_terms: Optional[str] = None,
                   metric_space: str = 't1',
                   output_space: str = 't1',
                   dwi_pattern: Optional[str] = None,
                   bval_pattern: Optional[str] = None,
                   skip_skullstrip_t1: bool = False,
                   skip_skullstrip_dwi: bool = False,
                   verbose: bool = False,
                   scheduler: Optional[str] = None,
                   cores: Optional[int] = None,
                   task_concurrency: Optional[int] = None,
                   config_yaml: Optional[str] = None,
                   force: bool = False) -> dict:
    """Run z-score workflow.

    Args:
        healthy_subjects_list: path to file listing healthy cohort subjects
        target_subjects_list: path to file listing target subjects
        metrics: dict mapping metric name -> pattern (must include {subject})
        output_dir: path for workflow outputs
        t1_pattern: pattern to locate subject T1s (use {subject})
        demographics_csv: CSV file with demographics used by GLM
        covariates: optional comma-separated covariates
        polynomial_terms: optional comma-separated polynomial terms with power notation (default power is 2, e.g., "age,age:2,bmi^3")
        metric_space: 't1' or 'dwi' (default: 't1')
        output_space: 't1' or 'dwi' (default: 't1')
        dwi_pattern: pattern to locate subject DWIs (use {subject}), required if metric_space='dwi'
        bval_pattern: pattern to locate bval files (use {subject}), optional
        skip_skullstrip_t1: if True, skip T1 skull stripping
        skip_skullstrip_dwi: if True, skip DWI b0 skull stripping
        verbose: enable verbose logging
        scheduler: 'local' or 'sge'
        cores: number of cores for Snakemake
        task_concurrency: max concurrently-running SGE array tasks (-tc).
            Defaults to 20 (see z_score_workflow.smk) to avoid overwhelming
            the shared NFS output tree when many array tasks start at once
            under --immediate-submit.

    Returns:
        dict with success, results_dir or error
    """
    try:
        if verbose:
            print("Running Z-score workflow")
            print(f"Healthy subjects list: {healthy_subjects_list}")
            print(f"Target subjects list: {target_subjects_list}")
            print(f"T1 pattern: {t1_pattern}")
            print(f"Metrics: {metrics}")
            print(f"Metric space: {metric_space}")
            print(f"Output space: {output_space}")
            print(f"DWI pattern: {dwi_pattern}")
            print(f"Skip T1 skullstrip: {skip_skullstrip_t1}")
            print(f"Skip DWI skullstrip: {skip_skullstrip_dwi}")
            print(f"Output dir: {output_dir}")

        results = apply_zscore(
            healthy_subjects_list=healthy_subjects_list,
            target_subjects_list=target_subjects_list,
            metrics=metrics,
            output_dir=output_dir,
            t1_pattern=t1_pattern,
            demographics_csv=demographics_csv,
            covariates=covariates,
            polynomial_terms=polynomial_terms,
            metric_space=metric_space,
            output_space=output_space,
            dwi_pattern=dwi_pattern,
            bval_pattern=bval_pattern,
            skip_skullstrip_t1=skip_skullstrip_t1,
            skip_skullstrip_dwi=skip_skullstrip_dwi,
            scheduler=scheduler,
            cores=cores,
            task_concurrency=task_concurrency,
            verbose=verbose,
            config_yaml=config_yaml,
            force_rules=_force_rules_for("zscore", force),
        )

        return results

    except Exception as e:
        return {
            "success": False,
            "error": str(e),
            "healthy_subjects_list": healthy_subjects_list,
            "target_subjects_list": target_subjects_list,
            "output_dir": output_dir,
        }


def process_tracula(subject: Optional[str] = None,
                    dwi: Optional[str] = None,
                    output_dir: Optional[str] = None,
                    freesurfer_recon_dir: Optional[str] = None,
                    t1: Optional[str] = None,
                    bvecs: Optional[str] = None,
                    bvals: Optional[str] = None,
                    verbose: bool = False,
                    scheduler: str = "local",
                    cores: int = 1,
                    parcellation: str = "freesurfer",
                    brain_parcellation: Optional[str] = None,
                    scratch: Optional[str] = None,
                    force: bool = False,
                    config_yaml: Optional[str] = None) -> dict:
    """Run the TRACULA workflow for one or more subjects.

    Args:
        subject:              Subject ID or path to a text file (one ID per line).
        dwi:                  DWI file or {subject} glob pattern.
        output_dir:           Directory for TRACULA outputs and config.
        freesurfer_recon_dir: Root containing per-subject recon-all outputs.
                              If omitted or missing for any subject, recon-all
                              is run automatically and t1 is required.
        t1:                   T1 file or {subject} pattern. Required when
                              recon-all has not been run yet.
        bvecs:                Optional bvecs file or {subject} pattern.
        bvals:                Optional bvals file or {subject} pattern.
    """
    try:
        if verbose:
            print(f"Running TRACULA for subject: {subject}")
            print(f"DWI input: {dwi}")
            if freesurfer_recon_dir:
                print(f"FreeSurfer recon dir: {freesurfer_recon_dir}")
            print(f"Output dir: {output_dir}")

        results = apply_tracula(
            subject_input=subject,
            output_dir=output_dir,
            dwi_pattern=dwi,
            freesurfer_recon_dir=freesurfer_recon_dir,
            t1_pattern=t1,
            bvecs_pattern=bvecs,
            bvals_pattern=bvals,
            scheduler=scheduler,
            cores=cores,
            verbose=verbose,
            parcellation=parcellation,
            brain_parcellation=brain_parcellation,
            force_rules=_force_rules_for_parcellation("trac_all", parcellation, force),
            config_yaml=config_yaml,
        )

        return results

    except Exception as e:
        return {
            "success": False,
            "error": str(e),
            "subject": subject,
            "dwi": dwi,
            "freesurfer_recon_dir": freesurfer_recon_dir,
            "output_dir": output_dir
        }


def process_metrics(subject: Optional[str] = None,
                    tractography_path: Optional[str] = None,
                    t1_path: Optional[str] = None,
                    dwi_path: Optional[str] = None,
                    lesion_path: Optional[str] = None,
                    metrics: Optional[dict] = None,
                    tract_mode: Optional[str] = None,
                    output_dir: Optional[str] = None,
                    scheduler: str = "local",
                    cores: int = 1,
                    parcellation: str = "freesurfer",
                    force: bool = False,
                    verbose: bool = False,
                    config_yaml: Optional[str] = None) -> dict:
    """Calculate metrics along tracts, lesions, and WMH regions.

    Args:
        subject: Subject ID or path to txt file (one subject per line)
        tractography_path: Tractography pattern in base:glob:space format
        t1_path: T1 pattern in base:glob:space format (optional)
        lesion_path: Lesion pattern in base:glob format (optional)
        metrics: Dict of metric_name -> base:glob:space patterns
        output_dir: Output directory for CSV files
        scheduler: Execution scheduler ("local" or "sge")
        cores: Number of cores for Snakemake
        parcellation: Comma-separated parcellation(s) to extract metrics for
            (e.g. "freesurfer" or "gif")
        force: Force-rerun metrics extraction even if outputs already exist
        verbose: Enable verbose output
    """
    try:
        if verbose:
            print("Running process-metrics")
            print(f"Subject input: {subject}")
            print(f"Tractography path: {tractography_path}")
            print(f"T1 path: {t1_path}")
            print(f"DWI path: {dwi_path}")
            print(f"Lesion path: {lesion_path}")
            print(f"Metrics: {metrics}")
            print(f"Output dir: {output_dir}")

        results = apply_metrics(
            subject=subject,
            tractography_path=tractography_path,
            t1_path=t1_path,
            dwi_path=dwi_path,
            lesion_path=lesion_path,
            metrics=metrics,
            tract_mode=tract_mode,
            output_dir=output_dir,
            scheduler=scheduler,
            cores=cores,
            parcellation=parcellation,
            force_rules=_force_rules_for_parcellation("extract_metrics", parcellation, force),
            verbose=verbose,
            config_yaml=config_yaml,
        )

        return results

    except Exception as e:
        return {
            "success": False,
            "error": str(e),
            "subject": subject,
            "tractography_path": tractography_path,
            "tract_mode": tract_mode,
            "output_dir": output_dir,
        }

def process_tract_qc(subject: Optional[str] = None,
                    tractography_path: Optional[str] = None,
                    output_dir: str = os.getcwd(),
                    scheduler: str = "local",
                    cores: int = 1,
                    force: bool = False,
                    verbose: bool = False,
                    config_yaml: Optional[str] = None) -> dict:
    """Run tractography QC workflow for a subject or batch.

    Args:
        subject: Subject ID or path to file with subject IDs
        tractography_path: Path to tractography root directory
        output_dir: Output directory for QC results
        scheduler: 'local' or 'sge'
        cores: Number of cores
        force: Force-rerun tract QC even if outputs already exist
        verbose: Verbose output

    Returns:
        Dictionary with success status and results
    """
    try:
        if verbose:
            print("Running Tractography QC workflow")
            print(f"Subject: {subject}")
            print(f"Tractography path: {tractography_path}")
            print(f"Output dir: {output_dir}")
            print(f"Scheduler: {scheduler}")
            print(f"Cores: {cores}")

        results = apply_tract_qc(
            subject=subject,
            tractography_path=tractography_path,
            output_dir=output_dir,
            scheduler=scheduler,
            cores=cores,
            force_rules=_force_rules_for("tract_qc", force),
            verbose=verbose,
            config_yaml=config_yaml,
        )

        return results

    except Exception as e:
        return {
            "success": False,
            "error": str(e),
            "subject": subject,
            "tractography_path": tractography_path,
            "output_dir": output_dir
        }