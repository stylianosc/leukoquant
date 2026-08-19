#!/usr/bin/env python3
"""
Command-line interface for leukoquant.

Atlas conversion toolkit for neuroimaging research.
"""

import click
from typing import Optional
from leukoquant.cli.commands import process_dti, process_gif, process_bamos, process_recon_all, process_noddi, process_atlas_conversion
from leukoquant.cli.commands import process_tracula, process_zscore
from leukoquant.cli.commands import process_metrics, process_tract_qc, process_all
from leukoquant.utils.external_utils import check_sge_plugin


@click.group()
@click.version_option(package_name="leukoquant")
def main():
    """LeukoQuant: Lesion-informed white matter damage metrics toolkit for cerebral small vessel disease."""
    pass


def _exit_on_failure(result: dict) -> dict:
    """Print the error and exit non-zero when a command's result dict reports failure.

    Every process-* command function catches its own exceptions internally
    and returns {"success": False, "error": ...} instead of raising -- so
    without this, Click has no return value handling of its own and a failed
    run (e.g. a required argument missing from both CLI and --config-yaml)
    would otherwise print nothing and exit 0, indistinguishable from success.
    """
    if isinstance(result, dict) and not result.get("success", True):
        click.echo(f"Error: {result.get('error', 'unknown error')}", err=True)
        raise SystemExit(1)
    return result


@main.command('process-gif')
@click.option('--subject', '-i', required=False,
              help='Subject ID or path to a text file with subject IDs (one per line) (optional if --config-yaml provides it)')
@click.option('--t1', '-t', required=False, default=None, type=str,
              help='T1 NIfTI file or {subject} glob pattern (optional if --flair is given)')
@click.option('--flair', '-f', required=False, default=None, type=str,
              help='FLAIR NIfTI file or {subject} glob pattern (optional if --t1 is given)')
@click.option('--output-dir', '-o', required=False, type=click.Path(),
              help='Output directory (optional if --config-yaml provides it)')
@click.option('--verbose', '-v', is_flag=True, help='Verbose output')
@click.option('--mask-file', type=str, help='Mask file or {subject} pattern for GIF segmentation')
@click.option('--scheduler', '-s', type=click.Choice(['local', 'sge'], case_sensitive=False),
              default=None, help="Scheduler to use: 'local' or 'sge' (default: 'local', or from --config-yaml)")
@click.option('--cores', '-c', default=None, type=int, help='Number of cores for Snakemake (default: 1, or from --config-yaml)')
@click.option('--config-yaml', required=False, type=click.Path(exists=True),
              help='Optional YAML config for GIF inputs (CLI options override config values)')
@click.option('--force', is_flag=True, default=False,
              help='Force-rerun GIF segmentation even if outputs already exist')
def process_gif_cmd(subject: Optional[str] = None, output_dir: Optional[str] = None,
                    t1: Optional[str] = None, flair: Optional[str] = None,
                    verbose: bool = False, mask_file: Optional[str] = None,
                    scheduler: Optional[str] = None, cores: Optional[int] = None,
                    config_yaml: Optional[str] = None, force: bool = False) -> dict:
    """Process image(s) with GIF segmentation. At least one of --t1 or --flair is required (via argument or --config-yaml)."""
    check_sge_plugin(scheduler or "local")
    return _exit_on_failure(process_gif(subject, output_dir, t1, flair, verbose, mask_file, scheduler, cores, force, config_yaml))


@main.command('process-bamos')
@click.option('--subject', '-i', required=False,
              help='Subject ID or path to a text file with subject IDs (one per line) (optional if --config-yaml provides it)')
@click.option('--flair', '-f', required=False, type=str,
              help='FLAIR NIfTI file or {subject} glob pattern (optional if --config-yaml provides it)')
@click.option('--t1', '-t', required=False, type=str,
              help='T1 NIfTI file or {subject} glob pattern (optional if --config-yaml provides it)')
@click.option('--gif-results-dir', '-g', required=False, type=str,
              help='GIF results directory or {subject} glob pattern (optional if GIF pipeline will run)')
@click.option('--output-dir', '-o', required=False, type=click.Path(),
              help='Output directory (optional if --config-yaml provides it)')
@click.option('--verbose', '-v', is_flag=True, help='Verbose output')
@click.option('--scheduler', '-s', type=click.Choice(['local', 'sge'], case_sensitive=False),
              default=None, help="Scheduler to use: 'local' or 'sge' (default: 'local', or from --config-yaml)")
@click.option('--cores', '-c', default=None, type=int, help='Number of cores for Snakemake (default: 1, or from --config-yaml)')
@click.option('--config-yaml', required=False, type=click.Path(exists=True),
              help='Optional YAML config for BaMoS inputs (CLI options override config values)')
@click.option('--force', is_flag=True, default=False,
              help='Force-rerun BaMoS lesion detection even if outputs already exist')
def process_bamos_cmd(subject: Optional[str] = None, flair: Optional[str] = None, t1: Optional[str] = None,
                      output_dir: Optional[str] = None, gif_results_dir: Optional[str] = None,
                      verbose: bool = False, scheduler: Optional[str] = None, cores: Optional[int] = None,
                      config_yaml: Optional[str] = None, force: bool = False) -> dict:
    """Process FLAIR and T1 images with BaMoS lesion detection and corrections."""
    check_sge_plugin(scheduler or "local")
    return _exit_on_failure(process_bamos(subject, flair, t1, output_dir, gif_results_dir, verbose, scheduler, cores, force, config_yaml))


@main.command('process-recon')
@click.option('--subject', '-i', required=False,
              help='Subject ID or path to a text file with subject IDs (one per line) (optional if --config-yaml provides it)')
@click.option('--t1', '-t', required=False, type=str,
              help='T1 NIfTI file or {subject} glob pattern (optional if --config-yaml provides it)')
@click.option('--output-dir', '-o', required=False, type=click.Path(),
              help='Output directory for workflow files (optional if --config-yaml provides it)')
@click.option('--verbose', '-v', is_flag=True, help='Verbose output')
@click.option('--scheduler', '-S', type=click.Choice(['local', 'sge'], case_sensitive=False),
              default=None, help="Scheduler to use: 'local' or 'sge' (default: 'local', or from --config-yaml)")
@click.option('--cores', '-c', default=None, type=int, help='Number of cores for Snakemake (default: 1, or from --config-yaml)')
@click.option('--engine', '-e', type=click.Choice(['snakemake', 'nextflow'], case_sensitive=False),
              default='snakemake', show_default=True,
              help='Workflow engine to use (snakemake or nextflow)')
@click.option('--config-yaml', required=False, type=click.Path(exists=True),
              help='Optional YAML config for recon-all inputs (CLI options override config values)')
@click.option('--force', is_flag=True, default=False,
              help='Force-rerun recon-all even if outputs already exist')
def process_recon_cmd(subject: Optional[str] = None, t1: Optional[str] = None, output_dir: Optional[str] = None,
                      verbose: bool = False, scheduler: Optional[str] = None, cores: Optional[int] = None,
                      engine: str = 'snakemake', config_yaml: Optional[str] = None, force: bool = False) -> dict:
    """Run FreeSurfer recon-all via Snakemake."""
    check_sge_plugin(scheduler or "local")
    return _exit_on_failure(process_recon_all(subject, t1, output_dir, verbose, scheduler, cores, engine, force, config_yaml))


@main.command('process-noddi')
@click.option('--subject', '-i', required=False,
              help='Subject ID or path to a text file with subject IDs (one per line) (optional if --config-yaml provides it)')
@click.option('--dwi', required=False, type=str, help='DWI file or {subject} glob pattern (optional if --config-yaml provides it)')
@click.option('--bvecs', required=False, type=str, help='Optional bvecs file or {subject} pattern')
@click.option('--bvals', required=False, type=str, help='Optional bvals file or {subject} pattern')
@click.option('--mask', required=False, default=None, type=str,
              help='Brain mask NIfTI or {subject} glob pattern. If omitted, NODDI is fit on the whole image.')
@click.option('--output-dir', '-o', required=False, type=click.Path(), help='Output directory for NODDI results (optional if --config-yaml provides it)')
@click.option('--skull-strip', is_flag=True,
              help='Auto-generate brain mask using mri_synthstrip (FreeSurfer) before fitting. '
                   'Ignored if --mask is provided.')
@click.option('--verbose', '-v', is_flag=True, help='Verbose output')
@click.option('--scheduler', '-s', type=click.Choice(['local', 'sge'], case_sensitive=False), default=None,
              help="Scheduler to use: 'local' or 'sge' (default: 'local', or from --config-yaml)")
@click.option('--cores', '-c', default=None, type=int, help='Number of cores for Snakemake (default: 1, or from --config-yaml)')
@click.option('--config-yaml', required=False, type=click.Path(exists=True),
              help='Optional YAML config for NODDI inputs (CLI options override config values)')
@click.option('--force', is_flag=True, default=False,
              help='Force-rerun NODDI fitting even if outputs already exist')
def process_noddi_cmd(subject: Optional[str] = None, dwi: Optional[str] = None, bvecs: Optional[str] = None,
                      bvals: Optional[str] = None, mask: Optional[str] = None, output_dir: Optional[str] = None,
                      skull_strip: bool = False, verbose: bool = False,
                      scheduler: Optional[str] = None, cores: Optional[int] = None,
                      config_yaml: Optional[str] = None, force: bool = False) -> dict:
    """Run NODDI fitting for one or more subjects."""
    check_sge_plugin(scheduler or "local")
    return _exit_on_failure(process_noddi(subject=subject, dwi=dwi, output_dir=output_dir, bvecs=bvecs, bvals=bvals,
                         mask_file=mask, skull_strip=skull_strip, verbose=verbose,
                         scheduler=scheduler, cores=cores, force=force, config_yaml=config_yaml))


@main.command('process-zscore')
@click.option('--healthy-list', '-H', required=False, type=click.Path(exists=True), help='Path to file with healthy cohort subject IDs (one per line)')
@click.option('--target-list', '-T', required=False, type=click.Path(exists=True), help='Path to file with target subject IDs (one per line)')
@click.option('--metric', '-m', required=False, multiple=True, help='Metric mapping in the form name=pattern. Can be repeated or comma-separated.')
@click.option('--t1-path', '--t1-pattern', required=False, help='Base:glob or {subject} pattern for T1 images')
@click.option('--demographics-csv', required=False, type=click.Path(exists=True), help='CSV file containing demographics for GLM')
@click.option('--output-dir', '-o', required=False, type=click.Path(), help='Output directory for z-score results')
@click.option('--covariates', required=False, help='Comma-separated covariates (optional)')
@click.option('--poly-terms', required=False, help='Comma-separated polynomial terms (optional)')
@click.option('--metric-space', type=click.Choice(['t1', 'dwi'], case_sensitive=False), default='t1', show_default=True, help='Space in which metrics are defined (t1 or dwi)')
@click.option('--output-space', type=click.Choice(['t1', 'dwi'], case_sensitive=False), default='t1', show_default=True, help='Space in which Z-score outputs are computed (t1 or dwi)')
@click.option('--dwi-pattern', required=False, help='Base:glob or {subject} pattern for DWI images (required if --metric-space=dwi or --output-space=dwi)')
@click.option('--bval-pattern', required=False, help='Base:glob or {subject} pattern for bval files (optional)')
@click.option('--skip-skullstrip-t1', is_flag=True, help='Skip skull stripping for target and healthy T1 images')
@click.option('--skip-skullstrip-dwi', is_flag=True, help='Skip skull stripping for DWI b0 images')
@click.option('--verbose', '-v', is_flag=True, help='Verbose output')
@click.option('--scheduler', '-s', type=click.Choice(['local', 'sge'], case_sensitive=False),
              default=None, help="Scheduler to use: 'local' or 'sge' (default: 'local', or from --config-yaml)")
@click.option('--cores', '-c', default=None, type=int,
              help='Number of cores for Snakemake (default: 1, or from --config-yaml)')
@click.option('--task-concurrency', default=None, type=int,
              help='Max concurrently-running SGE array tasks (-tc). Default: 20. '
                   'Each task writes to the shared NFS output tree on its first mkdir; '
                   'raising this too high on --scheduler sge can overwhelm the NFS server.')
@click.option('--config-yaml', required=False, type=click.Path(exists=True), help='Optional YAML config for z-score inputs (CLI options override config values)')
@click.option('--force', is_flag=True, default=False,
              help='Force-rerun z-score computation even if outputs already exist')
def process_zscore_cmd(healthy_list: Optional[str], target_list: Optional[str], metric: tuple, t1_path: Optional[str], demographics_csv: Optional[str], output_dir: Optional[str], covariates: Optional[str] = None, poly_terms: Optional[str] = None, metric_space: str = 't1', output_space: str = 't1', dwi_pattern: Optional[str] = None, bval_pattern: Optional[str] = None, skip_skullstrip_t1: bool = False, skip_skullstrip_dwi: bool = False, verbose: bool = False, scheduler: Optional[str] = None, cores: Optional[int] = None, task_concurrency: Optional[int] = None, config_yaml: Optional[str] = None, force: bool = False) -> dict:
    """Compute Z-scores for targets against a healthy cohort."""
    check_sge_plugin(scheduler or "local")
    # Parse repeated metric mappings into a dict
    metrics = {}
    metric_parts = []
    for m in metric:
        metric_parts.extend([item.strip() for item in m.split(",") if item.strip()])
    for item in metric_parts:
        if '=' not in item:
            raise click.BadParameter('Metric mappings must use name=pattern format')
        name, pattern = item.split('=', 1)
        metrics[name.strip()] = pattern.strip()

    if not metrics:
        metrics = None

    return _exit_on_failure(process_zscore(
        healthy_subjects_list=healthy_list,
        target_subjects_list=target_list,
        metrics=metrics,
        output_dir=output_dir,
        t1_pattern=t1_path,
        demographics_csv=demographics_csv,
        covariates=covariates,
        polynomial_terms=poly_terms,
        metric_space=metric_space,
        output_space=output_space,
        dwi_pattern=dwi_pattern,
        bval_pattern=bval_pattern,
        skip_skullstrip_t1=skip_skullstrip_t1,
        skip_skullstrip_dwi=skip_skullstrip_dwi,
        verbose=verbose,
        scheduler=scheduler,
        cores=cores,
        task_concurrency=task_concurrency,
        config_yaml=config_yaml,
        force=force,
    ))


@main.command('process-all')
@click.option('--subject', '-i', required=False,
              help='Subject ID or path to a text file with subject IDs (one per line) (optional if --config-yaml provides it)')
@click.option('--t1', '-t', required=False, type=str, help='T1 NIfTI file or {subject} glob pattern (optional if --config-yaml provides it)')
@click.option('--flair', '-f', required=False, type=str, help='FLAIR NIfTI file or {subject} glob pattern (optional if --config-yaml provides it)')
@click.option('--dwi', '-d', required=False, type=str, help='DWI file or {subject} glob pattern (optional if --config-yaml provides it)')
@click.option('--bvecs', required=False, type=str, help='Optional bvecs file or {subject} pattern')
@click.option('--bvals', required=False, type=str, help='Optional bvals file or {subject} pattern')
@click.option('--mask', required=False, type=str, help='Brain mask NIfTI (optional) or {subject} pattern')
@click.option('--output-dir', '-o', required=False, type=click.Path(), help='Output directory (optional if --config-yaml provides it)')
@click.option('--verbose', '-v', is_flag=True, help='Verbose output')
@click.option('--scheduler', '-S', type=click.Choice(['local', 'sge'], case_sensitive=False), default=None,
              help="Scheduler to use: 'local' or 'sge' (default: 'local', or from --config-yaml)")
@click.option('--cores', '-c', default=None, type=int, help='Number of cores for Snakemake (default: 1, or from --config-yaml)')
@click.option('--skip-zscore', is_flag=True, help='Skip Z-score pseudo-healthy computation')
@click.option('--healthy-subjects', required=False, type=click.Path(exists=True), help='Path to file with healthy cohort subject IDs (one per line)')
@click.option('--demographics-csv', required=False, type=click.Path(exists=True), help='CSV file containing demographics for GLM')
@click.option('--covariates', required=False, help='Comma-separated covariates (optional)')
@click.option('--poly-terms', required=False, help='Comma-separated polynomial terms (optional)')
@click.option('--parcellation', default='freesurfer', type=str,
              show_default=True,
              help='Parcellation(s) to use for TRACULA. Comma-separated for multiple '
                   '(e.g. "freesurfer,gif"). "freesurfer" uses recon-all; "gif" uses '
                   'the GIF output atlas converted to FreeSurfer label space.')
@click.option('--config-yaml', required=False, type=click.Path(exists=True),
              help='Optional YAML config for process-all inputs (CLI options override config values)')
@click.option('--force', is_flag=True, default=False,
              help='Force-rerun metrics extraction (for the given --parcellation) even if '
                   'outputs already exist. Anything upstream that is already up to date '
                   '(recon-all, GIF, BaMoS, DTI, NODDI, TRACULA) is left untouched.')
def process_all_cmd(subject: Optional[str] = None, t1: Optional[str] = None, flair: Optional[str] = None,
                    dwi: Optional[str] = None, output_dir: Optional[str] = None, bvecs: str = "",
                    bvals: str = "", mask: Optional[str] = None, verbose: bool = False,
                    scheduler: Optional[str] = None, cores: Optional[int] = None, skip_zscore: bool = False,
                    healthy_subjects: Optional[str] = None, demographics_csv: Optional[str] = None,
                    covariates: Optional[str] = None, poly_terms: Optional[str] = None,
                    parcellation: str = 'freesurfer', config_yaml: Optional[str] = None,
                    force: bool = False) -> dict:
    """Run full pipeline: recon-all, bamos, gif, tracula, dti, noddi, metrics."""
    check_sge_plugin(scheduler or "local")
    return _exit_on_failure(process_all(subject, t1, flair, dwi, bvecs, bvals, output_dir, mask, verbose, scheduler, cores, skip_zscore,
                       healthy_subjects, demographics_csv, covariates, poly_terms, parcellation, force, config_yaml))

@main.command('process-atlas-conversion')
@click.option('--subject', '-i', required=False,
              help='Subject ID (used as the output subdirectory name) (optional if --config-yaml provides it)')
@click.option('--input-parcellation', '-p', required=False, type=click.Path(exists=True),
              help='Input parcellation NIfTI file (optional if --config-yaml provides it)')
@click.option('--mapping-file', '-m', required=False, type=click.Path(exists=True),
              help='CSV label mapping file (optional if --config-yaml provides it)')
@click.option('--output-dir', '-o', type=click.Path(),
              help='Root output directory. Converted atlas is written to '
                   '{output_dir}/{subject}/outputs/converted_atlas.mgz. '
                   'Defaults to the parent directory of --input-parcellation.')
@click.option('--no-validate', is_flag=True, help='Skip post-conversion label validation')
@click.option('--verbose', '-v', is_flag=True, help='Verbose output')
@click.option('--scheduler', '-s', type=click.Choice(['local', 'sge'], case_sensitive=False),
              default=None, help="Scheduler to use: 'local' or 'sge' (default: 'local', or from --config-yaml)")
@click.option('--cores', '-c', default=None, type=int, help='Number of cores for Snakemake (default: 1, or from --config-yaml)')
@click.option('--config-yaml', required=False, type=click.Path(exists=True),
              help='Optional YAML config for atlas conversion inputs (CLI options override config values)')
@click.option('--force', is_flag=True, default=False,
              help='Force-rerun atlas conversion even if outputs already exist')
def convert_atlas_cmd(subject: Optional[str] = None, input_parcellation: Optional[str] = None,
                      mapping_file: Optional[str] = None,
                      output_dir: Optional[str] = None, no_validate: bool = False,
                      verbose: bool = False, scheduler: Optional[str] = None, cores: Optional[int] = None,
                      config_yaml: Optional[str] = None, force: bool = False) -> dict:
    """Convert brain parcellation to FreeSurfer label format via Snakemake."""
    check_sge_plugin(scheduler or "local")
    return _exit_on_failure(process_atlas_conversion(subject, input_parcellation, mapping_file,
                                    output_dir, not no_validate, verbose, scheduler, cores, force, config_yaml))


@main.command('process-dti')
@click.option('--subject', '-i', required=False,
              help='Subject ID or path to a text file with subject IDs (one per line) (optional if --config-yaml provides it)')
@click.option('--dwi', required=False, type=str, help='DWI file or {subject} glob pattern (optional if --config-yaml provides it)')
@click.option('--bvecs', required=False, type=str, help='Optional bvecs file or {subject} pattern')
@click.option('--bvals', required=False, type=str, help='Optional bvals file or {subject} pattern')
@click.option('--mask', required=False, default=None, type=str,
              help='Brain mask NIfTI or {subject} glob pattern. If omitted, DTI is fit on the whole image.')
@click.option('--output-dir', '-o', required=False, type=click.Path(), help='Output directory for DTI results (optional if --config-yaml provides it)')
@click.option('--skull-strip', is_flag=True,
              help='Auto-generate brain mask using mri_synthstrip (FreeSurfer) before fitting. '
                   'Ignored if --mask is provided.')
@click.option('--verbose', '-v', is_flag=True, help='Verbose output')
@click.option('--scheduler', '-s', type=click.Choice(['local', 'sge'], case_sensitive=False), default=None,
              help="Scheduler to use: 'local' or 'sge' (default: 'local', or from --config-yaml)")
@click.option('--cores', '-c', default=None, type=int, help='Number of cores for Snakemake (default: 1, or from --config-yaml)')
@click.option('--config-yaml', required=False, type=click.Path(exists=True),
              help='Optional YAML config for DTI inputs (CLI options override config values)')
@click.option('--force', is_flag=True, default=False,
              help='Force-rerun DTI fitting even if outputs already exist')
def process_dti_cmd(subject: Optional[str] = None, dwi: Optional[str] = None, bvecs: Optional[str] = None,
                    bvals: Optional[str] = None, mask: Optional[str] = None, output_dir: Optional[str] = None,
                    skull_strip: bool = False, verbose: bool = False,
                    scheduler: Optional[str] = None, cores: Optional[int] = None,
                    config_yaml: Optional[str] = None, force: bool = False) -> dict:
    """Run DTI fitting for one or more subjects."""
    check_sge_plugin(scheduler or "local")
    return _exit_on_failure(process_dti(subject=subject, dwi=dwi, output_dir=output_dir, bvecs=bvecs, bvals=bvals,
                       mask_file=mask, skull_strip=skull_strip, verbose=verbose,
                       scheduler=scheduler, cores=cores, force=force, config_yaml=config_yaml))


@main.command('process-tracula')
@click.option('--subject', '-s', required=False, help='Subject ID or path to a text file with subject IDs (one per line) (optional if --config-yaml provides it)')
@click.option('--dwi', required=False, type=str, help='DWI file or {subject} glob pattern (optional if --config-yaml provides it)')
@click.option('--bvecs', required=False, type=str, help='Optional bvecs file or {subject} pattern')
@click.option('--bvals', required=False, type=str, help='Optional bvals file or {subject} pattern')
@click.option('--t1', required=False, default=None, type=str,
              help='T1 NIfTI file or {subject} glob pattern. Required when recon-all has not been run yet.')
@click.option('--freesurfer-recon-dir', required=False, default=None, type=click.Path(),
              help='Root directory containing per-subject recon-all outputs '
                   '({dir}/{subject}/recon-all/outputs). If omitted or if outputs '
                   'are missing for any subject, recon-all is run automatically and --t1 is required.')
@click.option('--output-dir', '-o', required=False, type=click.Path(), help='Output directory for TRACULA results (optional if --config-yaml provides it)')
@click.option('--scratch', required=False, type=click.Path(), help='Optional scratch folder for bedpostx')
@click.option('--verbose', '-v', is_flag=True, help='Verbose output')
@click.option('--scheduler', type=click.Choice(['local', 'sge'], case_sensitive=False), default=None, help="Scheduler to use: 'local' or 'sge' (default: 'local', or from --config-yaml)")
@click.option('--cores', '-c', default=None, type=int, help='Number of cores for Snakemake (default: 1, or from --config-yaml)')
@click.option('--parcellation', default='freesurfer', type=str,
              help='Parcellation(s) to use. Comma-separated for multiple (e.g. "freesurfer,gif"). '
                   '"freesurfer" uses recon-all aparc+aseg; "gif" auto-resolves the GIF parcellation '
                   'from {output_dir}/{subject}/gif/outputs/ and requires process-gif to have run first. '
                   'Any other name requires --brain-parcellation and --mapping-file.')
@click.option('--brain-parcellation', type=click.Path(exists=True),
              help='Path to the brain parcellation NIfTI. Required for non-freesurfer, non-gif '
                   'parcellations. For gif, omit to auto-resolve from the GIF output directory.')
@click.option('--config-yaml', required=False, type=click.Path(exists=True),
              help='Optional YAML config for TRACULA inputs (CLI options override config values)')
@click.option('--force', is_flag=True, default=False,
              help='Force-rerun TRACULA tractography (for the given --parcellation) even if outputs already exist')
def process_tracula_cmd(subject: Optional[str] = None, dwi: Optional[str] = None, bvecs: Optional[str] = None,
                        bvals: Optional[str] = None, t1: Optional[str] = None,
                        freesurfer_recon_dir: Optional[str] = None,
                        output_dir: Optional[str] = None, scratch: Optional[str] = None,
                        verbose: bool = False, scheduler: Optional[str] = None, cores: Optional[int] = None,
                        parcellation: str = 'freesurfer', brain_parcellation: Optional[str] = None,
                        config_yaml: Optional[str] = None, force: bool = False) -> dict:
    """Run TRACULA tractography. Runs recon-all automatically if outputs are missing."""
    check_sge_plugin(scheduler or "local")
    return _exit_on_failure(process_tracula(subject=subject, dwi=dwi, t1=t1,
                           freesurfer_recon_dir=freesurfer_recon_dir,
                           output_dir=output_dir, bvecs=bvecs, bvals=bvals,
                           verbose=verbose, scheduler=scheduler, cores=cores,
                           parcellation=parcellation, brain_parcellation=brain_parcellation,
                           force=force, config_yaml=config_yaml))


@main.command('process-metrics')
@click.option('--subject', required=False,
              help='Subject ID or path to a txt file with subject IDs (one per line) (optional if --config-yaml provides it)')
@click.option('--tractography-path', required=False,
              help='Tractography pattern as base:glob:space. Example: ./examples/outputs/tracula_results:/dpath/*/path.pd.nii.gz:dwi (optional if --config-yaml provides it)')
@click.option('--t1-path', required=False,
              help='T1 pattern as base:glob:space. Example: ./examples/sample_files/SCANS/:/T1/I*.nii.gz:t1')
@click.option('--dwi-path', required=False,
              help='DWI pattern as base:glob:space. Example: ./examples/sample_files/SCANS/:/DWI/data.nii.gz:dwi')
@click.option('--lesion-path', required=False,
              help=(
                  'Lesion pattern(s) as [name=]base:glob[:space], comma-separated for multiple types. '
                  'Name (optional, default "lesion") is separated from the path spec by "=". '
                  'Example (single, unnamed): ./bamos_results:/bamos/CorrectLesion_*.nii.gz '
                  'Example (named): wmh=./bamos_results:/bamos/WMH_*.nii.gz:t1,lacunes=./lac_results:/lac_*.nii.gz:t1'
              ))
@click.option('--maps', '-m', required=False, multiple=True,
              help=(
                  'Metric mappings: metric_name=base:glob:space. Can be repeated or comma-separated. '
                  'Space is "dwi", "t1", or "atlas". '
                  'Example: --maps dti_fa=./data/dti:/fa.nii.gz:dwi dti_md=./data/dti:/md.nii.gz:dwi'
              ))
@click.option('--tract-mode', type=click.Choice(['tractography', 'tractography-atlas', 'atlas'], case_sensitive=False),
              default=None,
              help="Tract analysis mode: tractography, tractography-atlas, or atlas (default: tractography-atlas, or from --config-yaml)")
@click.option('--parcellation', default='freesurfer', type=str, show_default=True,
              help='Parcellation(s) to extract metrics for. Comma-separated for multiple '
                   '(e.g. "freesurfer,gif").')
@click.option('--output-dir', '-o', required=False, type=click.Path(),
              help='Output directory for metrics CSV files (optional if --config-yaml provides it)')
@click.option('--scheduler', type=click.Choice(['local', 'sge'], case_sensitive=False),
              default=None, help="Scheduler to use: 'local' or 'sge' (default: 'local', or from --config-yaml)")
@click.option('--cores', '-c', default=None, type=int,
              help='Number of cores for Snakemake (default: 1, or from --config-yaml)')
@click.option('--config-yaml', required=False, type=click.Path(exists=True),
              help='Optional YAML config for metrics inputs (CLI options override config values)')
@click.option('--force', is_flag=True, default=False,
              help='Force-rerun metrics extraction (for the given --parcellation) even if outputs already exist')
@click.option('--verbose', '-v', is_flag=True, help='Verbose output')
def process_metrics_cmd(subject: Optional[str] = None,
                        tractography_path: Optional[str] = None,
                        t1_path: Optional[str] = None,
                        dwi_path: Optional[str] = None,
                        lesion_path: Optional[str] = None,
                        maps: tuple = (),
                        tract_mode: Optional[str] = None,
                        parcellation: str = 'freesurfer',
                        output_dir: Optional[str] = None,
                        scheduler: Optional[str] = None,
                        cores: Optional[int] = None,
                        config_yaml: Optional[str] = None,
                        force: bool = False,
                        verbose: bool = False) -> dict:
    """Compute metrics along tracts, lesions, and WMH regions."""
    check_sge_plugin(scheduler or "local")
    # Parse repeated metric mappings into a dict
    metrics = {}
    metric_parts = []
    for m in maps:
        metric_parts.extend([item.strip() for item in m.split(",") if item.strip()])

    for item in metric_parts:
        if '=' not in item:
            raise click.BadParameter('Metric mappings must use name=pattern:space format')
        name, pattern = item.split('=', 1)
        metrics[name.strip()] = pattern.strip()

    if not metrics:
        metrics = None

    return _exit_on_failure(process_metrics(
        subject=subject,
        tractography_path=tractography_path,
        t1_path=t1_path,
        dwi_path=dwi_path,
        lesion_path=lesion_path,
        metrics=metrics,
        tract_mode=tract_mode,
        parcellation=parcellation,
        output_dir=output_dir,
        scheduler=scheduler,
        cores=cores,
        force=force,
        verbose=verbose,
        config_yaml=config_yaml,
    ))


@main.command('process-tract-qc')
@click.option('--subject', required=False, help='Subject ID or path to a txt file with subject IDs (one per line) (optional if --config-yaml provides it)')
@click.option('--tractography-path', required=False, help='Tractography pattern as base:glob:space. Example: ./outputs/tracula_results:/dpath/*/path.pd.nii.gz:dwi (optional if --config-yaml provides it)')
@click.option('--output-dir', '-o', required=False, type=click.Path(), help='Output directory for tract QC results (optional, defaults to current directory)')
@click.option('--scheduler', type=click.Choice(['local', 'sge'], case_sensitive=False), default=None, help="Scheduler to use: 'local' or 'sge' (default: 'local', or from --config-yaml)")
@click.option('--cores', '-c', default=None, type=int, help='Number of cores for Snakemake (default: 1, or from --config-yaml)')
@click.option('--config-yaml', required=False, type=click.Path(exists=True),
              help='Optional YAML config for tract QC inputs (CLI options override config values)')
@click.option('--force', is_flag=True, default=False,
              help='Force-rerun tract QC even if outputs already exist')
@click.option('--verbose', '-v', is_flag=True, help='Verbose output')
def process_tract_qc_cmd(subject: Optional[str] = None, tractography_path: Optional[str] = None,
                         output_dir: Optional[str] = None, scheduler: Optional[str] = None,
                         cores: Optional[int] = None, config_yaml: Optional[str] = None,
                         force: bool = False, verbose: bool = False) -> dict:
    """Run tractography QC workflow for a subject or batch."""
    check_sge_plugin(scheduler or "local")
    return _exit_on_failure(process_tract_qc(subject=subject, tractography_path=tractography_path, output_dir=output_dir, scheduler=scheduler, cores=cores, force=force, verbose=verbose, config_yaml=config_yaml))


if __name__ == '__main__':
    main()
