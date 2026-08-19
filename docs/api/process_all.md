# leukoquant process-all

Run the full LeukoQuant pipeline in the correct dependency order.

## Synopsis

```
leukoquant process-all [OPTIONS]
```

## Description

`process-all` is the top-level orchestration command.  It runs all pipeline
stages internally via Snakemake module composition, resolving inter-stage
dependencies automatically:

```
T1 + FLAIR + DWI
    │
    ├─ process-gif         (from T1 / FLAIR)
    ├─ process-recon       (from T1)          - freesurfer parcellation only
    ├─ process-bamos       (from T1 + FLAIR + GIF)
    ├─ process-dti         (from DWI)
    ├─ process-noddi       (from DWI)
    ├─ process-atlas-conversion  (from GIF)   - gif parcellation only
    │
    ├─ process-tracula     (from DWI + recon/gif)
    │
    ├─ process-metrics     (from tracula + bamos + dti/noddi)
    ├─ process-tract-qc    (from tracula)
    └─ process-zscore      (from dti/noddi + healthy cohort)   - optional
```

### Parcellation support

`--parcellation` controls which atlas is used for TRACULA.  Multiple
parcellations can be requested in a single run as a comma-separated list.

| `--parcellation` | Description |
|------------------|-------------|
| `freesurfer` (default) | FreeSurfer `aparc+aseg`-guided tractography |
| `gif` | GIF atlas-guided tractography |
| `freesurfer,gif` | Both; two independent tract sets and two metric outputs |

When multiple parcellations are used, metrics and tract QC outputs are written
to `metrics-{parcellation}/` and `tract_qc-{parcellation}/` subfolders,
respectively.

### Z-score

Z-score computation is enabled when `--healthy-subjects` is provided.
It can be explicitly suppressed with `--skip-zscore`.

Outputs are written to `{output_dir}/{subject}/` with one subfolder per stage.

## Parameters

| Flag | Short | Type | Default | Required | Description |
|------|-------|------|---------|----------|-------------|
| `--subject` | `-i` | `str` | - | **Yes**\* | Subject ID, or path to a `.txt` file listing one subject ID per line |
| `--t1` | `-t` | `str` | - | **Yes**\* | T1 NIfTI file or `{subject}` glob pattern |
| `--flair` | `-f` | `str` | - | **Yes**\* | FLAIR NIfTI file or `{subject}` glob pattern |
| `--dwi` | `-d` | `str` | - | **Yes**\* | DWI NIfTI (or `.zip` DICOM archive) file or `{subject}` glob pattern |
| `--bvecs` | - | `str` | `""` | No | bvecs file or `{subject}` glob pattern |
| `--bvals` | - | `str` | `""` | No | bvals file or `{subject}` glob pattern |
| `--mask` | - | `str` | `None` | No | Brain mask NIfTI or `{subject}` glob pattern (optional) |
| `--parcellation` | - | `str` | `freesurfer` | No | Parcellation(s) for TRACULA. Comma-separated for multiple (e.g. `freesurfer,gif`). |
| `--healthy-subjects` | - | `path` | `None` | No | Path to healthy cohort subject list. Enables Z-score computation when provided. |
| `--demographics-csv` | - | `path` | `None` | No | CSV with per-subject demographics for the Z-score GLM |
| `--covariates` | - | `str` | `None` | No | Comma-separated covariate columns for the GLM (e.g. `age,sex`) |
| `--poly-terms` | - | `str` | `None` | No | Comma-separated polynomial expansion terms (e.g. `age:2`) |
| `--skip-zscore` | - | flag | `False` | No | Explicitly skip Z-score computation even when `--healthy-subjects` is provided |
| `--output-dir` | `-o` | `path` | - | **Yes**\* | Root output directory |
| `--scheduler` | `-S` | `choice` | `local` | No | Execution scheduler: `local` or `sge` |
| `--cores` | `-c` | `int` | `1` | No | Number of Snakemake cores |
| `--config-yaml` | - | `path` | `None` | No | YAML file supplying defaults for any of the above - see [Config file](#config-file-config-yaml) below |
| `--force` | - | flag | `False` | No | Force-rerun metrics extraction (for `--parcellation`) even if outputs already exist. Anything upstream that's already up to date (recon-all, GIF, BaMoS, DTI, NODDI, TRACULA) is left untouched. |
| `--verbose` | `-v` | flag | `False` | No | Print Snakemake stdout to the console |

\* Required overall, but may be supplied via `--config-yaml` instead of the CLI flag - see below.

### Parcellation options

| Value | Description |
|-------|-------------|
| `freesurfer` | FreeSurfer `aparc+aseg` atlas (default) |
| `gif` | GIF multi-label parcellation |
| `freesurfer,gif` | Both parcellations; two parallel tract sets and metric outputs |

### Scheduler options

| Value | Description |
|-------|-------------|
| `local` | Run all jobs on the current machine |
| `sge` | Submit jobs to an SGE cluster |

## Config file (`--config-yaml`)

Instead of (or alongside) CLI flags, `--subject`, `--t1`, `--flair`, `--dwi`,
and `--output-dir` may be supplied by a YAML file. **CLI flags always take
priority** over the equivalent YAML key when both are given - the YAML file
only fills in whatever the CLI left unset.

Recognized keys:

| YAML key | Equivalent flag |
|----------|------------------|
| `subject` (or `subject_input`) | `--subject` |
| `t1` (or `t1_pattern`) | `--t1` |
| `flair` (or `flair_pattern`) | `--flair` |
| `dwi` (or `dwi_pattern`) | `--dwi` |
| `bvecs` (or `bvecs_pattern`) | `--bvecs` |
| `bvals` (or `bvals_pattern`) | `--bvals` |
| `mask` (or `mask_pattern`) | `--mask` |
| `parcellation` | `--parcellation` |
| `healthy_list` (or `healthy_subjects_list`) | `--healthy-subjects` |
| `demographics` (or `demographics_csv`) | `--demographics-csv` |
| `covariates` | `--covariates` |
| `poly_terms` (or `polynomial_terms`) | `--poly-terms` |
| `skip_zscore` | `--skip-zscore` |
| `output_dir` | `--output-dir` |
| `scheduler` | `--scheduler` |
| `cores` | `--cores` |

```yaml
# process_all_config.yaml
subject: ./subjects.txt
t1: "./data/{subject}/T1/I*.nii.gz"
flair: "./data/{subject}/FLAIR/I*.nii.gz"
dwi: "./data/{subject}/DWI/I*.zip"
output_dir: ./outputs
scheduler: sge
cores: 8
```

```bash
leukoquant process-all --config-yaml process_all_config.yaml

# override just the parcellation for this run; everything else still comes from the YAML
leukoquant process-all --config-yaml process_all_config.yaml --parcellation gif
```

## Returns

`dict` with keys:

| Key | Type | Description |
|-----|------|-------------|
| `success` | `bool` | `True` if the full pipeline completed without error |
| `results_dir` | `str` | Absolute path to the output directory |
| `error` | `str` | Error message (only present when `success` is `False`) |

## Notes

- `process-all` uses the Snakemake module system internally; it does not call
  `apply_metrics()` or `apply_tract_qc()` directly.  There is no need to run
  `process-metrics` or `process-tract-qc` separately after `process-all`.
- With `--parcellation freesurfer,gif`, `process-all` serialises the
  gif→atlas-conversion→tracula-gif chain while running freesurfer stages in
  parallel where possible.
- Outputs follow the parcellation-suffixed naming convention:
  `metrics-freesurfer/`, `metrics-gif/`, `tract_qc-freesurfer/`, `tract_qc-gif/`.

## Output structure

```
{output_dir}/
└── {subject}/
    ├── gif/outputs/
    ├── recon-all/outputs/          (freesurfer parcellation only)
    ├── bamos/outputs/
    ├── dti/outputs/
    ├── noddi/outputs/
    ├── tracula-freesurfer/outputs/ (freesurfer parcellation)
    ├── tracula-gif/outputs/        (gif parcellation)
    ├── metrics-freesurfer/outputs/
    │   └── whole_brain_metrics.csv
    ├── metrics-gif/outputs/        (gif parcellation)
    │   └── whole_brain_metrics.csv
    ├── tract_qc-freesurfer/outputs/
    │   └── qc_report.csv
    ├── tract_qc-gif/outputs/       (gif parcellation)
    │   └── qc_report.csv
    └── z_scores/outputs/           (when --healthy-subjects is provided)
```

## Examples

```bash
# Single subject, default freesurfer parcellation
leukoquant process-all \
  --subject sub-001 \
  --t1    ./data/sub-001/T1/scan.nii.gz \
  --flair ./data/sub-001/FLAIR/scan.nii.gz \
  --dwi   ./data/sub-001/DWI/data.nii.gz \
  --output-dir ./outputs

# Multiple subjects, SGE cluster, freesurfer parcellation
leukoquant process-all \
  --subject ./subjects.txt \
  --t1    "./data/{subject}/T1/I*.nii.gz" \
  --flair "./data/{subject}/FLAIR/I*.nii.gz" \
  --dwi   "./data/{subject}/DWI/I*.zip" \
  --output-dir ./outputs \
  --scheduler sge

# GIF parcellation only
leukoquant process-all \
  --subject ./subjects.txt \
  --t1    "./data/{subject}/T1/I*.nii.gz" \
  --flair "./data/{subject}/FLAIR/I*.nii.gz" \
  --dwi   "./data/{subject}/DWI/I*.zip" \
  --parcellation gif \
  --output-dir ./outputs_gif \
  --scheduler sge

# Both freesurfer and gif parcellations in a single run
leukoquant process-all \
  --subject ./subjects.txt \
  --t1    "./data/{subject}/T1/I*.nii.gz" \
  --flair "./data/{subject}/FLAIR/I*.nii.gz" \
  --dwi   "./data/{subject}/DWI/I*.zip" \
  --parcellation freesurfer,gif \
  --output-dir ./outputs_both \
  --scheduler sge

# Full pipeline with Z-score computation
leukoquant process-all \
  --subject ./patients.txt \
  --healthy-subjects ./healthy.txt \
  --t1    "./data/{subject}/T1/I*.nii.gz" \
  --flair "./data/{subject}/FLAIR/I*.nii.gz" \
  --dwi   "./data/{subject}/DWI/I*.zip" \
  --demographics-csv ./data/demographics.csv \
  --covariates age,sex \
  --poly-terms age:2 \
  --output-dir ./outputs_zscore \
  --scheduler sge
```

## See also

- [`process-gif`](process_gif.md) - GIF segmentation (run independently or via process-all)
- [`process-recon`](process_recon.md) - FreeSurfer recon-all
- [`process-bamos`](process_bamos.md) - WMH segmentation
- [`process-tracula`](process_tracula.md) - TRACULA tractography
- [`process-metrics`](process_metrics.md) - standalone metrics extraction
- [`process-tract-qc`](process_tract_qc.md) - standalone tract QC
- [`process-zscore`](process_zscore.md) - standalone Z-score computation

## Source

[`leukoquant/cli/main.py`](https://github.com/stylianosc/leukoquant/blob/main/leukoquant/cli/main.py) · [`leukoquant/core/process_all_processor.py`](https://github.com/stylianosc/leukoquant/blob/main/leukoquant/core/process_all_processor.py)
