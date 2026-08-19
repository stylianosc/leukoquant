# leukoquant process-zscore

Compute voxel-wise Z-score maps of target subjects relative to a healthy cohort.

## Synopsis

```
leukoquant process-zscore [OPTIONS]
```

## Description

`process-zscore` fits a general linear model (GLM) to a set of healthy control
subjects' metric maps (e.g. DTI FA), predicts the expected value and variance
at each voxel, and computes a Z-score for each target subject relative to
that normative distribution.

The GLM can optionally include demographic covariates (age, sex, etc.) and
polynomial expansion terms for non-linear relationships (e.g. quadratic age).

All images are co-registered to a common space before GLM fitting.
Registration is performed in T1 space by default; DWI space is also
supported.

Outputs are written to `{output_dir}/{subject}/z_scores/outputs/`.

!!! note "YAML configuration"
    For complex multi-metric setups, all options can be provided via
    `--config-yaml`.  CLI flags take precedence over YAML values when both
    are supplied.

## Parameters

| Flag | Short | Type | Default | Required | Description |
|------|-------|------|---------|----------|-------------|
| `--healthy-list` | `-H` | `path` | `None` | No* | Path to a text file listing healthy cohort subject IDs (one per line) |
| `--target-list` | `-T` | `path` | `None` | No* | Path to a text file listing target subject IDs (one per line) |
| `--metric` | `-m` | `str` | `None` | No (repeatable) | Metric map mapping: `name=base_dir/{subject}:glob:space`. Repeat for each metric, or comma-separate. |
| `--t1-path` | `--t1-pattern` | `str` | `None` | No | T1 image pattern: `base:glob` or `{subject}` glob |
| `--demographics-csv` | - | `path` | `None` | No | CSV file with per-subject demographics for the GLM |
| `--covariates` | - | `str` | `None` | No | Comma-separated covariate column names from `--demographics-csv` (e.g. `age,sex`) |
| `--poly-terms` | - | `str` | `None` | No | Comma-separated polynomial expansion terms (e.g. `age:2` for quadratic age) |
| `--metric-space` | - | `choice` | `t1` | No | Space in which the input metric maps are defined: `t1` or `dwi` |
| `--output-space` | - | `choice` | `t1` | No | Space in which Z-score outputs are computed: `t1` or `dwi` |
| `--dwi-pattern` | - | `str` | `None` | No | DWI image pattern. Required when `--metric-space=dwi` or `--output-space=dwi`. |
| `--bval-pattern` | - | `str` | `None` | No | bval file pattern (optional, for DWI skull-stripping) |
| `--skip-skullstrip-t1` | - | flag | `False` | No | Skip skull stripping for T1 images |
| `--skip-skullstrip-dwi` | - | flag | `False` | No | Skip skull stripping for DWI b0 images |
| `--output-dir` | `-o` | `path` | `None` | No* | Root output directory for Z-score results |
| `--config-yaml` | - | `path` | `None` | No | YAML config file providing any of the above options. CLI flags override YAML values. |
| `--scheduler` | `-s` | `choice` | `local` | No | Execution scheduler: `local` or `sge` |
| `--cores` | `-c` | `int` | `None` | No | Number of Snakemake cores (default: `1`, or from `--config-yaml`) |
| `--task-concurrency` | - | `int` | `20` | No | Max concurrently-running SGE array tasks (`-tc`). Each task writes to the shared NFS output tree on its first `mkdir`; raising this too high on `--scheduler sge` can overwhelm the NFS server. |
| `--force` | - | flag | `False` | No | Force-rerun Z-score computation even if outputs already exist |
| `--verbose` | `-v` | flag | `False` | No | Print Snakemake stdout to the console |

*`--healthy-list`, `--target-list`, and `--output-dir` may be provided via `--config-yaml` instead.

### Space options

| Value | Description |
|-------|-------------|
| `t1` | Images are in T1 (structural) space |
| `dwi` | Images are in DWI (diffusion) space |

### Scheduler options

| Value | Description |
|-------|-------------|
| `local` | Run all jobs on the current machine |
| `sge` | Submit jobs to an SGE cluster |

## Config file (`--config-yaml`)

Instead of (or alongside) CLI flags, any option above may be supplied by a
YAML file. **CLI flags always take priority** over the equivalent YAML key
when both are given - the YAML file only fills in whatever the CLI left
unset.

Recognized keys:

| YAML key | Equivalent flag |
|----------|------------------|
| `healthy_subjects_list` (or `healthy_list`) | `--healthy-list` |
| `target_subjects_list` (or `target_list`) | `--target-list` |
| `metrics` | `--metric` (as a `{name: pattern}` mapping, not the `name=pattern` string form) |
| `t1_pattern` (or `t1_path`) | `--t1-path` |
| `demographics_csv` | `--demographics-csv` |
| `covariates` | `--covariates` |
| `polynomial_terms` | `--poly-terms` |
| `metric_space` | `--metric-space` |
| `output_space` | `--output-space` |
| `dwi_pattern` (or `dwi_path`) | `--dwi-pattern` |
| `bval_pattern` (or `bval_path`) | `--bval-pattern` |
| `skip_skullstrip_t1` | `--skip-skullstrip-t1` |
| `skip_skullstrip_dwi` | `--skip-skullstrip-dwi` |
| `output_dir` | `--output-dir` |
| `scheduler` | `--scheduler` |
| `cores` | `--cores` |
| `task_concurrency` | `--task-concurrency` |

```yaml
# zscore_config.yaml
healthy_subjects_list: ./cohorts/healthy.txt
target_subjects_list: ./cohorts/patients.txt
t1_pattern: "./data/{subject}/T1/I*.nii.gz"
demographics_csv: ./data/demographics.csv
covariates: age,sex
poly_terms: "age:2"
metrics:
  dti_fa: "./dti_outputs/{subject}:fa.nii.gz:dwi"
  dti_md: "./dti_outputs/{subject}:md.nii.gz:dwi"
output_dir: ./zscores
scheduler: sge
cores: 4
```

```bash
leukoquant process-zscore --config-yaml zscore_config.yaml

# CLI flag overrides just the output directory for this run
leukoquant process-zscore --config-yaml zscore_config.yaml --output-dir ./zscores_v2
```

## Returns

`dict` with keys:

| Key | Type | Description |
|-----|------|-------------|
| `success` | `bool` | `True` if the workflow completed without error |
| `results_dir` | `str` | Absolute path to the output directory |
| `error` | `str` | Error message (only present when `success` is `False`) |

## Notes

- `--metric` mappings use the format `name=base_dir/{subject}:glob_pattern:space`.
  The `{subject}` placeholder is required in `base_dir` so the workflow can
  locate each subject's metric file independently.
- `--demographics-csv` must contain a column matching each subject ID; the
  column name can be configured via `--config-yaml`.
- Z-score outputs are NIfTI images in the specified `--output-space`.

## Output structure

```
{output_dir}/
└── {subject}/
    └── z_scores/
        └── outputs/
            ├── {metric_name}_zscore.nii.gz
            └── ...
```

## Examples

```bash
# DTI FA Z-scores with age + sex covariates (quadratic age term)
leukoquant process-zscore \
  --healthy-list ./cohorts/healthy.txt \
  --target-list  ./cohorts/patients.txt \
  --t1-path "./data/{subject}/T1/I*.nii.gz" \
  --demographics-csv ./data/demographics.csv \
  --covariates age,sex \
  --poly-terms age:2 \
  --metric "dti_fa=./dti_outputs/{subject}:fa.nii.gz:dwi" \
  --metric "dti_md=./dti_outputs/{subject}:md.nii.gz:dwi" \
  --output-dir ./zscores \
  --scheduler sge --cores 4

# Using a YAML config file
leukoquant process-zscore \
  --config-yaml ./configs/zscore_config.yaml \
  --output-dir ./zscores_override   # CLI flag overrides config
```

## See also

- [`process-dti`](process_dti.md) - produces the FA/MD maps used as metrics
- [`process-noddi`](process_noddi.md) - produces the ODI/ICVF maps used as metrics
- [`process-all`](process_all.md) - runs Z-score computation as part of the full pipeline

## Source

[`leukoquant/cli/main.py`](https://github.com/stylianosc/leukoquant/blob/main/leukoquant/cli/main.py) · [`leukoquant/core/zscore_processor.py`](https://github.com/stylianosc/leukoquant/blob/main/leukoquant/core/zscore_processor.py)
