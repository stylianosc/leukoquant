# leukoquant process-recon

Run FreeSurfer `recon-all` cortical reconstruction.

## Synopsis

```
leukoquant process-recon [OPTIONS]
```

## Description

`process-recon` runs the FreeSurfer `recon-all` pipeline on a T1-weighted
structural MRI image.  The output is a full FreeSurfer subject directory
including cortical parcellation (`aparc+aseg.mgz`), surface meshes, and
volumetric segmentations.

The resulting `recon-all` outputs are used as input to `process-tracula`
(via `--freesurfer-recon-dir`) when the FreeSurfer parcellation-based
tractography is desired.  Running `process-recon` separately allows you to
share a single reconstruction across multiple downstream analyses without
re-running the computationally expensive `recon-all`.

Outputs are written to `{output_dir}/{subject}/recon-all/outputs/`.

## Parameters

| Flag | Short | Type | Default | Required | Description |
|------|-------|------|---------|----------|-------------|
| `--subject` | `-i` | `str` | - | **Yes**\* | Subject ID, or path to a `.txt` file listing one subject ID per line |
| `--t1` | `-t` | `str` | - | **Yes**\* | T1 NIfTI file or `{subject}` glob pattern |
| `--output-dir` | `-o` | `path` | - | **Yes**\* | Root output directory for workflow files |
| `--engine` | `-e` | `choice` | `snakemake` | No | Workflow engine: `snakemake` or `nextflow` |
| `--scheduler` | `-S` | `choice` | `local` | No | Execution scheduler: `local` or `sge` |
| `--cores` | `-c` | `int` | `1` | No | Number of Snakemake cores |
| `--config-yaml` | - | `path` | `None` | No | YAML file supplying defaults for any of the above - see [Config file](#config-file-config-yaml) below |
| `--force` | - | flag | `False` | No | Force-rerun recon-all even if outputs already exist |
| `--verbose` | `-v` | flag | `False` | No | Print workflow stdout to the console |

\* Required overall, but may be supplied via `--config-yaml` instead of the CLI flag - see below.

### Reconstruction engine options

| Value | Description |
|-------|-------------|
| `snakemake` | Run via Snakemake (default) |
| `nextflow` | Run via Nextflow (experimental) |

### Scheduler options

| Value | Description |
|-------|-------------|
| `local` | Run all jobs on the current machine |
| `sge` | Submit jobs to an SGE cluster |

## Config file (`--config-yaml`)

Instead of (or alongside) CLI flags, `--subject`, `--t1`, and `--output-dir`
may be supplied by a YAML file. **CLI flags always take priority** over the
equivalent YAML key when both are given - the YAML file only fills in
whatever the CLI left unset.

Recognized keys:

| YAML key | Equivalent flag |
|----------|------------------|
| `subject` (or `subject_input`) | `--subject` |
| `t1` (or `t1_pattern`) | `--t1` |
| `output_dir` | `--output-dir` |
| `scheduler` | `--scheduler` |
| `cores` | `--cores` |
| `keep_intermediate` | *(no CLI equivalent - YAML-only)* |

```yaml
# recon_config.yaml
subject: ./subjects.txt
t1: "./data/{subject}/T1/I*.nii.gz"
output_dir: ./outputs
scheduler: sge
```

```bash
leukoquant process-recon --config-yaml recon_config.yaml
```

## Returns

`dict` with keys:

| Key | Type | Description |
|-----|------|-------------|
| `success` | `bool` | `True` if the workflow completed without error |
| `results_dir` | `str` | Absolute path to the output directory |
| `error` | `str` | Error message (only present when `success` is `False`) |

## Notes

- `recon-all` is computationally intensive (~8 hours per subject on a single
  core).  Use `--scheduler sge` to dispatch to a cluster.
- The FreeSurfer outputs in `{output_dir}/{subject}/recon-all/outputs/` are
  directly compatible with the `--freesurfer-recon-dir` option of
  `process-tracula`.

## Output structure

```
{output_dir}/
└── {subject}/
    └── recon-all/
        └── outputs/
            ├── mri/
            │   ├── aparc+aseg.mgz
            │   ├── T1.mgz
            │   └── ...
            ├── surf/
            └── ...
```

## Examples

```bash
# Single subject
leukoquant process-recon \
  --subject sub-001 \
  --t1 ./data/sub-001/T1/scan.nii.gz \
  --output-dir ./outputs

# Multiple subjects on SGE
leukoquant process-recon \
  --subject ./subjects.txt \
  --t1 "./data/{subject}/T1/I*.nii.gz" \
  --output-dir ./outputs \
  --scheduler sge
```

## See also

- [`process-tracula`](process_tracula.md) - uses recon-all output via `--freesurfer-recon-dir`
- [`process-all`](process_all.md) - runs recon-all internally as part of the full pipeline

## Source

[`leukoquant/cli/main.py`](https://github.com/stylianosc/leukoquant/blob/main/leukoquant/cli/main.py) · [`leukoquant/core/freesurfer_processor.py`](https://github.com/stylianosc/leukoquant/blob/main/leukoquant/core/freesurfer_processor.py)
