# leukoquant process-tract-qc

Run tractography quality control.

## Synopsis

```
leukoquant process-tract-qc [OPTIONS]
```

## Description

`process-tract-qc` inspects TRACULA probabilistic tract outputs and produces
a per-subject, per-tract QC report in CSV format.  The report includes
summary statistics (e.g. tract volume, mean path probability density) that
allow identification of failed or low-quality tract reconstructions before
running `process-metrics`.

The command accepts the same `base:glob:space` tractography path format as
`process-metrics`, including the parcellation-specific subfolder:

```
./tracula_outputs/{subject}/tracula-freesurfer/outputs:dpath/*/path.pd.nii.gz:dwi
```

Outputs are written to `{output_dir}/{subject}/tract_qc/outputs/`.

## Parameters

| Flag | Short | Type | Default | Required | Description |
|------|-------|------|---------|----------|-------------|
| `--subject` | - | `str` | - | **Yes**\* | Subject ID, or path to a `.txt` file listing one subject ID per line |
| `--tractography-path` | - | `str` | - | **Yes**\* | Tractography pattern: `base:glob:space`. Must point to the parcellation-specific subfolder (e.g. `tracula-freesurfer/outputs`). |
| `--output-dir` | `-o` | `path` | `./` | No | Root output directory for tract QC results (defaults to the current directory) |
| `--scheduler` | - | `choice` | `local` | No | Execution scheduler: `local` or `sge` |
| `--cores` | `-c` | `int` | `1` | No | Number of Snakemake cores |
| `--config-yaml` | - | `path` | `None` | No | YAML file supplying defaults for any of the above - see [Config file](#config-file-config-yaml) below |
| `--force` | - | flag | `False` | No | Force-rerun tract QC even if outputs already exist |
| `--verbose` | `-v` | flag | `False` | No | Print Snakemake stdout to the console |

\* Required overall, but may be supplied via `--config-yaml` instead of the CLI flag - see below.

### Scheduler options

| Value | Description |
|-------|-------------|
| `local` | Run all jobs on the current machine |
| `sge` | Submit jobs to an SGE cluster |

## Config file (`--config-yaml`)

Instead of (or alongside) CLI flags, `--subject`, `--tractography-path`, and
`--output-dir` may be supplied by a YAML file. **CLI flags always take
priority** over the equivalent YAML key when both are given - the YAML file
only fills in whatever the CLI left unset.

Recognized keys:

| YAML key | Equivalent flag |
|----------|------------------|
| `subject` | `--subject` |
| `tractography_path` | `--tractography-path` |
| `output_dir` | `--output-dir` |
| `scheduler` | `--scheduler` |
| `cores` | `--cores` |

```yaml
# tract_qc_config.yaml
subject: ./subjects.txt
tractography_path: "$TRACULA_OUT/{subject}/tracula-freesurfer/outputs:dpath/*/path.pd.nii.gz:dwi"
output_dir: ./tract_qc_outputs
scheduler: sge
cores: 4
```

```bash
leukoquant process-tract-qc --config-yaml tract_qc_config.yaml
```

## Returns

`dict` with keys:

| Key | Type | Description |
|-----|------|-------------|
| `success` | `bool` | `True` if the workflow completed without error |
| `results_dir` | `str` | Absolute path to the output directory |
| `error` | `str` | Error message (only present when `success` is `False`) |

## Notes

- `--tractography-path` must include the parcellation-specific output subfolder
  (`tracula-freesurfer/outputs` or `tracula-gif/outputs`).
- Standalone `process-tract-qc` always writes to `tract_qc/` (no parcellation
  suffix).  The suffixed folder `tract_qc-{parcellation}/` is used only when
  `process-all` orchestrates QC internally.
- Run tract QC before `process-metrics` to identify problematic subjects whose
  metrics should be excluded from downstream analysis.

## Output structure

```
{output_dir}/
└── {subject}/
    └── tract_qc/
        ├── outputs/
        │   └── qc_report.csv
        └── logs/
            ├── tract_qc_log.txt
            └── tract_qc_error.txt
```

## Examples

```bash
# Multi-subject, freesurfer tractography
leukoquant process-tract-qc \
  --subject ./subjects.txt \
  --tractography-path "$TRACULA_OUT/{subject}/tracula-freesurfer/outputs:dpath/*/path.pd.nii.gz:dwi" \
  --output-dir ./tract_qc_outputs \
  --scheduler sge --cores 4

# Multi-subject, gif tractography
leukoquant process-tract-qc \
  --subject ./subjects.txt \
  --tractography-path "$GIF_OUT/{subject}/tracula-gif/outputs:dpath/*/path.pd.nii.gz:dwi" \
  --output-dir ./tract_qc_gif \
  --scheduler sge --cores 4
```

## See also

- [`process-tracula`](process_tracula.md) - produces the tractography input
- [`process-metrics`](process_metrics.md) - computes tract metrics after QC passes
- [`process-all`](process_all.md) - runs tract QC as part of the full pipeline

## Source

[`leukoquant/cli/main.py`](https://github.com/stylianosc/leukoquant/blob/main/leukoquant/cli/main.py) · [`leukoquant/core/tract_qc_processor.py`](https://github.com/stylianosc/leukoquant/blob/main/leukoquant/core/tract_qc_processor.py)
