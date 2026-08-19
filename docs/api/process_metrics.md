# leukoquant process-metrics

Extract lesion-informed quantitative metrics along white matter tracts.

## Synopsis

```
leukoquant process-metrics [OPTIONS]
```

## Description

`process-metrics` samples voxel-wise microstructure maps (DTI, NODDI, custom)
along probabilistic white matter tracts, combined with lesion masks to produce
per-tract summary statistics.  The primary output is a CSV file containing
one row per subject with columns for each tract × metric combination.

The command operates in three **tract modes**:

| Mode | Description |
|------|-------------|
| `tractography` | Metrics are sampled along the raw probabilistic density maps from TRACULA |
| `tractography-atlas` | Tractography × atlas ROI intersection (default) |
| `atlas` | Metrics are sampled using atlas-defined ROIs only (no tractography required) |

### Tractography path format

`--tractography-path` uses the composite `base:glob:space` format:

```
base_dir/{subject}/tracula-freesurfer/outputs:dpath/*/path.pd.nii.gz:dwi
```

The path must include the parcellation subfolder (`tracula-freesurfer` or
`tracula-gif`) introduced in LeukoQuant ≥ 0.2.0.

Outputs are written to `{output_dir}/{subject}/metrics/outputs/`.

## Parameters

| Flag | Short | Type | Default | Required | Description |
|------|-------|------|---------|----------|-------------|
| `--subject` | - | `str` | - | **Yes**\* | Subject ID, or path to a `.txt` file listing one subject ID per line |
| `--tractography-path` | - | `str` | - | **Yes**\* | Tractography pattern: `base:glob:space`. Example: `./tracula/{subject}/tracula-freesurfer/outputs:dpath/*/path.pd.nii.gz:dwi` |
| `--t1-path` | - | `str` | `None` | No | T1 pattern as `base:glob` or `{subject}` glob. Used for registration. |
| `--dwi-path` | - | `str` | `None` | No | DWI pattern as `base:glob`. Used when registration to DWI space is needed. |
| `--lesion-path` | - | `str` | `None` | No | Lesion mask pattern(s). Format: `[name=]base:glob[:space]`, comma-separated for multiple lesion types. |
| `--maps` | `-m` | `str` | `None` | No (repeatable) | Metric map mappings: `metric_name=base:glob:space`. Repeat the flag for each map, or comma-separate multiple maps in one value. |
| `--tract-mode` | - | `choice` | `tractography-atlas` | No | Tract analysis mode: `tractography`, `tractography-atlas`, or `atlas` |
| `--parcellation` | - | `str` | `freesurfer` | No | Parcellation(s) to extract metrics for. Comma-separated for multiple (e.g. `freesurfer,gif`) |
| `--output-dir` | `-o` | `path` | - | **Yes**\* | Root output directory for metrics CSV files |
| `--scheduler` | - | `choice` | `local` | No | Execution scheduler: `local` or `sge` |
| `--cores` | `-c` | `int` | `1` | No | Number of Snakemake cores |
| `--config-yaml` | - | `path` | `None` | No | YAML file supplying defaults for any of the above - see [Config file](#config-file-config-yaml) below |
| `--force` | - | flag | `False` | No | Force-rerun metrics extraction (for `--parcellation`) even if outputs already exist |
| `--verbose` | `-v` | flag | `False` | No | Print Snakemake stdout to the console |

\* Required overall, but may be supplied via `--config-yaml` instead of the CLI flag - see below.

### Tract mode options

| Value | Description |
|-------|-------------|
| `tractography` | Sample metrics directly along probabilistic tract density maps |
| `tractography-atlas` | Intersection of tract density maps with atlas parcellation (default) |
| `atlas` | Atlas-only ROI-based sampling; no tractography files required |

### Lesion path format

```
[name=]base_dir:glob_pattern[:space]
```

Multiple lesion types can be passed as a comma-separated list in a single
`--lesion-path` value:

```
wmh=./bamos/{subject}:bamos/CorrectLesion_*.nii.gz:t1,lacunes=./lac/{subject}:lac_*.nii.gz:t1
```

When `name` is omitted, it defaults to `lesion`.

### Maps format

```
metric_name=base_dir:glob_pattern:space
```

```bash
--maps "dti_fa=./dti_outputs:fa.nii.gz:dwi" \
--maps "dti_md=./dti_outputs:md.nii.gz:dwi" \
--maps "noddi_odi=./noddi_outputs:odi.nii.gz:dwi"
```

### Space options

| Value | Description |
|-------|-------------|
| `t1` | Image is in T1 (structural) space |
| `dwi` | Image is in DWI (diffusion) space |
| `atlas` | Image is in atlas space |

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
| `t1_path` | `--t1-path` |
| `dwi_path` | `--dwi-path` |
| `lesion_path` | `--lesion-path` |
| `metrics` | `--maps` (as a `{name: pattern}` mapping, not the `name=pattern` string form) |
| `tract_mode` | `--tract-mode` |
| `parcellation` | `--parcellation` |
| `output_dir` | `--output-dir` |
| `scheduler` | `--scheduler` |
| `cores` | `--cores` |
| `keep_intermediate` | *(no CLI equivalent - YAML-only)* |

```yaml
# metrics_config.yaml
subject: ./subjects.txt
tractography_path: "$TRACULA_OUT/{subject}/tracula-freesurfer/outputs:dpath/*/path.pd.nii.gz:dwi"
t1_path: "./data/{subject}/T1/I*.nii.gz"
lesion_path: "wmh=$BAMOS_OUT:bamos/CorrectLesion_*.nii.gz:t1"
metrics:
  dti_fa: "$DTI_OUT:fa.nii.gz:dwi"
  dti_md: "$DTI_OUT:md.nii.gz:dwi"
output_dir: ./metrics_outputs
scheduler: sge
cores: 4
```

```bash
leukoquant process-metrics --config-yaml metrics_config.yaml

# override just the output directory for this run; everything else still comes from the YAML
leukoquant process-metrics --config-yaml metrics_config.yaml --output-dir ./metrics_outputs_v2
```

## Returns

`dict` with keys:

| Key | Type | Description |
|-----|------|-------------|
| `success` | `bool` | `True` if the workflow completed without error |
| `results_dir` | `str` | Absolute path to the output directory |
| `error` | `str` | Error message (only present when `success` is `False`) |

## Notes

- `--tractography-path` must point to the parcellation-specific subfolder
  produced by `process-tracula` (e.g. `tracula-freesurfer/outputs` or
  `tracula-gif/outputs`).  Using the bare `--output-dir` of an older run
  without the subfolder will fail to locate tract files.
- Standalone `process-metrics` always writes to `metrics/` (no parcellation
  suffix).  The suffixed folder `metrics-{parcellation}/` is used only when
  `process-all` orchestrates metrics internally.

## Output structure

```
{output_dir}/
└── {subject}/
    └── metrics/
        └── outputs/
            ├── whole_brain_metrics.csv
            └── ...
```

## Examples

```bash
# Multi-subject: freesurfer tractography + DTI maps + WMH lesion mask
leukoquant process-metrics \
  --subject ./subjects.txt \
  --tractography-path "$TRACULA_OUT/{subject}/tracula-freesurfer/outputs:dpath/*/path.pd.nii.gz:dwi" \
  --t1-path  "./data/{subject}/T1/I*.nii.gz" \
  --lesion-path "wmh=$BAMOS_OUT:bamos/CorrectLesion_*.nii.gz:t1" \
  --maps "dti_fa=$DTI_OUT:fa.nii.gz:dwi" \
  --maps "dti_md=$DTI_OUT:md.nii.gz:dwi" \
  --tract-mode tractography-atlas \
  --output-dir ./metrics_outputs \
  --scheduler sge --cores 4

# GIF tractography variant
leukoquant process-metrics \
  --subject ./subjects.txt \
  --tractography-path "$GIF_OUT/{subject}/tracula-gif/outputs:dpath/*/path.pd.nii.gz:dwi" \
  --t1-path  "./data/{subject}/T1/I*.nii.gz" \
  --lesion-path "wmh=$BAMOS_OUT:bamos/CorrectLesion_*.nii.gz:t1" \
  --maps "dti_fa=$DTI_OUT:fa.nii.gz:dwi" \
  --output-dir ./metrics_gif \
  --scheduler sge --cores 4
```

## See also

- [`process-tracula`](process_tracula.md) - produces the tractography input
- [`process-bamos`](process_bamos.md) - produces the lesion mask for `--lesion-path`
- [`process-dti`](process_dti.md) - produces FA/MD maps for `--maps`
- [`process-noddi`](process_noddi.md) - produces ODI/ICVF maps for `--maps`
- [`process-tract-qc`](process_tract_qc.md) - quality-control companion for tractography
- [`process-all`](process_all.md) - runs metrics extraction as part of the full pipeline

## Source

[`leukoquant/cli/main.py`](https://github.com/stylianosc/leukoquant/blob/main/leukoquant/cli/main.py) · [`leukoquant/core/metrics_processor.py`](https://github.com/stylianosc/leukoquant/blob/main/leukoquant/core/metrics_processor.py)
