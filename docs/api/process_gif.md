# leukoquant process-gif

Run GIF (Geodesic Information Flows) brain segmentation and multi-label parcellation.

## Synopsis

```
leukoquant process-gif [OPTIONS]
```

## Description

`process-gif` applies the GIF segmentation algorithm to structural MRI data,
producing a multi-label parcellation of the brain.  The output parcellation
NIfTI (`*_NeuroMorph_Parcellation.nii.gz`) is consumed downstream by
`process-bamos` (lesion segmentation) and `process-tracula` when
`--parcellation gif` is requested.

At least one of `--t1` or `--flair` must be provided.  When both are given,
GIF uses the T1 image as its primary input and the FLAIR is used for internal
co-registration; the parcellation is in T1 space.

Outputs are written to `{output_dir}/{subject}/gif/outputs/`.

## Parameters

| Flag | Short | Type | Default | Required | Description |
|------|-------|------|---------|----------|-------------|
| `--subject` | `-i` | `str` | - | **Yes**\*\* | Subject ID, or path to a `.txt` file listing one subject ID per line |
| `--t1` | `-t` | `str` | `None` | No* | T1 NIfTI file or `{subject}` glob pattern |
| `--flair` | `-f` | `str` | `None` | No* | FLAIR NIfTI file or `{subject}` glob pattern |
| `--mask-file` | - | `str` | `None` | No | Brain mask NIfTI or `{subject}` pattern for GIF segmentation |
| `--output-dir` | `-o` | `path` | - | **Yes**\*\* | Root output directory |
| `--scheduler` | `-s` | `choice` | `local` | No | Execution scheduler: `local` or `sge` |
| `--cores` | `-c` | `int` | `1` | No | Number of Snakemake cores |
| `--config-yaml` | - | `path` | `None` | No | YAML file supplying defaults for any of the above - see [Config file](#config-file-config-yaml) below |
| `--force` | - | flag | `False` | No | Force-rerun GIF segmentation even if outputs already exist |
| `--verbose` | `-v` | flag | `False` | No | Print Snakemake stdout to the console |

\* At least one of `--t1` or `--flair` is required.
\*\* Required overall, but may be supplied via `--config-yaml` instead of the CLI flag - see below.

### Scheduler options

| Value | Description |
|-------|-------------|
| `local` | Run all jobs on the current machine |
| `sge` | Submit jobs to an SGE cluster via the `snakemake-executor-plugin-sge` executor |

## Config file (`--config-yaml`)

Instead of (or alongside) CLI flags, `--subject`, `--t1`/`--flair`, and
`--output-dir` may be supplied by a YAML file. **CLI flags always take
priority** over the equivalent YAML key when both are given - the YAML file
only fills in whatever the CLI left unset.

Recognized keys:

| YAML key | Equivalent flag |
|----------|------------------|
| `subject` (or `subject_input`) | `--subject` |
| `t1` (or `t1_pattern`) | `--t1` |
| `flair` (or `flair_pattern`) | `--flair` |
| `mask_pattern` | `--mask-file` |
| `output_dir` | `--output-dir` |
| `scheduler` | `--scheduler` |
| `cores` | `--cores` |
| `flair_db` | *(no CLI equivalent - YAML-only)* |
| `keep_intermediate` | *(no CLI equivalent - YAML-only)* |

```yaml
# gif_config.yaml
subject: ./subjects.txt
t1: "./data/{subject}/T1/I*.nii.gz"
output_dir: ./outputs
scheduler: sge
cores: 4
```

```bash
leukoquant process-gif --config-yaml gif_config.yaml

# override the scheduler for this run only; subject/t1/output-dir still come from the YAML
leukoquant process-gif --config-yaml gif_config.yaml --scheduler local
```

## Returns

`dict` with keys:

| Key | Type | Description |
|-----|------|-------------|
| `success` | `bool` | `True` if the workflow completed without error |
| `results_dir` | `str` | Absolute path to the output directory |
| `error` | `str` | Error message (only present when `success` is `False`) |

## Notes

- Path patterns support `{subject}` as a per-subject placeholder: e.g.
  `./data/{subject}/T1/I*.nii.gz` is resolved independently for each subject.
- The parcellation filename stem matches the T1 image filename stem used as
  GIF's primary input.  Downstream commands that auto-resolve the GIF output
  use the glob `*_NeuroMorph_Parcellation.nii.gz`.
- When `--scheduler sge` is used, `leukoquant` blocks until all submitted SGE
  jobs complete before returning.

## Output structure

```
{output_dir}/
└── {subject}/
    └── gif/
        └── outputs/
            ├── {stem}_NeuroMorph_Parcellation.nii.gz
            ├── {stem}_NeuroMorph_Segmentation.nii.gz
            ├── {stem}_NeuroMorph_Brain.nii.gz
            └── ...
```

## Examples

```bash
# Single subject, T1 only
leukoquant process-gif \
  --subject sub-001 \
  --t1 ./data/sub-001/T1/scan.nii.gz \
  --output-dir ./outputs

# Multiple subjects, FLAIR only, SGE cluster
leukoquant process-gif \
  --subject ./subjects.txt \
  --flair ./data/{subject}/FLAIR/I*.nii.gz \
  --output-dir ./outputs \
  --scheduler sge

# Multiple subjects, T1 + FLAIR (T1 is primary input)
leukoquant process-gif \
  --subject ./subjects.txt \
  --t1    ./data/{subject}/T1/I*.nii.gz \
  --flair ./data/{subject}/FLAIR/I*.nii.gz \
  --output-dir ./outputs \
  --scheduler sge
```

## See also

- [`process-bamos`](process_bamos.md) - WMH segmentation that consumes GIF output
- [`process-tracula`](process_tracula.md) - tractography with optional GIF parcellation
- [`process-atlas-conversion`](process_atlas_conversion.md) - convert GIF labels to FreeSurfer space
- [`process-all`](process_all.md) - runs GIF as part of the full pipeline

## Source

[`leukoquant/cli/main.py`](https://github.com/stylianosc/leukoquant/blob/main/leukoquant/cli/main.py) · [`leukoquant/core/gif_processor.py`](https://github.com/stylianosc/leukoquant/blob/main/leukoquant/core/gif_processor.py)
